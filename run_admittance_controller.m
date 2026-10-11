function varargout = run_admittance_controller(mode, varargin)
%RUN_ADMITTANCE_CONTROLLER Unified entry for the four formal controller laws.
% mode:
%   'fixed'      : fixed nominal admittance with the shared supervisory shell
%   'zhang2018'  : Zhang-2018 literature-derived B/K scheduling comparator
%   'cao2024'    : Cao-2024-derived cycle-domain momentum B/K comparator
%   'proposed'   : proposed IE-AAC
%
% For backward compatibility, the legacy mode 'literature' maps to Zhang-2018.
    switch lower(char(mode))
        case 'fixed'
            [varargout{1:nargout}] = fixed_controller_core(varargin{:});
        case {'zhang2018','zhang','literature_zhang','literature'}
            [varargout{1:nargout}] = zhang2018_controller_core(varargin{:});
        case {'cao2024','cao','literature_cao'}
            [varargout{1:nargout}] = cao2024_controller_core(varargin{:});
        case 'proposed'
            [varargout{1:nargout}] = proposed_controller_core(varargin{:});
        otherwise
            error('Unknown controller mode: %s', char(mode));
    end
end
function [delta_q, delta_dq, delta_ddq, Ba_fixed, Ka_fixed, Ga_fixed, ...
          quiet_zone_flag, e_fb_used_all, e_fb_raw_all, interaction_diag] = ...
    fixed_controller_core(q_ref, dq_ref, params, pre_fb, t_all, ddq_ref)

    N = size(q_ref, 1);
    dt = params.dt;
    tic;

    %% ================= 输入验证 =================
    assert(size(q_ref,2)==3 && size(dq_ref,2)==3, ...
           '输入轨迹必须为 N x 3 矩阵');

    %% ================= 参数提取 =================
    Ma0 = params.Ma0;
    Ba0 = params.Ba0;
    Ka0 = params.Ka0;
    Ga0 = params.Ga0;

    delta_q_limit   = params.delta_q_limit;
    delta_dq_limit  = params.delta_dq_limit;
    delta_ddq_limit = params.delta_ddq_limit;
    if isscalar(delta_q_limit),   delta_q_limit   = repmat(delta_q_limit,  1, 3); end
    if isscalar(delta_dq_limit),  delta_dq_limit  = repmat(delta_dq_limit, 1, 3); end
    if isscalar(delta_ddq_limit), delta_ddq_limit = repmat(delta_ddq_limit,1, 3); end

    jerk_limit = params.jerk_limit;
    if isscalar(jerk_limit), jerk_limit = repmat(jerk_limit, 1, 3); end
    K_mis = params.K_mismatch_comp;
    if isscalar(K_mis), K_mis = repmat(K_mis, 1, 3); end

    % 关节极限（与自适应导纳保持一致）
    q1_max = params.q1_max; q2_max = params.q2_max;
    if isfield(params,'q3_max'), q3_max = params.q3_max; else q3_max = 0.5; end
    q_max = [q1_max, q2_max, q3_max];

    %% ================= 预分配 =================
    delta_q   = zeros(N, 3);
    delta_dq  = zeros(N, 3);
    delta_ddq = zeros(N, 3);
    quiet_zone_flag = zeros(N, 3);
    e_fb_used_all   = zeros(N, 3);
    e_fb_raw_all    = zeros(N, 3);

    % ================= 六分量交互力矩诊断（仅记录，不参与控制） =================
    % raw = mismatch + spasm + synergy + fatigue + event + noise
    interaction_diag = init_interaction_diag_fixed_controller_core(N);

    % 固定参数序列（全时段保持常量）
    Ba_fixed = repmat(Ba0(:)', N, 1);
    Ka_fixed = repmat(Ka0(:)', N, 1);
    Ga_fixed = repmat(Ga0(:)', N, 1);

    %% ================= 主循环 =================
    for i = 1:N-1
        i_prev = max(i-1, 1);   % 上一步索引（用于命令加速度估计）
        for j = 1:3
            % ---- 当前命令轨迹（标量） ----
            q_cmd_cur  = q_ref(i,j) + delta_q(i,j);
            dq_cmd_cur = dq_ref(i,j) + delta_dq(i,j);

            % ---- 实时反馈计算（与 adaptive_admittance 完全一致） ----
            B0 = params.B0(j);
            Kq = params.Kq_fb(j);
            Kv = params.Kv_fb(j);
            Ke = params.Ke_fb(j);

            % ① 参考基线（不随控制器行为漂移）
            phi_ref_j = min((abs(q_ref(i,j)) / q_max(j))^2, 1.0);
            fb_nominal_ref = B0 + Kq*abs(q_ref(i,j)) + Kv*abs(dq_ref(i,j)) + Ke*phi_ref_j;

            % ② 命令基线（实际命令状态下的正常阻抗）
            phi_cmd_j = min((abs(q_cmd_cur) / q_max(j))^2, 1.0);
            fb_nominal_cmd = B0 + Kq*abs(q_cmd_cur) + Kv*abs(dq_cmd_cur) + Ke*phi_cmd_j;

            % ③ 痉挛：门控使用实际命令角度(平滑sigmoid)，强度使用参考速度(避免正反馈)
            gate_st = params.spasm_gate_steepness;
            if j == 1  % 髋
                gate = 1 ./ (1 + exp(-gate_st * max(min(q_cmd_cur - params.spasm_q_hip_th, 5), -5)));
                spasm = (params.spasm_c1/100) * abs(dq_ref(i,j)) * gate;
            elseif j == 2  % 膝（耦合项改用髋参考速度）
                gate = 1 ./ (1 + exp(gate_st * max(min(q_cmd_cur + params.spasm_q_knee_th, 5), -5)));
                spasm = (params.spasm_c2/100) * abs(dq_ref(i,j) + 0.15*dq_ref(i,1)) * gate;
            else  % 踝（耦合项改用膝参考速度）
                spasm = params.spasm_c3 * abs(dq_ref(i,j) + 0.20*dq_ref(i,2));
            end

            % ④ 协同（其他关节使用命令状态 q_cmd = q_ref + delta_q，单步延迟耦合）
            if j == 1
                synergy = params.Csyn(1) * abs(q_ref(i,2) + delta_q(i,2));
            elseif j == 2
                synergy = params.Csyn(2) * abs(q_ref(i,1) + delta_q(i,1));
            else
                synergy = params.Csyn(3) * abs(q_ref(i,2) + delta_q(i,2));
            end

            % 疲劳、事件、噪声（预生成）
            fatigue = pre_fb.fatigue_slope(j) * params.cfatigue * pre_fb.fatigue_trend(i);
            event   = pre_fb.event_all(i,j);
            noise   = pre_fb.noise_all(i,j);

            % 真实反馈
            fb_true = fb_nominal_cmd + spasm + synergy + fatigue + event + noise;

            % ⑤ 原始误差：相对“参考基线”，显式包含阻抗失配项
            e_fb_raw = fb_true - fb_nominal_ref;

            % ⑥ 阻抗失配前馈补偿（防止 Kv*Δdq 自激正反馈）
            mismatch = fb_nominal_cmd - fb_nominal_ref;

            e_fb_comp = e_fb_raw - K_mis(j) * mismatch;

            % ---- 六分量诊断记录：此处只记录，不改变闭环计算 ----
            pathological = spasm + synergy + fatigue + event + noise;
            interaction_diag.mismatch(i,j) = mismatch;
            interaction_diag.spasm(i,j) = spasm;
            interaction_diag.synergy(i,j) = synergy;
            interaction_diag.fatigue(i,j) = fatigue;
            interaction_diag.event(i,j) = event;
            interaction_diag.noise(i,j) = noise;
            interaction_diag.pathological(i,j) = pathological;
            interaction_diag.compensated(i,j) = e_fb_comp;
            interaction_diag.raw(i,j) = e_fb_raw;
            interaction_diag.fb_nominal_ref(i,j) = fb_nominal_ref;
            interaction_diag.fb_nominal_cmd(i,j) = fb_nominal_cmd;
            interaction_diag.fb_true(i,j) = fb_true;

            % ⑦ 参考加速度（优先使用解析 ddq_ref）
            if nargin >= 6 && ~isempty(ddq_ref)
                ddq_ref_i = ddq_ref(i,j);
            else
                if i == 1
                    ddq_ref_i = (-3*dq_ref(i,j) + 4*dq_ref(i+1,j) - dq_ref(i+2,j)) / (2*dt);
                elseif i == N
                    ddq_ref_i = (3*dq_ref(i,j) - 4*dq_ref(i-1,j) + dq_ref(i-2,j)) / (2*dt);
                else
                    ddq_ref_i = (dq_ref(i+1,j) - dq_ref(i-1,j)) / (2*dt);
                end
            end
            ddq_cmd_est = ddq_ref_i + delta_ddq(i_prev, j);

            % ⑧ 死区 + 静止区处理（判定使用命令状态与命令加速度）
            [e_fb_use, quiet_flag] = apply_deadzone_and_quietzone(e_fb_comp, q_cmd_cur, dq_cmd_cur, ...
                ddq_cmd_est, j, params);
            quiet_zone_flag(i,j) = quiet_flag;
            e_fb_used_all(i,j)   = e_fb_use;
            e_fb_raw_all(i,j)    = e_fb_raw;
            interaction_diag.used(i,j) = e_fb_use;

            % ---- 固定参数（静止区缩放） ----
            Ba_tmp = Ba0(j);
            Ka_tmp = Ka0(j);
            Ga_tmp = Ga0(j);
            if quiet_flag
                Ga_tmp = params.quiet_gain_scale(j) * Ga_tmp;
                Ba_tmp = params.quiet_damp_scale(j)  * Ba_tmp;
            end
            % 二次限幅（防止缩放后越界，与自适应逻辑一致）
            Ba_tmp = min(max(Ba_tmp, 0.1*Ba0(j)), 4.0*Ba0(j));
            Ka_tmp = min(max(Ka_tmp, 0.01*Ka0(j)), 5.0*Ka0(j));
            Ga_tmp = min(max(Ga_tmp, 0.08*Ga0(j)), 2.5*Ga0(j));

            % ---- 固定导纳方程（后向欧拉） ----
            delta_ddq_tmp = (Ga_tmp * e_fb_use ...
                           - Ba_tmp * delta_dq(i,j) ...
                           - Ka_tmp * delta_q(i,j)) / Ma0(j);
            % 加速度限幅
            delta_ddq_tmp = max(min(delta_ddq_tmp, delta_ddq_limit(j)), -delta_ddq_limit(j));

            % ★ 急动度限幅：限制每步加速度变化率，消除速度/加速度尖峰
            if i > 1
                ddq_prev = delta_ddq(i-1, j);
            else
                ddq_prev = 0;
            end
            max_step = jerk_limit(j) * dt;
            delta_ddq_tmp = max(min(delta_ddq_tmp, ddq_prev + max_step), ddq_prev - max_step);
            delta_ddq(i,j) = delta_ddq_tmp;

            % ---- 积分与限幅 ----
            dq_pred = delta_dq(i,j) + delta_ddq_tmp * dt;
            dq_pred = max(min(dq_pred, delta_dq_limit(j)), -delta_dq_limit(j));

            q_pred = delta_q(i,j) + dq_pred * dt;
            q_pred = max(min(q_pred, delta_q_limit(j)), -delta_q_limit(j));

            delta_dq(i+1,j) = dq_pred;
            delta_q(i+1,j)  = q_pred;
        end
    end

    %% ================= 末尾补值 =================
    delta_q(end,:)   = delta_q(end-1,:);
    delta_dq(end,:)  = delta_dq(end-1,:);
    delta_ddq(end,:) = delta_ddq(end-1,:);
    quiet_zone_flag(end,:) = quiet_zone_flag(end-1,:);
    e_fb_used_all(end,:)   = e_fb_used_all(end-1,:);
    e_fb_raw_all(end,:)    = e_fb_raw_all(end-1,:);

    % 六分量诊断末尾补齐，并执行代数一致性检查
    interaction_diag = finalize_interaction_diag_fixed_controller_core(interaction_diag, params);

    exec_time = toc;
    fprintf('固定导纳（闭环反馈+失配补偿）运行完成：QuietZone 平均激活率 %.2f%% | 执行时间 %.3f s\n', ...
            mean(quiet_zone_flag, 'all')*100, exec_time);
end

function diag = init_interaction_diag_fixed_controller_core(N)
%INIT_INTERACTION_DIAG 预分配六分量诊断结构体。
    fields = {'mismatch','spasm','synergy','fatigue','event','noise', ...
              'pathological','compensated','used','raw', ...
              'fb_nominal_ref','fb_nominal_cmd','fb_true'};
    diag = struct();
    for k = 1:numel(fields)
        diag.(fields{k}) = zeros(N,3);
    end
end

function diag = finalize_interaction_diag_fixed_controller_core(diag, params)
%FINALIZE_INTERACTION_DIAG 补齐最后采样点，并检查六分量重构关系。
    fields = fieldnames(diag);
    for k = 1:numel(fields)
        x = diag.(fields{k});
        if size(x,1) >= 2
            x(end,:) = x(end-1,:);
        end
        diag.(fields{k}) = x;
    end

    raw_reconstructed = diag.mismatch + diag.spasm + diag.synergy + ...
                        diag.fatigue + diag.event + diag.noise;
    err_raw = max(abs(diag.raw(:) - raw_reconstructed(:)));
    assert(err_raw < 1e-10, ...
        'fixed_admittance:InteractionDecomposition', ...
        '六分量无法重构 raw interaction residual，最大误差 %.3e', err_raw);

    Kmis = params.K_mismatch_comp;
    if isscalar(Kmis), Kmis = repmat(Kmis,1,3); end
    if max(abs(Kmis(:)-1)) < 1e-12
        err_comp = max(abs(diag.compensated(:) - diag.pathological(:)));
        assert(err_comp < 1e-10, ...
            'fixed_admittance:MismatchCompensation', ...
            'K_mismatch_comp=1 时 compensated 应等于 pathological，最大误差 %.3e', err_comp);
    end
end

function [delta_q, delta_dq, delta_ddq, Ba_lit, Ka_lit, Ga_lit, quiet_zone_flag, e_fb_used_all, e_fb_raw_all, interaction_diag] = ...
    zhang2018_controller_core(q_ref, dq_ref, params, pre_fb, t_all, ddq_ref)
%ZHANG2018_CONTROLLER_CORE Zhang-2018 literature-derived variable-parameter admittance baseline.
%
% Source law adapted from:
%   M. Zhang, S. Q. Xie, X. Li, G. Zhu, W. Meng, X. Huang, and A. Veale,
%   "Adaptive Patient-Cooperative Control of a Compliant Ankle
%   Rehabilitation Robot (CARR) with Enhanced Training Safety,"
%   IEEE Trans. Ind. Electron., vol. 65, no. 2, pp. 1398-1407, 2018,
%   doi:10.1109/TIE.2017.2733425.
%
% The cited controller keeps the admittance inertia fixed and schedules
% damping B and stiffness K from joint posture and the ratio between a
% model-predicted passive torque and the measured patient-robot torque.
% This file ports that parameter-scheduling law into the present joint-space
% three-DOF simulation while preserving the project's common closed-loop
% patient model and supervisory shell.
%
% Fair-comparison rule used here:
%   - Same reference, patient model, pathology, event/noise realization,
%     mismatch compensation, dead-zone, quiet-zone, projection, motion and
%     jerk limits as Fixed and Proposed controllers.
%   - Only the admittance-parameter scheduling law is changed.
%   - M_a and G_a remain nominal. B_a and K_a follow the Zhang et al. law.
%
% Signal correspondence for the port:
%   theta        -> controller-specific q_cmd
%   T_modeling   -> nominal passive interaction at q_cmd (fb_nominal_cmd)
%   T_measured   -> total simulated interaction torque (fb_true)
%
% The source paper states that its weighting coefficients were selected
% empirically/trial-and-error. Therefore the dimensionless weights in
% params.literature_zhang2018 are calibrated only from deterministic nominal
% behavior and then frozen before Monte Carlo evaluation.

    N = size(q_ref, 1);
    dt = params.dt;
    tic;

    %% ================= 参数提取 =================
    Ma0 = params.Ma0; Ba0 = params.Ba0; Ka0 = params.Ka0; Ga0 = params.Ga0;

    assert(isfield(params,'literature_zhang2018'), ...
        'zhang2018_variable_admittance:MissingConfig', ...
        'params.literature_zhang2018 未定义；请在 parameters.m 中配置文献型 VA 参数。');
    lit = params.literature_zhang2018;
    a1 = lit.a1; a2 = lit.a2; a3 = lit.a3;
    b1 = lit.b1; b2 = lit.b2; b3 = lit.b3;
    torque_floor = lit.torque_floor;
    names = {'a1','a2','a3','b1','b2','b3'};
    vals = {a1,a2,a3,b1,b2,b3};
    for kk = 1:numel(vals)
        if isscalar(vals{kk}), vals{kk} = repmat(vals{kk},1,3); end
        assert(numel(vals{kk}) == 3 && all(isfinite(vals{kk})) && all(vals{kk} > 0), ...
            'literature_zhang2018.%s 必须为正标量或 3 元素正向量', names{kk});
    end
    [a1,a2,a3,b1,b2,b3] = vals{:};
    assert(isscalar(torque_floor) && isfinite(torque_floor) && torque_floor > 0, ...
        'literature_zhang2018.torque_floor 必须为正标量');

    delta_q_limit   = params.delta_q_limit;
    delta_dq_limit  = params.delta_dq_limit;
    delta_ddq_limit = params.delta_ddq_limit;
    if isscalar(delta_q_limit),   delta_q_limit   = repmat(delta_q_limit,  1, 3); end
    if isscalar(delta_dq_limit),  delta_dq_limit  = repmat(delta_dq_limit, 1, 3); end
    if isscalar(delta_ddq_limit), delta_ddq_limit = repmat(delta_ddq_limit,1, 3); end

    jerk_limit = params.jerk_limit;
    if isscalar(jerk_limit), jerk_limit = repmat(jerk_limit, 1, 3); end
    K_mis = params.K_mismatch_comp;
    if isscalar(K_mis), K_mis = repmat(K_mis, 1, 3); end

    % 关节极限（用于边界接近度 φ(q)）
    q1_max = params.q1_max; q2_max = params.q2_max;
    if isfield(params,'q3_max'), q3_max = params.q3_max; else q3_max = 0.5; end
    q_max = [q1_max, q2_max, q3_max];

    %% ================= 预分配 =================
    delta_q   = zeros(N, 3);
    delta_dq  = zeros(N, 3);
    delta_ddq = zeros(N, 3);
    Ba_lit = zeros(N, 3); Ka_lit = zeros(N, 3); Ga_lit = zeros(N, 3);

    quiet_zone_flag = zeros(N, 3);
    e_fb_used_all   = zeros(N, 3);
    e_fb_raw_all    = zeros(N, 3);

    % ================= 六分量交互力矩诊断（仅记录，不参与控制） =================
    % raw = mismatch + spasm + synergy + fatigue + event + noise
    interaction_diag = init_interaction_diag_zhang2018_controller_core(N);

    %% ================= 主循环 =================
    for i = 1:N-1
        i_prev = max(i-1, 1);   % 上一步索引（用于上一步修正量）
        for j = 1:3
            % ---- 当前命令轨迹（标量） ----
            q_cmd_cur  = q_ref(i,j)   + delta_q(i,j);
            dq_cmd_cur = dq_ref(i,j)  + delta_dq(i,j);

            % ---- 实时反馈计算（基于 q_cmd 闭环） ----
            B0 = params.B0(j);
            Kq = params.Kq_fb(j);
            Kv = params.Kv_fb(j);
            Ke = params.Ke_fb(j);

            % ① 参考基线（健康参考状态下的正常阻抗，不随控制器行为漂移）
            phi_ref_j = min((abs(q_ref(i,j)) / q_max(j))^2, 1.0);
            fb_nominal_ref = B0 + Kq*abs(q_ref(i,j)) + Kv*abs(dq_ref(i,j)) + Ke*phi_ref_j;

            % ② 命令基线（机器人实际命令状态下的正常阻抗）
            phi_cmd_j = min((abs(q_cmd_cur) / q_max(j))^2, 1.0);
            fb_nominal_cmd = B0 + Kq*abs(q_cmd_cur) + Kv*abs(dq_cmd_cur) + Ke*phi_cmd_j;

            % ③ 异常扰动
            % 痉挛：门控使用实际命令角度(平滑sigmoid)，强度使用参考速度(避免正反馈)
            gate_st = params.spasm_gate_steepness;
            if j == 1  % 髋
                gate = 1 ./ (1 + exp(-gate_st * max(min(q_cmd_cur - params.spasm_q_hip_th, 5), -5)));
                spasm = (params.spasm_c1/100) * abs(dq_ref(i,j)) * gate;
            elseif j == 2  % 膝（耦合项改用髋参考速度）
                gate = 1 ./ (1 + exp(gate_st * max(min(q_cmd_cur + params.spasm_q_knee_th, 5), -5)));
                spasm = (params.spasm_c2/100) * abs(dq_ref(i,j) + 0.15*dq_ref(i,1)) * gate;
            else  % 踝（耦合项改用膝参考速度）
                spasm = params.spasm_c3 * abs(dq_ref(i,j) + 0.20*dq_ref(i,2));
            end

            % 协同（其他关节使用命令状态 q_cmd = q_ref + delta_q，单步延迟耦合）
            if j == 1
                synergy = params.Csyn(1) * abs(q_ref(i,2) + delta_q(i,2));
            elseif j == 2
                synergy = params.Csyn(2) * abs(q_ref(i,1) + delta_q(i,1));
            else
                synergy = params.Csyn(3) * abs(q_ref(i,2) + delta_q(i,2));
            end

            % 疲劳、事件、噪声（预生成）
            fatigue = pre_fb.fatigue_slope(j) * params.cfatigue * pre_fb.fatigue_trend(i);
            event   = pre_fb.event_all(i,j);
            noise   = pre_fb.noise_all(i,j);

            % 真实反馈（暂不考虑延迟）
            fb_true = fb_nominal_cmd + spasm + synergy + fatigue + event + noise;

            % 原始误差：相对“参考基线”，显式包含阻抗失配项
            e_fb_raw = fb_true - fb_nominal_ref;

            % ④ 阻抗失配前馈补偿（防止 Kv*Δdq 自激正反馈）
            mismatch = fb_nominal_cmd - fb_nominal_ref;

            e_fb_comp = e_fb_raw - K_mis(j) * mismatch;

            % ---- 六分量诊断记录：此处只记录，不改变闭环计算 ----
            pathological = spasm + synergy + fatigue + event + noise;
            interaction_diag.mismatch(i,j) = mismatch;
            interaction_diag.spasm(i,j) = spasm;
            interaction_diag.synergy(i,j) = synergy;
            interaction_diag.fatigue(i,j) = fatigue;
            interaction_diag.event(i,j) = event;
            interaction_diag.noise(i,j) = noise;
            interaction_diag.pathological(i,j) = pathological;
            interaction_diag.compensated(i,j) = e_fb_comp;
            interaction_diag.raw(i,j) = e_fb_raw;
            interaction_diag.fb_nominal_ref(i,j) = fb_nominal_ref;
            interaction_diag.fb_nominal_cmd(i,j) = fb_nominal_cmd;
            interaction_diag.fb_true(i,j) = fb_true;

            % ⑤ 参考加速度（优先使用解析 ddq_ref）
            if nargin >= 6 && ~isempty(ddq_ref)
                ddq_ref_i = ddq_ref(i,j);
            else
                if i == 1
                    ddq_ref_i = (-3*dq_ref(i,j) + 4*dq_ref(i+1,j) - dq_ref(i+2,j)) / (2*dt);
                elseif i == N
                    ddq_ref_i = (3*dq_ref(i,j) - 4*dq_ref(i-1,j) + dq_ref(i-2,j)) / (2*dt);
                else
                    ddq_ref_i = (dq_ref(i+1,j) - dq_ref(i-1,j)) / (2*dt);
                end
            end
            ddq_cmd_est = ddq_ref_i + delta_ddq(i_prev, j);

            % ⑥ 死区 + 静止区处理（对补偿后的误差进行）
            [e_fb_use, quiet_flag] = apply_deadzone_and_quietzone(e_fb_comp, q_cmd_cur, dq_cmd_cur, ...
                ddq_cmd_est, j, params);
            quiet_zone_flag(i,j) = quiet_flag;
            e_fb_used_all(i,j)   = e_fb_use;
            e_fb_raw_all(i,j)    = e_fb_raw;
            interaction_diag.used(i,j) = e_fb_use;

            % ---- Zhang et al. (2018) variable-parameter B/K law ----
            % Source structure (Eqs. 17-18):
            %   B = B0/(a1*exp(|a2*theta|)) * |a3*T_modeling/T_measured|
            %   K = K0/(b1*exp(|b2*theta|)) * |b3*T_modeling/T_measured|
            % followed by lower/upper saturation.
            %
            % In this simulation, fb_nominal_cmd is the model-predicted
            % passive interaction at the controller-specific command state,
            % while fb_true is the corresponding total simulated interaction.
            % This preserves the source law's posture/interaction-ratio logic
            % without adding a new estimator or inner-loop controller.
            T_modeling = abs(fb_nominal_cmd);
            T_measured = max(abs(fb_true), torque_floor);
            torque_ratio = T_modeling / T_measured;

            Ba_raw = Ba0(j) / (a1(j) * exp(abs(a2(j) * q_cmd_cur))) * ...
                     abs(a3(j) * torque_ratio);
            Ka_raw = Ka0(j) / (b1(j) * exp(abs(b2(j) * q_cmd_cur))) * ...
                     abs(b3(j) * torque_ratio);

            % G_a is deliberately NOT adapted in this literature baseline.
            Ga_raw = Ga0(j);

            % ---- 第一次限幅：与 Fixed/Proposed 完全共享 ----
            Ba_tmp = min(max(Ba_raw, 0.1 * Ba0(j)), 4.0 * Ba0(j));
            Ka_tmp = min(max(Ka_raw, 0.01 * Ka0(j)), 5.0 * Ka0(j));
            Ga_tmp = min(max(Ga_raw, 0.08 * Ga0(j)), 2.5 * Ga0(j));

            % ---- 静止区缩放：与 Fixed/Proposed 完全共享 ----
            if quiet_flag
                Ga_tmp = params.quiet_gain_scale(j) * Ga_tmp;
                Ba_tmp = params.quiet_damp_scale(j) * Ba_tmp;
            end

            % ---- 第二次限幅：与 Fixed/Proposed 完全共享 ----
            Ba_tmp = min(max(Ba_tmp, 0.1 * Ba0(j)), 4.0 * Ba0(j));
            Ka_tmp = min(max(Ka_tmp, 0.01 * Ka0(j)), 5.0 * Ka0(j));
            Ga_tmp = min(max(Ga_tmp, 0.08 * Ga0(j)), 2.5 * Ga0(j));

            % 保存参数
            Ba_lit(i,j) = Ba_tmp;
            Ka_lit(i,j) = Ka_tmp;
            Ga_lit(i,j) = Ga_tmp;

            % ---- 导纳离散（后向欧拉） ----
            delta_ddq_tmp = (Ga_tmp * e_fb_use ...
                           - Ba_tmp * delta_dq(i,j) ...
                           - Ka_tmp * delta_q(i,j)) / Ma0(j);
            % 加速度限幅
            delta_ddq_tmp = max(min(delta_ddq_tmp, delta_ddq_limit(j)), -delta_ddq_limit(j));

            % ★ 急动度限幅：限制每步加速度变化率，消除速度/加速度尖峰
            if i > 1
                ddq_prev = delta_ddq(i-1, j);
            else
                ddq_prev = 0;
            end
            max_step = jerk_limit(j) * dt;
            delta_ddq_tmp = max(min(delta_ddq_tmp, ddq_prev + max_step), ddq_prev - max_step);
            delta_ddq(i,j) = delta_ddq_tmp;

            % 积分
            dq_pred = delta_dq(i,j) + delta_ddq_tmp * dt;
            dq_pred = max(min(dq_pred, delta_dq_limit(j)), -delta_dq_limit(j));

            q_pred = delta_q(i,j) + dq_pred * dt;
            q_pred = max(min(q_pred, delta_q_limit(j)), -delta_q_limit(j));

            delta_dq(i+1,j) = dq_pred;
            delta_q(i+1,j)  = q_pred;
        end
    end

    %% ================= 末尾补齐 =================
    delta_q(end,:)   = delta_q(end-1,:);
    delta_dq(end,:)  = delta_dq(end-1,:);
    delta_ddq(end,:) = delta_ddq(end-1,:);
    Ba_lit(end,:) = Ba_lit(end-1,:);
    Ka_lit(end,:) = Ka_lit(end-1,:);
    Ga_lit(end,:) = Ga_lit(end-1,:);
    quiet_zone_flag(end,:) = quiet_zone_flag(end-1,:);
    e_fb_used_all(end,:)   = e_fb_used_all(end-1,:);
    e_fb_raw_all(end,:)    = e_fb_raw_all(end-1,:);

    % 六分量诊断末尾补齐，并执行代数一致性检查
    interaction_diag = finalize_interaction_diag_zhang2018_controller_core(interaction_diag, params);

    exec_time = toc;
    fprintf('✅ Zhang-2018-derived VA（posture/interaction-ratio B/K + shared shell）运行完成：QuietZone 平均激活率 %.2f%% | 执行时间 %.3f s\n', ...
            mean(quiet_zone_flag,'all')*100, exec_time);
end

function diag = init_interaction_diag_zhang2018_controller_core(N)
%INIT_INTERACTION_DIAG 预分配六分量诊断结构体。
    fields = {'mismatch','spasm','synergy','fatigue','event','noise', ...
              'pathological','compensated','used','raw', ...
              'fb_nominal_ref','fb_nominal_cmd','fb_true'};
    diag = struct();
    for k = 1:numel(fields)
        diag.(fields{k}) = zeros(N,3);
    end
end

function diag = finalize_interaction_diag_zhang2018_controller_core(diag, params)
%FINALIZE_INTERACTION_DIAG 补齐最后采样点，并检查六分量重构关系。
    fields = fieldnames(diag);
    for k = 1:numel(fields)
        x = diag.(fields{k});
        if size(x,1) >= 2
            x(end,:) = x(end-1,:);
        end
        diag.(fields{k}) = x;
    end

    raw_reconstructed = diag.mismatch + diag.spasm + diag.synergy + ...
                        diag.fatigue + diag.event + diag.noise;
    err_raw = max(abs(diag.raw(:) - raw_reconstructed(:)));
    assert(err_raw < 1e-10, ...
        'zhang2018_variable_admittance:InteractionDecomposition', ...
        '六分量无法重构 raw interaction residual，最大误差 %.3e', err_raw);

    Kmis = params.K_mismatch_comp;
    if isscalar(Kmis), Kmis = repmat(Kmis,1,3); end
    if max(abs(Kmis(:)-1)) < 1e-12
        err_comp = max(abs(diag.compensated(:) - diag.pathological(:)));
        assert(err_comp < 1e-10, ...
            'zhang2018_variable_admittance:MismatchCompensation', ...
            'K_mismatch_comp=1 时 compensated 应等于 pathological，最大误差 %.3e', err_comp);
    end
end

function [delta_q, delta_dq, delta_ddq, Ba_lit, Ka_lit, Ga_lit, quiet_zone_flag, e_fb_used_all, e_fb_raw_all, interaction_diag] = ...
    cao2024_controller_core(q_ref, dq_ref, params, pre_fb, t_all, ddq_ref)
%CAO2024_CONTROLLER_CORE Cao-2024-derived frozen momentum-based variable admittance comparator.
%
% Source:
%   J. Cao et al., "Variable Admittance Control of High Compatibility
%   Exoskeleton Based on Human-Robotic Interaction Force," Chinese Journal
%   of Mechanical Engineering, 37:119, 2024, doi:10.1186/s10033-024-01113-6.
%
% Source-supported structure retained here:
%   M_d*Delta_ddq + B_d*Delta_dq + K_d*Delta_q = tau_HRI
%   M_d fixed; B_d/K_d adapted by momentum updates across repetitions.
%
% Declared project-specific reduction (NOT exact source reproduction):
%   - Source HRI force is already represented here in joint-torque domain.
%   - HRI deviation/cost signal -> e_fb_comp.
%   - Source repetition index n -> rehabilitation cycle index. The same
%     within-cycle phase is matched across cycles.
%   - The source ILC environmental-gradient layer is replaced by the directly
%     available compensated joint-torque deviation. Effective eta_B/eta_K
%     absorb source beta*alpha/R scaling, because the source does not report
%     all numerical learning parameters needed for exact reproduction.
%   - The original active controller concerns hip/knee; the same joint-wise
%     law is explicitly extended to the ankle for the present 3-DOF study.
%   - Inner PID/robot rigid-body dynamics are not added, consistent with the
%     present project's command-level admittance scope.
%
% IMPORTANT: formal Monte Carlo requires params.literature_cao2024.frozen=true.
% The effective port gains must not be retuned from nominal/MC performance.

    N = size(q_ref, 1);
    dt = params.dt;
    tic;

    %% ================= 参数提取 =================
    Ma0 = params.Ma0; Ba0 = params.Ba0; Ka0 = params.Ka0; Ga0 = params.Ga0;

    assert(isfield(params,'literature_cao2024'), ...
        'cao2024_variable_admittance:MissingConfig', ...
        'params.literature_cao2024 未定义；请在 parameters.m 中配置文献型 VA 参数。');
    lit = params.literature_cao2024;
    zeta_B = expand3_cao(lit.zeta_B, 'zeta_B');
    zeta_K = expand3_cao(lit.zeta_K, 'zeta_K');
    eta_B  = expand3_cao(lit.eta_B,  'eta_B');
    eta_K  = expand3_cao(lit.eta_K,  'eta_K');
    assert(all(zeta_B >= 0 & zeta_B < 1), 'literature_cao2024.zeta_B 必须满足 0<=zeta<1');
    assert(all(zeta_K >= 0 & zeta_K < 1), 'literature_cao2024.zeta_K 必须满足 0<=zeta<1');
    assert(all(eta_B > 0) && all(eta_K > 0), 'literature_cao2024.eta_B/eta_K 必须为正');

    verbose = true;
    if isfield(lit,'verbose'), verbose = logical(lit.verbose); end

    delta_q_limit   = params.delta_q_limit;
    delta_dq_limit  = params.delta_dq_limit;
    delta_ddq_limit = params.delta_ddq_limit;
    if isscalar(delta_q_limit),   delta_q_limit   = repmat(delta_q_limit,  1, 3); end
    if isscalar(delta_dq_limit),  delta_dq_limit  = repmat(delta_dq_limit, 1, 3); end
    if isscalar(delta_ddq_limit), delta_ddq_limit = repmat(delta_ddq_limit,1, 3); end

    jerk_limit = params.jerk_limit;
    if isscalar(jerk_limit), jerk_limit = repmat(jerk_limit, 1, 3); end
    K_mis = params.K_mismatch_comp;
    if isscalar(K_mis), K_mis = repmat(K_mis, 1, 3); end

    q1_max = params.q1_max; q2_max = params.q2_max;
    if isfield(params,'q3_max'), q3_max = params.q3_max; else q3_max = 0.5; end
    q_max = [q1_max, q2_max, q3_max];

    %% ================= Cao cycle/phase mapping =================
    % Current trajectory template has (#key states-1) equal-duration segments.
    n_segments = size(params.traj1_deg,1) - 1;
    samples_per_cycle = round((n_segments * params.T_seg) / dt);
    assert(samples_per_cycle >= 1, 'Cao port: samples_per_cycle 无效');
    inferred_cycles = (N-1) / samples_per_cycle;
    assert(abs(inferred_cycles - round(inferred_cycles)) < 1e-10, ...
        'Cao port 要求 N-1 能被每周期采样间隔数整除；当前 N=%d, samples/cycle=%d', N, samples_per_cycle);
    inferred_cycles = round(inferred_cycles);
    assert(inferred_cycles == params.loop_num, ...
        'Cao port 周期映射不一致：轨迹推断 %d 周期，params.loop_num=%d', inferred_cycles, params.loop_num);

    % Each phase stores the base (pre-quiet-zone) B/K used when that phase
    % appears next time. Cycle 1 starts from the common nominal B0/K0.
    B_phase = repmat(Ba0(:)', samples_per_cycle, 1);
    K_phase = repmat(Ka0(:)', samples_per_cycle, 1);
    VB_phase = zeros(samples_per_cycle, 3);
    VK_phase = zeros(samples_per_cycle, 3);

    %% ================= 预分配 =================
    delta_q   = zeros(N, 3);
    delta_dq  = zeros(N, 3);
    delta_ddq = zeros(N, 3);
    Ba_lit = zeros(N, 3); Ka_lit = zeros(N, 3); Ga_lit = zeros(N, 3);
    quiet_zone_flag = zeros(N, 3);
    e_fb_used_all   = zeros(N, 3);
    e_fb_raw_all    = zeros(N, 3);
    interaction_diag = init_interaction_diag_cao2024_controller_core(N);

    %% ================= 主循环 =================
    for i = 1:N-1
        i_prev = max(i-1, 1);
        cycle_idx = floor((i-1)/samples_per_cycle) + 1;
        phase_idx = mod(i-1, samples_per_cycle) + 1;

        for j = 1:3
            q_cmd_cur  = q_ref(i,j)  + delta_q(i,j);
            dq_cmd_cur = dq_ref(i,j) + delta_dq(i,j);

            % ---- shared closed-loop patient interaction model ----
            B0 = params.B0(j);
            Kq = params.Kq_fb(j);
            Kv = params.Kv_fb(j);
            Ke = params.Ke_fb(j);

            phi_ref_j = min((abs(q_ref(i,j)) / q_max(j))^2, 1.0);
            fb_nominal_ref = B0 + Kq*abs(q_ref(i,j)) + Kv*abs(dq_ref(i,j)) + Ke*phi_ref_j;

            phi_cmd_j = min((abs(q_cmd_cur) / q_max(j))^2, 1.0);
            fb_nominal_cmd = B0 + Kq*abs(q_cmd_cur) + Kv*abs(dq_cmd_cur) + Ke*phi_cmd_j;

            gate_st = params.spasm_gate_steepness;
            if j == 1
                gate = 1 ./ (1 + exp(-gate_st * max(min(q_cmd_cur - params.spasm_q_hip_th, 5), -5)));
                spasm = (params.spasm_c1/100) * abs(dq_ref(i,j)) * gate;
            elseif j == 2
                gate = 1 ./ (1 + exp(gate_st * max(min(q_cmd_cur + params.spasm_q_knee_th, 5), -5)));
                spasm = (params.spasm_c2/100) * abs(dq_ref(i,j) + 0.15*dq_ref(i,1)) * gate;
            else
                spasm = params.spasm_c3 * abs(dq_ref(i,j) + 0.20*dq_ref(i,2));
            end

            if j == 1
                synergy = params.Csyn(1) * abs(q_ref(i,2) + delta_q(i,2));
            elseif j == 2
                synergy = params.Csyn(2) * abs(q_ref(i,1) + delta_q(i,1));
            else
                synergy = params.Csyn(3) * abs(q_ref(i,2) + delta_q(i,2));
            end

            fatigue = pre_fb.fatigue_slope(j) * params.cfatigue * pre_fb.fatigue_trend(i);
            event   = pre_fb.event_all(i,j);
            noise   = pre_fb.noise_all(i,j);

            fb_true = fb_nominal_cmd + spasm + synergy + fatigue + event + noise;
            e_fb_raw = fb_true - fb_nominal_ref;
            mismatch = fb_nominal_cmd - fb_nominal_ref;
            e_fb_comp = e_fb_raw - K_mis(j) * mismatch;

            pathological = spasm + synergy + fatigue + event + noise;
            interaction_diag.mismatch(i,j) = mismatch;
            interaction_diag.spasm(i,j) = spasm;
            interaction_diag.synergy(i,j) = synergy;
            interaction_diag.fatigue(i,j) = fatigue;
            interaction_diag.event(i,j) = event;
            interaction_diag.noise(i,j) = noise;
            interaction_diag.pathological(i,j) = pathological;
            interaction_diag.compensated(i,j) = e_fb_comp;
            interaction_diag.raw(i,j) = e_fb_raw;
            interaction_diag.fb_nominal_ref(i,j) = fb_nominal_ref;
            interaction_diag.fb_nominal_cmd(i,j) = fb_nominal_cmd;
            interaction_diag.fb_true(i,j) = fb_true;

            if nargin >= 6 && ~isempty(ddq_ref)
                ddq_ref_i = ddq_ref(i,j);
            else
                if i == 1
                    ddq_ref_i = (-3*dq_ref(i,j) + 4*dq_ref(i+1,j) - dq_ref(i+2,j)) / (2*dt);
                elseif i == N
                    ddq_ref_i = (3*dq_ref(i,j) - 4*dq_ref(i-1,j) + dq_ref(i-2,j)) / (2*dt);
                else
                    ddq_ref_i = (dq_ref(i+1,j) - dq_ref(i-1,j)) / (2*dt);
                end
            end
            ddq_cmd_est = ddq_ref_i + delta_ddq(i_prev, j);

            [e_fb_use, quiet_flag] = apply_deadzone_and_quietzone(e_fb_comp, q_cmd_cur, dq_cmd_cur, ...
                ddq_cmd_est, j, params);
            quiet_zone_flag(i,j) = quiet_flag;
            e_fb_used_all(i,j)   = e_fb_use;
            e_fb_raw_all(i,j)    = e_fb_raw;
            interaction_diag.used(i,j) = e_fb_use;

            % ---- Cao-2024-derived cycle-domain momentum B/K law ----
            % Current phase's pre-supervision B/K are inherited from the
            % previous repetition. Cycle 1 begins at the common nominal scale.
            Ba_base = B_phase(phase_idx,j);
            Ka_base = K_phase(phase_idx,j);
            Ga_raw = Ga0(j);  % source keeps inertia fixed; present comparator keeps G common/nominal

            % Common projection + quiet-zone supervision for the APPLIED law.
            Ba_tmp = min(max(Ba_base, 0.1 * Ba0(j)), 4.0 * Ba0(j));
            Ka_tmp = min(max(Ka_base, 0.01 * Ka0(j)), 5.0 * Ka0(j));
            Ga_tmp = min(max(Ga_raw, 0.08 * Ga0(j)), 2.5 * Ga0(j));
            if quiet_flag
                Ga_tmp = params.quiet_gain_scale(j) * Ga_tmp;
                Ba_tmp = params.quiet_damp_scale(j) * Ba_tmp;
            end
            Ba_tmp = min(max(Ba_tmp, 0.1 * Ba0(j)), 4.0 * Ba0(j));
            Ka_tmp = min(max(Ka_tmp, 0.01 * Ka0(j)), 5.0 * Ka0(j));
            Ga_tmp = min(max(Ga_tmp, 0.08 * Ga0(j)), 2.5 * Ga0(j));

            Ba_lit(i,j) = Ba_tmp;
            Ka_lit(i,j) = Ka_tmp;
            Ga_lit(i,j) = Ga_tmp;

            % ---- Common admittance discrete dynamics ----
            delta_ddq_tmp = (Ga_tmp * e_fb_use ...
                           - Ba_tmp * delta_dq(i,j) ...
                           - Ka_tmp * delta_q(i,j)) / Ma0(j);
            delta_ddq_tmp = max(min(delta_ddq_tmp, delta_ddq_limit(j)), -delta_ddq_limit(j));
            if i > 1
                ddq_prev = delta_ddq(i-1, j);
            else
                ddq_prev = 0;
            end
            max_step = jerk_limit(j) * dt;
            delta_ddq_tmp = max(min(delta_ddq_tmp, ddq_prev + max_step), ddq_prev - max_step);
            delta_ddq(i,j) = delta_ddq_tmp;

            dq_pred = delta_dq(i,j) + delta_ddq_tmp * dt;
            dq_pred = max(min(dq_pred, delta_dq_limit(j)), -delta_dq_limit(j));
            q_pred = delta_q(i,j) + dq_pred * dt;
            q_pred = max(min(q_pred, delta_q_limit(j)), -delta_q_limit(j));
            delta_dq(i+1,j) = dq_pred;
            delta_q(i+1,j)  = q_pred;

            % ---- Causal update for the SAME phase in the NEXT repetition ----
            % The source obtains dW/dF through an ILC environmental model.
            % In this reduced-order port the directly available compensated
            % joint-torque deviation is used as the gradient surrogate.
            g_hri = e_fb_comp;
            VB_new = zeta_B(j) * VB_phase(phase_idx,j) + ...
                     eta_B(j) * delta_dq(i,j) * g_hri;
            VK_new = zeta_K(j) * VK_phase(phase_idx,j) + ...
                     eta_K(j) * delta_q(i,j) * g_hri;

            B_next = Ba_base - VB_new;
            K_next = Ka_base - VK_new;
            % Project the learned base values to the SAME common admissible set.
            B_next = min(max(B_next, 0.1 * Ba0(j)), 4.0 * Ba0(j));
            K_next = min(max(K_next, 0.01 * Ka0(j)), 5.0 * Ka0(j));

            interaction_diag.cao_cycle_index(i,j) = cycle_idx;
            interaction_diag.cao_phase_index(i,j) = phase_idx;
            interaction_diag.cao_gradient_surrogate(i,j) = g_hri;
            interaction_diag.cao_base_B(i,j) = Ba_base;
            interaction_diag.cao_base_K(i,j) = Ka_base;
            interaction_diag.cao_VB(i,j) = VB_new;
            interaction_diag.cao_VK(i,j) = VK_new;
            interaction_diag.cao_next_B(i,j) = B_next;
            interaction_diag.cao_next_K(i,j) = K_next;

            VB_phase(phase_idx,j) = VB_new;
            VK_phase(phase_idx,j) = VK_new;
            B_phase(phase_idx,j) = B_next;
            K_phase(phase_idx,j) = K_next;
        end
    end

    %% ================= 末尾补齐 =================
    delta_q(end,:)   = delta_q(end-1,:);
    delta_dq(end,:)  = delta_dq(end-1,:);
    delta_ddq(end,:) = delta_ddq(end-1,:);
    Ba_lit(end,:) = Ba_lit(end-1,:);
    Ka_lit(end,:) = Ka_lit(end-1,:);
    Ga_lit(end,:) = Ga_lit(end-1,:);
    quiet_zone_flag(end,:) = quiet_zone_flag(end-1,:);
    e_fb_used_all(end,:)   = e_fb_used_all(end-1,:);
    e_fb_raw_all(end,:)    = e_fb_raw_all(end-1,:);

    interaction_diag = finalize_interaction_diag_cao2024_controller_core(interaction_diag, params);

    exec_time = toc;
    if verbose
        fprintf(['✅ Cao-2024-derived VA（cycle-domain momentum B/K + shared shell）运行完成：', ...
                 'cycles=%d | QuietZone %.2f%% | %.3f s\n'], ...
                inferred_cycles, mean(quiet_zone_flag,'all')*100, exec_time);
    end
end

function v = expand3_cao(v, name)
    if isscalar(v), v = repmat(v,1,3); end
    assert(numel(v)==3 && all(isfinite(v)), 'literature_cao2024.%s 必须为有限标量或3元素向量', name);
    v = v(:)';
end

function diag = init_interaction_diag_cao2024_controller_core(N)
%INIT_INTERACTION_DIAG Preallocate shared + Cao-port diagnostic fields.
    fields = {'mismatch','spasm','synergy','fatigue','event','noise', ...
              'pathological','compensated','used','raw', ...
              'fb_nominal_ref','fb_nominal_cmd','fb_true', ...
              'cao_cycle_index','cao_phase_index','cao_gradient_surrogate', ...
              'cao_base_B','cao_base_K','cao_VB','cao_VK','cao_next_B','cao_next_K'};
    diag = struct();
    for k = 1:numel(fields)
        diag.(fields{k}) = zeros(N,3);
    end
end

function diag = finalize_interaction_diag_cao2024_controller_core(diag, params)
%FINALIZE_INTERACTION_DIAG Fill final sample and check interaction algebra.
    fields = fieldnames(diag);
    for k = 1:numel(fields)
        x = diag.(fields{k});
        if size(x,1) >= 2
            x(end,:) = x(end-1,:);
        end
        diag.(fields{k}) = x;
    end

    raw_reconstructed = diag.mismatch + diag.spasm + diag.synergy + ...
                        diag.fatigue + diag.event + diag.noise;
    err_raw = max(abs(diag.raw(:) - raw_reconstructed(:)));
    assert(err_raw < 1e-10, ...
        'cao2024_variable_admittance:InteractionDecomposition', ...
        '六分量无法重构 raw interaction residual，最大误差 %.3e', err_raw);

    Kmis = params.K_mismatch_comp;
    if isscalar(Kmis), Kmis = repmat(Kmis,1,3); end
    if max(abs(Kmis(:)-1)) < 1e-12
        err_comp = max(abs(diag.compensated(:) - diag.pathological(:)));
        assert(err_comp < 1e-10, ...
            'cao2024_variable_admittance:MismatchCompensation', ...
            'K_mismatch_comp=1 时 compensated 应等于 pathological，最大误差 %.3e', err_comp);
    end
end

function [delta_q, delta_dq, delta_ddq, Ba_adp, Ka_adp, Ga_adp, quiet_zone_flag, e_fb_used_all, e_fb_raw_all, interaction_diag] = ...
    proposed_controller_core(q_ref, dq_ref, params, pre_fb, t_all, ddq_ref)

    N = size(q_ref, 1);
    dt = params.dt;
    tic;

    %% ================= 参数提取 =================
    Ma0 = params.Ma0; Ba0 = params.Ba0; Ka0 = params.Ka0; Ga0 = params.Ga0;
    alpha_B = params.alpha_B; beta_B = params.beta_B;
    alpha_C = params.alpha_C;
    alpha_K = params.alpha_K;
    alpha_G = params.alpha_G;

    % ================= Core-ablation switches =================
    % Default behavior is exactly the original full Proposed controller.
    % The switches are read only when params.ablation is supplied.
    adapt_B = true;
    adapt_K = true;
    adapt_G = true;
    ablation_name = 'Full Proposed';
    if isfield(params, 'ablation') && ~isempty(params.ablation)
        if isfield(params.ablation, 'adapt_B'), adapt_B = logical(params.ablation.adapt_B); end
        if isfield(params.ablation, 'adapt_K'), adapt_K = logical(params.ablation.adapt_K); end
        if isfield(params.ablation, 'adapt_G'), adapt_G = logical(params.ablation.adapt_G); end
        if isfield(params.ablation, 'name') && ~isempty(params.ablation.name)
            ablation_name = char(params.ablation.name);
        end
    end

    delta_q_limit   = params.delta_q_limit;
    delta_dq_limit  = params.delta_dq_limit;
    delta_ddq_limit = params.delta_ddq_limit;
    if isscalar(delta_q_limit),   delta_q_limit   = repmat(delta_q_limit,  1, 3); end
    if isscalar(delta_dq_limit),  delta_dq_limit  = repmat(delta_dq_limit, 1, 3); end
    if isscalar(delta_ddq_limit), delta_ddq_limit = repmat(delta_ddq_limit,1, 3); end

    jerk_limit = params.jerk_limit;
    if isscalar(jerk_limit), jerk_limit = repmat(jerk_limit, 1, 3); end
    K_mis = params.K_mismatch_comp;
    if isscalar(K_mis), K_mis = repmat(K_mis, 1, 3); end

    % 关节极限（用于边界接近度 φ(q)）
    q1_max = params.q1_max; q2_max = params.q2_max;
    if isfield(params,'q3_max'), q3_max = params.q3_max; else q3_max = 0.5; end
    q_max = [q1_max, q2_max, q3_max];

    %% ================= 预分配 =================
    delta_q   = zeros(N, 3);
    delta_dq  = zeros(N, 3);
    delta_ddq = zeros(N, 3);
    Ba_adp = zeros(N, 3); Ka_adp = zeros(N, 3); Ga_adp = zeros(N, 3);
    quiet_zone_flag = zeros(N, 3);
    e_fb_used_all   = zeros(N, 3);
    e_fb_raw_all    = zeros(N, 3);

    % ================= 六分量交互力矩诊断（仅记录，不参与控制） =================
    % raw = mismatch + spasm + synergy + fatigue + event + noise
    interaction_diag = init_interaction_diag_proposed_controller_core(N);

    % 为差分计算保留上一步的 e_fb_use
    e_fb_prev = zeros(1, 3);

    %% ================= 主循环 =================
    for i = 1:N-1
        i_prev = max(i-1, 1);   % 上一步索引（用于上一步修正量）
        for j = 1:3
            % ---- 当前命令轨迹（标量） ----
            q_cmd_cur  = q_ref(i,j)   + delta_q(i,j);
            dq_cmd_cur = dq_ref(i,j)  + delta_dq(i,j);

            % ---- 实时反馈计算（基于 q_cmd 闭环） ----
            B0 = params.B0(j);
            Kq = params.Kq_fb(j);
            Kv = params.Kv_fb(j);
            Ke = params.Ke_fb(j);

            % ① 参考基线（健康参考状态下的正常阻抗，不随控制器行为漂移）
            phi_ref_j = min((abs(q_ref(i,j)) / q_max(j))^2, 1.0);
            fb_nominal_ref = B0 + Kq*abs(q_ref(i,j)) + Kv*abs(dq_ref(i,j)) + Ke*phi_ref_j;

            % ② 命令基线（机器人实际命令状态下的正常阻抗）
            phi_cmd_j = min((abs(q_cmd_cur) / q_max(j))^2, 1.0);
            fb_nominal_cmd = B0 + Kq*abs(q_cmd_cur) + Kv*abs(dq_cmd_cur) + Ke*phi_cmd_j;

            % ③ 异常扰动
            % 痉挛：门控使用实际命令角度(平滑sigmoid)，强度使用参考速度(避免正反馈)
            gate_st = params.spasm_gate_steepness;
            if j == 1  % 髋
                gate = 1 ./ (1 + exp(-gate_st * max(min(q_cmd_cur - params.spasm_q_hip_th, 5), -5)));
                spasm = (params.spasm_c1/100) * abs(dq_ref(i,j)) * gate;
            elseif j == 2  % 膝（耦合项改用髋参考速度）
                gate = 1 ./ (1 + exp(gate_st * max(min(q_cmd_cur + params.spasm_q_knee_th, 5), -5)));
                spasm = (params.spasm_c2/100) * abs(dq_ref(i,j) + 0.15*dq_ref(i,1)) * gate;
            else  % 踝（耦合项改用膝参考速度）
                spasm = params.spasm_c3 * abs(dq_ref(i,j) + 0.20*dq_ref(i,2));
            end

            % 协同（其他关节使用命令状态 q_cmd = q_ref + delta_q，单步延迟耦合）
            if j == 1
                synergy = params.Csyn(1) * abs(q_ref(i,2) + delta_q(i,2));
            elseif j == 2
                synergy = params.Csyn(2) * abs(q_ref(i,1) + delta_q(i,1));
            else
                synergy = params.Csyn(3) * abs(q_ref(i,2) + delta_q(i,2));
            end

            % 疲劳、事件、噪声（预生成）
            fatigue = pre_fb.fatigue_slope(j) * params.cfatigue * pre_fb.fatigue_trend(i);
            event   = pre_fb.event_all(i,j);
            noise   = pre_fb.noise_all(i,j);

            % 真实反馈（暂不考虑延迟）
            fb_true = fb_nominal_cmd + spasm + synergy + fatigue + event + noise;

            % 原始误差：相对“参考基线”，显式包含阻抗失配项
            e_fb_raw = fb_true - fb_nominal_ref;

            % ④ 阻抗失配前馈补偿（防止 Kv*Δdq 自激正反馈）
            mismatch = fb_nominal_cmd - fb_nominal_ref;

            e_fb_comp = e_fb_raw - K_mis(j) * mismatch;

            % ---- 六分量诊断记录：此处只记录，不改变闭环计算 ----
            pathological = spasm + synergy + fatigue + event + noise;
            interaction_diag.mismatch(i,j) = mismatch;
            interaction_diag.spasm(i,j) = spasm;
            interaction_diag.synergy(i,j) = synergy;
            interaction_diag.fatigue(i,j) = fatigue;
            interaction_diag.event(i,j) = event;
            interaction_diag.noise(i,j) = noise;
            interaction_diag.pathological(i,j) = pathological;
            interaction_diag.compensated(i,j) = e_fb_comp;
            interaction_diag.raw(i,j) = e_fb_raw;
            interaction_diag.fb_nominal_ref(i,j) = fb_nominal_ref;
            interaction_diag.fb_nominal_cmd(i,j) = fb_nominal_cmd;
            interaction_diag.fb_true(i,j) = fb_true;

            % ⑤ 参考加速度（优先使用解析 ddq_ref）
            if nargin >= 6 && ~isempty(ddq_ref)
                ddq_ref_i = ddq_ref(i,j);
            else
                if i == 1
                    ddq_ref_i = (-3*dq_ref(i,j) + 4*dq_ref(i+1,j) - dq_ref(i+2,j)) / (2*dt);
                elseif i == N
                    ddq_ref_i = (3*dq_ref(i,j) - 4*dq_ref(i-1,j) + dq_ref(i-2,j)) / (2*dt);
                else
                    ddq_ref_i = (dq_ref(i+1,j) - dq_ref(i-1,j)) / (2*dt);
                end
            end
            ddq_cmd_est = ddq_ref_i + delta_ddq(i_prev, j);

            % ⑥ 死区 + 静止区处理（对补偿后的误差进行）
            [e_fb_use, quiet_flag] = apply_deadzone_and_quietzone(e_fb_comp, q_cmd_cur, dq_cmd_cur, ...
                ddq_cmd_est, j, params);
            quiet_zone_flag(i,j) = quiet_flag;
            e_fb_used_all(i,j)   = e_fb_use;
            e_fb_raw_all(i,j)    = e_fb_raw;
            interaction_diag.used(i,j) = e_fb_use;

            % ---- 用于参数更新的误差幅值（使用处理后的 e_fb_use） ----
            e_abs  = abs(e_fb_use);
            if i == 1
                de_abs = 0;
                e_fb_prev(j) = e_fb_use;
            else
                de_abs = abs(e_fb_use - e_fb_prev(j)) / dt;
                e_fb_prev(j) = e_fb_use;
            end

            % ---- 边界接近度 φ (使用命令状态 q_cmd) ----
            phi = phi_cmd_j;

            % ---- 计算自适应参数（未饱和） ----
            % Ablation rule: a disabled channel is held at its nominal value
            % before the common projection/quiet-zone supervisory shell.
            if adapt_B
                Ba_raw = Ba0(j) - alpha_B(j) * e_abs + beta_B(j) * de_abs;
            else
                Ba_raw = Ba0(j);
            end

            if adapt_K
                Ka_raw = Ka0(j) - alpha_C(j) * e_abs + alpha_K(j) * phi;
            else
                Ka_raw = Ka0(j);
            end

            if adapt_G
                Ga_raw = Ga0(j) - alpha_G(j) * e_abs;
            else
                Ga_raw = Ga0(j);
            end

            % ---- 第一次限幅（宽松下限） ----
            Ba_tmp = min(max(Ba_raw, 0.1 * Ba0(j)), 4.0 * Ba0(j));
            Ka_tmp = min(max(Ka_raw, 0.01 * Ka0(j)), 5.0 * Ka0(j));
            Ga_tmp = min(max(Ga_raw, 0.08 * Ga0(j)), 2.5 * Ga0(j));

            % ---- 静止区缩放 ----
            if quiet_flag
                Ga_tmp = params.quiet_gain_scale(j) * Ga_tmp;
                Ba_tmp = params.quiet_damp_scale(j) * Ba_tmp;
            end

            % ---- 第二次限幅（防止缩放后越界） ----
            Ba_tmp = min(max(Ba_tmp, 0.1 * Ba0(j)), 4.0 * Ba0(j));
            Ka_tmp = min(max(Ka_tmp, 0.01 * Ka0(j)), 5.0 * Ka0(j));
            Ga_tmp = min(max(Ga_tmp, 0.08 * Ga0(j)), 2.5 * Ga0(j));

            % 保存参数
            Ba_adp(i,j) = Ba_tmp;
            Ka_adp(i,j) = Ka_tmp;
            Ga_adp(i,j) = Ga_tmp;

            % ---- 导纳离散（后向欧拉） ----
            delta_ddq_tmp = (Ga_tmp * e_fb_use ...
                           - Ba_tmp * delta_dq(i,j) ...
                           - Ka_tmp * delta_q(i,j)) / Ma0(j);
            % 加速度限幅
            delta_ddq_tmp = max(min(delta_ddq_tmp, delta_ddq_limit(j)), -delta_ddq_limit(j));

            % ★ 急动度限幅：限制每步加速度变化率，消除速度/加速度尖峰
            if i > 1
                ddq_prev = delta_ddq(i-1, j);
            else
                ddq_prev = 0;
            end
            max_step = jerk_limit(j) * dt;
            delta_ddq_tmp = max(min(delta_ddq_tmp, ddq_prev + max_step), ddq_prev - max_step);
            delta_ddq(i,j) = delta_ddq_tmp;

            % 积分
            dq_pred = delta_dq(i,j) + delta_ddq_tmp * dt;
            dq_pred = max(min(dq_pred, delta_dq_limit(j)), -delta_dq_limit(j));

            q_pred = delta_q(i,j) + dq_pred * dt;
            q_pred = max(min(q_pred, delta_q_limit(j)), -delta_q_limit(j));

            delta_dq(i+1,j) = dq_pred;
            delta_q(i+1,j)  = q_pred;
        end
    end

    %% ================= 末尾补齐 =================
    delta_q(end,:)   = delta_q(end-1,:);
    delta_dq(end,:)  = delta_dq(end-1,:);
    delta_ddq(end,:) = delta_ddq(end-1,:);
    Ba_adp(end,:) = Ba_adp(end-1,:);
    Ka_adp(end,:) = Ka_adp(end-1,:);
    Ga_adp(end,:) = Ga_adp(end-1,:);
    quiet_zone_flag(end,:) = quiet_zone_flag(end-1,:);
    e_fb_used_all(end,:)   = e_fb_used_all(end-1,:);
    e_fb_raw_all(end,:)    = e_fb_raw_all(end-1,:);

    % 六分量诊断末尾补齐，并执行代数一致性检查
    interaction_diag = finalize_interaction_diag_proposed_controller_core(interaction_diag, params);

    exec_time = toc;
    fprintf('✅ Proposed variant [%s] 运行完成：K_mis=[%.2f %.2f %.2f] | Adapt(B,K,G)=[%d %d %d] | QuietZone %.2f%% | %.3f s\n', ...
            ablation_name, K_mis(1), K_mis(2), K_mis(3), adapt_B, adapt_K, adapt_G, ...
            mean(quiet_zone_flag,'all')*100, exec_time);
end

function diag = init_interaction_diag_proposed_controller_core(N)
%INIT_INTERACTION_DIAG 预分配六分量诊断结构体。
    fields = {'mismatch','spasm','synergy','fatigue','event','noise', ...
              'pathological','compensated','used','raw', ...
              'fb_nominal_ref','fb_nominal_cmd','fb_true'};
    diag = struct();
    for k = 1:numel(fields)
        diag.(fields{k}) = zeros(N,3);
    end
end

function diag = finalize_interaction_diag_proposed_controller_core(diag, params)
%FINALIZE_INTERACTION_DIAG 补齐最后采样点，并检查六分量重构关系。
    fields = fieldnames(diag);
    for k = 1:numel(fields)
        x = diag.(fields{k});
        if size(x,1) >= 2
            x(end,:) = x(end-1,:);
        end
        diag.(fields{k}) = x;
    end

    raw_reconstructed = diag.mismatch + diag.spasm + diag.synergy + ...
                        diag.fatigue + diag.event + diag.noise;
    err_raw = max(abs(diag.raw(:) - raw_reconstructed(:)));
    assert(err_raw < 1e-10, ...
        'adaptive_admittance:InteractionDecomposition', ...
        '六分量无法重构 raw interaction residual，最大误差 %.3e', err_raw);

    Kmis = params.K_mismatch_comp;
    if isscalar(Kmis), Kmis = repmat(Kmis,1,3); end
    if max(abs(Kmis(:)-1)) < 1e-12
        err_comp = max(abs(diag.compensated(:) - diag.pathological(:)));
        assert(err_comp < 1e-10, ...
            'adaptive_admittance:MismatchCompensation', ...
            'K_mismatch_comp=1 时 compensated 应等于 pathological，最大误差 %.3e', err_comp);
    end
end

function [e_fb_used, quiet_zone_flag] = apply_deadzone_and_quietzone(e_fb_raw, q_ref_i, dq_ref_i, ddq_ref_i, j, params)
    % 统一死区与静止区处理函数（v2：换向邻域静止区判定）
    % 输入：
    %   e_fb_raw   : 当前原始反馈误差 (标量)
    %   q_ref_i    : 当前参考位置 (标量) – 保留以备未来使用
    %   dq_ref_i   : 当前参考速度 (标量)
    %   ddq_ref_i  : 当前参考加速度 (标量, 新增)
    %   j          : 关节索引 (1-3)
    %   params     : 参数结构体（需包含 quiet_ddq_th）
    % 输出：
    %   e_fb_used       : 处理后的误差（用于导纳计算）
    %   quiet_zone_flag : 静止区标志 (0或1)
    
    %% ================= 参数提取 =================
    quiet_dq_th     = params.quiet_dq_th(j);
    quiet_gain_scale = params.quiet_gain_scale(j);
    quiet_damp_scale = params.quiet_damp_scale(j);
    e_deadzone      = params.e_deadzone(j);
    
    % 换向邻域加速度阈值（若未定义则默认取速度阈值的 10 倍量级）
    if isfield(params, 'quiet_ddq_th')
        quiet_ddq_th = params.quiet_ddq_th(j);
    else
        quiet_ddq_th = 10.0 * pi/180;   % 默认 10 deg/s² ≈ 0.1745 rad/s²
    end
    
    %% ================= 静止区判断（换向邻域） =================
    % 条件：速度接近零 且 加速度接近零 → 换向点附近
    in_quiet_zone = (abs(dq_ref_i) < quiet_dq_th) && (abs(ddq_ref_i) < quiet_ddq_th);
    quiet_zone_flag = double(in_quiet_zone);
    
    %% ================= 误差死区处理 =================
    if abs(e_fb_raw) < e_deadzone
        e_fb_used = 0;
    else
        e_fb_used = sign(e_fb_raw) * (abs(e_fb_raw) - e_deadzone);
    end
    
    %% ================= 可选：日志记录（调试用） =================
    % 若需在论文Supplementary中展示静止区行为，可在此添加持久化记录
    % persistent log_data;
end
