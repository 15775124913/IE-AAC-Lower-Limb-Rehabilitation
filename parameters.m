function params = parameters()


    params = struct();
    
    % ========================= 仿真基本设置 =========================
    params.dt        = 0.01;                % 采样时间 (s)
    params.T_seg     = 2.0;                 % 每段轨迹时间 (s)
    params.loop_num  = 14;                  % 循环次数（14-cycle evaluation horizon）
    params.seed      = 42;                  % 固定随机种子
    params.config_date = datetime('now','Format','yyyyMMdd_HHmmss');
    
    % ===================== 连杆几何尺寸 (表5.1) =====================
    params.model.L1              = 0.45;     % 大腿长度 (m)
    params.model.L2              = 0.40;     % 小腿长度 (m)
    params.model.L_ankle_to_heel = 0.07;     % 踝至脚跟 (m)
    params.model.L_heel_to_toe   = 0.15;     % 脚跟至脚尖 (m)
    
    % =================== 生理极点与安全限幅 (表5.1) =================
    params.q1_max          = 2.09;           % 髋关节控制归一化边界 (rad)
    params.q2_max          = 2.356;           % 膝关节控制归一化边界 (rad)
    params.q3_max          = 0.785;           % 踝关节控制归一化边界 (rad)
    params.eps_singular    = 0.01;           % 末端防奇异小正数
    
    % ========================= 安全运动限幅 =========================
    params.delta_q_limit   = 0.08;           % 最大位置修正 (rad)
    params.delta_dq_limit  = 1.5;            % 最大速度修正 (rad/s)
    params.delta_ddq_limit = 5.0;            % 最大加速度修正 (rad/s²)
    params.jerk_limit      = [20.0 20.0 15.0]; % 导纳修正加速度变化率上限 (rad/s³)
    params.motion_limit_saturation_ratio = 0.999; % 认为触及共享运动限幅的归一化阈值；用于诊断，不参与控制
    
    % ================ 正常生理阻抗 (表5.2) ========================
    params.B0     = [0.5  0.5  0.5];         % 关节固有基线摩擦力矩 (Nm)
    params.Kq_fb  = [2.0  2.0  2.0];         % 角度刚度增益
    params.Kv_fb  = [1.2  1.2  1.2];         % 速度黏滞增益
    params.Ke_fb  = [50.0 50.0 50.0];        % 末端保护刚度增益 (Nm/rad)
    
    % ================ 异常病理扰动参数 (表5.2) ======================
    params.spasm_q_hip_th  = 0.8;            % 髋痉挛阈值 (rad)
    params.spasm_q_knee_th = 0.5;            % 膝痉挛阈值 (rad)
    params.spasm_gate_steepness = 20.0;      % 痉挛门控 Sigmoid 陡度 (1/rad)
    % 痉挛速度灵敏度高于正常黏滞项，用于形成异常速度相关交互
    % 实际系数 = c/100: 髋3.5, 膝4.5, 踝1.0
    params.spasm_c1        = 350.0;          % 髋痉挛等效阻尼增益1
    params.spasm_c2        = 450.0;          % 膝痉挛等效阻尼增益2
    params.spasm_c3        = 1.0;            % 踝痉挛等效增益
    params.event_amp     = [0.50, 0.60, 0.70];  % 事件扰动峰值 (Nm)

    % ★ Repeated-cycle event configuration
    % Event is injected after adaptation has evolved for several cycles.
    % Absolute event time is generated in event_model.m from cycle index + phase.
    params.event_cycle = 8;                       % event occurs in the 8th cycle
    params.event_phase_base = [2.4, 3.8, 5.2];    % within-cycle event phase (s)
    params.event_time_base = params.event_phase_base; % compatibility placeholder
    params.event_sigma       = [0.25, 0.28, 0.30];
    params.event_time_offset = [0, 0, 0];
    params.event_window_nsigma = 3.0;

    % 仅作旧接口兼容；正式指标/绘图均通过 event_model() 动态计算
    event_centers_nominal = params.event_time_base + params.event_time_offset;
    params.t_abnormal = [ ...
        min(event_centers_nominal - params.event_window_nsigma*params.event_sigma), ...
        max(event_centers_nominal + params.event_window_nsigma*params.event_sigma)];
    params.event_win = params.t_abnormal;

    params.fatigue_slope = [0.015, 0.018, 0.022];  % 疲劳斜率系数
    
    params.Csyn            = [0.2 0.4 0.1];  % 跨关节运动耦合系数
    params.cfatigue        = 10.0;           % 最大疲劳阻力增益 (Nm)
    params.feedback_noise_level = 0.02;      % 反馈高斯噪声标准差/幅值参数
    params.sigma_N = params.feedback_noise_level; % 论文符号兼容别名，不单独参与计算
    
    % ================ 自适应导纳控制参数 (表5.3) ======================
    % 基础虚拟导纳参数
    params.Ma0 = [2.0 2.0 1.5];              % 虚拟惯量对角 (kg·m²)
    params.Ba0 = [10  10   6];               % 基础虚拟阻尼 (Nms/rad)
    params.Ka0 = [50  40  30];               % 基础虚拟刚度 (Nm/rad)
    params.Ga0 = [2.0 2.0 2.5];              % 初始力矩误差映射增益
    
    % 交互误差软死区
    params.e_deadzone = [0.05 0.05 0.05];    % 交互误差死区阈值 Dth (Nm)
    params.K_mismatch_comp = [1 1 1];  % 阻抗失配前馈补偿系数
                                             % 0=不补偿(原行为), 1=完全补偿(只保留病理力矩进入)
    
    % IE-AAC 自适应系数
    params.alpha_B = [3.0  3.0  3.0];        % 阻尼比例系数
    params.beta_B  = [0.5  0.5  0.5];        % 阻尼微分系数
    params.alpha_K = [20.0 20.0 20.0];       % 刚度比例系数
    params.alpha_G = [1.0  0.5  1.0];        % 增益自适应率；保持当前已执行参数值
    params.alpha_C = [8.0  8.0  8.0];        % 中心顺应系数

    % ================ 文献型 Variable Admittance 对照 =================
    % 本项目保留两种相互独立的文献派生变导纳比较器：
    %   1) Zhang-2018：姿态 + model/measured interaction ratio 驱动 B/K 调度；
    %   2) Cao-2024-derived：cycle-domain momentum 机制驱动 B/K 更新。
    % 两种比较器均保持 M_a 和 G_a 为项目公共 nominal 值，并共享患者模型、
    % mismatch compensation、soft dead-zone、quiet-zone、parameter projection
    % 以及 Delta-q / Delta-dq / Delta-ddq / jerk 约束。

    % ---------- Zhang et al. 2018 literature-derived comparator ----------
    % M. Zhang et al., IEEE Trans. Ind. Electron., 65(2):1398-1407, 2018,
    % doi:10.1109/TIE.2017.2733425.
    %
    % Source structure:
    %   B = B0/(a1*exp(|a2*theta|)) * |a3*T_modeling/T_measured|
    %   K = K0/(b1*exp(|b2*theta|)) * |b3*T_modeling/T_measured|
    % with fixed admittance inertia. The project port uses q_cmd as theta,
    % command-state nominal passive interaction as T_modeling, and total
    % simulated interaction torque as T_measured.
    params.literature_zhang2018 = struct();
    params.literature_zhang2018.source = 'Zhang2018_TIE_variable_parameter_admittance';
    params.literature_zhang2018.doi = '10.1109/TIE.2017.2733425';
    params.literature_zhang2018.a1 = [1.0 1.0 1.0];
    params.literature_zhang2018.a2 = [0.15 0.15 0.15];
    params.literature_zhang2018.a3 = [1.0 1.0 1.0];
    params.literature_zhang2018.b1 = [1.0 1.0 1.0];
    params.literature_zhang2018.b2 = [0.15 0.15 0.15];
    params.literature_zhang2018.b3 = [1.0 1.10 1.0];
    params.literature_zhang2018.torque_floor = 0.05;
    params.literature_zhang2018.frozen = true;
    params.literature_zhang2018.freeze_basis = ...
        'Nominal audit 20260923_184921: zero B/K bound hits; Knee DeltaQ saturation 7.498125%; all DeltaDQ/DeltaDDQ saturation 0%.';
    params.literature_zhang2018.port_note = ...
        ['Zhang2018 Eqs17-18 B/K scheduling; theta=q_cmd, Tmodel=command-state nominal passive torque, ', ...
         'Tmeasured=total interaction torque; M/G nominal; shared supervisory shell; frozen before formal Monte Carlo.'];

    % ---------- Cao et al. 2024-derived comparator ----------
    % J. Cao et al., Chinese Journal of Mechanical Engineering, 37:119, 2024,
    % doi:10.1186/s10033-024-01113-6.
    %
    % Source-supported structure:
    %   M_d*Delta_ddq + B_d*Delta_dq + K_d*Delta_q = tau_HRI
    %   M_d fixed; B_d/K_d updated by momentum terms over repeated operations.
    %
    % Present reduced-order port:
    %   q-q0, dq-dq0        -> delta_q, delta_dq
    %   HRI deviation       -> e_fb_comp (joint-torque residual, Nm)
    %   repetition index n  -> rehabilitation-cycle index with matched phase
    %   beta*alpha/R terms  -> effective eta_B / eta_K
    % The source ILC environmental-gradient layer is not recreated; the
    % directly available compensated joint-torque deviation is used as the
    % gradient surrogate. The same joint-wise law is explicitly extended to
    % the ankle because the source active controller is centered on hip/knee.
    params.literature_cao2024 = struct();
    params.literature_cao2024.source = 'Cao2024_CJME_momentum_variable_admittance_derived';
    params.literature_cao2024.doi = '10.1186/s10033-024-01113-6';
    params.literature_cao2024.original_scope = 'Hip/knee active joints; HRI force converted to joint torque; fixed M; momentum+ILC B/K adaptation.';
    params.literature_cao2024.port_type = 'reduced_order_cycle_domain_momentum_port';
    params.literature_cao2024.input_mapping = 'HRI deviation -> e_fb_comp (joint torque residual, Nm)';
    params.literature_cao2024.iteration_mapping = 'source repetition n -> rehabilitation cycle n; same within-cycle phase is matched';
    params.literature_cao2024.ankle_extension = true;
    params.literature_cao2024.zeta_B = 0.50;
    params.literature_cao2024.zeta_K = 0.50;
    params.literature_cao2024.eta_B  = 2.5;
    params.literature_cao2024.eta_K  = 10.0;
    params.literature_cao2024.verbose = true;
    params.literature_cao2024.frozen = true;
    params.literature_cao2024.freeze_basis = ...
        ['Comparator-only structural calibration CaoCalibration_20261002_150348; no Proposed performance metric was used for selection. ', ...
         'Frozen effective gains: zeta_B=zeta_K=0.5, eta_B=2.5, eta_K=10. ', ...
         'Selected candidate: B/K bound-hit 0%; max DeltaQ saturation 5.949%; DeltaDQ/DeltaDDQ saturation 0%; ', ...
         'mean cycle-1-to-last B/K RMS change 5.309%/2.890% of nominal; ', ...
         'joint B changes [3.123, 11.442, 1.361]% and K changes [2.098, 5.998, 0.575]% for hip/knee/ankle.'];
    params.literature_cao2024.port_note = ...
        ['Cao2024-derived cycle-domain momentum B/K adaptation. e_fb_comp is the joint-torque HRI-deviation surrogate; ', ...
         'four rehabilitation repetitions are treated as finite iterations; M/G remain nominal; ankle uses a declared ', ...
         'joint-wise extension; the project does not reproduce the source environmental ILC state-space layer. ', ...
         'Effective eta_B/eta_K were frozen by comparator-only structural calibration before the renewed formal Monte Carlo run.'];

    % ================ 静止区 & 死区 (控制辅助参数, 非表格项) =========
    params.quiet_dq_th = deg2rad([0.1 0.1 0.1]);   % 速度静止阈值 (rad/s)
    params.quiet_ddq_th = deg2rad([0.5, 0.5, 0.5]);   % 换向邻域加速度阈值 (rad/s²)
    params.quiet_gain_scale = [1.0 1.0 1.0];     % 静止区增益缩放
    params.quiet_damp_scale = [1.35 1.35 1.60];     % 静止区阻尼缩放
    
    
    % ==================== 示教重建统一配置 ====================
    % nominal 与 Monte Carlo 独立运行时均由 TeachingDataPipeline 读取此处，
    % 避免额外 prepare_teaching_reference.m 包装脚本。
    params.teaching = struct();
    params.teaching.noise_level = 0.03;
    params.teaching.noise_vel = 0.5;
    params.teaching.noise_acc = 5.0;
    params.teaching.warp_amp = 0.05;
    params.teaching.warp_freq = 2*pi/0.8;
    params.teaching.sample_dt = 0.10;         % 示教采样周期 0.10 s = 10 Hz（控制器仍为 100 Hz）
    params.teaching.eps_tol = 1.0;
    params.teaching.max_points = 10;
    params.teaching.min_points = 6;           % 当前10 s周期 + 2 s最大间隔约束下的实际最小值
    params.teaching.max_key_gap = params.T_seg;
    params.teaching.min_key_spacing = 0.70;  % 抑制0.4~0.5 s局部关键点扎堆，降低分段五次多项式高曲率

    % ==================== Monte Carlo 不确定性配置 ====================
    % 第5项修复：所有随机扰动均围绕 nominal 参数，采用有界均匀分布。
    % 只扰动当前闭环控制/统计指标中真正生效的变量。
    params.mc = struct();
    params.mc.distribution = 'bounded_uniform';

    params.mc.rel_admittance = 0.10;          % Ma0/Ba0/Ka0/Ga0: ±10%
    params.mc.rel_normal_feedback = 0.10;     % Kq_fb/Kv_fb/Ke_fb: ±10%
    params.mc.rel_spasm_threshold = 0.10;     % 痉挛门控阈值: ±10%
    params.mc.rel_pathology = 0.15;           % 病理强度参数: ±15%
    params.mc.rel_noise = 0.20;               % feedback_noise_level: ±20%
    params.mc.rel_quiet_threshold = 0.10;     % quiet_dq_th/quiet_ddq_th: ±10%
    params.mc.event_time_offset_abs = 0.10;   % 三关节事件中心偏移: ±0.10 s

    % 参数抽样与时序噪声使用独立局部随机流，避免随机序列相互耦合
    params.mc.parameter_seed_offset = 10000;
    params.mc.noise_seed_offset = 20000;

    % 论文审计阈值：低于此成功率时给出强警告
    params.mc.min_success_rate = 0.95;

    % 明确不随机化：
    % model.L1/L2 - 当前 MC 性能指标均为关节空间/交互指标，连杆长度不进入控制闭环；
    % B0          - 常数基线项在当前反馈误差构造中抵消。

    % ==================== 方案A：基础髋膝联动被动屈伸（末段动力学协调版） ====================
    % 原末段 Hip [40,-10] -> [0,0]、Knee [-25,8] -> [0,0] 在 2 s 内
    % 需要较强的中途再加速/再减速。这里仅调整倒数第二关键状态，使返回中立位时
    % 的位移需求与速度方向/幅值更一致；首末静止边界保持 [0,0,0]。
    params.traj1_deg = [0   0   0; 20  8   0; 45  12  0; 65  0   0; 35 -18  0; 0   0   0];
    params.traj2_deg = [0   0   0; -15  -6  0; -35 -10  0; -55  -8  0; -22  11  0; 0   0   0];
    params.traj3_deg = [0   0   0; 5   2   0; 10   0   0; 8   -2   0; -3  -2   0; 0   0   0];
    % ==================== 方案B：多周期温和步态模拟 ====================
    % params.traj1_deg = [0   0   0; 25  10  0; 55  15  0; 70   8  0; 35 -12  0; 0   0   0];
    % params.traj2_deg = [0   0   0; -20  -8  0; -45 -12  0; -58  -6  0; -22  10  0; 0   0   0];
    % params.traj3_deg = [0   0   0; 8   3   0; 14   4   0; 6   2   0; -10  -4   0; 0   0   0];
    % ==================== 方案C：强调踝关节+小范围膝髋 ====================
    % params.traj1_deg = [0    0   0; 15    6   0; 28    7   0; 35    8   0; 25   -6   0; 0    0   0];
    % params.traj2_deg = [0    0   0; -10   -4   0; -25   -7   0; -35   -9   0; -20    6   0; 0    0   0];
    % params.traj3_deg = [0    0   0; 12    4   0; -8   -5   0; 10    4   0; -10   -4   0; 0    0   0];
    
    % ================ 版本与验证 ====================================
    time_stamp = datestr(now,'yyyymmdd_HHMMSS');
    params.version = ['Final_v12_four_controller_two_literature_plus_ablation_', time_stamp];
    
    validate_parameters(params); % 参数合理性检查（本文件局部函数）
    
    if ~exist('data','dir'), mkdir('data'); end
    save(fullfile('data',['config_',params.version,'.mat']),'params');
    export_paper_tables(params);

end

% ==================================================================
% 函数：export_paper_tables
% 功能：基于当前 params 结构体生成论文表5.1,5.2,5.3的LaTeX三线表格及CSV
% ==================================================================
function export_paper_tables(params)
    
    % 使用与 save_results 完全相同的子文件夹路径
    if isfield(params, 'config_date')
        timestamp = char(params.config_date);
        out_dir = fullfile('results', timestamp);
    else
        % 兼容无 config_date 的情况（兜底）
        out_dir = 'results';
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    
    % ------------- 表5.1：连杆几何尺寸与安全限幅 -------------
    % ★ 同步更新：限幅与代码一致
    tab1 = cell(10,4);
    tab1(1,:)  = {'连杆长度','$L_1$ / $L_2$','大腿/小腿长度','0.45 m / 0.40 m'};
    tab1(2,:)  = {'','$L_{\mathrm{heel}}$ / $L_{\mathrm{toe}}$','踝至脚跟 / 脚跟至脚尖','0.07 m / 0.15 m'};
    tab1(3,:)  = {'控制边界','$q_{1,\max}$ / $q_{2,\max}$ / $q_{3,\max}$','髋/膝/踝归一化边界','2.09 / 2.356 / 0.785 rad'};
    tab1(4,:)  = {'','$\varepsilon$','防奇异小正数','0.01'};
    tab1(5,:)  = {'安全限幅','$\Delta q_{\max}$','单周期最大位置修正','0.08 rad'};
    tab1(6,:)  = {'','$\Delta\dot{q}_{\max}$','最大速度修正','1.5 rad/s'};
    tab1(7,:)  = {'','$\Delta\ddot{q}_{\max}$','最大加速度修正','5.0 rad/s$^2$'};
    tab1(8,:)  = {'','','',''};
    tab1(9,:)  = {'','','',''};
    tab1(10,:) = {'','','',''};
    
    % ------------- 表5.2：正常生理阻抗与异常病理扰动 -------------
    % ★ 同步更新：c1/c2/c3、cfatigue 与代码一致
    tab2 = cell(13,4);
    tab2(1,:)  = {'正常基线','$B_{0,j}$','关节固有基线摩擦力矩','0.5 Nm'};
    tab2(2,:)  = {'','$K_{q,j}$ / $K_{v,j}$','角度刚度增益 / 速度黏滞增益','2.0 / 1.2'};
    tab2(3,:)  = {'','$K_{e,j}$','末端保护刚度增益','50.0 Nm/rad'};
    tab2(4,:)  = {'痉挛特征','$q_{\mathrm{hip,th}}$ / $q_{\mathrm{knee,th}}$','痉挛角度阈值','0.8 rad / 0.5 rad'};
    tab2(5,:)  = {'','$c_1$ / $c_2$ / $c_3$','痉挛等效阻尼增益','350.0 / 450.0 / 1.0'};
    tab2(6,:)  = {'异常协同','$C_{\mathrm{syn}}$','跨关节耦合系数','[0.2, 0.4, 0.1]'};
    tab2(7,:)  = {'疲劳','$c_{\mathrm{fatigue}}$','最大疲劳阻力','10.0 Nm'};
    tab2(8,:)  = {'噪声','$\sigma_N$','反馈高斯噪声水平','0.02'};
    tab2(9,:)  = {'','','',''};
    tab2(10,:) = {'','','',''};
    tab2(11,:) = {'','','',''};
    tab2(12,:) = {'','','',''};
    tab2(13,:) = {'','','',''};
    
    % ------------- 表5.3：自适应导纳控制动态演化参数 -------------
    % ★ 同步更新：Ma0/Ba0/Ka0/Ga0/Dth/自适应率 与代码一致
    tab3 = cell(15,4);
    tab3(1,:)  = {'基础导纳','$M_{a,j}$','虚拟惯量对角','[2.0, 2.0, 1.5] kg$\cdot$m$^2$'};
    tab3(2,:)  = {'','$B_{a0,j}$ / $K_{a0,j}$','基础虚拟阻尼/刚度','[10,10,6] Nms/rad / [50,40,30] Nm/rad'};
    tab3(3,:)  = {'','$G_{a0,j}$','初始力矩映射增益','[2.0, 2.0, 2.5]'};
    tab3(4,:)  = {'滤波机制','$D_{\mathrm{th}}$','交互误差死区阈值','0.05 Nm'};
    tab3(5,:)  = {'自适应率','$\alpha_{B,j}$ / $\beta_{B,j}$','阻尼比例/微分系数','3.0 / 0.5'};
    tab3(6,:)  = {'','$\alpha_{K,j}$','刚度比例增益','20.0'};
    tab3(7,:)  = {'','$\alpha_{G,j}$','映射增益自适应率（髋/膝/踝）','[1.0, 0.5, 1.0]'};
    tab3(8,:)  = {'','$\alpha_{C,j}$','中心顺应系数','8.0'};
    tab3(9,:)  = {'','','',''};
    tab3(10,:) = {'','','',''};
    tab3(11,:) = {'','','',''};
    tab3(12,:) = {'','','',''};
    tab3(13,:) = {'','','',''};
    tab3(14,:) = {'','','',''};
    tab3(15,:) = {'','','',''};
    
    % ---------- 生成 .tex 文件（传入 out_dir） ----------
    write_latex_table(tab1, 'Table_5_1_Link_Safety.tex', '连杆几何尺寸与安全限幅参数', out_dir);
    write_latex_table(tab2, 'Table_5_2_Impedance_Pathology.tex', '正常生理阻抗与异常病理扰动模型参数', out_dir);
    write_latex_table(tab3, 'Table_5_3_Admittance_Evolution.tex', '自适应导纳控制动态演化参数', out_dir);
    
    % ---------- 生成 CSV 汇总 ----------
    fid = fopen(fullfile(out_dir, 'all_params.csv'), 'w');
    for i = 1:size(tab1,1)
        if ~all(cellfun(@isempty, tab1(i,1:2)))
            fprintf(fid, '%s,%s,%s,%s\n', tab1{i,1}, tab1{i,2}, tab1{i,3}, tab1{i,4});
        else
            fprintf(fid, '\n');
        end
    end
    fprintf(fid, '\n');
    for i = 1:size(tab2,1)
        if ~all(cellfun(@isempty, tab2(i,1:2)))
            fprintf(fid, '%s,%s,%s,%s\n', tab2{i,1}, tab2{i,2}, tab2{i,3}, tab2{i,4});
        else
            fprintf(fid, '\n');
        end
    end
    fprintf(fid, '\n');
    for i = 1:size(tab3,1)
        if ~all(cellfun(@isempty, tab3(i,1:2)))
            fprintf(fid, '%s,%s,%s,%s\n', tab3{i,1}, tab3{i,2}, tab3{i,3}, tab3{i,4});
        else
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
    
    fprintf('✅ 论文表格已导出至子文件夹：%s\n', out_dir);
end

% ---------- 辅助函数：写出LaTeX三线表（增加 out_dir 参数） ----------
function write_latex_table(data_rows, filename, caption, out_dir)
    fid = fopen(fullfile(out_dir, filename), 'w');
    fprintf(fid, '%% 自动生成 - 请勿手动修改\n');
    fprintf(fid, '\\begin{table}[htbp]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{%s}\n', caption);
    fprintf(fid, '\\begin{tabular}{lllp{4cm}}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, '参数类别 & 物理符号 & 描述说明 & 设定数值 \\\\\n');
    fprintf(fid, '\\midrule\n');
    for i = 1:size(data_rows,1)
        row = data_rows(i,:);
        fprintf(fid, '%s & %s & %s & %s \\\\\n', row{1}, row{2}, row{3}, row{4});
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);
end

% ---------- 辅助函数：写出CSV（保留但不再直接调用） ----------
function csvwrite_with_header(filename, data)
    fid = fopen(filename, 'w');
    for i = 1:size(data,1)
        fprintf(fid, '%s,%s,%s,%s\n', data{i,1}, data{i,2}, data{i,3}, data{i,4});
    end
    fclose(fid);
end

% ==================================================================
% 局部函数：参数合理性验证（原 config_validation.m 已合并到本文件）
% ==================================================================
function validate_parameters(params)
    % 参数合理性验证函数（作为 parameters.m 的子函数）
    % 在 parameters.m 末尾调用，确保仿真参数物理合理、维度一致
    % 发现问题时给出清晰警告或错误，便于调试和Supplementary说明
    
    fprintf('=== 开始参数配置验证 ===\n');
    errors = {};
    warnings = {};
    
    %% ================= 基本维度检查 =================
    n_joints = 3;  % 髋、膝、踝
    
    % 检查关节相关向量长度
    fields_to_check = {'Ba0', 'Ka0', 'Ga0', 'Ma0', 'alpha_B', 'beta_B', ...
                       'alpha_K', 'alpha_G', 'B0', 'Kq_fb', 'Kv_fb', 'Ke_fb', ...
                       'quiet_gain_scale', 'quiet_damp_scale'};
    
    for i = 1:length(fields_to_check)
        fld = fields_to_check{i};
        if isfield(params, fld)
            val = params.(fld);
            if length(val) ~= n_joints
                errors{end+1} = sprintf('参数 %s 长度应为 %d，当前为 %d', ...
                    fld, n_joints, length(val));
            end
        end
    end
    
    %% ================= 物理参数检查 =================
    model = params.model;
    
    if model.L1 <= 0 || model.L2 <= 0
        errors{end+1} = '连杆长度 L1、L2 必须大于0';
    end
    if model.L_ankle_to_heel <= 0 || model.L_heel_to_toe <= 0
        errors{end+1} = '脚部尺寸参数必须大于0';
    end
    
    %% ================= 控制增益合理性检查 =================
    if any(params.Ba0 <= 0) || any(params.Ka0 <= 0)
        warnings{end+1} = 'Ba0 或 Ka0 存在非正值，可能导致系统不稳定';
    end
    if any(params.Ga0 <= 0)
        errors{end+1} = 'Ga0 必须为正值';
    end
    
    % 静止区阈值
    if any(params.quiet_dq_th < 0) || any(params.quiet_ddq_th < 0)
        warnings{end+1} = '静止区速度/加速度阈值不应为负值';
    end

    % 两种文献型 variable-admittance comparator 参数
    if ~isfield(params,'literature_zhang2018') || ~isstruct(params.literature_zhang2018)
        errors{end+1} = '缺少 params.literature_zhang2018 配置';
    else
        lit = params.literature_zhang2018;
        req_lit = {'a1','a2','a3','b1','b2','b3'};
        for ii = 1:numel(req_lit)
            fld = req_lit{ii};
            if ~isfield(lit, fld)
                errors{end+1} = sprintf('params.literature_zhang2018 缺少字段 %s', fld);
            elseif ~(isscalar(lit.(fld)) || numel(lit.(fld)) == n_joints)
                errors{end+1} = sprintf('params.literature_zhang2018.%s 必须为标量或3元素向量', fld);
            elseif any(~isfinite(lit.(fld))) || any(lit.(fld) <= 0)
                errors{end+1} = sprintf('params.literature_zhang2018.%s 必须为正有限值', fld);
            end
        end
        if ~isfield(lit,'torque_floor') || ~isscalar(lit.torque_floor) || ...
                ~isfinite(lit.torque_floor) || lit.torque_floor <= 0
            errors{end+1} = 'params.literature_zhang2018.torque_floor 必须为正有限标量';
        end
        if ~isfield(lit,'frozen') || ~isscalar(lit.frozen)
            errors{end+1} = 'params.literature_zhang2018.frozen 必须存在且为标量逻辑量';
        end
    end

    if ~isfield(params,'literature_cao2024') || ~isstruct(params.literature_cao2024)
        errors{end+1} = '缺少 params.literature_cao2024 配置';
    else
        lit = params.literature_cao2024;
        req_lit = {'zeta_B','zeta_K','eta_B','eta_K'};
        for ii = 1:numel(req_lit)
            fld = req_lit{ii};
            if ~isfield(lit, fld)
                errors{end+1} = sprintf('params.literature_cao2024 缺少字段 %s', fld);
            elseif ~(isscalar(lit.(fld)) || numel(lit.(fld)) == n_joints)
                errors{end+1} = sprintf('params.literature_cao2024.%s 必须为标量或3元素向量', fld);
            elseif any(~isfinite(lit.(fld)))
                errors{end+1} = sprintf('params.literature_cao2024.%s 必须为有限值', fld);
            end
        end
        if isfield(lit,'zeta_B') && (any(lit.zeta_B < 0) || any(lit.zeta_B >= 1))
            errors{end+1} = 'params.literature_cao2024.zeta_B 必须满足 0 <= zeta_B < 1';
        end
        if isfield(lit,'zeta_K') && (any(lit.zeta_K < 0) || any(lit.zeta_K >= 1))
            errors{end+1} = 'params.literature_cao2024.zeta_K 必须满足 0 <= zeta_K < 1';
        end
        if isfield(lit,'eta_B') && any(lit.eta_B <= 0)
            errors{end+1} = 'params.literature_cao2024.eta_B 必须为正值';
        end
        if isfield(lit,'eta_K') && any(lit.eta_K <= 0)
            errors{end+1} = 'params.literature_cao2024.eta_K 必须为正值';
        end
        if ~isfield(lit,'frozen') || ~isscalar(lit.frozen)
            errors{end+1} = 'params.literature_cao2024.frozen 必须存在且为标量逻辑量';
        end
    end

    %% ================= 轨迹参数检查 =================
    traj_fields = {'traj1_deg', 'traj2_deg', 'traj3_deg'};
    for i = 1:length(traj_fields)
        fld = traj_fields{i};
        traj = params.(fld);
        if size(traj, 2) ~= 3
            errors{end+1} = sprintf('%s 每行应为 [q, dq, ddq] 三列，当前列数 %d', ...
                fld, size(traj,2));
        end
        % 检查角度是否在合理生理范围内（松散检查）
        if any(abs(traj(:,1)) > 180)  % 度
            warnings{end+1} = sprintf('%s 中角度超过 ±180°，请确认生理合理性', fld);
        end
        % 当前项目采用循环被动训练模板，首尾 q/dq/ddq 应闭合。
        if max(abs(traj(1,:) - traj(end,:))) > 1e-12
            warnings{end+1} = sprintf('%s 首尾状态不一致；多周期训练可能产生连接不连续', fld);
        end
    end
    
    %% ================= 时间与仿真参数 =================
    if params.dt <= 0 || params.T_seg <= 0 || params.loop_num < 1
        errors{end+1} = 'dt、T_seg、loop_num 必须为正值';
    end

    %% ================= 事件扰动时间参数检查 =================
    event_fields = {'event_time_base','event_sigma','event_time_offset'};
    for i = 1:length(event_fields)
        fld = event_fields{i};
        if ~isfield(params, fld) || numel(params.(fld)) ~= 3
            errors{end+1} = sprintf('%s 必须存在且长度为3（Hip/Knee/Ankle）', fld);
        end
    end
    if isfield(params,'event_sigma') && any(params.event_sigma <= 0)
        errors{end+1} = 'event_sigma 必须全部大于0';
    end
    if ~isfield(params,'event_window_nsigma') || ...
            ~isscalar(params.event_window_nsigma) || params.event_window_nsigma <= 0
        errors{end+1} = 'event_window_nsigma 必须为正标量';
    end
    
    %% ================= Monte Carlo 配置检查 =================
    if ~isfield(params, 'mc') || ~isstruct(params.mc)
        errors{end+1} = '缺少 params.mc Monte Carlo 配置';
    else
        rel_fields = {'rel_admittance','rel_normal_feedback', ...
                      'rel_spasm_threshold','rel_pathology', ...
                      'rel_noise','rel_quiet_threshold'};
        for ii = 1:numel(rel_fields)
            fld_mc = rel_fields{ii};
            if ~isfield(params.mc, fld_mc) || ~isscalar(params.mc.(fld_mc)) || ...
                    params.mc.(fld_mc) < 0 || params.mc.(fld_mc) >= 1
                errors{end+1} = sprintf('params.mc.%s 必须在 [0,1) 内', fld_mc);
            end
        end

        if ~isfield(params.mc,'event_time_offset_abs') || ...
                ~isscalar(params.mc.event_time_offset_abs) || ...
                params.mc.event_time_offset_abs < 0
            errors{end+1} = 'params.mc.event_time_offset_abs 必须为非负标量';
        end

        if ~isfield(params.mc,'min_success_rate') || ...
                ~isscalar(params.mc.min_success_rate) || ...
                params.mc.min_success_rate <= 0 || params.mc.min_success_rate > 1
            errors{end+1} = 'params.mc.min_success_rate 必须在 (0,1] 内';
        end
    end


    %% ================= 示教重建配置检查 =================
    if ~isfield(params,'teaching') || ~isstruct(params.teaching)
        errors{end+1} = '缺少 params.teaching 示教重建配置';
    else
        if params.teaching.noise_level < 0 || ...
           params.teaching.noise_vel < 0 || ...
           params.teaching.noise_acc < 0
            errors{end+1} = 'params.teaching 噪声幅值必须为非负';
        end
        if params.teaching.max_points < params.teaching.min_points
            errors{end+1} = 'teaching.max_points 必须 >= teaching.min_points';
        end
        if params.teaching.max_key_gap <= 0 || ...
           params.teaching.min_key_spacing < 0
            errors{end+1} = '示教关键点时间约束必须合理';
        end
    end

    %% ================= 输出验证结果 =================
    if ~isempty(errors)
        fprintf('❌ 发现 %d 个错误：\n', length(errors));
        for e = errors
            disp(['   - ' e{1}]);
        end
        error('参数验证失败，请修正以上错误后重新运行。');
    end
    
    if ~isempty(warnings)
        fprintf('⚠️  发现 %d 个警告：\n', length(warnings));
        for w = warnings
            disp(['   - ' w{1}]);
        end
    else
        fprintf('✅ 参数配置验证通过，所有参数合理。\n');
    end
end
