function [params_updated, metrics, q_controller_ref, dq_controller_ref, ddq_controller_ref, t_controller_ref] = TeachingDataPipeline(params, varargin)
% TeachingDataPipeline  稀疏模拟编码器示教数据 → 自适应关键点抽取 → 五次多项式重建 → 误差统计
% 改进：物理特征点（速度/加速度/Jerk过零点）+覆盖约束+基于示教观测的贪婪误差驱动
%
% 输入:
%   params   - 原始参数结构体
%   varargin - 'NoiseLevel','WarpAmp','TeachDt','PlotFlag','SaveFlag','Verbose','EpsTol','MaxPoints','MinPoints','Seed'
% 输出:
%   params_updated - 更新后的参数结构体（关键点数自适应）
%   metrics        - 结构体，包含拟合误差、重建轨迹和示教数据等
%   q_controller_ref / dq_controller_ref / ddq_controller_ref / t_controller_ref
%                  - 最终实际送入控制器的参考轨迹（由本函数唯一生成并直接返回）

% ---------- 0. 设置保存路径（与 save_results.m 一致） ----------
timestampStr = '';
saveDir = '';
if isfield(params, 'config_date')
    try
        timestampStr = datestr(params.config_date, 'yyyymmdd_HHMMSS');
        saveDir = fullfile('results', timestampStr);
    catch
        timestampStr = '';
    end
end
if isempty(timestampStr)
    timestampStr = datestr(now, 'yyyymmdd_HHMMSS');
    saveDir = fullfile('results', timestampStr);
end
% 保存目录是否真正创建，由 SaveFlag 在参数解析后决定

% ---------- 1. 解析可选参数 ----------
% 数值默认值统一从 parameters.m -> params.teaching 读取。
teach_cfg = teaching_defaults(params);
p = inputParser;
addParameter(p, 'NoiseLevel', teach_cfg.noise_level);
addParameter(p, 'NoiseVel', teach_cfg.noise_vel);
addParameter(p, 'NoiseAcc', teach_cfg.noise_acc);
addParameter(p, 'WarpAmp', teach_cfg.warp_amp);
addParameter(p, 'WarpFreq', teach_cfg.warp_freq);
if isfield(teach_cfg, 'sample_dt') && ~isempty(teach_cfg.sample_dt)
    teach_dt_default = teach_cfg.sample_dt;
else
    teach_dt_default = 0.10;  % 默认 10 Hz 示教采样，控制器仍使用 params.dt
end
addParameter(p, 'TeachDt', teach_dt_default); % 示教采样周期(s)，与控制周期解耦
addParameter(p, 'PlotFlag', true);          % 是否显示示教管线图
addParameter(p, 'SaveFlag', true);             % ★ 是否保存示教中间 MAT/CSV/PDF
addParameter(p, 'Verbose', true);              % ★ 是否打印关键点详细信息
addParameter(p, 'EpsTol', teach_cfg.eps_tol);            % 重建RMSE容差（°）
addParameter(p, 'MaxPoints', teach_cfg.max_points);           % 最大关键点数
addParameter(p, 'MinPoints', teach_cfg.min_points);            % 最少关键点数（含起止点）
addParameter(p, 'Seed', []);                 % ★ 示教管线独立随机种子；为空时继承 params.seed
addParameter(p, 'MaxKeyGap', teach_cfg.max_key_gap);            % ★ 最大关键点时间间隔；为空时取原始最大轨迹段时长
addParameter(p, 'MinKeySpacing', teach_cfg.min_key_spacing);      % ★ 最小关键点间隔(s)，防止关键点局部扎堆
parse(p, varargin{:});
opts = p.Results;
MinPoints = opts.MinPoints;

if opts.SaveFlag && ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% ---------- 1.1 示教随机数流：固定种子 + 不污染全局 RNG ----------
if isempty(opts.Seed)
    if ~isfield(params, 'seed') || isempty(params.seed)
        error('TeachingDataPipeline:MissingSeed', ...
              '未提供 Seed，且 params.seed 不存在，无法保证示教数据可重复性。');
    end
    teaching_seed = params.seed;
else
    teaching_seed = opts.Seed;
end

validateattributes(teaching_seed, {'numeric'}, ...
    {'scalar','real','finite','nonnegative'}, ...
    mfilename, 'Seed');

% 使用局部 RandStream：
% 1) 相同 Seed 得到严格一致的示教噪声
% 2) 不调用 rng(...)，避免重置主仿真 / Monte Carlo 的全局随机状态
teaching_stream = RandStream('mt19937ar', 'Seed', double(teaching_seed));

% ---------- 2. 生成完美离散轨迹（真值） ----------
% ★ 第2项修复：示教评估与控制器统一使用 params.dt 的离散时间网格
% 原代码使用固定 4000 点 linspace，会导致示教评估网格与控制器网格不一致。
if isfield(params, 'traj_times')
    seg_dur = params.traj_times(:);
else
    nseg = size(params.traj1_deg,1)-1;
    if isfield(params, 'T_seg')
        seg_dur = params.T_seg * ones(nseg,1);
    else
        seg_dur = 2 * ones(nseg,1);
    end
end

dt_ctrl = params.dt;
T_total = sum(seg_dur);
n_steps = round(T_total / dt_ctrl);

if abs(n_steps * dt_ctrl - T_total) > 1e-9
    error('TeachingDataPipeline:TimeGridMismatch', ...
        ['原始示教周期总时长 %.12g s 不能被控制周期 dt=%.12g s 整除。' ...
         '请先使轨迹时长与控制采样周期一致。'], T_total, dt_ctrl);
end

t_all = (0:n_steps)' * dt_ctrl;
NUM_CTRL_SAMPLES = length(t_all);

% ---------- 2.1 示教采样网格：低于控制频率，减少示教数据点 ----------
% 控制器保持 dt_ctrl=0.01 s（100 Hz）；示教采样默认 TeachDt=0.10 s（10 Hz）。
% 为保证关键点时刻可无误差落到控制器网格，TeachDt 被吸附为 dt_ctrl 的整数倍。
validateattributes(opts.TeachDt, {'numeric'}, {'scalar','real','finite','positive'}, ...
                   mfilename, 'TeachDt');
teach_stride = max(1, round(opts.TeachDt / dt_ctrl));
dt_teach = teach_stride * dt_ctrl;
if abs(dt_teach - opts.TeachDt) > 1e-12 && opts.Verbose
    fprintf('提示：TeachDt=%.6f s 已调整为控制周期的整数倍 %.6f s。\n', ...
            opts.TeachDt, dt_teach);
end

t_teach = t_all(1:teach_stride:end);
if abs(t_teach(end) - t_all(end)) > 1e-12
    t_teach = [t_teach; t_all(end)];
end
NUM_TEACH_SAMPLES = length(t_teach);

% ★ 关键点覆盖约束
if isempty(opts.MaxKeyGap)
    max_key_gap = max(seg_dur);   % 当前 nominal 参数下为 2 s
else
    max_key_gap = opts.MaxKeyGap;
end
min_key_spacing = max(opts.MinKeySpacing, 2*dt_teach);

validateattributes(max_key_gap, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(min_key_spacing, {'numeric'}, {'scalar','real','finite','nonnegative'});

% 检查 MaxPoints 是否足以覆盖整个周期
min_points_for_coverage = ceil(T_total / max_key_gap) + 1;
if opts.MaxPoints < min_points_for_coverage
    error('TeachingDataPipeline:InsufficientMaxPoints', ...
        ['MaxPoints=%d 不足以满足 MaxKeyGap=%.3f s 的全周期覆盖约束；' ...
         '至少需要 %d 个关键点。'], ...
        opts.MaxPoints, max_key_gap, min_points_for_coverage);
end

joint_bc = {params.traj1_deg, params.traj2_deg, params.traj3_deg};
q_ref = zeros(length(t_all),3);
dq_ref = zeros(length(t_all),3);
ddq_ref = zeros(length(t_all),3);
for j=1:3
    [q_ref(:,j), dq_ref(:,j), ddq_ref(:,j)] = ...
        multi_segment_quintic(joint_bc{j}, seg_dur, t_all);
end

% ---------- 3. 模拟编码器采集（低频示教采样 + 时间扭曲 + 独立有色噪声） ----------
% 示教测量只在 t_teach 上采样，而不是在每个控制周期采样。
% t_warp_source(k) 表示名义示教采样时刻 t_teach(k) 对应的源轨迹时刻。
t_warp_source = t_teach + opts.WarpAmp*sin(opts.WarpFreq*t_teach);
t_warp_source(1)   = t_all(1);
t_warp_source(end) = t_all(end);

% 防止参数改变后出现越界或时间逆序。
t_warp_source = min(max(t_warp_source, t_all(1)), t_all(end));
if any(diff(t_warp_source) <= 0)
    error('TeachingDataPipeline:NonMonotonicTimeWarp', ...
          ['时间扭曲后源时间轴不是严格单调递增。请减小 WarpAmp 或 WarpFreq，' ...
           '以避免不物理的时间折返。']);
end

% 先在时间扭曲后的源时刻读取潜在真实轨迹，再加入编码器测量噪声。
q_teach_nom   = interp1(t_all, q_ref,   t_warp_source, 'pchip');
dq_teach_nom  = interp1(t_all, dq_ref,  t_warp_source, 'pchip');
ddq_teach_nom = interp1(t_all, ddq_ref, t_warp_source, 'pchip');

% 使用示教管线局部随机流；噪声长度按示教采样点数生成。
b=0.5; a=[1,-0.5];
noise_q = filter(b,a,randn(teaching_stream,NUM_TEACH_SAMPLES,3));
noise_q = opts.NoiseLevel * (noise_q./std(noise_q,0,1));

noise_dq = filter(b,a,randn(teaching_stream,NUM_TEACH_SAMPLES,3));
noise_dq = opts.NoiseVel * (noise_dq./std(noise_dq,0,1));

noise_ddq = filter(b,a,randn(teaching_stream,NUM_TEACH_SAMPLES,3));
noise_ddq = opts.NoiseAcc * (noise_ddq./std(noise_ddq,0,1));

% 同一噪声下的未时间扭曲示教序列，仅用于审计/诊断。
q_teach_unwarped   = interp1(t_all, q_ref,   t_teach, 'pchip') + noise_q;
dq_teach_unwarped  = interp1(t_all, dq_ref,  t_teach, 'pchip') + noise_dq;
ddq_teach_unwarped = interp1(t_all, ddq_ref, t_teach, 'pchip') + noise_ddq;

q_teaching   = q_teach_nom   + noise_q;
dq_teaching  = dq_teach_nom  + noise_dq;
ddq_teaching = ddq_teach_nom + noise_ddq;

% ---------- 3.1 运动学一致的导数代理 ----------
% 位置是编码器的主观测量；若把独立加噪的 dq/ddq 直接作为 quintic 的硬边界，
% 则 q、dq、ddq 三者并不满足导数一致性，短分段会把这种矛盾放大成窄峰。
% 因此：关键点位置仍来自真实示教位置 q_teaching；用于重建的 dq/ddq 则
% 从轻度平滑后的同一位置序列求导，保证三者来自同一个运动学信号。
kin_win = 2*floor(0.40/dt_teach) + 1;  % 默认约0.9 s窗口(10 Hz时为9点)
kin_win = max(kin_win, 5);
if kin_win >= NUM_TEACH_SAMPLES
    kin_win = NUM_TEACH_SAMPLES - 1 + mod(NUM_TEACH_SAMPLES,2);
    if mod(kin_win,2)==0, kin_win = kin_win - 1; end
end
if kin_win <= 3
    kin_order = 2;
else
    kin_order = 3;
end

q_kinematic   = q_teaching;
dq_kinematic  = zeros(size(q_teaching));
ddq_kinematic = zeros(size(q_teaching));
for j = 1:3
    if kin_win >= 5
        q_kinematic(:,j) = sgolayfilt(q_teaching(:,j), kin_order, kin_win);
    end
    dq_kinematic(:,j) = gradient(q_kinematic(:,j), t_teach);
    ddq_kinematic(:,j) = gradient(dq_kinematic(:,j), t_teach);
end

% 已知周期静止端点在导数代理中也显式固定。
for j = 1:3
    base_j = params.(['traj' num2str(j) '_deg']);
    q_kinematic(1,j) = base_j(1,1);
    q_kinematic(end,j) = base_j(end,1);
    dq_kinematic(1,j) = base_j(1,2);
    dq_kinematic(end,j) = base_j(end,2);
    ddq_kinematic(1,j) = base_j(1,3);
    ddq_kinematic(end,j) = base_j(end,3);
end

% 用实际示教位置数据构造控制网格上的稠密观测代理，仅用于关键点误差驱动。
% 这样关键点选择不再依赖仿真中不可观测的 latent ground truth q_ref。
q_teaching_dense = interp1(t_teach, q_teaching, t_all, 'pchip');

if any(~isfinite(q_teaching(:))) || any(~isfinite(dq_teaching(:))) || ...
        any(~isfinite(ddq_teaching(:)))
    error('TeachingDataPipeline:InvalidWarpedTeachingData', ...
          '时间扭曲示教数据包含 NaN/Inf，请检查 WarpAmp/WarpFreq 和时间范围。');
end

% ---------- 4. 关键点自适应抽取：物理特征 + 全周期覆盖 + 误差驱动 ----------
% 第2项补丁：防止候选点因噪声过零而在局部密集，导致关键点扎堆和长时间空段。

% 4.1 仅用于“特征检测”的轻度平滑（不改变 q_teaching/dq_teaching/ddq_teaching 本身）
feature_win = 2*floor(0.10/dt_teach) + 1;  % 短时特征窗口；10 Hz 下受最小5点约束
feature_win = max(feature_win, 5);
if feature_win >= NUM_TEACH_SAMPLES
    feature_win = NUM_TEACH_SAMPLES - 1 + mod(NUM_TEACH_SAMPLES,2);
    if mod(feature_win,2)==0, feature_win = feature_win - 1; end
end

dq_feature  = dq_kinematic;
ddq_feature = ddq_kinematic;
if feature_win >= 5
    for j = 1:3
        dq_feature(:,j)  = sgolayfilt(dq_kinematic(:,j),  3, feature_win);
        ddq_feature(:,j) = sgolayfilt(ddq_kinematic(:,j), 3, feature_win);
    end
end

% 4.2 物理候选点：速度/加速度/Jerk过零点
candidate_t = [];
for j=1:3
    sig_v = sign(dq_feature(:,j));
    idx_v = find(diff(sig_v)~=0);
    candidate_t = [candidate_t; (t_teach(idx_v)+t_teach(idx_v+1))/2];

    sig_a = sign(ddq_feature(:,j));
    idx_a = find(diff(sig_a)~=0);
    candidate_t = [candidate_t; (t_teach(idx_a)+t_teach(idx_a+1))/2];

    dddq_feature = gradient(ddq_feature(:,j), t_teach);
    sig_j = sign(dddq_feature);
    idx_j = find(diff(sig_j)~=0);
    candidate_t = [candidate_t; (t_teach(idx_j)+t_teach(idx_j+1))/2];
end

candidate_t = [t_all(1); candidate_t(:); t_all(end)];
candidate_t = unique(snap_times_to_grid(candidate_t, t_all, dt_ctrl));

% 候选点先做最小间距去重，降低噪声造成的局部候选密度偏置
candidate_t = enforce_min_spacing(candidate_t, min_key_spacing, t_all(1), t_all(end));

% 4.3 初始关键点：按“时间轴”覆盖，而不是按 candidate_t 的数组索引覆盖
% 每个均匀时间目标优先吸附到附近物理候选点；附近无候选时直接使用网格目标点。
n_init = max(MinPoints, min_points_for_coverage);
uniform_targets = linspace(t_all(1), t_all(end), n_init)';
uniform_targets = unique(snap_times_to_grid(uniform_targets, t_all, dt_ctrl));

key_t = [t_all(1); t_all(end)];
coverage_radius = 0.45 * (T_total / max(n_init-1,1));

for k = 2:length(uniform_targets)-1
    target = uniform_targets(k);

    cand_local = candidate_t( ...
        candidate_t > key_t(1) + min_key_spacing & ...
        candidate_t < key_t(end) - min_key_spacing & ...
        abs(candidate_t - target) <= coverage_radius);

    if ~isempty(cand_local)
        [~, idx_near] = min(abs(cand_local - target));
        tc = cand_local(idx_near);
    else
        tc = target;
    end

    key_t = sort(unique([key_t; tc]));
end

% 再次保证最小间距；端点永久保留
key_t = enforce_min_spacing(key_t, min_key_spacing, t_all(1), t_all(end));

% 4.4 强制最大时间间隔：只要存在 > MaxKeyGap 的空段，就优先在该段补点
while max(diff(key_t)) > max_key_gap + 1e-12
    if length(key_t) >= opts.MaxPoints
        error('TeachingDataPipeline:CoverageConstraintFailed', ...
              ['关键点已达到 MaxPoints=%d，但最大空段仍为 %.3f s，' ...
               '超过 MaxKeyGap=%.3f s。'], ...
              opts.MaxPoints, max(diff(key_t)), max_key_gap);
    end

    gaps = diff(key_t);
    [~, igap] = max(gaps);
    tL = key_t(igap);
    tR = key_t(igap+1);
    target = 0.5*(tL+tR);

    cand_seg = candidate_t(candidate_t > tL + min_key_spacing & ...
                           candidate_t < tR - min_key_spacing);

    if ~isempty(cand_seg)
        [~, idx_near] = min(abs(cand_seg-target));
        tc = cand_seg(idx_near);
    else
        tc = snap_times_to_grid(target, t_all, dt_ctrl);
    end

    if any(abs(key_t-tc) < dt_ctrl/2)
        error('TeachingDataPipeline:CannotInsertCoveragePoint', ...
              '无法在 %.3f~%.3f s 区间继续插入有效覆盖关键点。', tL, tR);
    end

    key_t = sort([key_t; tc]);
end

% 4.5 误差驱动贪婪加点
% 与原算法不同：先确定“当前误差最大的时间段”，再在该段选一个物理候选点，
% 从根本上避免候选点数量多的局部区域连续抢占多个关键点。
N_max = opts.MaxPoints;
eps_tol = opts.EpsTol;

while length(key_t) < N_max
    seg_dur_cur = diff(key_t);
    if any(seg_dur_cur <= 1e-12), break; end

    q_key   = interp1(t_teach, q_teaching,   key_t, 'pchip');
    dq_key  = interp1(t_teach, dq_kinematic,  key_t, 'pchip');
    ddq_key = interp1(t_teach, ddq_kinematic, key_t, 'pchip');

    q_hat = zeros(length(t_all),3);
    for j=1:3
        bc_cur = [q_key(:,j), dq_key(:,j), ddq_key(:,j)];
        [q_hat(:,j),~,~] = multi_segment_quintic(bc_cur, seg_dur_cur, t_all);
    end

    error_q_global = q_hat - q_teaching_dense;
    rmse_fit_global = sqrt(mean(error_q_global(:).^2));

    % 达到对示教观测代理的拟合精度，同时覆盖约束已满足，则停止
    if rmse_fit_global <= eps_tol && max(diff(key_t)) <= max_key_gap + 1e-12
        break;
    end

    % 每个现有时间段计算积分误差
    nseg_key = length(key_t)-1;
    seg_error = zeros(nseg_key,1);
    for s = 1:nseg_key
        idx_seg = (t_all >= key_t(s)) & (t_all <= key_t(s+1));
        if sum(idx_seg) >= 2
            e_seg = sum(abs(error_q_global(idx_seg,:)),2);
            seg_error(s) = trapz(t_all(idx_seg), e_seg);
        end
    end

    % 优先选择积分误差最大的段
    [~, s_best] = max(seg_error);
    tL = key_t(s_best);
    tR = key_t(s_best+1);

    idx_seg = find(t_all >= tL & t_all <= tR);
    point_error = sum(abs(error_q_global(idx_seg,:)),2);
    [~, ip] = max(point_error);
    target = t_all(idx_seg(ip));

    % 在该段内、且离已有关键点足够远的物理候选点中，选最接近最大误差位置者
    cand_seg = candidate_t(candidate_t > tL + min_key_spacing & ...
                           candidate_t < tR - min_key_spacing);

    if ~isempty(cand_seg)
        [~, idx_near] = min(abs(cand_seg-target));
        tc = cand_seg(idx_near);
    else
        tc = snap_times_to_grid(target, t_all, dt_ctrl);
        % 最大误差点太靠近端点时退回区间中点
        if tc <= tL + min_key_spacing || tc >= tR - min_key_spacing
            tc = snap_times_to_grid(0.5*(tL+tR), t_all, dt_ctrl);
        end
    end

    if any(abs(key_t-tc) < min_key_spacing - 1e-12)
        % 当前高误差段无法继续安全插点，屏蔽该段后尝试次高误差段
        inserted = false;
        [~, order_seg] = sort(seg_error, 'descend');
        for ss = order_seg(:)'
            tLs = key_t(ss);
            tRs = key_t(ss+1);
            if (tRs-tLs) <= 2*min_key_spacing
                continue;
            end
            cand_s = candidate_t(candidate_t > tLs + min_key_spacing & ...
                                 candidate_t < tRs - min_key_spacing);
            if ~isempty(cand_s)
                tc2 = cand_s(round((length(cand_s)+1)/2));
            else
                tc2 = snap_times_to_grid(0.5*(tLs+tRs), t_all, dt_ctrl);
            end
            if all(abs(key_t-tc2) >= min_key_spacing - 1e-12)
                tc = tc2;
                inserted = true;
                break;
            end
        end
        if ~inserted
            break;
        end
    end

    key_t = sort(unique([key_t; tc]));
end

% 4.6 最终覆盖检查
key_t = unique(snap_times_to_grid(key_t, t_all, dt_ctrl));
key_t(1)=t_all(1);
key_t(end)=t_all(end);

if max(diff(key_t)) > max_key_gap + 1e-12
    error('TeachingDataPipeline:FinalCoverageViolation', ...
          '最终关键点最大间隔 %.3f s 超过 MaxKeyGap %.3f s。', ...
          max(diff(key_t)), max_key_gap);
end

seg_dur_new = diff(key_t);

% 插值得到关键点状态。位置来自示教位置；速度/加速度来自同一位置序列的
% 运动学一致导数代理，避免独立噪声导数被作为 quintic 硬约束。
key_q_deg   = interp1(t_teach, q_teaching,    key_t, 'pchip');
key_dq_deg  = interp1(t_teach, dq_kinematic,  key_t, 'pchip');
key_ddq_deg = interp1(t_teach, ddq_kinematic, key_t, 'pchip');

% ---------- 4.7 周期边界物理约束 ----------
% 该项目的基础被动训练模板是已知的周期运动：首、末关键状态均为
% 规定的静止/闭合状态。示教噪声和 SG 端点效应不应改变这一先验。
% 因此仅锁定周期首尾 q/dq/ddq；内部关键点仍完全来自示教数据。
base_bc = {params.traj1_deg, params.traj2_deg, params.traj3_deg};
for j = 1:3
    key_q_deg(1,j)    = base_bc{j}(1,1);
    key_dq_deg(1,j)   = base_bc{j}(1,2);
    key_ddq_deg(1,j)  = base_bc{j}(1,3);

    key_q_deg(end,j)   = base_bc{j}(end,1);
    key_dq_deg(end,j)  = base_bc{j}(end,2);
    key_ddq_deg(end,j) = base_bc{j}(end,3);
end

% 周期闭合诊断：对本项目循环训练，首尾状态应严格一致。
boundary_state_gap = zeros(3,3); % rows: q,dq,ddq; cols: joints
for j = 1:3
    boundary_state_gap(:,j) = [ ...
        key_q_deg(end,j)   - key_q_deg(1,j); ...
        key_dq_deg(end,j)  - key_dq_deg(1,j); ...
        key_ddq_deg(end,j) - key_ddq_deg(1,j)];
end
if max(abs(boundary_state_gap(:))) > 1e-10
    warning('TeachingDataPipeline:NonPeriodicBoundary', ...
        ['示教模板首尾状态并非严格周期闭合；当前仍保留各自已知边界。' ...
         '若用于循环训练，请检查 parameters.m 的首末 q/dq/ddq。']);
end

% ---------- 5. 形成最终关键点参数 ----------
% ★ 第2项修复：最终重建不再使用 TeachingDataPipeline 内部的另一套轨迹生成器。
% 先把提取出的关键点写入 params_updated，再统一调用 trajectory_generator。
params_updated = params;
params_updated.traj1_deg = [key_q_deg(:,1), key_dq_deg(:,1), key_ddq_deg(:,1)];
params_updated.traj2_deg = [key_q_deg(:,2), key_dq_deg(:,2), key_ddq_deg(:,2)];
params_updated.traj3_deg = [key_q_deg(:,3), key_dq_deg(:,3), key_ddq_deg(:,3)];
params_updated.teaching_pipeline_flag = true;
params_updated.traj_times = seg_dur_new;
params_updated.teaching_seed = teaching_seed;

if isfield(params_updated, 'T_seg')
    params_updated.T_seg = mean(seg_dur_new);
end

% ---------- 6. 用“控制器实际轨迹生成器”生成唯一权威参考轨迹 ----------
% 这一份数据稍后由 main_simulation 直接送入 fixed/adaptive controller，
% 不再在 main_simulation 中二次调用 trajectory_generator。
[q_controller_ref, dq_controller_ref, ddq_controller_ref, t_controller_ref] = ...
    trajectory_generator(params_updated);

% 单周期重建轨迹 = 控制器实际输入轨迹的第一个周期（严格同一数组的切片）
N_cycle = length(t_all);
if size(q_controller_ref,1) < N_cycle
    error('TeachingDataPipeline:ControllerReferenceTooShort', ...
          '控制器参考轨迹长度小于单周期示教评估长度。');
end

t_rebuild   = t_controller_ref(1:N_cycle);
q_rebuilt   = rad2deg(q_controller_ref(1:N_cycle,:));
dq_rebuilt  = rad2deg(dq_controller_ref(1:N_cycle,:));
ddq_rebuilt = rad2deg(ddq_controller_ref(1:N_cycle,:));

% 严格检查：控制器首周期时间网格必须与示教评估网格完全一致
time_grid_error = max(abs(t_rebuild - t_all));
if time_grid_error > 1e-12
    error('TeachingDataPipeline:ReferenceGridNotIdentical', ...
          '示教评估与控制器参考轨迹时间网格不一致，最大差值 %.3e s。', ...
          time_grid_error);
end

% ---------- 7. 基于“真正送入控制器”的轨迹计算重建误差 ----------
% q_ref 是原始完美示教真值（deg），q_rebuilt 是控制器真实使用的首周期参考（deg）
error_q = q_rebuilt - q_ref;
rmse_global = sqrt(mean(error_q(:).^2));
rmse = sqrt(mean(error_q.^2,1));
mae  = mean(abs(error_q),1);

% 另存对示教观测代理的拟合误差，用于审计关键点选择目标；不替代真值评估指标。
error_q_fit = q_rebuilt - q_teaching_dense;
rmse_fit_global = sqrt(mean(error_q_fit(:).^2));

if opts.Verbose
    fprintf('\n========== 自适应关键点抽取完成 ==========\n');
    fprintf('最终关键点数: %d\n', length(key_t));
    fprintf('重建全局RMSE（相对latent truth，仅用于仿真评估）: %.4f°\n', rmse_global);
    fprintf('关键点拟合RMSE（相对示教观测代理）: %.4f°\n', rmse_fit_global);
    fprintf('参考轨迹一致性检查: PASS（评估轨迹 = 控制器实际输入首周期）\n');
end

% ---------- 7.1 打包指标 ----------
metrics = struct();
metrics.teaching_seed = teaching_seed;
metrics.reference_pipeline = 'TeachingDataPipeline';
metrics.reference_pipeline_version = 'v6_periodic_boundary_kinematic_derivatives';
metrics.teaching_config = teach_cfg;
metrics.num_keypoints = length(key_t);
metrics.rmse_global   = rmse_global;
metrics.rmse_fit_teaching_global = rmse_fit_global;
metrics.selection_target = 'PCHIP-dense teaching position observations';
metrics.rmse_hip_deg   = rmse(1);
metrics.rmse_knee_deg  = rmse(2);
metrics.rmse_ankle_deg = rmse(3);
metrics.mae_hip_deg    = mae(1);
metrics.mae_knee_deg   = mae(2);
metrics.mae_ankle_deg  = mae(3);

% 保持原字段语义：单周期、单位为 deg / deg/s / deg/s^2
metrics.t_rebuild      = t_rebuild;
metrics.q_rebuilt      = q_rebuilt;
metrics.dq_rebuilt     = dq_rebuilt;
metrics.ddq_rebuilt    = ddq_rebuilt;

metrics.q_true         = q_ref;
metrics.dq_true        = dq_ref;
metrics.ddq_true       = ddq_ref;

metrics.t_teaching     = t_teach;
metrics.t_warp_source  = t_warp_source;
metrics.teaching_sample_dt = dt_teach;
metrics.teaching_sample_hz = 1/dt_teach;
metrics.num_teaching_samples = NUM_TEACH_SAMPLES;
metrics.q_noisy_unwarped = q_teach_unwarped;
metrics.dq_noisy_unwarped = dq_teach_unwarped;
metrics.ddq_noisy_unwarped = ddq_teach_unwarped;
metrics.q_teaching     = q_teaching;
metrics.q_teaching_dense = q_teaching_dense;
metrics.dq_teaching    = dq_teaching;
metrics.ddq_teaching   = ddq_teaching;
metrics.q_kinematic_derivative_source = q_kinematic;
metrics.dq_kinematic = dq_kinematic;
metrics.ddq_kinematic = ddq_kinematic;
metrics.kinematic_derivative_window_samples = kin_win;
metrics.derivative_boundary_mode = 'derived_from_smoothed_teaching_position';

metrics.keypoints_t       = key_t;
metrics.keypoints_q_deg   = key_q_deg;
metrics.keypoints_dq_deg  = key_dq_deg;
metrics.keypoints_ddq_deg = key_ddq_deg;
metrics.boundary_state_gap_deg = boundary_state_gap;
metrics.periodic_boundary_constrained = true;

% 新增：保存完整的、真正送入控制器的多周期参考轨迹（SI单位）
metrics.t_controller_ref   = t_controller_ref;
metrics.q_controller_ref   = q_controller_ref;      % rad
metrics.dq_controller_ref  = dq_controller_ref;     % rad/s
metrics.ddq_controller_ref = ddq_controller_ref;    % rad/s^2
metrics.time_grid_identity_error = time_grid_error;
metrics.keypoint_intervals_s = diff(key_t);
metrics.max_keypoint_gap_s = max(diff(key_t));
metrics.min_keypoint_gap_s = min(diff(key_t));
metrics.max_key_gap_limit_s = max_key_gap;
metrics.min_key_spacing_limit_s = min_key_spacing;

% ---------- 8. 打印关键点 ----------
if opts.Verbose
    fprintf('\n========== 示教数据处理：自适应关键点（单位：°，°/s，°/s²）==========\n');
    joint_labels = {'髋关节(Hip)','膝关节(Knee)','踝关节(Ankle)'};
    for j=1:3
        fprintf('\n▸ %s (%d个关键点)\n', joint_labels{j}, size(key_q_deg,1));
        fprintf('   %-10s %-12s %-12s %-12s\n','序号','角度(°)','角速度(°/s)','角加速度(°/s²)');
        for k=1:size(key_q_deg,1)
            fprintf('   %-10d %-12.4f %-12.4f %-12.4f\n',...
                k,key_q_deg(k,j),key_dq_deg(k,j),key_ddq_deg(k,j));
        end
    end
    fprintf('   (关键点时间: '); fprintf('%.3f ',key_t); fprintf(')\n');
    fprintf('   (关键点间隔: '); fprintf('%.3f ',diff(key_t)); fprintf(')\n');
    fprintf('   (最大间隔 %.3f s / 限制 %.3f s；最小间隔 %.3f s)\n', ...
            max(diff(key_t)), max_key_gap, min(diff(key_t)));
    fprintf('   (全局RMSE=%.4f°, MinPoints=%d, MaxPoints=%d, EpsTol=%.2f)\n',...
            rmse_global, MinPoints, N_max, eps_tol);
    fprintf('   (示教采样: %.1f Hz, %d points；控制采样: %.1f Hz, %d points/cycle)\n', ...
            1/dt_teach, NUM_TEACH_SAMPLES, 1/dt_ctrl, NUM_CTRL_SAMPLES);

end

% ---------- 9. 可视化（显示与保存解耦） ----------
if opts.PlotFlag
    % 与 generate_paper_figures.m 使用同一论文级布局/导出策略：
    % 显式 axes 位置、白色背景、分图统一对齐、PDF 页面真实留白。
    [fig_teaching, lay_teaching] = paper_figure_utils('create', ...
        'Teaching Data Pipeline',[100 45 1050 1180],3,1, ...
        sprintf('Adaptive Key-State Reconstruction (N=%d, RMSE=%.4f deg)', ...
                length(key_t),rmse_global));
    joint_names={'Hip','Knee','Ankle'};
    axs = gobjects(3,1); lgds = gobjects(3,1);
    for j=1:3
        ax = paper_figure_utils('tile',lay_teaching,j); axs(j)=ax;
        plot(ax,t_all,q_ref(:,j),'k--','LineWidth',1.5,'DisplayName','Perfect Ref'); hold(ax,'on');
        scatter(ax,t_teach,q_teaching(:,j),10,'b','filled','MarkerFaceAlpha',0.45, ...
            'DisplayName','Teaching Data');
        plot(ax,t_rebuild,q_rebuilt(:,j),'r-','LineWidth',1.2,'DisplayName','Rebuilt');
        scatter(ax,key_t,key_q_deg(:,j),40,'m','v','filled','MarkerFaceAlpha',0.6, ...
            'DisplayName','Keypoints');
        ylabel(ax,'Angle (deg)');
        if j==3, xlabel(ax,'Time (s)'); end
        lgds(j)=legend(ax,'Location','northeast','AutoUpdate','off');
        paper_figure_utils('panel_title',ax,sprintf('(%c) %s Joint',char('a'+j-1),joint_names{j}));
    end
    paper_figure_utils('align_legends',axs,lgds,'NumColumns',2);

    if opts.SaveFlag
        paper_figure_utils('export',fig_teaching,saveDir, ...
            ['teaching_pipeline_' timestampStr],'SavePNG',true);
        if opts.Verbose
            fprintf('   ✅ 可视化 PDF/PNG 保存至: %s\n', saveDir);
        end
    end
end

% ---------- 10. 保存结果（MAT + CSV 四份） ----------
if opts.SaveFlag
    % 10a. 保存 metrics 结构体（MAT）
    matFile = fullfile(saveDir, ['teaching_metrics_' timestampStr '.mat']);
    save(matFile, 'metrics');
    fprintf('   ✅ 度量数据保存至: %s\n', matFile);

    % 10b. 示教数据 CSV
    T_teaching = table(t_teach, t_warp_source, ...
        q_teaching(:,1), dq_teaching(:,1), ddq_teaching(:,1), ...
        q_teaching(:,2), dq_teaching(:,2), ddq_teaching(:,2), ...
        q_teaching(:,3), dq_teaching(:,3), ddq_teaching(:,3), ...
        'VariableNames', {'Time_s','WarpedSourceTime_s', ...
            'Hip_Angle_deg', 'Hip_Vel_deg_s', 'Hip_Acc_deg_s2', ...
            'Knee_Angle_deg','Knee_Vel_deg_s','Knee_Acc_deg_s2', ...
            'Ankle_Angle_deg','Ankle_Vel_deg_s','Ankle_Acc_deg_s2'});
    csvFile = fullfile(saveDir, ['teaching_data_' timestampStr '.csv']);
    writetable(T_teaching, csvFile);
    fprintf('   ✅ 示教数据保存至: %s\n', csvFile);

    % 10c. 关键点 CSV（已有）
    T_keypoints = table(key_t, ...
        key_q_deg(:,1), key_dq_deg(:,1), key_ddq_deg(:,1), ...
        key_q_deg(:,2), key_dq_deg(:,2), key_ddq_deg(:,2), ...
        key_q_deg(:,3), key_dq_deg(:,3), key_ddq_deg(:,3), ...
        'VariableNames', {'Time_s', ...
            'Hip_Angle_deg','Hip_Vel_deg_s','Hip_Acc_deg_s2', ...
            'Knee_Angle_deg','Knee_Vel_deg_s','Knee_Acc_deg_s2', ...
            'Ankle_Angle_deg','Ankle_Vel_deg_s','Ankle_Acc_deg_s2'});
    csvFile = fullfile(saveDir, ['teaching_keypoints_' timestampStr '.csv']);
    writetable(T_keypoints, csvFile);
    fprintf('   ✅ 关键点数据保存至: %s\n', csvFile);

    % 10d. 重建拟合数据 CSV（包含重建的速度、加速度）
    T_fitted = table(t_rebuild, ...
        q_rebuilt(:,1), dq_rebuilt(:,1), ddq_rebuilt(:,1), ...
        q_rebuilt(:,2), dq_rebuilt(:,2), ddq_rebuilt(:,2), ...
        q_rebuilt(:,3), dq_rebuilt(:,3), ddq_rebuilt(:,3), ...
        'VariableNames', {'Time_s', ...
            'Hip_Angle_deg','Hip_Vel_deg_s','Hip_Acc_deg_s2', ...
            'Knee_Angle_deg','Knee_Vel_deg_s','Knee_Acc_deg_s2', ...
            'Ankle_Angle_deg','Ankle_Vel_deg_s','Ankle_Acc_deg_s2'});
    csvFile = fullfile(saveDir, ['teaching_fitted_' timestampStr '.csv']);
    writetable(T_fitted, csvFile);
    fprintf('   ✅ 重建拟合数据保存至: %s\n', csvFile);

    % 10e. 关键点拟合数据（重建轨迹在关键点处的值）
    % 直接从重建结果中插值
    q_key_fit   = interp1(t_all, q_rebuilt,   key_t, 'pchip');
    dq_key_fit  = interp1(t_all, dq_rebuilt,  key_t, 'pchip');
    ddq_key_fit = interp1(t_all, ddq_rebuilt, key_t, 'pchip');

    T_key_fit = table(key_t, ...
        q_key_fit(:,1), dq_key_fit(:,1), ddq_key_fit(:,1), ...
        q_key_fit(:,2), dq_key_fit(:,2), ddq_key_fit(:,2), ...
        q_key_fit(:,3), dq_key_fit(:,3), ddq_key_fit(:,3), ...
        'VariableNames', {'Time_s', ...
            'Hip_Angle_deg','Hip_Vel_deg_s','Hip_Acc_deg_s2', ...
            'Knee_Angle_deg','Knee_Vel_deg_s','Knee_Acc_deg_s2', ...
            'Ankle_Angle_deg','Ankle_Vel_deg_s','Ankle_Acc_deg_s2'});
    csvFile = fullfile(saveDir, ['teaching_keypoints_fitted_' timestampStr '.csv']);
    writetable(T_key_fit, csvFile);
    fprintf('   ✅ 关键点拟合数据保存至: %s\n', csvFile);

end

end % 主函数

%% ================= 辅助函数：示教参数默认值 =================
function cfg = teaching_defaults(params)
    % 优先读取 parameters.m 中的 params.teaching，并补齐缺失字段。
    if isfield(params,'teaching') && isstruct(params.teaching)
        cfg = params.teaching;
    else
        cfg = struct();
    end

    defaults = struct();
    defaults.noise_level = 0.03;
    defaults.noise_vel = 0.5;
    defaults.noise_acc = 5.0;
    defaults.warp_amp = 0.05;
    defaults.warp_freq = 2*pi/0.8;
    defaults.sample_dt = 0.10;      % 10 Hz 示教采样
    defaults.eps_tol = 1.0;
    defaults.max_points = 10;
    defaults.min_points = 6;
    defaults.max_key_gap = params.T_seg;
    defaults.min_key_spacing = 0.20;

    names = fieldnames(defaults);
    for ii = 1:numel(names)
        f = names{ii};
        if ~isfield(cfg, f) || isempty(cfg.(f))
            cfg.(f) = defaults.(f);
        end
    end
end

%% ================= 辅助函数：关键点最小时间间距 =================
function t_out = enforce_min_spacing(t_in, min_spacing, t_start, t_end)
    t = sort(unique(t_in(:)));
    if isempty(t)
        t_out = [t_start; t_end];
        return;
    end

    % 强制端点存在
    t = sort(unique([t_start; t; t_end]));

    keep = false(size(t));
    keep(1) = true;
    last_kept = t(1);

    % 中间点按时间顺序保留
    for i = 2:length(t)-1
        if (t(i)-last_kept) >= min_spacing && ...
           (t_end-t(i)) >= min_spacing
            keep(i) = true;
            last_kept = t(i);
        end
    end
    keep(end) = true;

    t_out = t(keep);

    % 若最后一个中间点离终点过近，则去掉该中间点
    if length(t_out) >= 3 && (t_out(end)-t_out(end-1)) < min_spacing
        t_out(end-1) = [];
    end
end

%% ================= 辅助函数：时间吸附到控制器网格 =================
function t_snap = snap_times_to_grid(t_in, t_grid, dt)
    idx = round(t_in(:) / dt) + 1;
    idx = max(1, min(length(t_grid), idx));
    t_snap = t_grid(idx);
end

%% ================= 辅助函数 multi_segment_quintic =================
function [q_seg,dq_seg,ddq_seg]=multi_segment_quintic(bc,seg_dur,t_eval)
    N=size(bc,1); t_end=[0;cumsum(seg_dur(:))]; num_seg=N-1;
    q_seg=zeros(size(t_eval)); dq_seg=zeros(size(t_eval)); ddq_seg=zeros(size(t_eval));
    seg_idx=discretize(t_eval,t_end);
    seg_idx(t_eval==t_end(end))=num_seg;
    valid=(seg_idx>=1)&(seg_idx<=num_seg);
    if ~any(valid); return; end
    t_rel=t_eval(valid)-t_end(seg_idx(valid));
    T_vec=seg_dur(seg_idx(valid)); tau=t_rel./T_vec;
    H0=10*tau.^3-15*tau.^4+6*tau.^5;
    H1=tau-6*tau.^3+8*tau.^4-3*tau.^5;
    H2=-4*tau.^3+7*tau.^4-3*tau.^5;
    H3=0.5*(tau.^2-3*tau.^3+3*tau.^4-tau.^5);
    H4=0.5*(tau.^3-2*tau.^4+tau.^5);
    q0=bc(seg_idx(valid),1); v0=bc(seg_idx(valid),2); a0=bc(seg_idx(valid),3);
    qf=bc(seg_idx(valid)+1,1); vf=bc(seg_idx(valid)+1,2); af=bc(seg_idx(valid)+1,3);
    q_loc=q0+(qf-q0).*H0+v0.*T_vec.*H1+vf.*T_vec.*H2+a0.*T_vec.^2.*H3+af.*T_vec.^2.*H4;
    dH0=(30*tau.^2-60*tau.^3+30*tau.^4)./T_vec;
    dH1=1-18*tau.^2+32*tau.^3-15*tau.^4;
    dH2=-12*tau.^2+28*tau.^3-15*tau.^4;
    dH3=(tau-4.5*tau.^2+6*tau.^3-2.5*tau.^4).*T_vec;
    dH4=(1.5*tau.^2-4*tau.^3+2.5*tau.^4).*T_vec;
    dq_loc=(qf-q0).*dH0+v0.*dH1+vf.*dH2+a0.*dH3+af.*dH4;
    ddH0=(60*tau-180*tau.^2+120*tau.^3)./T_vec.^2;
    ddH1=(-36*tau+96*tau.^2-60*tau.^3)./T_vec;
    ddH2=(-24*tau+84*tau.^2-60*tau.^3)./T_vec;
    ddH3=1-9*tau+18*tau.^2-10*tau.^3;
    ddH4=3*tau-12*tau.^2+10*tau.^3;
    ddq_loc=(qf-q0).*ddH0+v0.*ddH1+vf.*ddH2+a0.*ddH3+af.*ddH4;
    q_seg(valid)=q_loc; dq_seg(valid)=dq_loc; ddq_seg(valid)=ddq_loc;
    q_seg(t_eval==0)=bc(1,1); dq_seg(t_eval==0)=bc(1,2); ddq_seg(t_eval==0)=bc(1,3);
    q_seg(t_eval==t_eval(end))=bc(end,1); dq_seg(t_eval==t_eval(end))=bc(end,2); ddq_seg(t_eval==t_eval(end))=bc(end,3);
end

% ===== Internal quintic trajectory generator (formerly separate file) =====
function [q_ref, dq_ref, ddq_ref, t_all] = trajectory_generator(params)
    % 生成多关节参考轨迹（对应论文 Results 中 Fig.1 参考轨迹）
    % 已优化：支持示教管线提取的不均匀关键点时长 traj_times

    dt = params.dt;
    T_seg = params.T_seg;

    % 转换为弧度
    traj1 = deg2rad(params.traj1_deg);  % Hip
    traj2 = deg2rad(params.traj2_deg);  % Knee
    traj3 = deg2rad(params.traj3_deg);  % Ankle

    % ★ 关键修改：优先使用示教管线提取的 traj_times ★
    if isfield(params, 'traj_times') && ~isempty(params.traj_times)
        seg_durs = params.traj_times(:)';   % 1×n_seg
        if length(seg_durs) ~= size(traj1,1)-1
            warning('traj_times 长度(%d)与轨迹段数(%d)不匹配，回退到 T_seg', ...
                    length(seg_durs), size(traj1,1)-1);
            seg_durs = T_seg * ones(1, size(traj1,1)-1);
        end
    else
        seg_durs = T_seg * ones(1, size(traj1,1)-1);
    end

    % ★ 第2项修复：离散控制要求每个关键点时刻落在 dt 网格上
    grid_ratio = seg_durs / dt;
    if any(abs(grid_ratio - round(grid_ratio)) > 1e-9)
        error('trajectory_generator:SegmentNotOnControlGrid', ...
              ['存在轨迹段时长不能被 dt=%.12g s 整除。' ...
               '示教关键点必须与控制器时间网格一致。'], dt);
    end

    % 生成单个周期（每个关节独立）
    [q1_cycle, dq1_cycle, ddq1_cycle, t_cycle] = generate_multi_segment(traj1, seg_durs, dt);
    [q2_cycle, dq2_cycle, ddq2_cycle, ~]       = generate_multi_segment(traj2, seg_durs, dt);
    [q3_cycle, dq3_cycle, ddq3_cycle, ~]       = generate_multi_segment(traj3, seg_durs, dt);

    % 多周期拼接：共享周期连接处的同一个时间样本
    % 单周期含 [0,T] 两个端点。第2周期起只追加其第2个样本到末端，
    % 避免把每个周期的 t=0 作为一个额外 dt 样本重复插入。
    % 因而总采样点数严格为 loop_num*(N_cycle-1)+1，
    % 总时长严格为 loop_num*T_cycle。
    N_cycle = length(t_cycle);
    N_step_cycle = N_cycle - 1;
    N_total = params.loop_num * N_step_cycle + 1;

    q_ref   = zeros(N_total, 3);
    dq_ref  = zeros(N_total, 3);
    ddq_ref = zeros(N_total, 3);

    cycle_q   = [q1_cycle, q2_cycle, q3_cycle];
    cycle_dq  = [dq1_cycle, dq2_cycle, dq3_cycle];
    cycle_ddq = [ddq1_cycle, ddq2_cycle, ddq3_cycle];

    % 第一周期完整保留，保证 TeachingDataPipeline 可直接切取 0~T 周期。
    q_ref(1:N_cycle,:)   = cycle_q;
    dq_ref(1:N_cycle,:)  = cycle_dq;
    ddq_ref(1:N_cycle,:) = cycle_ddq;

    for k = 2:params.loop_num
        idx_start = (k-1)*N_step_cycle + 2;
        idx_end   = k*N_step_cycle + 1;
        q_ref(idx_start:idx_end,:)   = cycle_q(2:end,:);
        dq_ref(idx_start:idx_end,:)  = cycle_dq(2:end,:);
        ddq_ref(idx_start:idx_end,:) = cycle_ddq(2:end,:);
    end

    t_all = (0:N_total-1)' * dt;

    % 轨迹连续性与合理性检查
    check_trajectory_continuity(q_ref, dq_ref, ddq_ref, params);
    check_acceleration_profile(ddq_ref, params);

    fprintf('轨迹生成完成：%d 个周期，总时长 %.2f s，采样点数 %d\n', ...
            params.loop_num, t_all(end), N_total);
    fprintf('角加速度峰值（deg/s²）：Hip %.1f, Knee %.1f, Ankle %.1f\n', ...
            max(abs(ddq_ref(:,1)))*180/pi, max(abs(ddq_ref(:,2)))*180/pi, ...
            max(abs(ddq_ref(:,3)))*180/pi);
end

% ====================== 辅助函数 ======================
function [q_all, dq_all, ddq_all, t_all] = generate_multi_segment(traj, seg_durs, dt)
    % 多段轨迹拼接（支持不等段时长；边界速度/加速度严格使用示教提取值，不做额外缩放）
    q_all = [];
    dq_all = [];
    ddq_all = [];
    t_all = [];

    n_seg = size(traj,1) - 1;
    if isscalar(seg_durs)
        seg_durs = seg_durs * ones(1, n_seg);
    end
    assert(length(seg_durs) == n_seg, '段时长数量必须等于轨迹段数');

    t_offset = 0;

    for s = 1:n_seg
        q0   = traj(s,1);
        dq0  = traj(s,2);
        ddq0 = traj(s,3);

        qf   = traj(s+1,1);
        dqf  = traj(s+1,2);
        ddqf = traj(s+1,3);

        T_seg = seg_durs(s);

        % 调用 quintic_poly
        [q, dq, ddq, t] = local_quintic_poly(q0, dq0, ddq0, qf, dqf, ddqf, T_seg, dt);

        % 避免首点重复
        if s > 1
            start_idx = 2;
        else
            start_idx = 1;
        end

        q_all   = [q_all;   q(start_idx:end)];
        dq_all  = [dq_all;  dq(start_idx:end)];
        ddq_all = [ddq_all; ddq(start_idx:end)];
        t_all   = [t_all;   t(start_idx:end) + t_offset];

        t_offset = t_all(end);
    end
end

function check_trajectory_continuity(q_ref, dq_ref, ddq_ref, params)
    % 连续性检查：数值导数一致性 + 周期连接邻域跳变量
    max_dq_jump = max(abs(diff(q_ref,1,1)/params.dt - dq_ref(1:end-1,:)), [], 'all');
    if max_dq_jump > 0.02
        warning('轨迹速度与位置差分存在较大偏差（%.4f rad/s），建议检查关键点设置', max_dq_jump);
    end

    if params.loop_num > 1
        % 每个周期理论长度由总步数反推，检查周期连接前后的一阶离散变化。
        n_step_cycle = (size(q_ref,1)-1) / params.loop_num;
        if abs(n_step_cycle-round(n_step_cycle)) < 1e-12
            n_step_cycle = round(n_step_cycle);
            join_idx = (1:params.loop_num-1)*n_step_cycle + 1;
            q_join_step = q_ref(join_idx+1,:) - q_ref(join_idx,:);
            dq_join_step = dq_ref(join_idx+1,:) - dq_ref(join_idx,:);
            ddq_join_step = ddq_ref(join_idx+1,:) - ddq_ref(join_idx,:);
            fprintf('周期连接最大单步变化：q %.4g rad, dq %.4g rad/s, ddq %.4g rad/s²\n', ...
                max(abs(q_join_step),[],'all'), max(abs(dq_join_step),[],'all'), ...
                max(abs(ddq_join_step),[],'all'));
        end
    end
end

function check_acceleration_profile(ddq_ref, params)
    % 加速度轮廓合理性检查
    max_ddq = max(abs(ddq_ref), [], 1) * 180/pi;
    fprintf('当前角加速度峰值（deg/s²）：Hip %.1f, Knee %.1f, Ankle %.1f\n', max_ddq);

    if any(max_ddq < 10)
        warning('角加速度峰值仍偏低（<10 °/s²）。建议：缩短 parameters.m 中的 T_seg（至1.5~1.8s）或直接增大 traj*_deg 中的 dq/ddq 分量。');
    else
        fprintf('✅ 角加速度轮廓正常，峰值合理（满足动态要求）。\n');
    end
end

% ===== Local helper: quintic segment generator =====
function [q, dq, ddq, t] = local_quintic_poly(q0, dq0, ddq0, qf, dqf, ddqf, T, dt)
    % 五次多项式轨迹规划（最小jerk经典实现）
    % 输入：
    %   q0, dq0, ddq0   : 起始位置、速度、加速度 (rad, rad/s, rad/s²)
    %   qf, dqf, ddqf   : 终止位置、速度、加速度
    %   T               : 段持续时间 (s)
    %   dt              : 采样时间 (s)
    % 输出：
    %   q, dq, ddq      : 位置、速度、加速度序列 (列向量)
    %   t               : 时间序列
    
    % 参数验证
    if T <= 0 || dt <= 0
        error('T 和 dt 必须为正值');
    end
    if nargin < 8
        dt = 0.01;  % 默认采样时间
    end
    
    % 严格按整数控制步数构造时间网格，避免 0:dt:T 的浮点端点丢失。
    n_steps = round(T / dt);
    if abs(n_steps*dt - T) > 1e-10
        error('trajectory_generator:LocalSegmentGridMismatch', ...
              '局部轨迹段 T=%.12g s 不能被 dt=%.12g s 整除。', T, dt);
    end
    t = (0:n_steps)' * dt;
    t(end) = T;
    N = length(t);
    
    % 五次多项式系数计算
    a0 = q0;
    a1 = dq0;
    a2 = ddq0 / 2;
    
    % 构造线性方程组求解 a3,a4,a5
    M = [T^3,     T^4,      T^5;
         3*T^2,   4*T^3,    5*T^4;
         6*T,    12*T^2,   20*T^3];
    
    b = [qf  - (a0 + a1*T + a2*T^2);
         dqf - (a1 + 2*a2*T);
         ddqf - (2*a2)];
    
    % 求解系数
    x = M \ b;
    a3 = x(1);
    a4 = x(2);
    a5 = x(3);
    
    % 计算轨迹
    q   = a0 + a1*t + a2*t.^2 + a3*t.^3 + a4*t.^4 + a5*t.^5;
    dq  = a1 + 2*a2*t + 3*a3*t.^2 + 4*a4*t.^3 + 5*a5*t.^4;
    ddq = 2*a2 + 6*a3*t + 12*a4*t.^2 + 20*a5*t.^3;
    
    % 可选：轻微数值平滑（提升论文图表质量）
    % q = smoothdata(q, 1, 'movmean', 3);
end
