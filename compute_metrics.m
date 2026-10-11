function metrics = compute_metrics(q_ref, q_actual, dq_actual, ddq_actual, ...
                                   e_fb_all, e_fb_used, params, interaction_diag)
% 计算康复机器人控制性能指标（统一指标定义，支持自适应 & 固定导纳）
% 输入：
%   q_ref / q_actual : 参考 / 实际关节角度 (N x 3)
%   dq_actual        : 实际关节速度 (N x 3)
%   ddq_actual       : 实际关节加速度 (N x 3)
%   e_fb_all         : 原始反馈误差 (N x 3)，用于反映真实交互力矩
%   e_fb_used        : 经死区/静止区处理后的反馈误差 (N x 3)，用于反映控制器实际输入
%   params           : 参数结构体（含 dt 等）

    if nargin < 8, interaction_diag = []; end
    dt = params.dt;
    N  = size(q_ref, 1);
    joints = {'Hip', 'Knee', 'Ankle'};

    % 兼容：若 e_fb_all 为空则回退到 e_fb_used
    if isempty(e_fb_all)
        e_fb_all = e_fb_used;
    end

    metrics = struct();

    %% ================= 1. 轨迹跟踪误差指标 =================
    q_err = q_actual - q_ref;                       % 角度误差 (rad)

    metrics.RMSE = sqrt(mean(q_err.^2, 1)) * 180/pi;             % (deg)
    metrics.MAE  = mean(abs(q_err), 1) * 180/pi;                 % (deg)
    metrics.IAE  = sum(abs(q_err), 1) * dt;                      % (rad·s)
    % 最大跟踪误差保留为追溯/诊断量。由于 q_actual-q_ref=Delta q 且
    % Delta q 存在共享硬限幅，该统计量在 Monte Carlo 中可能被右删失，
    % 因此不作为主性能指标。
    metrics.MaxError = max(abs(q_err), [], 1) * 180/pi;          % (deg), diagnostic only

    % 共享位置修正限幅诊断：比 Peak Error 更能区分“偶发触限”与
    % “长时间由限幅塑形”。该量不参与控制，仅用于 supervisory audit。
    q_lim = params.delta_q_limit;
    if isscalar(q_lim), q_lim = repmat(q_lim,1,3); end
    assert(numel(q_lim)==3 && all(q_lim>0), 'delta_q_limit 必须为正标量或 1x3');
    sat_ratio_th = 0.999;
    if isfield(params,'motion_limit_saturation_ratio') && ~isempty(params.motion_limit_saturation_ratio)
        sat_ratio_th = params.motion_limit_saturation_ratio;
    end
    q_corr_ratio = abs(q_err) ./ q_lim(:)';
    metrics.PositionCorrectionMaxRatio = max(q_corr_ratio, [], 1);
    metrics.PositionCorrectionSaturationPct = 100*mean(q_corr_ratio >= sat_ratio_th, 1);

    %% ================= 2. 速度与平滑度指标 =================
    metrics.RMS_velocity = sqrt(mean(dq_actual.^2, 1));          % (rad/s)

    % ★ Jerk 必须沿时间维（第1维）逐关节计算。
    % 对 N×3 矩阵直接调用 gradient(ddq_actual, dt) 会把矩阵维度解释为
    % 空间维度，不能保证得到每个关节沿时间轴的导数。这里显式逐列求导。
    assert(size(ddq_actual,1) == N && size(ddq_actual,2) == 3, ...
           'ddq_actual 必须为与 q_ref 同长度的 N×3 矩阵');
    jerk = zeros(size(ddq_actual));
    for jj = 1:size(ddq_actual,2)
        jerk(:,jj) = gradient(ddq_actual(:,jj), dt);             % (rad/s^3)
    end
    metrics.RMS_jerk = sqrt(mean(jerk.^2, 1));                   % (rad/s³)
    metrics.RMS_jerk_deg = metrics.RMS_jerk * 180/pi;             % (deg/s³), paper reporting

    %% ================= 3. 交互力矩指标 =================
    % ★ 关键修改：主要对比指标基于原始误差 e_fb_all，真实反映人机交互力矩
    if ~isempty(e_fb_all)
        metrics.RMS_InteractionTorque = sqrt(mean(e_fb_all.^2, 1));      % (Nm)
        metrics.Mean_InteractionTorque = mean(abs(e_fb_all), 1);         % (Nm)
        metrics.PeakInteractionTorque = max(abs(e_fb_all), [], 1);       % (Nm)
    else
        metrics.RMS_InteractionTorque = zeros(1,3);
        metrics.Mean_InteractionTorque = zeros(1,3);
        metrics.PeakInteractionTorque = zeros(1,3);
    end

    % 控制器实际使用的补偿/死区后信号
    if ~isempty(e_fb_used)
        metrics.RMS_ControllerTorque = sqrt(mean(e_fb_used.^2, 1));      % (Nm)
        metrics.ControllerInputEnergy = sum(abs(e_fb_used).*abs(dq_actual),1)*dt; % (J)
    else
        metrics.RMS_ControllerTorque = zeros(1,3);
        metrics.ControllerInputEnergy = zeros(1,3);
    end

    % 病理交互分量：来自控制器六分量分解中的 pathological。
    % 这是模型内定义的病理扰动残差，不等同于传感器实测患者力矩。
    if ~isempty(interaction_diag) && isfield(interaction_diag,'pathological')
        tau_path = interaction_diag.pathological;
        metrics.RMS_PathologicalInteraction = sqrt(mean(tau_path.^2,1));
        metrics.Peak_PathologicalInteraction = max(abs(tau_path),[],1);
        metrics.PathologicalInteractionEnergy = sum(abs(tau_path).*abs(dq_actual),1)*dt;
    else
        metrics.RMS_PathologicalInteraction = nan(1,3);
        metrics.Peak_PathologicalInteraction = nan(1,3);
        metrics.PathologicalInteractionEnergy = nan(1,3);
    end

    %% ================= 4. 交互能量 =================
    % Raw residual energy 仅作为诊断量；它受 mismatch-pathology cancellation 影响。
    metrics.InteractionEnergy = sum(abs(e_fb_all) .* abs(dq_actual), 1) * dt;
    metrics.RawInteractionResidualEnergy = metrics.InteractionEnergy;
    metrics.RMS_RawInteractionResidual = metrics.RMS_InteractionTorque;

    %% ================= 5. 追溯/诊断指标 =================
    % 与 MaxError 等价，保留旧字段以兼容既有 MAT/脚本；不作为 MC 主指标。
    metrics.PeakTrackingError = metrics.MaxError;                  % (deg), diagnostic only

    % 建立时间 SettlingTime (s)
    steady_start = round(0.8 * N) + 1;
    if steady_start > N, steady_start = 1; end
    steady_err = mean(abs(q_err(steady_start:end, :)), 1);
    threshold  = 1.1 * steady_err;

    settling_time = zeros(1,3);
    for j = 1:3
        found = false;
        for i = 1:N-1
            if mean(abs(q_err(i:end, j))) <= threshold(j)
                settling_time(j) = (i-1) * dt;
                found = true;
                break;
            end
        end
        if ~found
            settling_time(j) = NaN;
        end
    end
    metrics.SettlingTime = settling_time;

    % 超调量 Overshoot (%)
    max_err_abs = max(abs(q_err), [], 1);
    overshoot = zeros(1,3);
    for j = 1:3
        if steady_err(j) > 1e-8
            overshoot(j) = (max_err_abs(j) - steady_err(j)) / steady_err(j) * 100;
        else
            overshoot(j) = 0;
        end
    end
    metrics.Overshoot = overshoot;

    %% ================= 事件扰动时段分段指标 =================
    % ★ 第4项修复：禁止使用固定 [4,8] s。
    % 根据本次 params 的真实 event_time_base + event_time_offset +
    % event_sigma 自动生成事件有效区间；Monte Carlo 每个 trial 会自动移动。
    t_all = (0:N-1)' * dt;
    [~, event_timing] = event_model(t_all, params);
    idx_ab = event_timing.event_mask;

    metrics.EventCenters_s = event_timing.centers;
    metrics.EventSigma_s = event_timing.sigma;
    metrics.EventWindow_s = event_timing.abnormal_window;
    metrics.EventWindowNSigma = event_timing.n_sigma;
    metrics.EventSampleCount = sum(idx_ab);

    for j = 1:3
        % 对每个关节使用该关节自身的事件窗口，而不是笼统固定窗口
        idx_event_j = event_timing.event_mask_joint(:,j);

        if sum(idx_event_j) > 1
            metrics.RMS_IT_abnormal(j) = ...
                sqrt(mean(e_fb_all(idx_event_j, j).^2));
            metrics.Jerk_RMS_abnormal(j) = ...
                sqrt(mean(jerk(idx_event_j, j).^2));
            metrics.Peak_IT(j) = ...
                max(abs(e_fb_all(idx_event_j, j)));
        else
            metrics.RMS_IT_abnormal(j) = NaN;
            metrics.Jerk_RMS_abnormal(j) = NaN;
            metrics.Peak_IT(j) = NaN;
        end
    end

    %% ================= 6. 汇总统计 =================
    metrics.joints = joints;
    metrics.total_RMSE = mean(metrics.RMSE);
    metrics.total_RMS_InteractionTorque = mean(metrics.RMS_InteractionTorque);
    metrics.total_Mean_InteractionTorque = mean(metrics.Mean_InteractionTorque);
    metrics.total_PositionCorrectionSaturationPct = mean(metrics.PositionCorrectionSaturationPct);

    fprintf('   RMSE (deg): %.3f  %.3f  %.3f | 平均: %.3f\n', ...
            metrics.RMSE, metrics.total_RMSE);
    fprintf('   RMS Raw Interaction Residual (Nm): %.4f  %.4f  %.4f | 平均: %.4f\n', ...
            metrics.RMS_InteractionTorque, metrics.total_RMS_InteractionTorque);
    fprintf('   Raw Residual Energy (J): %.3f  %.3f  %.3f\n', ...
            metrics.InteractionEnergy);
    fprintf('   Position-correction saturation (%%): %.2f  %.2f  %.2f\n', ...
            metrics.PositionCorrectionSaturationPct);
end