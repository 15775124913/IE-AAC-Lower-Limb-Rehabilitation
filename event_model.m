function [event_all, timing] = event_model(t_all, params)
%EVENT_MODEL 事件扰动生成 + 事件时间窗口计算的统一入口
%
% 统一替代原来的：
%   get_event_timing.m
%   build_event_disturbance.m
%
% 调用：
%   [event_all, timing] = event_model(t_all, params);
%   [~, timing]         = event_model(t_all, params);
%
% 事件模型：
%   e_j(t) = A_j * exp(-((t - t0_j)/sigma_j)^2)

    if nargin < 1
        t_all = [];
    end
    if nargin < 2
        error('event_model:MissingInput', '必须提供 t_all 和 params');
    end

    %% ================= 1. 读取统一事件时间参数 =================
    % Event timing: cycle index + within-cycle phase
    if isfield(params,'event_cycle')
        event_cycle = params.event_cycle;
    else
        event_cycle = 1;
    end

    if isfield(params,'event_phase_base')
        event_phase_base = params.event_phase_base;
    else
        event_phase_base = [2.4, 3.8, 5.2];
    end

    if isfield(params,'traj_times') && ~isempty(params.traj_times)
        cycle_duration = sum(params.traj_times);
    else
        cycle_duration = (size(params.traj1_deg,1)-1)*params.T_seg;
    end

    event_time_base = (event_cycle-1)*cycle_duration + event_phase_base;

    if isfield(params, 'event_time_offset') && ~isempty(params.event_time_offset)
        event_time_offset = params.event_time_offset;
    else
        event_time_offset = [0, 0, 0];
    end

    if isfield(params, 'event_sigma') && ~isempty(params.event_sigma)
        event_sigma = params.event_sigma;
    else
        event_sigma = [0.25, 0.28, 0.30];
    end

    if isfield(params, 'event_window_nsigma') && ~isempty(params.event_window_nsigma)
        n_sigma = params.event_window_nsigma;
    else
        n_sigma = 3.0;
    end

    if isscalar(event_time_base),   event_time_base = repmat(event_time_base,1,3); end
    if isscalar(event_time_offset), event_time_offset = repmat(event_time_offset,1,3); end
    if isscalar(event_sigma),       event_sigma = repmat(event_sigma,1,3); end

    event_time_base   = event_time_base(:)';
    event_time_offset = event_time_offset(:)';
    event_sigma       = event_sigma(:)';

    assert(numel(event_time_base)==3,   'event_time_base 必须为 1x3');
    assert(numel(event_time_offset)==3, 'event_time_offset 必须为 1x3');
    assert(numel(event_sigma)==3,       'event_sigma 必须为 1x3');
    assert(all(event_sigma > 0),        'event_sigma 必须全部为正');
    assert(isscalar(n_sigma) && n_sigma > 0, ...
           'event_window_nsigma 必须为正标量');

    centers = event_time_base + event_time_offset;

    individual_windows = [ ...
        centers(:) - n_sigma*event_sigma(:), ...
        centers(:) + n_sigma*event_sigma(:)];

    abnormal_window = [ ...
        min(individual_windows(:,1)), ...
        max(individual_windows(:,2))];

    timing = struct();
    timing.base_centers = event_time_base;
    timing.offsets = event_time_offset;
    timing.centers = centers;
    timing.sigma = event_sigma;
    timing.n_sigma = n_sigma;
    timing.individual_windows = individual_windows;
    timing.abnormal_window_unclipped = abnormal_window;
    timing.envelope_threshold_ratio = exp(-(n_sigma^2));

    %% ================= 2. 时间掩码 =================
    if ~isempty(t_all)
        t_all = t_all(:);
        sim_start = t_all(1);
        sim_end = t_all(end);

        clipped = [ ...
            max(abnormal_window(1), sim_start), ...
            min(abnormal_window(2), sim_end)];

        if clipped(2) < clipped(1)
            clipped = [sim_start, sim_start];
        end

        timing.abnormal_window = clipped;

        event_mask_joint = false(length(t_all),3);
        for j = 1:3
            event_mask_joint(:,j) = ...
                (t_all >= individual_windows(j,1)) & ...
                (t_all <= individual_windows(j,2));
        end

        timing.event_mask_joint = event_mask_joint;
        timing.event_mask = any(event_mask_joint,2);

        pre_end = min(max(clipped(1), sim_start), sim_end);
        timing.pre_event_window = [sim_start, pre_end];
    else
        timing.abnormal_window = abnormal_window;
        timing.event_mask_joint = [];
        timing.event_mask = [];
        timing.pre_event_window = [];
    end

    %% ================= 3. 生成高斯事件扰动 =================
    if isfield(params, 'event_amp') && ~isempty(params.event_amp)
        event_amp = params.event_amp;
    else
        event_amp = [0.50, 0.60, 0.70];
    end

    if isscalar(event_amp)
        event_amp = repmat(event_amp,1,3);
    end
    event_amp = event_amp(:)';
    assert(numel(event_amp)==3, 'event_amp 必须为 1x3');

    if isempty(t_all)
        event_all = [];
        return;
    end

    event_all = zeros(length(t_all),3);
    for j = 1:3
        event_all(:,j) = event_amp(j) .* ...
            exp(-((t_all - centers(j)) ./ event_sigma(j)).^2);
    end
end
