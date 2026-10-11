function T = nominal_diagnostics(mode, varargin)
%NOMINAL_DIAGNOSTICS Central nominal-only audits for the four-controller project.
% These diagnostics are structural/supervisory checks, not performance-tuning
% objectives. They must not be used to retune either literature comparator
% against the Proposed controller or Monte Carlo outcomes.
    switch lower(char(mode))
        case {'zhang2018','zhang'}
            T = assess_literature_va_nominal_core('Zhang-2018', varargin{:});
        case {'cao2024','cao'}
            T = assess_literature_va_nominal_core('Cao-2024-derived', varargin{:});
        case 'motion_limits'
            T = assess_four_controller_nominal_core(varargin{:});
        otherwise
            error('Unknown nominal diagnostic: %s', char(mode));
    end
end

function T = assess_literature_va_nominal_core(comparator_name, Ba_lit, Ka_lit, Ga_lit, ...
    delta_q_lit, delta_dq_lit, delta_ddq_lit, quiet_zone_flag_lit, params, varargin)
%ASSESS_LITERATURE_VA_NOMINAL Structural audit for one literature comparator.

    interaction_diag = [];
    if ~isempty(varargin)
        interaction_diag = varargin{1};
    end

    joints = {'Hip';'Knee';'Ankle'};
    Comparator = repmat({comparator_name},3,1);
    Ba0 = params.Ba0(:)'; Ka0 = params.Ka0(:)'; Ga0 = params.Ga0(:)';

    q_lim = params.delta_q_limit; dq_lim = params.delta_dq_limit; ddq_lim = params.delta_ddq_limit;
    if isscalar(q_lim), q_lim = repmat(q_lim,1,3); end
    if isscalar(dq_lim), dq_lim = repmat(dq_lim,1,3); end
    if isscalar(ddq_lim), ddq_lim = repmat(ddq_lim,1,3); end
    sat_ratio_th = 0.999;
    if isfield(params,'motion_limit_saturation_ratio') && ~isempty(params.motion_limit_saturation_ratio)
        sat_ratio_th = params.motion_limit_saturation_ratio;
    end

    B_std_pct = zeros(3,1); K_std_pct = zeros(3,1);
    B_min_ratio = zeros(3,1); B_max_ratio = zeros(3,1);
    K_min_ratio = zeros(3,1); K_max_ratio = zeros(3,1);
    B_bound_pct = zeros(3,1); K_bound_pct = zeros(3,1);
    q_sat_pct = zeros(3,1); dq_sat_pct = zeros(3,1); ddq_sat_pct = zeros(3,1);
    quiet_pct = zeros(3,1); Ga_max_dev_pct = zeros(3,1);

    for j = 1:3
        nonquiet = quiet_zone_flag_lit(:,j) < 0.5;
        if ~any(nonquiet), nonquiet = true(size(nonquiet)); end
        B = Ba_lit(nonquiet,j); K = Ka_lit(nonquiet,j);

        B_std_pct(j) = 100*std(B)/max(Ba0(j),eps);
        K_std_pct(j) = 100*std(K)/max(Ka0(j),eps);
        B_min_ratio(j) = min(B)/Ba0(j); B_max_ratio(j) = max(B)/Ba0(j);
        K_min_ratio(j) = min(K)/Ka0(j); K_max_ratio(j) = max(K)/Ka0(j);

        B_lo = 0.1*Ba0(j); B_hi = 4.0*Ba0(j);
        K_lo = 0.01*Ka0(j); K_hi = 5.0*Ka0(j);
        btol = 1e-8*max(1,Ba0(j)); ktol = 1e-8*max(1,Ka0(j));
        B_bound_pct(j) = 100*mean(abs(B-B_lo)<=btol | abs(B-B_hi)<=btol);
        K_bound_pct(j) = 100*mean(abs(K-K_lo)<=ktol | abs(K-K_hi)<=ktol);

        q_sat_pct(j) = 100*mean(abs(delta_q_lit(:,j)) >= sat_ratio_th*q_lim(j));
        dq_sat_pct(j) = 100*mean(abs(delta_dq_lit(:,j)) >= sat_ratio_th*dq_lim(j));
        ddq_sat_pct(j) = 100*mean(abs(delta_ddq_lit(:,j)) >= sat_ratio_th*ddq_lim(j));
        quiet_pct(j) = 100*mean(quiet_zone_flag_lit(:,j) > 0.5);
        Ga_max_dev_pct(j) = 100*max(abs(Ga_lit(:,j)-Ga0(j)))/max(Ga0(j),eps);
    end

    % Cao-only cycle-domain diagnostics. Zhang rows remain NaN by design.
    B_cycle_delta_pct = nan(3,1); K_cycle_delta_pct = nan(3,1);
    VB_rms_pct = nan(3,1); VK_rms_pct = nan(3,1);
    GradAbsMedian = nan(3,1); CycleCount = nan(3,1);
    if ~isempty(interaction_diag) && isfield(interaction_diag,'cao_cycle_index') && ...
            isfield(interaction_diag,'cao_base_B') && isfield(interaction_diag,'cao_base_K')
        cyc = interaction_diag.cao_cycle_index(:,1);
        valid_cyc = cyc > 0;
        if any(valid_cyc)
            first_c = min(cyc(valid_cyc)); last_c = max(cyc(valid_cyc));
            for j = 1:3
                m1 = interaction_diag.cao_cycle_index(:,j) == first_c;
                mL = interaction_diag.cao_cycle_index(:,j) == last_c;
                B1 = interaction_diag.cao_base_B(m1,j); BL = interaction_diag.cao_base_B(mL,j);
                K1 = interaction_diag.cao_base_K(m1,j); KL = interaction_diag.cao_base_K(mL,j);
                nmatch = min([numel(B1),numel(BL),numel(K1),numel(KL)]);
                if nmatch > 0
                    B_cycle_delta_pct(j) = 100*sqrt(mean((BL(1:nmatch)-B1(1:nmatch)).^2))/max(Ba0(j),eps);
                    K_cycle_delta_pct(j) = 100*sqrt(mean((KL(1:nmatch)-K1(1:nmatch)).^2))/max(Ka0(j),eps);
                end
                if isfield(interaction_diag,'cao_VB')
                    VB_rms_pct(j) = 100*sqrt(mean(interaction_diag.cao_VB(valid_cyc,j).^2))/max(Ba0(j),eps);
                end
                if isfield(interaction_diag,'cao_VK')
                    VK_rms_pct(j) = 100*sqrt(mean(interaction_diag.cao_VK(valid_cyc,j).^2))/max(Ka0(j),eps);
                end
                if isfield(interaction_diag,'cao_gradient_surrogate')
                    GradAbsMedian(j) = median(abs(interaction_diag.cao_gradient_surrogate(valid_cyc,j)));
                end
                CycleCount(j) = last_c;
            end
        end
    end

    T = table(Comparator, joints, B_std_pct, K_std_pct, ...
        B_min_ratio, B_max_ratio, K_min_ratio, K_max_ratio, ...
        B_bound_pct, K_bound_pct, q_sat_pct, dq_sat_pct, ddq_sat_pct, ...
        quiet_pct, Ga_max_dev_pct, ...
        B_cycle_delta_pct, K_cycle_delta_pct, VB_rms_pct, VK_rms_pct, GradAbsMedian, CycleCount, ...
        'VariableNames', {'Comparator','Joint','B_NonQuietStd_pct_B0','K_Std_pct_K0', ...
        'B_MinRatio_B0','B_MaxRatio_B0','K_MinRatio_K0','K_MaxRatio_K0', ...
        'B_NonQuietBoundHit_pct','K_BoundHit_pct', ...
        'DeltaQ_Saturation_pct','DeltaDQ_Saturation_pct','DeltaDDQ_Saturation_pct', ...
        'QuietZone_pct','Ga_MaxDeviation_pct_G0', ...
        'B_Cycle1ToLast_RMSDelta_pct_B0','K_Cycle1ToLast_RMSDelta_pct_K0', ...
        'VB_RMS_pct_B0','VK_RMS_pct_K0','HRIGradientAbsMedian_Nm','CycleCount'});
end

function T = assess_four_controller_nominal_core( ...
    delta_q_fixed, delta_dq_fixed, delta_ddq_fixed, ...
    delta_q_zhang, delta_dq_zhang, delta_ddq_zhang, ...
    delta_q_cao, delta_dq_cao, delta_ddq_cao, ...
    delta_q_adp, delta_dq_adp, delta_ddq_adp, params)
%ASSESS_FOUR_CONTROLLER_NOMINAL Common motion-limit audit for all formal controllers.

    joints = {'Hip';'Knee';'Ankle'};
    controllers = {'Fixed','Zhang2018VA','Cao2024DerivedVA','Proposed'};

    q_lim = params.delta_q_limit;
    dq_lim = params.delta_dq_limit;
    ddq_lim = params.delta_ddq_limit;
    if isscalar(q_lim), q_lim = repmat(q_lim,1,3); end
    if isscalar(dq_lim), dq_lim = repmat(dq_lim,1,3); end
    if isscalar(ddq_lim), ddq_lim = repmat(ddq_lim,1,3); end

    sat_ratio_th = 0.999;
    if isfield(params,'motion_limit_saturation_ratio') && ~isempty(params.motion_limit_saturation_ratio)
        sat_ratio_th = params.motion_limit_saturation_ratio;
    end

    qsets = {delta_q_fixed, delta_q_zhang, delta_q_cao, delta_q_adp};
    dqsets = {delta_dq_fixed, delta_dq_zhang, delta_dq_cao, delta_dq_adp};
    ddqsets = {delta_ddq_fixed, delta_ddq_zhang, delta_ddq_cao, delta_ddq_adp};

    nrows = numel(controllers)*3;
    Controller = cell(nrows,1); Joint = cell(nrows,1);
    DeltaQ_MaxRatio = zeros(nrows,1); DeltaQ_Saturation_pct = zeros(nrows,1);
    DeltaDQ_MaxRatio = zeros(nrows,1); DeltaDQ_Saturation_pct = zeros(nrows,1);
    DeltaDDQ_MaxRatio = zeros(nrows,1); DeltaDDQ_Saturation_pct = zeros(nrows,1);

    r = 0;
    for c = 1:numel(controllers)
        for j = 1:3
            r = r + 1;
            Controller{r} = controllers{c};
            Joint{r} = joints{j};
            qr = abs(qsets{c}(:,j)) / q_lim(j);
            dqr = abs(dqsets{c}(:,j)) / dq_lim(j);
            ddqr = abs(ddqsets{c}(:,j)) / ddq_lim(j);
            DeltaQ_MaxRatio(r) = max(qr);
            DeltaQ_Saturation_pct(r) = 100*mean(qr >= sat_ratio_th);
            DeltaDQ_MaxRatio(r) = max(dqr);
            DeltaDQ_Saturation_pct(r) = 100*mean(dqr >= sat_ratio_th);
            DeltaDDQ_MaxRatio(r) = max(ddqr);
            DeltaDDQ_Saturation_pct(r) = 100*mean(ddqr >= sat_ratio_th);
        end
    end

    T = table(Controller, Joint, DeltaQ_MaxRatio, DeltaQ_Saturation_pct, ...
        DeltaDQ_MaxRatio, DeltaDQ_Saturation_pct, ...
        DeltaDDQ_MaxRatio, DeltaDDQ_Saturation_pct);
end
