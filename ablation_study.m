% ABLATION_STUDY.M - consolidated from V1/V2/V3 experimental scripts.
% Standalone test wrappers and duplicate figure/runner files were merged here.

function results = ablation_study(mode, varargin)
%ABLATION_STUDY Unified core-ablation, mismatch-sensitivity and evidence module.
%
% Modes
%   'core'      - Full / No-Mismatch / Simplified B-K ablation.
%   'mismatch'  - K_mis sensitivity and mechanism diagnostics.
%   'finalize'  - Aggregate completed 100-run outputs and generate final tables/figures.
%   'all'       - Run core + mismatch + finalization sequentially.
%
% Examples
%   R = ablation_study('core','NumRuns',100);
%   R = ablation_study('core','NumRuns',5,'MakeFigures',false);  % smoke test
%   R = ablation_study('mismatch','NumRuns',100);
%   R = ablation_study('finalize');
%   R = ablation_study('all','NumRuns',100);
%
% The full Proposed controller remains unchanged by default. Ablation variants
% are applied through params.ablation and K_mismatch_comp only.

    if nargin < 1 || isempty(mode), mode = 'core'; end
    mode = lower(char(mode));
    switch mode
        case 'core'
            results = ab_core_runner(varargin{:});
        case 'mismatch'
            results = ab_mismatch_runner(varargin{:});
        case 'finalize'
            results = ab_finalize_evidence(varargin{:});
        case 'all'
            ip = inputParser;
            addParameter(ip,'NumRuns',100,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
            addParameter(ip,'MakeFigures',true,@(x)islogical(x)||isnumeric(x));
            addParameter(ip,'Seed',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=0));
            addParameter(ip,'OutputRoot','results',@(x)ischar(x)||isstring(x));
            addParameter(ip,'KmisGrid',[0 0.25 0.5 0.75 1],@(x)isnumeric(x)&&isvector(x));
            parse(ip,varargin{:});
            o=ip.Results;
            core = ab_core_runner('NumRuns',o.NumRuns,'MakeFigures',o.MakeFigures, ...
                'Seed',o.Seed,'OutputRoot',o.OutputRoot);
            mismatch = ab_mismatch_runner('NumRuns',o.NumRuns,'MakeFigures',o.MakeFigures, ...
                'Seed',o.Seed,'OutputRoot',o.OutputRoot,'KmisGrid',o.KmisGrid);
            final = ab_finalize_evidence('AblationDir',core.output_dir, ...
                'MismatchDir',mismatch.output_dir,'ExpectedRuns',o.NumRuns, ...
                'RequireExactRuns',true,'MakeFigures',o.MakeFigures,'OutputRoot',o.OutputRoot);
            results = struct('core',core,'mismatch',mismatch,'final',final);
        otherwise
            error('ablation_study:UnknownMode','Unknown mode: %s',mode);
    end
end


function results = ab_core_runner(varargin)
%CORE_ABLATION_RUNNER Core-module ablation experiment for Proposed IE-AAC.
%
% Primary variants:
%   1) Full Proposed
%   2) No Mismatch Compensation: K_mismatch_comp = 0
%   3) Simplified Adaptation (B/K): adapt B_a and K_a, hold G_a at G_a0
%
% Optional:
%   4) B-only: adapt B_a, hold K_a and G_a at their nominal values
%
% The experiment keeps the reference, sampled patient/controller parameters,
% pregenerated event/noise realizations, dead-zone, quiet-zone supervision,
% parameter projection, and motion/jerk limits common within each paired trial.
% State-dependent interaction components are recomputed from each variant's
% own commanded trajectory, as required by the closed-loop model.
%
% Examples:
%   ab_core_runner();
%   ab_core_runner('NumRuns', 10, 'MakeFigures', false); % smoke test
%   ab_core_runner('IncludeBOnly', true);
%
% Outputs are written to:
%   results/ABLATION_yyyymmdd_HHMMSS/
%
% This file is intentionally independent of monte_carlo_runner.m so the
% existing main paper experiment remains untouched.

    %% ---------------- Options ----------------
    ip = inputParser;
    addParameter(ip, 'NumRuns', 100, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(ip, 'RunMonteCarlo', true, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'IncludeBOnly', false, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'MakeFigures', true, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'SaveNominalTimeseries', true, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'Seed', [], @(x)isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));
    addParameter(ip, 'OutputRoot', 'results', @(x)ischar(x) || isstring(x));
    parse(ip, varargin{:});
    opt = ip.Results;

    num_runs = round(opt.NumRuns);
    run_mc = logical(opt.RunMonteCarlo);
    include_b_only = logical(opt.IncludeBOnly);
    make_figures = logical(opt.MakeFigures);
    save_ts = logical(opt.SaveNominalTimeseries);

    %% ---------------- Base configuration and output ----------------
    base_params = parameters();
    if ~isempty(opt.Seed)
        base_params.seed = double(opt.Seed);
    end

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    out_root = char(opt.OutputRoot);
    if ~exist(out_root, 'dir'), mkdir(out_root); end
    out_dir = fullfile(out_root, ['ABLATION_' timestamp]);
    mkdir(out_dir);

    fprintf('\n============================================================\n');
    fprintf(' Core ablation experiment: Proposed IE-AAC\n');
    fprintf(' Output: %s\n', out_dir);
    fprintf(' Seed: %d | Monte Carlo runs: %d\n', base_params.seed, num_runs);
    fprintf('============================================================\n\n');

    %% ---------------- One authoritative reference ----------------
    [base_params, metrics_teaching, q_ref, dq_ref, ddq_ref, t_all] = ...
        TeachingDataPipeline(base_params, ...
            'Seed', base_params.seed, ...
            'PlotFlag', false, ...
            'SaveFlag', false, ...
            'Verbose', false);

    assert(isequal(size(q_ref), size(dq_ref), size(ddq_ref)), ...
        'Reference position/velocity/acceleration sizes are inconsistent.');

    specs = ab_core_build_variant_specs(base_params, include_b_only);
    nV = numel(specs);

    variant_manifest = ab_core_build_variant_manifest(specs);
    writetable(variant_manifest, fullfile(out_dir, 'ablation_variant_manifest.csv'));

    %% ---------------- Nominal paired ablation ----------------
    pre_fb_nominal = mc_trial_utils('pre_fb', base_params, q_ref, t_all, base_params.seed);

    nominal = struct();
    for v = 1:nV
        pv = ab_core_apply_variant(base_params, specs(v));
        fprintf('\n--- Nominal: %s ---\n', specs(v).label);
        nominal.(specs(v).key) = ab_core_run_one_variant( ...
            q_ref, dq_ref, ddq_ref, t_all, pv, pre_fb_nominal);
    end

    T_nominal = ab_core_build_nominal_table(nominal, specs);
    writetable(T_nominal, fullfile(out_dir, 'ablation_nominal_metrics.csv'));

    if save_ts
        T_ts = ab_core_build_nominal_timeseries(t_all, q_ref, nominal, specs);
        writetable(T_ts, fullfile(out_dir, 'ablation_nominal_timeseries.csv'));
    end

    %% ---------------- Paired Monte Carlo ----------------
    stats = struct();
    mc_summary = table();
    paired_comparison = table();
    if run_mc
        stats = ab_core_init_stats(specs, num_runs);
        stats.success_mask = false(num_runs,1);
        stats.parameter_seed = nan(num_runs,1);
        stats.noise_seed = nan(num_runs,1);
        stats.failure_message = strings(num_runs,1);
        stats.trial_audit = cell(num_runs,1);

        fprintf('\n============================================================\n');
        fprintf(' Paired Monte Carlo ablation: %d trials x %d variants\n', num_runs, nV);
        fprintf('============================================================\n');

        for r = 1:num_runs
            fprintf('Ablation MC %d/%d ... ', r, num_runs);
            try
                [trial_params, audit] = mc_trial_utils('sample', base_params, r);
                pre_fb = mc_trial_utils('pre_fb', trial_params, q_ref, t_all, audit.noise_seed);

                tmp = struct();
                all_ok = true;
                fail_msg = "";

                for v = 1:nV
                    pv = ab_core_apply_variant(trial_params, specs(v));
                    rv = ab_core_run_one_variant(q_ref, dq_ref, ddq_ref, t_all, pv, pre_fb);
                    if ~ab_core_variant_is_finite(rv)
                        all_ok = false;
                        fail_msg = "Non-finite output in variant " + string(specs(v).label);
                        break;
                    end
                    tmp.(specs(v).key) = ab_core_extract_trial_metrics(rv);
                end

                if all_ok
                    for v = 1:nV
                        key = specs(v).key;
                        stats.(key) = ab_core_commit_trial(stats.(key), tmp.(key), r);
                    end
                    stats.success_mask(r) = true;
                    stats.parameter_seed(r) = audit.parameter_seed;
                    stats.noise_seed(r) = audit.noise_seed;
                    stats.trial_audit{r} = audit;
                    fprintf('OK\n');
                else
                    stats.failure_message(r) = fail_msg;
                    fprintf('FAILED: %s\n', fail_msg);
                end
            catch ME
                stats.failure_message(r) = string(ME.message);
                fprintf('FAILED: %s\n', ME.message);
            end
        end

        success_rate = mean(stats.success_mask);
        fprintf('\nPaired success rate: %.1f%% (%d/%d)\n', ...
            100*success_rate, sum(stats.success_mask), num_runs);

        if isfield(base_params, 'mc') && isfield(base_params.mc, 'min_success_rate') ...
                && success_rate < base_params.mc.min_success_rate
            warning('ab_core_runner:LowSuccessRate', ...
                'Paired success rate %.1f%% is below configured threshold %.1f%%.', ...
                100*success_rate, 100*base_params.mc.min_success_rate);
        end

        T_mc_long = ab_core_build_mc_long_table(stats, specs);
        writetable(T_mc_long, fullfile(out_dir, 'ablation_mc_trial_metrics_long.csv'));

        mc_summary = ab_core_build_mc_summary(stats, specs);
        writetable(mc_summary, fullfile(out_dir, 'ablation_mc_summary_mean_std.csv'));

        paired_comparison = ab_core_build_paired_comparison(stats, specs);
        writetable(paired_comparison, fullfile(out_dir, 'ablation_mc_paired_full_vs_ablations.csv'));
    end

    %% ---------------- Save all traceable data ----------------
    save(fullfile(out_dir, 'ablation_full_results.mat'), ...
        'base_params', 'metrics_teaching', 'q_ref', 'dq_ref', 'ddq_ref', 't_all', ...
        'specs', 'nominal', 'T_nominal', 'stats', 'mc_summary', 'paired_comparison', ...
        'timestamp', '-v7.3');

    if make_figures
        try
            ab_core_generate_figures(out_dir, t_all, q_ref, nominal, stats, specs);
        catch ME
            warning('ab_core_runner:FigureGenerationFailed', ...
                'Ablation data were saved, but figure generation failed: %s', ME.message);
        end
    end

    %% ---------------- Return compact result ----------------
    results = struct();
    results.output_dir = out_dir;
    results.timestamp = timestamp;
    results.specs = specs;
    results.nominal_table = T_nominal;
    results.mc_summary = mc_summary;
    results.paired_comparison = paired_comparison;
    if run_mc
        results.success_rate = mean(stats.success_mask);
    else
        results.success_rate = NaN;
    end

    fprintf('\nCore ablation experiment complete.\n');
    fprintf('Results saved to: %s\n', out_dir);
end


%% ========================================================================

function specs = ab_core_build_variant_specs(base_params, include_b_only)
    specs = struct([]);

    specs(1).key = 'full';
    specs(1).label = 'Full Proposed';
    specs(1).K_mis = base_params.K_mismatch_comp;
    specs(1).adapt_B = true;
    specs(1).adapt_K = true;
    specs(1).adapt_G = true;

    specs(2).key = 'no_mismatch';
    specs(2).label = 'No Mismatch Compensation';
    specs(2).K_mis = [0 0 0];
    specs(2).adapt_B = true;
    specs(2).adapt_K = true;
    specs(2).adapt_G = true;

    specs(3).key = 'simplified_bk';
    specs(3).label = 'Simplified Adaptation (B/K)';
    specs(3).K_mis = base_params.K_mismatch_comp;
    specs(3).adapt_B = true;
    specs(3).adapt_K = true;
    specs(3).adapt_G = false;

    if include_b_only
        specs(4).key = 'simplified_b';
        specs(4).label = 'Simplified Adaptation (B only)';
        specs(4).K_mis = base_params.K_mismatch_comp;
        specs(4).adapt_B = true;
        specs(4).adapt_K = false;
        specs(4).adapt_G = false;
    end
end



function T = ab_core_build_variant_manifest(specs)
    nV = numel(specs);
    Variant = strings(nV,1);
    Key = strings(nV,1);
    Kmis_Hip = zeros(nV,1);
    Kmis_Knee = zeros(nV,1);
    Kmis_Ankle = zeros(nV,1);
    Adapt_B = false(nV,1);
    Adapt_K = false(nV,1);
    Adapt_G = false(nV,1);

    for v = 1:nV
        km = specs(v).K_mis;
        if isscalar(km), km = repmat(km,1,3); end
        Variant(v) = string(specs(v).label);
        Key(v) = string(specs(v).key);
        Kmis_Hip(v) = km(1);
        Kmis_Knee(v) = km(2);
        Kmis_Ankle(v) = km(3);
        Adapt_B(v) = specs(v).adapt_B;
        Adapt_K(v) = specs(v).adapt_K;
        Adapt_G(v) = specs(v).adapt_G;
    end

    T = table(Variant, Key, Kmis_Hip, Kmis_Knee, Kmis_Ankle, ...
        Adapt_B, Adapt_K, Adapt_G);
end



function p = ab_core_apply_variant(p, spec)
    p.K_mismatch_comp = spec.K_mis;
    p.ablation = struct();
    p.ablation.name = spec.label;
    p.ablation.adapt_B = spec.adapt_B;
    p.ablation.adapt_K = spec.adapt_K;
    p.ablation.adapt_G = spec.adapt_G;
end



function R = ab_core_run_one_variant(q_ref, dq_ref, ddq_ref, t_all, params, pre_fb)
    [delta_q, delta_dq, delta_ddq, ...
     Ba, Ka, Ga, quiet_zone_flag, e_used, e_raw, diag] = ...
        run_admittance_controller('proposed', ...
            q_ref, dq_ref, params, pre_fb, t_all, ddq_ref);

    R = struct();
    R.delta_q = delta_q;
    R.delta_dq = delta_dq;
    R.delta_ddq = delta_ddq;
    R.q = q_ref + delta_q;
    R.dq = dq_ref + delta_dq;
    R.ddq = ddq_ref + delta_ddq;
    R.Ba = Ba;
    R.Ka = Ka;
    R.Ga = Ga;
    R.quiet_zone_flag = quiet_zone_flag;
    R.e_used = e_used;
    R.e_raw = e_raw;
    R.diag = diag;
    R.metrics = compute_metrics(q_ref, R.q, R.dq, R.ddq, ...
        e_raw, e_used, params, diag);

    % Mechanism diagnostics; these are not primary performance metrics.
    R.mechanism = struct();
    R.mechanism.MismatchRMS = sqrt(mean(diag.mismatch.^2, 1));
    R.mechanism.CompensatedResidualRMS = sqrt(mean(diag.compensated.^2, 1));
    R.mechanism.MismatchLeakageRMS = ...
        sqrt(mean((diag.compensated - diag.pathological).^2, 1));
end



function ok = ab_core_variant_is_finite(R)
    ok = all(isfinite(R.q(:))) && all(isfinite(R.dq(:))) && ...
         all(isfinite(R.ddq(:))) && all(isfinite(R.Ba(:))) && ...
         all(isfinite(R.Ka(:))) && all(isfinite(R.Ga(:))) && ...
         all(isfinite(R.e_raw(:))) && all(isfinite(R.e_used(:))) && ...
         all(isfinite(R.metrics.RMSE(:))) && ...
         all(isfinite(R.metrics.RMS_jerk_deg(:))) && ...
         all(isfinite(R.metrics.RMS_InteractionTorque(:))) && ...
         all(isfinite(R.mechanism.MismatchLeakageRMS(:)));
end



function T = ab_core_build_nominal_table(nominal, specs)
    joints = ["Hip"; "Knee"; "Ankle"];
    nV = numel(specs);
    nRows = nV*3;

    Variant = strings(nRows,1);
    Joint = strings(nRows,1);
    K_mis = zeros(nRows,1);
    Adapt_B = false(nRows,1);
    Adapt_K = false(nRows,1);
    Adapt_G = false(nRows,1);
    RMSE_deg = zeros(nRows,1);
    JerkRMS_deg_s3 = zeros(nRows,1);
    RawResidualRMS_Nm = zeros(nRows,1);
    ControllerInputRMS_Nm = zeros(nRows,1);
    PositionSaturation_pct = zeros(nRows,1);
    PathologicalRMS_Nm = zeros(nRows,1);
    MismatchRMS_Nm = zeros(nRows,1);
    CompensatedResidualRMS_Nm = zeros(nRows,1);
    MismatchLeakageRMS_Nm = zeros(nRows,1);

    row = 0;
    for v = 1:nV
        R = nominal.(specs(v).key);
        km = specs(v).K_mis;
        if isscalar(km), km = repmat(km,1,3); end
        for j = 1:3
            row = row + 1;
            Variant(row) = string(specs(v).label);
            Joint(row) = joints(j);
            K_mis(row) = km(j);
            Adapt_B(row) = specs(v).adapt_B;
            Adapt_K(row) = specs(v).adapt_K;
            Adapt_G(row) = specs(v).adapt_G;
            RMSE_deg(row) = R.metrics.RMSE(j);
            JerkRMS_deg_s3(row) = R.metrics.RMS_jerk_deg(j);
            RawResidualRMS_Nm(row) = R.metrics.RMS_InteractionTorque(j);
            ControllerInputRMS_Nm(row) = R.metrics.RMS_ControllerTorque(j);
            PositionSaturation_pct(row) = R.metrics.PositionCorrectionSaturationPct(j);
            PathologicalRMS_Nm(row) = R.metrics.RMS_PathologicalInteraction(j);
            MismatchRMS_Nm(row) = R.mechanism.MismatchRMS(j);
            CompensatedResidualRMS_Nm(row) = R.mechanism.CompensatedResidualRMS(j);
            MismatchLeakageRMS_Nm(row) = R.mechanism.MismatchLeakageRMS(j);
        end
    end

    T = table(Variant, Joint, K_mis, Adapt_B, Adapt_K, Adapt_G, ...
        RMSE_deg, JerkRMS_deg_s3, RawResidualRMS_Nm, ...
        ControllerInputRMS_Nm, PositionSaturation_pct, PathologicalRMS_Nm, ...
        MismatchRMS_Nm, CompensatedResidualRMS_Nm, MismatchLeakageRMS_Nm);
end



function T = ab_core_build_nominal_timeseries(t_all, q_ref, nominal, specs)
    joints = {'Hip','Knee','Ankle'};
    T = table(t_all(:), 'VariableNames', {'Time_s'});
    for j = 1:3
        T.(sprintf('Reference_%s_deg', joints{j})) = q_ref(:,j)*180/pi;
    end

    for v = 1:numel(specs)
        key = specs(v).key;
        R = nominal.(key);
        for j = 1:3
            T.(sprintf('%s_qcmd_%s_deg', key, joints{j})) = R.q(:,j)*180/pi;
            T.(sprintf('%s_error_%s_deg', key, joints{j})) = R.delta_q(:,j)*180/pi;
            T.(sprintf('%s_Ba_%s', key, joints{j})) = R.Ba(:,j);
            T.(sprintf('%s_Ka_%s', key, joints{j})) = R.Ka(:,j);
            T.(sprintf('%s_Ga_%s', key, joints{j})) = R.Ga(:,j);
            T.(sprintf('%s_raw_%s_Nm', key, joints{j})) = R.e_raw(:,j);
            T.(sprintf('%s_used_%s_Nm', key, joints{j})) = R.e_used(:,j);
            T.(sprintf('%s_mismatch_%s_Nm', key, joints{j})) = R.diag.mismatch(:,j);
        end
    end
end



function stats = ab_core_init_stats(specs, num_runs)
    stats = struct();
    for v = 1:numel(specs)
        key = specs(v).key;
        S = struct();
        S.RMSE = nan(num_runs,3);
        S.JerkRMS = nan(num_runs,3);
        S.RawResidualRMS = nan(num_runs,3);
        S.ControllerInputRMS = nan(num_runs,3);
        S.PositionSaturationPct = nan(num_runs,3);
        S.PathologicalRMS = nan(num_runs,3);
        S.MismatchRMS = nan(num_runs,3);
        S.CompensatedResidualRMS = nan(num_runs,3);
        S.MismatchLeakageRMS = nan(num_runs,3);
        stats.(key) = S;
    end
end



function M = ab_core_extract_trial_metrics(R)
    M = struct();
    M.RMSE = R.metrics.RMSE;
    M.JerkRMS = R.metrics.RMS_jerk_deg;
    M.RawResidualRMS = R.metrics.RMS_InteractionTorque;
    M.ControllerInputRMS = R.metrics.RMS_ControllerTorque;
    M.PositionSaturationPct = R.metrics.PositionCorrectionSaturationPct;
    M.PathologicalRMS = R.metrics.RMS_PathologicalInteraction;
    M.MismatchRMS = R.mechanism.MismatchRMS;
    M.CompensatedResidualRMS = R.mechanism.CompensatedResidualRMS;
    M.MismatchLeakageRMS = R.mechanism.MismatchLeakageRMS;
end



function S = ab_core_commit_trial(S, M, r)
    fields = fieldnames(M);
    for k = 1:numel(fields)
        S.(fields{k})(r,:) = M.(fields{k});
    end
end



function T = ab_core_build_mc_long_table(stats, specs)
    joints = ["Hip","Knee","Ankle"];
    num_runs = numel(stats.success_mask);
    nV = numel(specs);
    nRows = num_runs*nV*3;

    Trial = zeros(nRows,1);
    Success = false(nRows,1);
    Variant = strings(nRows,1);
    Joint = strings(nRows,1);
    RMSE_deg = nan(nRows,1);
    JerkRMS_deg_s3 = nan(nRows,1);
    RawResidualRMS_Nm = nan(nRows,1);
    ControllerInputRMS_Nm = nan(nRows,1);
    PositionSaturation_pct = nan(nRows,1);
    PathologicalRMS_Nm = nan(nRows,1);
    MismatchRMS_Nm = nan(nRows,1);
    CompensatedResidualRMS_Nm = nan(nRows,1);
    MismatchLeakageRMS_Nm = nan(nRows,1);

    row = 0;
    for r = 1:num_runs
        for v = 1:nV
            key = specs(v).key;
            for j = 1:3
                row = row + 1;
                Trial(row) = r;
                Success(row) = stats.success_mask(r);
                Variant(row) = string(specs(v).label);
                Joint(row) = joints(j);
                RMSE_deg(row) = stats.(key).RMSE(r,j);
                JerkRMS_deg_s3(row) = stats.(key).JerkRMS(r,j);
                RawResidualRMS_Nm(row) = stats.(key).RawResidualRMS(r,j);
                ControllerInputRMS_Nm(row) = stats.(key).ControllerInputRMS(r,j);
                PositionSaturation_pct(row) = stats.(key).PositionSaturationPct(r,j);
                PathologicalRMS_Nm(row) = stats.(key).PathologicalRMS(r,j);
                MismatchRMS_Nm(row) = stats.(key).MismatchRMS(r,j);
                CompensatedResidualRMS_Nm(row) = stats.(key).CompensatedResidualRMS(r,j);
                MismatchLeakageRMS_Nm(row) = stats.(key).MismatchLeakageRMS(r,j);
            end
        end
    end

    T = table(Trial, Success, Variant, Joint, RMSE_deg, JerkRMS_deg_s3, ...
        RawResidualRMS_Nm, ControllerInputRMS_Nm, PositionSaturation_pct, ...
        PathologicalRMS_Nm, MismatchRMS_Nm, CompensatedResidualRMS_Nm, ...
        MismatchLeakageRMS_Nm);
end



function T = ab_core_build_mc_summary(stats, specs)
    joints = ["Hip","Knee","Ankle"];
    metric_names = { ...
        'RMSE','JerkRMS','RawResidualRMS','ControllerInputRMS', ...
        'PositionSaturationPct','PathologicalRMS','MismatchRMS', ...
        'CompensatedResidualRMS','MismatchLeakageRMS'};

    nRows = numel(specs)*3*numel(metric_names);
    Variant = strings(nRows,1);
    Joint = strings(nRows,1);
    Metric = strings(nRows,1);
    N = zeros(nRows,1);
    Mean = nan(nRows,1);
    Std = nan(nRows,1);
    Median = nan(nRows,1);

    row = 0;
    for v = 1:numel(specs)
        key = specs(v).key;
        for m = 1:numel(metric_names)
            field = metric_names{m};
            X = stats.(key).(field);
            for j = 1:3
                row = row + 1;
                x = X(:,j);
                x = x(isfinite(x));
                Variant(row) = string(specs(v).label);
                Joint(row) = joints(j);
                Metric(row) = string(field);
                N(row) = numel(x);
                if ~isempty(x)
                    Mean(row) = mean(x);
                    Std(row) = std(x);
                    Median(row) = median(x);
                end
            end
        end
    end

    T = table(Variant, Joint, Metric, N, Mean, Std, Median);
end



function T = ab_core_build_paired_comparison(stats, specs)
    % Primary paper metrics only. Positive ComparatorMinusFull means
    % the Full Proposed controller has a lower (better) metric.
    joints = ["Hip","Knee","Ankle"];
    metric_names = {'RMSE','JerkRMS','RawResidualRMS'};

    if numel(specs) < 2
        T = table();
        return;
    end

    nRows = (numel(specs)-1)*numel(metric_names)*3;
    Comparator = strings(nRows,1);
    Joint = strings(nRows,1);
    Metric = strings(nRows,1);
    N = zeros(nRows,1);
    FullMean = nan(nRows,1);
    ComparatorMean = nan(nRows,1);
    ComparatorMinusFull_Mean = nan(nRows,1);
    Difference95CI_Low = nan(nRows,1);
    Difference95CI_High = nan(nRows,1);
    FullImprovementPct_Mean = nan(nRows,1);

    row = 0;
    full_key = specs(1).key;

    for v = 2:numel(specs)
        comp_key = specs(v).key;
        for m = 1:numel(metric_names)
            field = metric_names{m};
            XF = stats.(full_key).(field);
            XC = stats.(comp_key).(field);

            for j = 1:3
                row = row + 1;
                f = XF(:,j);
                c = XC(:,j);
                idx = isfinite(f) & isfinite(c);
                f = f(idx);
                c = c(idx);

                Comparator(row) = string(specs(v).label);
                Joint(row) = joints(j);
                Metric(row) = string(field);
                N(row) = numel(f);

                if ~isempty(f)
                    d = c - f;
                    FullMean(row) = mean(f);
                    ComparatorMean(row) = mean(c);
                    ComparatorMinusFull_Mean(row) = mean(d);
                    [lo, hi] = ab_core_normal95_ci(d);
                    Difference95CI_Low(row) = lo;
                    Difference95CI_High(row) = hi;

                    safe = abs(c) > eps;
                    if any(safe)
                        FullImprovementPct_Mean(row) = mean(100*(c(safe)-f(safe))./c(safe));
                    end
                end
            end
        end
    end

    T = table(Comparator, Joint, Metric, N, FullMean, ComparatorMean, ...
        ComparatorMinusFull_Mean, Difference95CI_Low, Difference95CI_High, ...
        FullImprovementPct_Mean);
end



function [lo, hi] = ab_core_normal95_ci(x)
    x = x(isfinite(x));
    if isempty(x)
        lo = NaN; hi = NaN;
    elseif numel(x) == 1
        lo = x; hi = x;
    else
        mu = mean(x);
        se = std(x)/sqrt(numel(x));
        lo = mu - 1.96*se;
        hi = mu + 1.96*se;
    end
end


%% ========================================================================
% Monte Carlo sampler copied from the current project logic so this runner
% does not modify or depend on monte_carlo_runner.m internals.


function ab_core_generate_figures(out_dir, t_all, q_ref, nominal, stats, specs)
%GENERATE_ABLATION_FIGURES Paper-oriented figures for the core ablation study.
% Presentation is standardized through paper_figure_utils.m. Every exported
% PDF/PNG also receives an identically named CSV and MAT in figure_data/.

    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    fig_dir = fullfile(out_dir, 'figures');
    dat_dir = fullfile(out_dir, 'figure_data');
    if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
    if ~exist(dat_dir, 'dir'), mkdir(dat_dir); end

    labels = {specs.label};
    nV = numel(specs);
    joint_names = {'Hip','Knee','Ankle'};

    %% Figure A1: Nominal primary metrics
    metric_fields = {'RMSE','RMS_jerk_deg','RMS_InteractionTorque'};
    metric_titles = {'Tracking RMSE','Jerk RMS','Raw interaction-residual RMS'};
    ylabels = {'RMSE (deg)','Jerk RMS (deg/s^3)','Raw residual RMS (N m)'};

    base_name = 'FigA1_nominal_core_ablation_primary_metrics';
    [f1,lay1] = paper_figure_utils('create','Core Ablation - Nominal Metrics', ...
        [100 60 1080 1160],3,1,'Core Ablation: Nominal Primary Metrics');
    ax1 = gobjects(3,1);
    nr=3*3*nV; Metric=strings(nr,1); Joint=strings(nr,1); Variant=strings(nr,1); Value=nan(nr,1); row=0;
    for m = 1:3
        ax = paper_figure_utils('tile',lay1,m); ax1(m)=ax;
        M = zeros(3,nV);
        for v = 1:nV
            R = nominal.(specs(v).key);
            M(:,v) = R.metrics.(metric_fields{m})(:);
            for j=1:3
                row=row+1; Metric(row)=string(metric_titles{m}); Joint(row)=string(joint_names{j});
                Variant(row)=string(labels{v}); Value(row)=M(j,v);
            end
        end
        bar(ax,M,'grouped');
        xticks(ax,1:3); xticklabels(ax,joint_names);
        ylabel(ax,ylabels{m});
        if m==3, xlabel(ax,'Joint'); end
        paper_figure_utils('panel_title',ax,sprintf('(%c) %s',char('a'+m-1),metric_titles{m}));
    end
    lgd = gobjects(numel(ax1),1);
    for ii = 1:numel(ax1)
        lgd(ii) = legend(ax1(ii),labels,'Location','northeast','Interpreter','none','AutoUpdate','off');
    end
    paper_figure_utils('align_legends',ax1,lgd,'NumColumns',min(2,nV));
    Tfig=table(Metric,Joint,Variant,Value);
    meta=struct('figure',base_name,'description','Nominal core-ablation primary metrics', ...
        'metric_fields',{metric_fields},'units',{ylabels},'variants',{labels},'joints',{joint_names});
    ab_save_figure_data(dat_dir,base_name,Tfig,meta);
    ab_core_export_figure_pair(f1, fig_dir, base_name);

    %% Figure A2: Nominal mechanism diagnostics
    base_name = 'FigA2_nominal_ablation_mechanism_diagnostics';
    [f2,lay2] = paper_figure_utils('create','Core Ablation - Mechanism Diagnostics', ...
        [130 80 1080 850],2,1,'Core Ablation: Mechanism Diagnostics');
    Mleak = zeros(3,nV); Minput = zeros(3,nV);
    for v = 1:nV
        R = nominal.(specs(v).key);
        Mleak(:,v) = R.mechanism.MismatchLeakageRMS(:);
        Minput(:,v) = R.metrics.RMS_ControllerTorque(:);
    end
    ax21=paper_figure_utils('tile',lay2,1);
    bar(ax21,Mleak,'grouped'); xticks(ax21,1:3); xticklabels(ax21,joint_names);
    ylabel(ax21,'RMS (N m)');
    paper_figure_utils('panel_title',ax21,'(a) Mismatch Leakage into Compensated Residual');
    ax22=paper_figure_utils('tile',lay2,2);
    bar(ax22,Minput,'grouped'); xticks(ax22,1:3); xticklabels(ax22,joint_names);
    ylabel(ax22,'RMS (N m)'); xlabel(ax22,'Joint');
    paper_figure_utils('panel_title',ax22,'(b) Processed Controller-Input RMS');
    lgd = gobjects(2,1);
    lgd(1)=legend(ax21,labels,'Location','northeast','Interpreter','none','AutoUpdate','off');
    lgd(2)=legend(ax22,labels,'Location','northeast','Interpreter','none','AutoUpdate','off');
    paper_figure_utils('align_legends',[ax21 ax22],lgd,'NumColumns',min(2,nV));
    nr=2*3*nV; Diagnostic=strings(nr,1); Joint=strings(nr,1); Variant=strings(nr,1); Value=nan(nr,1); row=0;
    diag_names={'Mismatch Leakage RMS','Processed Controller-Input RMS'}; mats={Mleak,Minput};
    for d=1:2
        for v=1:nV
            for j=1:3
                row=row+1; Diagnostic(row)=string(diag_names{d}); Joint(row)=string(joint_names{j});
                Variant(row)=string(labels{v}); Value(row)=mats{d}(j,v);
            end
        end
    end
    Tfig=table(Diagnostic,Joint,Variant,Value);
    meta=struct('figure',base_name,'description','Nominal core-ablation mechanism diagnostics', ...
        'unit','N m RMS','variants',{labels},'joints',{joint_names});
    ab_save_figure_data(dat_dir,base_name,Tfig,meta);
    ab_core_export_figure_pair(f2, fig_dir, base_name);

    %% Figure A3: Nominal adaptive-parameter trajectories
    base_name = 'FigA3_nominal_ablation_parameter_evolution';
    [f3,lay3] = paper_figure_utils('create','Core Ablation - Parameters', ...
        [160 45 1150 1180],3,1,'Core Ablation: Effective Admittance Parameters (Nominal)');
    param_names = {'Ba','Ka','Ga'};
    y_param = {'B_a','K_a','G_a'};
    ax3 = gobjects(3,1);
    nt=numel(t_all); nr=nt*3*nV; Time_s=nan(nr,1); Parameter=strings(nr,1); Variant=strings(nr,1); Value=nan(nr,1); row=0;
    for pp = 1:3
        ax = paper_figure_utils('tile',lay3,pp); ax3(pp)=ax; hold(ax,'on');
        for v = 1:nV
            R = nominal.(specs(v).key);
            X = R.(param_names{pp}); y=mean(X,2);
            plot(ax,t_all,y,'LineWidth',1.25,'DisplayName',labels{v});
            idx=row+(1:nt); Time_s(idx)=t_all(:); Parameter(idx)=string(y_param{pp});
            Variant(idx)=string(labels{v}); Value(idx)=y(:); row=row+nt;
        end
        ylabel(ax,y_param{pp}); if pp==3, xlabel(ax,'Time (s)'); end
        paper_figure_utils('panel_title',ax,sprintf('(%c) %s (Three-Joint Mean)',char('a'+pp-1),y_param{pp}));
    end
    lgd = gobjects(numel(ax3),1);
    for ii = 1:numel(ax3)
        lgd(ii)=legend(ax3(ii),'Location','northeast','Interpreter','none','AutoUpdate','off');
    end
    paper_figure_utils('align_legends',ax3,lgd,'NumColumns',min(2,nV));
    Tfig=table(Time_s,Parameter,Variant,Value);
    meta=struct('figure',base_name,'description','Three-joint mean nominal effective admittance-parameter trajectories', ...
        'parameters',{y_param},'variants',{labels});
    ab_save_figure_data(dat_dir,base_name,Tfig,meta);
    ab_core_export_figure_pair(f3, fig_dir, base_name);

    %% Figure A4: Monte Carlo mean +/- SD for primary metrics
    if ~isempty(stats) && isstruct(stats) && isfield(stats,'success_mask') && any(stats.success_mask)
        base_name = 'FigA4_mc_core_ablation_primary_metrics';
        mc_fields = {'RMSE','JerkRMS','RawResidualRMS'};
        mc_titles = {'Tracking RMSE','Jerk RMS','Raw interaction-residual RMS'};
        mc_ylabels = {'RMSE (deg)','Jerk RMS (deg/s^3)','Raw residual RMS (N m)'};
        [f4,lay4] = paper_figure_utils('create','Core Ablation - Monte Carlo', ...
            [190 45 1120 1180],3,1, ...
            sprintf('Paired Monte Carlo Core Ablation (%d Successful Trials)',sum(stats.success_mask)));
        ax4 = gobjects(3,1);
        nr=3*3*nV; Metric=strings(nr,1); Joint=strings(nr,1); Variant=strings(nr,1); Mean=nan(nr,1); SD=nan(nr,1); N=zeros(nr,1); row=0;
        for m = 1:3
            ax = paper_figure_utils('tile',lay4,m); ax4(m)=ax;
            Means=zeros(3,nV); Stds=zeros(3,nV);
            for v=1:nV
                X=stats.(specs(v).key).(mc_fields{m});
                Means(:,v)=mean(X,1,'omitnan')';
                Stds(:,v)=std(X,0,1,'omitnan')';
                for j=1:3
                    x=X(:,j); row=row+1; Metric(row)=string(mc_titles{m}); Joint(row)=string(joint_names{j});
                    Variant(row)=string(labels{v}); Mean(row)=Means(j,v); SD(row)=Stds(j,v); N(row)=sum(isfinite(x));
                end
            end
            b=bar(ax,Means,'grouped'); hold(ax,'on');
            for v=1:numel(b)
                try
                    x=b(v).XEndPoints;
                catch
                    groupwidth=min(0.8,nV/(nV+1.5));
                    x=(1:3)-groupwidth/2+(2*v-1)*groupwidth/(2*nV);
                end
                errorbar(ax,x,Means(:,v),Stds(:,v),'k.','LineWidth',1.0,'CapSize',6,'HandleVisibility','off');
            end
            xticks(ax,1:3); xticklabels(ax,joint_names);
            ylabel(ax,mc_ylabels{m}); if m==3, xlabel(ax,'Joint'); end
            paper_figure_utils('panel_title',ax,sprintf('(%c) %s (Mean +/- SD)',char('a'+m-1),mc_titles{m}));
        end
        lgd = gobjects(numel(ax4),1);
        for ii = 1:numel(ax4)
            lgd(ii)=legend(ax4(ii),labels,'Location','northeast','Interpreter','none','AutoUpdate','off');
        end
        paper_figure_utils('align_legends',ax4,lgd,'NumColumns',min(2,nV));
        Tfig=table(Metric,Joint,Variant,Mean,SD,N);
        meta=struct('figure',base_name,'description','Paired Monte Carlo core-ablation mean and standard deviation', ...
            'successful_trials',sum(stats.success_mask),'metric_fields',{mc_fields},'units',{mc_ylabels});
        ab_save_figure_data(dat_dir,base_name,Tfig,meta);
        ab_core_export_figure_pair(f4, fig_dir, base_name);
    end

    %% Figure A5: Nominal tracking-error traces during full time sequence
    base_name = 'FigA5_nominal_ablation_tracking_error';
    [f5,lay5] = paper_figure_utils('create','Core Ablation - Tracking Error', ...
        [220 45 1150 1180],3,1,'Core Ablation: Nominal Tracking Error');
    ax5=gobjects(3,1);
    nt=numel(t_all); nr=nt*3*nV; Time_s=nan(nr,1); Joint=strings(nr,1); Variant=strings(nr,1); TrackingError_deg=nan(nr,1); row=0;
    for j=1:3
        ax=paper_figure_utils('tile',lay5,j); ax5(j)=ax; hold(ax,'on');
        for v=1:nV
            R=nominal.(specs(v).key); y=(R.q(:,j)-q_ref(:,j))*180/pi;
            plot(ax,t_all,y,'LineWidth',1.15,'DisplayName',labels{v});
            idx=row+(1:nt); Time_s(idx)=t_all(:); Joint(idx)=string(joint_names{j});
            Variant(idx)=string(labels{v}); TrackingError_deg(idx)=y(:); row=row+nt;
        end
        ylabel(ax,'Error (deg)'); if j==3, xlabel(ax,'Time (s)'); end
        paper_figure_utils('panel_title',ax,sprintf('(%c) %s',char('a'+j-1),joint_names{j}));
    end
    lgd = gobjects(numel(ax5),1);
    for ii = 1:numel(ax5)
        lgd(ii)=legend(ax5(ii),'Location','northeast','Interpreter','none','AutoUpdate','off');
    end
    paper_figure_utils('align_legends',ax5,lgd,'NumColumns',min(2,nV));
    Tfig=table(Time_s,Joint,Variant,TrackingError_deg);
    meta=struct('figure',base_name,'description','Nominal full-sequence tracking-error traces','unit','deg');
    ab_save_figure_data(dat_dir,base_name,Tfig,meta);
    ab_core_export_figure_pair(f5, fig_dir, base_name);

    % Keep figures visible after export, consistent with generate_paper_figures.m.
end

function ab_core_export_figure_pair(fig, fig_dir, base_name)
    paper_figure_utils('export',fig,fig_dir,base_name, ...
        'SavePNG',true,'SaveJPEG',true,'JPEGResolution',1000);
end


function results = ab_mismatch_runner(varargin)
%MISMATCH_MECHANISM_RUNNER Mechanism diagnostics and K_mis sensitivity study.
%
% Purpose
%   Explain the closed-loop effect of nominal-mismatch compensation without
%   changing the Proposed B_a/K_a/G_a adaptation law. The only swept factor
%   is K_mismatch_comp. All sensitivity variants retain full B/K/G adaptation.
%
% Default sensitivity levels
%   K_mis = 0, 0.25, 0.50, 0.75, 1.00, applied equally to all three joints.
%
% Strict paired design
%   Within every Monte Carlo trial, all K_mis variants share the same:
%     - authoritative reconstructed reference;
%     - sampled controller/patient-model parameter realization;
%     - event timing and event realization;
%     - stochastic-noise realization;
%     - fatigue trend;
%     - dead-zone, quiet-zone, projection, and motion/jerk constraints.
%   State-dependent interaction terms are recomputed from each variant's own
%   commanded trajectory, as required by the closed-loop model.
%
% Main outputs
%   mismatch_nominal_sensitivity.csv
%   mismatch_nominal_timeseries.csv
%   mismatch_mc_trial_metrics_long.csv
%   mismatch_mc_summary_long.csv
%   mismatch_mc_sensitivity_compact.csv
%   mismatch_mc_endpoint_full_vs_no_mismatch.csv
%   mismatch_sensitivity_results.mat
%
% Examples
%   ab_mismatch_runner('NumRuns',5,'MakeFigures',false); % smoke test
%   ab_mismatch_runner();                               % 100 runs
%
% Notes
%   This runner is intentionally independent of the main paper MC runner.
%   It does not tune K_mis. The grid is a sensitivity/diagnostic experiment.

    %% ---------------- Options ----------------
    ip = inputParser;
    addParameter(ip, 'NumRuns', 100, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(ip, 'KmisGrid', [0 0.25 0.50 0.75 1.00], ...
        @(x)isnumeric(x) && isvector(x) && ~isempty(x) && all(isfinite(x)) && all(x>=0) && all(x<=1));
    addParameter(ip, 'RunMonteCarlo', true, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'MakeFigures', true, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'SaveNominalTimeseries', true, @(x)islogical(x) || isnumeric(x));
    addParameter(ip, 'Seed', [], @(x)isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));
    addParameter(ip, 'OutputRoot', 'results', @(x)ischar(x) || isstring(x));
    parse(ip, varargin{:});
    opt = ip.Results;

    num_runs = round(opt.NumRuns);
    run_mc = logical(opt.RunMonteCarlo);
    make_figures = logical(opt.MakeFigures);
    save_ts = logical(opt.SaveNominalTimeseries);
    kmis_grid = unique(double(opt.KmisGrid(:)'), 'stable');

    % Ensure the two endpoints needed by the core ablation are present.
    if ~any(abs(kmis_grid-0) < 1e-12) || ~any(abs(kmis_grid-1) < 1e-12)
        error('ab_mismatch_runner:MissingEndpoints', ...
            'KmisGrid must contain both 0 and 1 for the endpoint ablation comparison.');
    end

    %% ---------------- Base configuration and output ----------------
    base_params = parameters();
    if ~isempty(opt.Seed)
        base_params.seed = double(opt.Seed);
    end

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    out_root = char(opt.OutputRoot);
    if ~exist(out_root, 'dir'), mkdir(out_root); end
    out_dir = fullfile(out_root, ['MISMATCH_DIAG_' timestamp]);
    mkdir(out_dir);

    fprintf('\n============================================================\n');
    fprintf(' Nominal-mismatch mechanism diagnostics and sensitivity\n');
    fprintf(' Output: %s\n', out_dir);
    fprintf(' Seed: %d | MC runs: %d\n', base_params.seed, num_runs);
    fprintf(' K_mis grid: %s\n', mat2str(kmis_grid,3));
    fprintf('============================================================\n\n');

    %% ---------------- One authoritative reference ----------------
    [base_params, metrics_teaching, q_ref, dq_ref, ddq_ref, t_all] = ...
        TeachingDataPipeline(base_params, ...
            'Seed', base_params.seed, ...
            'PlotFlag', false, ...
            'SaveFlag', false, ...
            'Verbose', false);

    assert(isequal(size(q_ref), size(dq_ref), size(ddq_ref)), ...
        'Reference position/velocity/acceleration sizes are inconsistent.');

    specs = ab_mis_build_kmis_specs(kmis_grid);
    nV = numel(specs);
    T_manifest = ab_mis_build_manifest(specs);
    writetable(T_manifest, fullfile(out_dir, 'mismatch_sensitivity_variant_manifest.csv'));

    %% ---------------- Nominal sensitivity ----------------
    pre_fb_nominal = mc_trial_utils('pre_fb', base_params, q_ref, t_all, base_params.seed);
    nominal = struct();

    for v = 1:nV
        pv = ab_mis_apply_kmis_variant(base_params, specs(v));
        fprintf('\n--- Nominal K_mis = %.2f ---\n', specs(v).K_mis);
        nominal.(specs(v).key) = ab_mis_run_one_variant( ...
            q_ref, dq_ref, ddq_ref, t_all, pv, pre_fb_nominal);
    end

    T_nominal = ab_mis_build_nominal_table(nominal, specs);
    writetable(T_nominal, fullfile(out_dir, 'mismatch_nominal_sensitivity.csv'));

    if save_ts
        T_ts = ab_mis_build_nominal_timeseries(t_all, q_ref, nominal, specs);
        writetable(T_ts, fullfile(out_dir, 'mismatch_nominal_timeseries.csv'));
    end

    %% ---------------- Paired Monte Carlo sensitivity ----------------
    stats = struct();
    T_mc_long = table();
    T_mc_summary = table();
    T_mc_compact = table();
    T_endpoint = table();

    if run_mc
        stats = ab_mis_init_stats(specs, num_runs);
        stats.success_mask = false(num_runs,1);
        stats.parameter_seed = nan(num_runs,1);
        stats.noise_seed = nan(num_runs,1);
        stats.failure_message = strings(num_runs,1);
        stats.trial_audit = cell(num_runs,1);

        fprintf('\n============================================================\n');
        fprintf(' Paired K_mis sensitivity MC: %d trials x %d levels\n', num_runs, nV);
        fprintf('============================================================\n');

        for r = 1:num_runs
            fprintf('Mismatch MC %d/%d ... ', r, num_runs);
            try
                [trial_params, audit] = mc_trial_utils('sample', base_params, r);
                pre_fb = mc_trial_utils('pre_fb', trial_params, q_ref, t_all, audit.noise_seed);

                tmp = struct();
                all_ok = true;
                fail_msg = "";

                for v = 1:nV
                    pv = ab_mis_apply_kmis_variant(trial_params, specs(v));
                    Rv = ab_mis_run_one_variant(q_ref, dq_ref, ddq_ref, t_all, pv, pre_fb);
                    if ~ab_mis_variant_is_finite(Rv)
                        all_ok = false;
                        fail_msg = "Non-finite output at K_mis=" + string(specs(v).K_mis);
                        break;
                    end
                    tmp.(specs(v).key) = ab_mis_extract_trial_metrics(Rv);
                end

                if all_ok
                    for v = 1:nV
                        key = specs(v).key;
                        stats.(key) = ab_mis_commit_trial(stats.(key), tmp.(key), r);
                    end
                    stats.success_mask(r) = true;
                    stats.parameter_seed(r) = audit.parameter_seed;
                    stats.noise_seed(r) = audit.noise_seed;
                    stats.trial_audit{r} = audit;
                    fprintf('OK\n');
                else
                    stats.failure_message(r) = fail_msg;
                    fprintf('FAILED: %s\n', fail_msg);
                end
            catch ME
                stats.failure_message(r) = string(ME.message);
                fprintf('FAILED: %s\n', ME.message);
            end
        end

        success_rate = mean(stats.success_mask);
        fprintf('\nPaired success rate: %.1f%% (%d/%d)\n', ...
            100*success_rate, sum(stats.success_mask), num_runs);

        if isfield(base_params,'mc') && isfield(base_params.mc,'min_success_rate') ...
                && success_rate < base_params.mc.min_success_rate
            warning('ab_mismatch_runner:LowSuccessRate', ...
                'Paired success rate %.1f%% is below configured threshold %.1f%%.', ...
                100*success_rate, 100*base_params.mc.min_success_rate);
        end

        T_mc_long = ab_mis_build_mc_long_table(stats, specs);
        writetable(T_mc_long, fullfile(out_dir, 'mismatch_mc_trial_metrics_long.csv'));

        T_mc_summary = ab_mis_build_mc_summary_long(stats, specs);
        writetable(T_mc_summary, fullfile(out_dir, 'mismatch_mc_summary_long.csv'));

        T_mc_compact = ab_mis_build_mc_compact_summary(stats, specs);
        writetable(T_mc_compact, fullfile(out_dir, 'mismatch_mc_sensitivity_compact.csv'));

        T_endpoint = ab_mis_build_endpoint_comparison(stats, specs);
        writetable(T_endpoint, fullfile(out_dir, 'mismatch_mc_endpoint_full_vs_no_mismatch.csv'));
    end

    %% ---------------- Save traceable data ----------------
    save(fullfile(out_dir, 'mismatch_sensitivity_results.mat'), ...
        'base_params','metrics_teaching','q_ref','dq_ref','ddq_ref','t_all', ...
        'specs','nominal','T_nominal','stats','T_mc_long','T_mc_summary', ...
        'T_mc_compact','T_endpoint','timestamp','-v7.3');

    if make_figures
        try
            ab_mismatch_generate_figures( ...
                out_dir, t_all, q_ref, nominal, stats, specs);
        catch ME
            warning('ab_mismatch_runner:FigureGenerationFailed', ...
                'Data were saved, but diagnostic figure generation failed: %s', ME.message);
        end
    end

    %% ---------------- Return compact result ----------------
    results = struct();
    results.output_dir = out_dir;
    results.timestamp = timestamp;
    results.specs = specs;
    results.nominal_table = T_nominal;
    results.mc_summary_long = T_mc_summary;
    results.mc_compact_summary = T_mc_compact;
    results.endpoint_comparison = T_endpoint;
    if run_mc
        results.success_rate = mean(stats.success_mask);
    else
        results.success_rate = NaN;
    end

    fprintf('\nMismatch diagnostic experiment complete.\n');
    fprintf('Results saved to: %s\n', out_dir);
end


%% ========================================================================

function specs = ab_mis_build_kmis_specs(kmis_grid)
    specs = struct([]);
    for v = 1:numel(kmis_grid)
        k = kmis_grid(v);
        specs(v).K_mis = k;
        specs(v).key = sprintf('kmis_%03d', round(100*k));
        if abs(k-1) < 1e-12
            specs(v).label = 'Full Proposed (Kmis=1)';
        elseif abs(k) < 1e-12
            specs(v).label = 'No Mismatch Compensation (Kmis=0)';
        else
            specs(v).label = sprintf('Kmis=%.2f', k);
        end
    end
end



function T = ab_mis_build_manifest(specs)
    nV = numel(specs);
    Variant = strings(nV,1);
    Key = strings(nV,1);
    K_mis = zeros(nV,1);
    Adapt_B = true(nV,1);
    Adapt_K = true(nV,1);
    Adapt_G = true(nV,1);
    for v = 1:nV
        Variant(v) = string(specs(v).label);
        Key(v) = string(specs(v).key);
        K_mis(v) = specs(v).K_mis;
    end
    T = table(Variant,Key,K_mis,Adapt_B,Adapt_K,Adapt_G);
end



function p = ab_mis_apply_kmis_variant(p, spec)
    p.K_mismatch_comp = repmat(spec.K_mis,1,3);
    p.ablation = struct();
    p.ablation.name = spec.label;
    p.ablation.adapt_B = true;
    p.ablation.adapt_K = true;
    p.ablation.adapt_G = true;
end



function R = ab_mis_run_one_variant(q_ref, dq_ref, ddq_ref, t_all, params, pre_fb)
    [delta_q, delta_dq, delta_ddq, ...
     Ba, Ka, Ga, quiet_zone_flag, e_used, e_raw, diag] = ...
        run_admittance_controller('proposed', ...
            q_ref, dq_ref, params, pre_fb, t_all, ddq_ref);

    R = struct();
    R.delta_q = delta_q;
    R.delta_dq = delta_dq;
    R.delta_ddq = delta_ddq;
    R.q = q_ref + delta_q;
    R.dq = dq_ref + delta_dq;
    R.ddq = ddq_ref + delta_ddq;
    R.Ba = Ba;
    R.Ka = Ka;
    R.Ga = Ga;
    R.quiet_zone_flag = quiet_zone_flag;
    R.e_used = e_used;
    R.e_raw = e_raw;
    R.diag = diag;
    R.metrics = compute_metrics(q_ref,R.q,R.dq,R.ddq,e_raw,e_used,params,diag);
    R.mechanism = ab_mis_compute_mechanism_metrics(R,q_ref,t_all,params);
end



function M = ab_mis_compute_mechanism_metrics(R, q_ref, t_all, params)
    M = struct();
    mismatch = R.diag.mismatch;
    pathological = R.diag.pathological;
    compensated = R.diag.compensated;
    raw = R.diag.raw;
    used = R.diag.used;

    M.MismatchRMS = ab_mis_rms_cols(mismatch);
    M.PathologicalRMS = ab_mis_rms_cols(pathological);
    M.CompensatedResidualRMS = ab_mis_rms_cols(compensated);
    M.RawResidualRMS = ab_mis_rms_cols(raw);
    M.ControllerInputRMS = ab_mis_rms_cols(used);
    M.MismatchLeakageRMS = ab_mis_rms_cols(compensated - pathological);

    M.CorrMismatchPathological = zeros(1,3);
    M.MismatchPathologicalMeanProduct = zeros(1,3);
    M.CancellationRatio = zeros(1,3);
    for j = 1:3
        M.CorrMismatchPathological(j) = ab_mis_safe_corr(mismatch(:,j), pathological(:,j));
        M.MismatchPathologicalMeanProduct(j) = mean(mismatch(:,j).*pathological(:,j));
        M.CancellationRatio(j) = M.RawResidualRMS(j) / ...
            max(M.MismatchRMS(j) + M.PathologicalRMS(j), eps);
    end

    % Event-window metrics. The mask is regenerated from the exact params of
    % this trial so event-time offsets are respected.
    [~, timing] = event_model(t_all, params);
    dt = params.dt;
    jerk = zeros(size(R.ddq));
    for j = 1:3
        jerk(:,j) = gradient(R.ddq(:,j),dt);
    end
    qerr = R.q - q_ref;

    M.EventTrackingRMSE_deg = nan(1,3);
    M.EventJerkRMS_deg_s3 = nan(1,3);
    M.EventRawResidualRMS = nan(1,3);
    M.EventControllerInputRMS = nan(1,3);
    M.EventMismatchRMS = nan(1,3);
    M.EventPathologicalRMS = nan(1,3);
    M.EventCorrMismatchPathological = nan(1,3);
    M.EventMismatchPathologicalMeanProduct = nan(1,3);

    for j = 1:3
        idx = timing.event_mask_joint(:,j);
        if sum(idx) > 1
            M.EventTrackingRMSE_deg(j) = sqrt(mean(qerr(idx,j).^2))*180/pi;
            M.EventJerkRMS_deg_s3(j) = sqrt(mean(jerk(idx,j).^2))*180/pi;
            M.EventRawResidualRMS(j) = sqrt(mean(raw(idx,j).^2));
            M.EventControllerInputRMS(j) = sqrt(mean(used(idx,j).^2));
            M.EventMismatchRMS(j) = sqrt(mean(mismatch(idx,j).^2));
            M.EventPathologicalRMS(j) = sqrt(mean(pathological(idx,j).^2));
            M.EventCorrMismatchPathological(j) = ab_mis_safe_corr(mismatch(idx,j),pathological(idx,j));
            M.EventMismatchPathologicalMeanProduct(j) = mean(mismatch(idx,j).*pathological(idx,j));
        end
    end
end



function y = ab_mis_rms_cols(X)
    y = sqrt(mean(X.^2,1));
end



function r = ab_mis_safe_corr(x,y)
    x = x(:); y = y(:);
    idx = isfinite(x) & isfinite(y);
    x = x(idx); y = y(idx);
    if numel(x) < 3 || std(x) < 1e-14 || std(y) < 1e-14
        r = NaN;
        return;
    end
    C = corrcoef(x,y);
    r = C(1,2);
end



function ok = ab_mis_variant_is_finite(R)
    % Correlations are allowed to be NaN only in the degenerate case where a
    % signal has essentially zero variance. All physical/performance outputs
    % must remain finite.
    required = [R.q(:);R.dq(:);R.ddq(:);R.Ba(:);R.Ka(:);R.Ga(:); ...
        R.e_raw(:);R.e_used(:);R.metrics.RMSE(:);R.metrics.RMS_jerk_deg(:); ...
        R.metrics.RMS_InteractionTorque(:);R.mechanism.MismatchRMS(:); ...
        R.mechanism.PathologicalRMS(:);R.mechanism.MismatchLeakageRMS(:); ...
        R.mechanism.EventTrackingRMSE_deg(:);R.mechanism.EventJerkRMS_deg_s3(:); ...
        R.mechanism.EventRawResidualRMS(:)];
    ok = all(isfinite(required));
end



function T = ab_mis_build_nominal_table(nominal, specs)
    joints = ["Hip";"Knee";"Ankle"];
    nRows = numel(specs)*3;

    K_mis = zeros(nRows,1);
    Joint = strings(nRows,1);
    RMSE_deg = nan(nRows,1);
    JerkRMS_deg_s3 = nan(nRows,1);
    RawResidualRMS_Nm = nan(nRows,1);
    ControllerInputRMS_Nm = nan(nRows,1);
    PositionSaturation_pct = nan(nRows,1);
    MismatchRMS_Nm = nan(nRows,1);
    PathologicalRMS_Nm = nan(nRows,1);
    CompensatedResidualRMS_Nm = nan(nRows,1);
    MismatchLeakageRMS_Nm = nan(nRows,1);
    Corr_Mismatch_Pathological = nan(nRows,1);
    MeanProduct_Mismatch_Pathological_Nm2 = nan(nRows,1);
    CancellationRatio = nan(nRows,1);
    Event_RMSE_deg = nan(nRows,1);
    Event_JerkRMS_deg_s3 = nan(nRows,1);
    Event_RawResidualRMS_Nm = nan(nRows,1);
    Event_ControllerInputRMS_Nm = nan(nRows,1);
    Event_MismatchRMS_Nm = nan(nRows,1);
    Event_PathologicalRMS_Nm = nan(nRows,1);
    Event_Corr_Mismatch_Pathological = nan(nRows,1);

    row = 0;
    for v = 1:numel(specs)
        R = nominal.(specs(v).key);
        for j = 1:3
            row = row + 1;
            K_mis(row) = specs(v).K_mis;
            Joint(row) = joints(j);
            RMSE_deg(row) = R.metrics.RMSE(j);
            JerkRMS_deg_s3(row) = R.metrics.RMS_jerk_deg(j);
            RawResidualRMS_Nm(row) = R.metrics.RMS_InteractionTorque(j);
            ControllerInputRMS_Nm(row) = R.metrics.RMS_ControllerTorque(j);
            PositionSaturation_pct(row) = R.metrics.PositionCorrectionSaturationPct(j);
            MismatchRMS_Nm(row) = R.mechanism.MismatchRMS(j);
            PathologicalRMS_Nm(row) = R.mechanism.PathologicalRMS(j);
            CompensatedResidualRMS_Nm(row) = R.mechanism.CompensatedResidualRMS(j);
            MismatchLeakageRMS_Nm(row) = R.mechanism.MismatchLeakageRMS(j);
            Corr_Mismatch_Pathological(row) = R.mechanism.CorrMismatchPathological(j);
            MeanProduct_Mismatch_Pathological_Nm2(row) = R.mechanism.MismatchPathologicalMeanProduct(j);
            CancellationRatio(row) = R.mechanism.CancellationRatio(j);
            Event_RMSE_deg(row) = R.mechanism.EventTrackingRMSE_deg(j);
            Event_JerkRMS_deg_s3(row) = R.mechanism.EventJerkRMS_deg_s3(j);
            Event_RawResidualRMS_Nm(row) = R.mechanism.EventRawResidualRMS(j);
            Event_ControllerInputRMS_Nm(row) = R.mechanism.EventControllerInputRMS(j);
            Event_MismatchRMS_Nm(row) = R.mechanism.EventMismatchRMS(j);
            Event_PathologicalRMS_Nm(row) = R.mechanism.EventPathologicalRMS(j);
            Event_Corr_Mismatch_Pathological(row) = R.mechanism.EventCorrMismatchPathological(j);
        end
    end

    T = table(K_mis,Joint,RMSE_deg,JerkRMS_deg_s3,RawResidualRMS_Nm, ...
        ControllerInputRMS_Nm,PositionSaturation_pct,MismatchRMS_Nm, ...
        PathologicalRMS_Nm,CompensatedResidualRMS_Nm,MismatchLeakageRMS_Nm, ...
        Corr_Mismatch_Pathological,MeanProduct_Mismatch_Pathological_Nm2, ...
        CancellationRatio,Event_RMSE_deg,Event_JerkRMS_deg_s3, ...
        Event_RawResidualRMS_Nm,Event_ControllerInputRMS_Nm, ...
        Event_MismatchRMS_Nm,Event_PathologicalRMS_Nm, ...
        Event_Corr_Mismatch_Pathological);
end



function T = ab_mis_build_nominal_timeseries(t_all,q_ref,nominal,specs)
    joints = {'Hip','Knee','Ankle'};
    T = table(t_all(:),'VariableNames',{'Time_s'});
    for j = 1:3
        T.(sprintf('Reference_%s_deg',joints{j})) = q_ref(:,j)*180/pi;
    end

    for v = 1:numel(specs)
        key = specs(v).key;
        R = nominal.(key);
        for j = 1:3
            T.(sprintf('%s_error_%s_deg',key,joints{j})) = (R.q(:,j)-q_ref(:,j))*180/pi;
            T.(sprintf('%s_raw_%s_Nm',key,joints{j})) = R.diag.raw(:,j);
            T.(sprintf('%s_comp_%s_Nm',key,joints{j})) = R.diag.compensated(:,j);
            T.(sprintf('%s_used_%s_Nm',key,joints{j})) = R.diag.used(:,j);
            T.(sprintf('%s_mismatch_%s_Nm',key,joints{j})) = R.diag.mismatch(:,j);
            T.(sprintf('%s_pathological_%s_Nm',key,joints{j})) = R.diag.pathological(:,j);
            T.(sprintf('%s_Ba_%s',key,joints{j})) = R.Ba(:,j);
            T.(sprintf('%s_Ka_%s',key,joints{j})) = R.Ka(:,j);
            T.(sprintf('%s_Ga_%s',key,joints{j})) = R.Ga(:,j);
        end
    end
end



function stats = ab_mis_init_stats(specs,num_runs)
    stats = struct();
    metric_fields = ab_mis_metric_field_names();
    for v = 1:numel(specs)
        key = specs(v).key;
        S = struct();
        for m = 1:numel(metric_fields)
            S.(metric_fields{m}) = nan(num_runs,3);
        end
        stats.(key) = S;
    end
end



function names = ab_mis_metric_field_names()
    names = { ...
        'RMSE','JerkRMS','RawResidualRMS','ControllerInputRMS','PositionSaturationPct', ...
        'MismatchRMS','PathologicalRMS','CompensatedResidualRMS','MismatchLeakageRMS', ...
        'CorrMismatchPathological','MismatchPathologicalMeanProduct','CancellationRatio', ...
        'EventRMSE','EventJerkRMS','EventRawResidualRMS','EventControllerInputRMS', ...
        'EventMismatchRMS','EventPathologicalRMS','EventCorrMismatchPathological', ...
        'EventMismatchPathologicalMeanProduct'};
end



function M = ab_mis_extract_trial_metrics(R)
    M = struct();
    M.RMSE = R.metrics.RMSE;
    M.JerkRMS = R.metrics.RMS_jerk_deg;
    M.RawResidualRMS = R.metrics.RMS_InteractionTorque;
    M.ControllerInputRMS = R.metrics.RMS_ControllerTorque;
    M.PositionSaturationPct = R.metrics.PositionCorrectionSaturationPct;
    M.MismatchRMS = R.mechanism.MismatchRMS;
    M.PathologicalRMS = R.mechanism.PathologicalRMS;
    M.CompensatedResidualRMS = R.mechanism.CompensatedResidualRMS;
    M.MismatchLeakageRMS = R.mechanism.MismatchLeakageRMS;
    M.CorrMismatchPathological = R.mechanism.CorrMismatchPathological;
    M.MismatchPathologicalMeanProduct = R.mechanism.MismatchPathologicalMeanProduct;
    M.CancellationRatio = R.mechanism.CancellationRatio;
    M.EventRMSE = R.mechanism.EventTrackingRMSE_deg;
    M.EventJerkRMS = R.mechanism.EventJerkRMS_deg_s3;
    M.EventRawResidualRMS = R.mechanism.EventRawResidualRMS;
    M.EventControllerInputRMS = R.mechanism.EventControllerInputRMS;
    M.EventMismatchRMS = R.mechanism.EventMismatchRMS;
    M.EventPathologicalRMS = R.mechanism.EventPathologicalRMS;
    M.EventCorrMismatchPathological = R.mechanism.EventCorrMismatchPathological;
    M.EventMismatchPathologicalMeanProduct = R.mechanism.EventMismatchPathologicalMeanProduct;
end



function S = ab_mis_commit_trial(S,M,r)
    fields = fieldnames(M);
    for k = 1:numel(fields)
        S.(fields{k})(r,:) = M.(fields{k});
    end
end



function T = ab_mis_build_mc_long_table(stats,specs)
    joints = ["Hip","Knee","Ankle"];
    num_runs = numel(stats.success_mask);
    nRows = num_runs*numel(specs)*3;
    names = ab_mis_metric_field_names();

    Trial = zeros(nRows,1);
    Success = false(nRows,1);
    K_mis = nan(nRows,1);
    Joint = strings(nRows,1);
    data = struct();
    for m = 1:numel(names)
        data.(names{m}) = nan(nRows,1);
    end

    row = 0;
    for r = 1:num_runs
        for v = 1:numel(specs)
            key = specs(v).key;
            for j = 1:3
                row = row + 1;
                Trial(row) = r;
                Success(row) = stats.success_mask(r);
                K_mis(row) = specs(v).K_mis;
                Joint(row) = joints(j);
                for m = 1:numel(names)
                    data.(names{m})(row) = stats.(key).(names{m})(r,j);
                end
            end
        end
    end

    T = table(Trial,Success,K_mis,Joint);
    for m = 1:numel(names)
        T.(names{m}) = data.(names{m});
    end
end



function T = ab_mis_build_mc_summary_long(stats,specs)
    joints = ["Hip","Knee","Ankle"];
    names = ab_mis_metric_field_names();
    nRows = numel(specs)*3*numel(names);

    K_mis = nan(nRows,1);
    Joint = strings(nRows,1);
    Metric = strings(nRows,1);
    N = zeros(nRows,1);
    Mean = nan(nRows,1);
    Std = nan(nRows,1);
    Median = nan(nRows,1);
    Mean95CI_Low = nan(nRows,1);
    Mean95CI_High = nan(nRows,1);

    row = 0;
    for v = 1:numel(specs)
        key = specs(v).key;
        for m = 1:numel(names)
            X = stats.(key).(names{m});
            for j = 1:3
                row = row + 1;
                x = X(:,j);
                x = x(isfinite(x));
                K_mis(row) = specs(v).K_mis;
                Joint(row) = joints(j);
                Metric(row) = string(names{m});
                N(row) = numel(x);
                if ~isempty(x)
                    Mean(row) = mean(x);
                    Std(row) = std(x);
                    Median(row) = median(x);
                    [Mean95CI_Low(row),Mean95CI_High(row)] = ab_mis_normal95_ci(x);
                end
            end
        end
    end

    T = table(K_mis,Joint,Metric,N,Mean,Std,Median,Mean95CI_Low,Mean95CI_High);
end



function T = ab_mis_build_mc_compact_summary(stats,specs)
    joints = ["Hip";"Knee";"Ankle"];
    nRows = numel(specs)*3;

    K_mis = nan(nRows,1);
    Joint = strings(nRows,1);
    N = zeros(nRows,1);
    RMSE_Mean = nan(nRows,1); RMSE_Std = nan(nRows,1);
    JerkRMS_Mean = nan(nRows,1); JerkRMS_Std = nan(nRows,1);
    RawResidualRMS_Mean = nan(nRows,1); RawResidualRMS_Std = nan(nRows,1);
    ControllerInputRMS_Mean = nan(nRows,1);
    MismatchRMS_Mean = nan(nRows,1);
    PathologicalRMS_Mean = nan(nRows,1);
    MismatchLeakageRMS_Mean = nan(nRows,1);
    CorrMismatchPathological_Mean = nan(nRows,1);
    MeanProductMismatchPathological_Mean = nan(nRows,1);
    CancellationRatio_Mean = nan(nRows,1);
    EventRMSE_Mean = nan(nRows,1);
    EventJerkRMS_Mean = nan(nRows,1);
    EventRawResidualRMS_Mean = nan(nRows,1);
    EventCorrMismatchPathological_Mean = nan(nRows,1);
    PositionSaturationPct_Mean = nan(nRows,1);

    row = 0;
    for v = 1:numel(specs)
        S = stats.(specs(v).key);
        for j = 1:3
            row = row + 1;
            K_mis(row) = specs(v).K_mis;
            Joint(row) = joints(j);
            good = isfinite(S.RMSE(:,j));
            N(row) = sum(good);
            RMSE_Mean(row) = mean(S.RMSE(good,j)); RMSE_Std(row) = std(S.RMSE(good,j));
            JerkRMS_Mean(row) = mean(S.JerkRMS(good,j)); JerkRMS_Std(row) = std(S.JerkRMS(good,j));
            RawResidualRMS_Mean(row) = mean(S.RawResidualRMS(good,j)); RawResidualRMS_Std(row) = std(S.RawResidualRMS(good,j));
            ControllerInputRMS_Mean(row) = mean(S.ControllerInputRMS(good,j));
            MismatchRMS_Mean(row) = mean(S.MismatchRMS(good,j));
            PathologicalRMS_Mean(row) = mean(S.PathologicalRMS(good,j));
            MismatchLeakageRMS_Mean(row) = mean(S.MismatchLeakageRMS(good,j));
            CorrMismatchPathological_Mean(row) = mean(S.CorrMismatchPathological(good,j),'omitnan');
            MeanProductMismatchPathological_Mean(row) = mean(S.MismatchPathologicalMeanProduct(good,j),'omitnan');
            CancellationRatio_Mean(row) = mean(S.CancellationRatio(good,j));
            EventRMSE_Mean(row) = mean(S.EventRMSE(good,j));
            EventJerkRMS_Mean(row) = mean(S.EventJerkRMS(good,j));
            EventRawResidualRMS_Mean(row) = mean(S.EventRawResidualRMS(good,j));
            EventCorrMismatchPathological_Mean(row) = mean(S.EventCorrMismatchPathological(good,j),'omitnan');
            PositionSaturationPct_Mean(row) = mean(S.PositionSaturationPct(good,j));
        end
    end

    T = table(K_mis,Joint,N,RMSE_Mean,RMSE_Std,JerkRMS_Mean,JerkRMS_Std, ...
        RawResidualRMS_Mean,RawResidualRMS_Std,ControllerInputRMS_Mean, ...
        MismatchRMS_Mean,PathologicalRMS_Mean,MismatchLeakageRMS_Mean, ...
        CorrMismatchPathological_Mean,MeanProductMismatchPathological_Mean, ...
        CancellationRatio_Mean,EventRMSE_Mean,EventJerkRMS_Mean, ...
        EventRawResidualRMS_Mean,EventCorrMismatchPathological_Mean, ...
        PositionSaturationPct_Mean);
end



function T = ab_mis_build_endpoint_comparison(stats,specs)
    idx0 = find(abs([specs.K_mis]-0)<1e-12,1);
    idx1 = find(abs([specs.K_mis]-1)<1e-12,1);
    if isempty(idx0) || isempty(idx1)
        T = table();
        return;
    end

    joints = ["Hip","Knee","Ankle"];
    metrics = {'RMSE','JerkRMS','RawResidualRMS','ControllerInputRMS', ...
        'MismatchLeakageRMS','MismatchRMS','PathologicalRMS','CancellationRatio', ...
        'EventRMSE','EventJerkRMS','EventRawResidualRMS'};
    nRows = numel(metrics)*3;

    Joint = strings(nRows,1);
    Metric = strings(nRows,1);
    N = zeros(nRows,1);
    Full_Kmis1_Mean = nan(nRows,1);
    NoComp_Kmis0_Mean = nan(nRows,1);
    NoCompMinusFull_Mean = nan(nRows,1);
    Difference95CI_Low = nan(nRows,1);
    Difference95CI_High = nan(nRows,1);
    FullImprovementPct_Mean = nan(nRows,1);

    S0 = stats.(specs(idx0).key);
    S1 = stats.(specs(idx1).key);
    row = 0;
    for m = 1:numel(metrics)
        field = metrics{m};
        for j = 1:3
            row = row + 1;
            f = S1.(field)(:,j);
            n = S0.(field)(:,j);
            good = isfinite(f) & isfinite(n);
            f = f(good); n = n(good);
            Joint(row) = joints(j);
            Metric(row) = string(field);
            N(row) = numel(f);
            if ~isempty(f)
                d = n-f;
                Full_Kmis1_Mean(row) = mean(f);
                NoComp_Kmis0_Mean(row) = mean(n);
                NoCompMinusFull_Mean(row) = mean(d);
                [Difference95CI_Low(row),Difference95CI_High(row)] = ab_mis_normal95_ci(d);
                safe = abs(n)>eps;
                if any(safe)
                    FullImprovementPct_Mean(row) = mean(100*(n(safe)-f(safe))./n(safe));
                end
            end
        end
    end

    T = table(Joint,Metric,N,Full_Kmis1_Mean,NoComp_Kmis0_Mean, ...
        NoCompMinusFull_Mean,Difference95CI_Low,Difference95CI_High, ...
        FullImprovementPct_Mean);
end



function [lo,hi] = ab_mis_normal95_ci(x)
    x = x(isfinite(x));
    if isempty(x)
        lo = NaN; hi = NaN;
    elseif numel(x)==1
        lo = x; hi = x;
    else
        mu = mean(x);
        se = std(x)/sqrt(numel(x));
        lo = mu-1.96*se;
        hi = mu+1.96*se;
    end
end


%% ========================================================================
% Monte Carlo sampler mirrors the current project sampling design.


function ab_mismatch_generate_figures(out_dir,t_all,q_ref,nominal,stats,specs)
%GENERATE_MISMATCH_DIAGNOSTIC_FIGURES Figures for K_mis mechanism/sensitivity.
% Each exported PDF/PNG receives an identically named CSV and MAT in
% figure_data/ containing the values actually plotted.

    fig_dir=fullfile(out_dir,'figures');
    dat_dir=fullfile(out_dir,'figure_data');
    if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
    if ~exist(dat_dir,'dir'), mkdir(dat_dir); end
    kgrid=[specs.K_mis]; [kgrid,order]=sort(kgrid); specs=specs(order);
    joints={'Hip','Knee','Ankle'};

    %% M1: Monte Carlo primary-metric sensitivity
    if ~isempty(stats) && isfield(stats,'success_mask') && any(stats.success_mask)
        base_name='FigM1_mc_kmis_primary_metric_sensitivity';
        fields={'RMSE','JerkRMS','RawResidualRMS'};
        titles={'Tracking RMSE','Jerk RMS','Raw interaction-residual RMS'};
        ylabs={'RMSE (deg)','Jerk RMS (deg/s^3)','Raw residual RMS (N m)'};
        [f1,lay1]=paper_figure_utils('create','Kmis Primary Sensitivity', ...
            [80 45 1120 1180],3,1, ...
            sprintf('Kmis Sensitivity (%d Paired Monte Carlo Trials)',sum(stats.success_mask)));
        ax1=gobjects(3,1);
        nr=3*3*numel(kgrid); Metric=strings(nr,1); Joint=strings(nr,1); K_mis=nan(nr,1); Mean=nan(nr,1); CI95_Low=nan(nr,1); CI95_High=nan(nr,1); N=zeros(nr,1); row=0;
        for m=1:3
            ax=paper_figure_utils('tile',lay1,m); ax1(m)=ax; hold(ax,'on');
            for j=1:3
                mu=nan(size(kgrid)); lo=mu; hi=mu; nn=zeros(size(kgrid));
                for v=1:numel(specs)
                    x=stats.(specs(v).key).(fields{m})(:,j); x=x(isfinite(x)); nn(v)=numel(x);
                    if ~isempty(x)
                        mu(v)=mean(x); se=std(x)/sqrt(numel(x)); lo(v)=mu(v)-1.96*se; hi(v)=mu(v)+1.96*se;
                    end
                    row=row+1; Metric(row)=string(titles{m}); Joint(row)=string(joints{j}); K_mis(row)=kgrid(v);
                    Mean(row)=mu(v); CI95_Low(row)=lo(v); CI95_High(row)=hi(v); N(row)=nn(v);
                end
                errorbar(ax,kgrid,mu,mu-lo,hi-mu,'-o','LineWidth',1.2,'MarkerSize',5,'DisplayName',joints{j});
            end
            xlim(ax,[min(kgrid) max(kgrid)]); xticks(ax,kgrid); ylabel(ax,ylabs{m});
            if m==3, xlabel(ax,'Mismatch-Compensation Coefficient K_{mis}'); end
            paper_figure_utils('panel_title',ax,sprintf('(%c) %s (Mean and 95%% CI)',char('a'+m-1),titles{m}));
        end
        lgd = gobjects(numel(ax1),1);
        for ii = 1:numel(ax1)
            lgd(ii)=legend(ax1(ii),'Location','northeast','Interpreter','none','AutoUpdate','off');
        end
        paper_figure_utils('align_legends',ax1,lgd);
        Tfig=table(Metric,Joint,K_mis,Mean,CI95_Low,CI95_High,N);
        meta=struct('figure',base_name,'description','K_mis Monte Carlo primary-metric sensitivity with plotted normal-approximation 95% CI', ...
            'successful_trials',sum(stats.success_mask),'units',{ylabs});
        ab_save_figure_data(dat_dir,base_name,Tfig,meta);
        ab_mis_export_pair(f1,fig_dir,base_name);
    end

    %% M2: Monte Carlo mechanism sensitivity
    if ~isempty(stats) && isfield(stats,'success_mask') && any(stats.success_mask)
        base_name='FigM2_mc_kmis_mechanism_sensitivity';
        fields={'MismatchLeakageRMS','ControllerInputRMS','CorrMismatchPathological','CancellationRatio'};
        titles={'Mismatch Leakage RMS','Processed Controller-Input RMS', ...
            'Corr(Mismatch, Pathological)','Raw-Residual Cancellation Ratio'};
        ylabs={'RMS (N m)','RMS (N m)','Correlation','RMS(raw)/(RMS(mis)+RMS(path))'};
        [f2,lay2]=paper_figure_utils('create','Kmis Mechanism Diagnostics', ...
            [110 65 1250 930],2,2,'Kmis Mechanism Diagnostics');
        ax2=gobjects(4,1);
        nr=4*3*numel(kgrid); Metric=strings(nr,1); Joint=strings(nr,1); K_mis=nan(nr,1); Mean=nan(nr,1); N=zeros(nr,1); row=0;
        for m=1:4
            ax=paper_figure_utils('tile',lay2,m); ax2(m)=ax; hold(ax,'on');
            for j=1:3
                mu=nan(size(kgrid));
                for v=1:numel(specs)
                    x=stats.(specs(v).key).(fields{m})(:,j); x=x(isfinite(x));
                    if ~isempty(x), mu(v)=mean(x); end
                    row=row+1; Metric(row)=string(titles{m}); Joint(row)=string(joints{j}); K_mis(row)=kgrid(v); Mean(row)=mu(v); N(row)=numel(x);
                end
                plot(ax,kgrid,mu,'-o','LineWidth',1.2,'MarkerSize',5,'DisplayName',joints{j});
            end
            xlim(ax,[min(kgrid) max(kgrid)]); xticks(ax,kgrid); ylabel(ax,ylabs{m}); xlabel(ax,'K_{mis}');
            paper_figure_utils('panel_title',ax,sprintf('(%c) %s',char('a'+m-1),titles{m}));
        end
        lgd = gobjects(numel(ax2),1);
        for ii = 1:numel(ax2)
            lgd(ii)=legend(ax2(ii),'Location','northeast','Interpreter','none','AutoUpdate','off');
        end
        paper_figure_utils('align_legends',ax2,lgd);
        Tfig=table(Metric,Joint,K_mis,Mean,N);
        meta=struct('figure',base_name,'description','K_mis Monte Carlo mechanism sensitivity; values are the plotted cross-trial means', ...
            'successful_trials',sum(stats.success_mask),'units',{ylabs});
        ab_save_figure_data(dat_dir,base_name,Tfig,meta);
        ab_mis_export_pair(f2,fig_dir,base_name);
    end

    %% M3: Nominal knee mechanism traces at the two endpoints
    idx0=find(abs([specs.K_mis]-0)<1e-12,1); idx1=find(abs([specs.K_mis]-1)<1e-12,1);
    if ~isempty(idx0) && ~isempty(idx1)
        base_name='FigM3_nominal_knee_endpoint_mechanism';
        R0=nominal.(specs(idx0).key); R1=nominal.(specs(idx1).key); j=2;
        [f3,lay3]=paper_figure_utils('create','Knee Compensation Mechanism', ...
            [140 35 1150 1320],4,1,'Knee: Full Compensation vs No Compensation (Nominal)');
        ax3=gobjects(4,1); lgd3=gobjects(4,1);

        track_full=(R1.q(:,j)-q_ref(:,j))*180/pi;
        track_none=(R0.q(:,j)-q_ref(:,j))*180/pi;
        ax=paper_figure_utils('tile',lay3,1); ax3(1)=ax; hold(ax,'on');
        plot(ax,t_all,track_full,'LineWidth',1.2,'DisplayName','Full Kmis=1');
        plot(ax,t_all,track_none,'LineWidth',1.2,'DisplayName','No comp. Kmis=0');
        ylabel(ax,'Error (deg)'); lgd3(1)=legend(ax,'Location','northeast','AutoUpdate','off');
        paper_figure_utils('panel_title',ax,'(a) Tracking Error');

        ax=paper_figure_utils('tile',lay3,2); ax3(2)=ax; hold(ax,'on');
        plot(ax,t_all,R0.diag.mismatch(:,j),'LineWidth',1.1,'DisplayName','Mismatch');
        plot(ax,t_all,R0.diag.pathological(:,j),'LineWidth',1.1,'DisplayName','Pathological');
        plot(ax,t_all,R0.diag.raw(:,j),'LineWidth',1.1,'DisplayName','Raw = mismatch + pathological');
        ylabel(ax,'N m'); lgd3(2)=legend(ax,'Location','northeast','AutoUpdate','off');
        paper_figure_utils('panel_title',ax,'(b) No-Compensation Interaction Decomposition');

        ax=paper_figure_utils('tile',lay3,3); ax3(3)=ax; hold(ax,'on');
        plot(ax,t_all,R1.diag.used(:,j),'LineWidth',1.1,'DisplayName','Full Kmis=1');
        plot(ax,t_all,R0.diag.used(:,j),'LineWidth',1.1,'DisplayName','No comp. Kmis=0');
        ylabel(ax,'N m'); lgd3(3)=legend(ax,'Location','northeast','AutoUpdate','off');
        paper_figure_utils('panel_title',ax,'(c) Processed Controller Input e_use');

        ax=paper_figure_utils('tile',lay3,4); ax3(4)=ax; hold(ax,'on');
        plot(ax,t_all,R1.Ga(:,j),'LineWidth',1.1,'DisplayName','Full Kmis=1');
        plot(ax,t_all,R0.Ga(:,j),'LineWidth',1.1,'DisplayName','No comp. Kmis=0');
        ylabel(ax,'G_a'); xlabel(ax,'Time (s)'); lgd3(4)=legend(ax,'Location','northeast','AutoUpdate','off');
        paper_figure_utils('panel_title',ax,'(d) Interaction-Mapping Gain Response');
        paper_figure_utils('align_legends',ax3,lgd3);

        Time_s=t_all(:); TrackingError_Full_deg=track_full(:); TrackingError_NoComp_deg=track_none(:);
        NoComp_Mismatch_Nm=R0.diag.mismatch(:,j); NoComp_Pathological_Nm=R0.diag.pathological(:,j); NoComp_Raw_Nm=R0.diag.raw(:,j);
        ControllerInput_Full_Nm=R1.diag.used(:,j); ControllerInput_NoComp_Nm=R0.diag.used(:,j);
        Ga_Full=R1.Ga(:,j); Ga_NoComp=R0.Ga(:,j);
        Tfig=table(Time_s,TrackingError_Full_deg,TrackingError_NoComp_deg,NoComp_Mismatch_Nm, ...
            NoComp_Pathological_Nm,NoComp_Raw_Nm,ControllerInput_Full_Nm,ControllerInput_NoComp_Nm,Ga_Full,Ga_NoComp);
        meta=struct('figure',base_name,'description','Nominal knee endpoint mechanism traces for K_mis=1 and K_mis=0','joint','Knee');
        ab_save_figure_data(dat_dir,base_name,Tfig,meta);
        ab_mis_export_pair(f3,fig_dir,base_name);
    end

    %% M4: Event-window sensitivity
    if ~isempty(stats) && isfield(stats,'success_mask') && any(stats.success_mask)
        base_name='FigM4_mc_event_window_kmis_sensitivity';
        fields={'EventRMSE','EventJerkRMS','EventRawResidualRMS'};
        titles={'Event-Window Tracking RMSE','Event-Window Jerk RMS','Event-Window Raw Residual RMS'};
        ylabs={'RMSE (deg)','Jerk RMS (deg/s^3)','Raw residual RMS (N m)'};
        [f4,lay4]=paper_figure_utils('create','Kmis Event-Window Sensitivity', ...
            [170 45 1120 1180],3,1,'Kmis Sensitivity During Joint-Specific Event Windows');
        ax4=gobjects(3,1);
        nr=3*3*numel(kgrid); Metric=strings(nr,1); Joint=strings(nr,1); K_mis=nan(nr,1); Mean=nan(nr,1); N=zeros(nr,1); row=0;
        for m=1:3
            ax=paper_figure_utils('tile',lay4,m); ax4(m)=ax; hold(ax,'on');
            for j=1:3
                mu=nan(size(kgrid));
                for v=1:numel(specs)
                    x=stats.(specs(v).key).(fields{m})(:,j); x=x(isfinite(x));
                    if ~isempty(x), mu(v)=mean(x); end
                    row=row+1; Metric(row)=string(titles{m}); Joint(row)=string(joints{j}); K_mis(row)=kgrid(v); Mean(row)=mu(v); N(row)=numel(x);
                end
                plot(ax,kgrid,mu,'-o','LineWidth',1.2,'MarkerSize',5,'DisplayName',joints{j});
            end
            xlim(ax,[min(kgrid) max(kgrid)]); xticks(ax,kgrid); ylabel(ax,ylabs{m});
            if m==3, xlabel(ax,'Mismatch-Compensation Coefficient K_{mis}'); end
            paper_figure_utils('panel_title',ax,sprintf('(%c) %s',char('a'+m-1),titles{m}));
        end
        lgd = gobjects(numel(ax4),1);
        for ii = 1:numel(ax4)
            lgd(ii)=legend(ax4(ii),'Location','northeast','AutoUpdate','off');
        end
        paper_figure_utils('align_legends',ax4,lgd);
        Tfig=table(Metric,Joint,K_mis,Mean,N);
        meta=struct('figure',base_name,'description','K_mis sensitivity of metrics computed within joint-specific event windows', ...
            'successful_trials',sum(stats.success_mask),'units',{ylabs});
        ab_save_figure_data(dat_dir,base_name,Tfig,meta);
        ab_mis_export_pair(f4,fig_dir,base_name);
    end

    % Keep figures visible after export, consistent with generate_paper_figures.m.
end

function ab_mis_export_pair(fig,fig_dir,base_name)
    paper_figure_utils('export',fig,fig_dir,base_name, ...
        'SavePNG',true,'SaveJPEG',true,'JPEGResolution',1000);
end


function R = ab_finalize_evidence(varargin)
%FINALIZE_ABLATION_EVIDENCE Aggregate the completed core-ablation evidence.
%
% This function does NOT rerun simulations. It consumes the CSV outputs from:
%   1) core_ablation_runner.m
%   2) mismatch_mechanism_runner.m
%
% Primary core comparison:
%   Full Proposed
%   No Mismatch Compensation
%   Simplified Adaptation (B/K)
%
% Primary performance metrics (all lower is better):
%   Tracking RMSE (deg)
%   Jerk RMS (deg/s^3)
%   Raw interaction-residual RMS (N m)
%
% Statistical analysis:
%   - trial-paired comparator-minus-Full differences;
%   - paired Student t test (two-sided), implemented without Statistics Toolbox;
%   - t-based 95% CI of the paired mean difference;
%   - paired standardized mean difference d_z;
%   - percentage of paired trials in which Full has the lower metric;
%   - Holm correction across all 18 primary endpoint tests
%       (2 comparators x 3 joints x 3 metrics).
%
% Mismatch sensitivity is retained as a mechanism/interpretation experiment;
% it is not used to retune K_mis or redefine the Full Proposed controller.
%
% Examples
% --------
% Automatic selection of the completed 100-run result folders under ./results:
%   R = ab_finalize_evidence();
%
% Explicit folders:
%   R = ab_finalize_evidence( ...
%       'AblationDir','results/ABLATION_20260924_195811', ...
%       'MismatchDir','results/MISMATCH_DIAG_20260924_202239');
%
% Smoke-test finalization (normally not for paper output):
%   R = ab_finalize_evidence('ExpectedRuns',5,'RequireExactRuns',true);

    %% ---------------- Options ----------------
    ip = inputParser;
    addParameter(ip,'ResultsRoot','results',@(x)ischar(x)||isstring(x));
    addParameter(ip,'AblationDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'MismatchDir','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'OutputRoot','',@(x)ischar(x)||isstring(x));
    addParameter(ip,'ExpectedRuns',100,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(ip,'RequireExactRuns',true,@(x)islogical(x)||isnumeric(x));
    addParameter(ip,'MakeFigures',true,@(x)islogical(x)||isnumeric(x));
    parse(ip,varargin{:});
    opt = ip.Results;

    expected_runs = round(opt.ExpectedRuns);
    require_exact = logical(opt.RequireExactRuns);
    make_figures = logical(opt.MakeFigures);
    results_root = char(opt.ResultsRoot);

    %% ---------------- Locate completed data ----------------
    ablation_dir = char(opt.AblationDir);
    mismatch_dir = char(opt.MismatchDir);

    if isempty(ablation_dir)
        ablation_dir = ab_fin_resolve_completed_result_dir( ...
            results_root,'ABLATION_','ablation_mc_trial_metrics_long.csv', ...
            expected_runs,require_exact,'Variant');
    end
    if isempty(mismatch_dir)
        mismatch_dir = ab_fin_resolve_completed_result_dir( ...
            results_root,'MISMATCH_DIAG_','mismatch_mc_trial_metrics_long.csv', ...
            expected_runs,require_exact,'K_mis');
    end

    ab_fin_validate_required_file(ablation_dir,'ablation_mc_trial_metrics_long.csv');
    ab_fin_validate_required_file(ablation_dir,'ablation_nominal_metrics.csv');
    ab_fin_validate_required_file(mismatch_dir,'mismatch_mc_sensitivity_compact.csv');
    ab_fin_validate_required_file(mismatch_dir,'mismatch_nominal_sensitivity.csv');

    %% ---------------- Output folder ----------------
    output_root = char(opt.OutputRoot);
    if isempty(output_root), output_root = results_root; end
    if ~exist(output_root,'dir'), mkdir(output_root); end
    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    out_dir = fullfile(output_root,['FINAL_ABLATION_EVIDENCE_' timestamp]);
    mkdir(out_dir);

    fprintf('\n============================================================\n');
    fprintf(' Final core-ablation evidence aggregation\n');
    fprintf(' Ablation source : %s\n',ablation_dir);
    fprintf(' Mismatch source : %s\n',mismatch_dir);
    fprintf(' Output          : %s\n',out_dir);
    fprintf('============================================================\n\n');

    %% ---------------- Read source tables ----------------
    Ttrial = readtable(fullfile(ablation_dir,'ablation_mc_trial_metrics_long.csv'), ...
        'TextType','string');
    Tnom = readtable(fullfile(ablation_dir,'ablation_nominal_metrics.csv'), ...
        'TextType','string');
    Tmis = readtable(fullfile(mismatch_dir,'mismatch_mc_sensitivity_compact.csv'), ...
        'TextType','string');
    TmisNom = readtable(fullfile(mismatch_dir,'mismatch_nominal_sensitivity.csv'), ...
        'TextType','string');

    required_trial_vars = {'Trial','Success','Variant','Joint','RMSE_deg', ...
        'JerkRMS_deg_s3','RawResidualRMS_Nm'};
    ab_fin_assert_table_vars(Ttrial,required_trial_vars,'ablation trial table');

    success_trial_ids = unique(Ttrial.Trial(Ttrial.Success==1));
    n_success_trials = numel(success_trial_ids);
    if require_exact && n_success_trials ~= expected_runs
        error('ab_finalize_evidence:UnexpectedRunCount', ...
            'Expected exactly %d successful paired trials, but found %d in %s.', ...
            expected_runs,n_success_trials,ablation_dir);
    elseif n_success_trials < expected_runs
        warning('ab_finalize_evidence:FewerRuns', ...
            'Only %d successful paired trials found (expected %d).', ...
            n_success_trials,expected_runs);
    end

    %% ---------------- Final MC mean +/- std table ----------------
    Tmc = ab_fin_build_core_mc_table(Ttrial);
    writetable(Tmc,fullfile(out_dir,'TableA_core_ablation_mc_mean_std.csv'));

    %% ---------------- Paired inferential statistics ----------------
    Tstats = ab_fin_build_pairwise_statistics(Ttrial);
    Tstats.P_Holm = ab_fin_holm_adjust(Tstats.P_Raw);
    Tstats.Significant_Holm_0p05 = Tstats.P_Holm < 0.05;
    writetable(Tstats,fullfile(out_dir,'TableB_core_ablation_pairwise_statistics.csv'));

    %% ---------------- Nominal compact table ----------------
    Tnominal = ab_fin_build_nominal_compact(Tnom);
    writetable(Tnominal,fullfile(out_dir,'TableC_core_ablation_nominal.csv'));

    %% ---------------- Mismatch mechanism tables ----------------
    Tmechanism = ab_fin_build_mismatch_mechanism_table(Tmis);
    writetable(Tmechanism,fullfile(out_dir,'TableD_mismatch_sensitivity_mechanism.csv'));

    Tendpoint = ab_fin_build_mismatch_endpoint_table(Tmis);
    writetable(Tendpoint,fullfile(out_dir,'TableE_mismatch_endpoint_mechanism.csv'));

    TmechanismNom = ab_fin_build_mismatch_nominal_endpoint_table(TmisNom);
    writetable(TmechanismNom,fullfile(out_dir,'TableF_mismatch_nominal_endpoint_mechanism.csv'));

    %% ---------------- Machine-readable manifest ----------------
    key = ["AblationSource";"MismatchSource";"SuccessfulPairedTrials"; ...
           "ExpectedTrials";"PrimaryMetricDirection";"DifferenceConvention"; ...
           "MultiplicityCorrection";"MismatchSensitivityPurpose"];
    value = [string(ablation_dir);string(mismatch_dir);string(n_success_trials); ...
             string(expected_runs);"Lower is better"; ...
             "Comparator minus Full; positive means Full has lower metric"; ...
             "Holm across 18 core endpoint tests"; ...
             "Mechanism/sensitivity only; no K_mis retuning"];
    Tmanifest = table(key,value,'VariableNames',{'Item','Value'});
    writetable(Tmanifest,fullfile(out_dir,'evidence_manifest.csv'));

    %% ---------------- Figures ----------------
    if make_figures
        ab_final_generate_figures(out_dir,Ttrial,Tstats,Tmis);
    end

    %% ---------------- Save traceable MAT ----------------
    save(fullfile(out_dir,'final_ablation_evidence.mat'), ...
        'Ttrial','Tnom','Tmis','TmisNom','Tmc','Tstats','Tnominal', ...
        'Tmechanism','Tendpoint','TmechanismNom','Tmanifest', ...
        'ablation_dir','mismatch_dir','timestamp','-v7.3');

    ab_fin_write_readme(out_dir,n_success_trials);

    %% ---------------- Return ----------------
    R = struct();
    R.output_dir = out_dir;
    R.ablation_source = ablation_dir;
    R.mismatch_source = mismatch_dir;
    R.num_successful_paired_trials = n_success_trials;
    R.mc_table = Tmc;
    R.pairwise_statistics = Tstats;
    R.nominal_table = Tnominal;
    R.mismatch_mechanism = Tmechanism;

    fprintf('\nFinal evidence aggregation complete.\n');
    fprintf('Paired trials used: %d\n',n_success_trials);
    fprintf('Saved to: %s\n',out_dir);
end

%% ========================================================================

function out_dir = ab_fin_resolve_completed_result_dir(root_dir,prefix,required_csv,expected_runs,require_exact,group_var)
    if ~exist(root_dir,'dir')
        error('ab_finalize_evidence:MissingResultsRoot', ...
            'Results root does not exist: %s',root_dir);
    end

    D = dir(fullfile(root_dir,[prefix '*']));
    D = D([D.isdir]);
    if isempty(D)
        error('ab_finalize_evidence:NoResultFolder', ...
            'No %s* folders found under %s.',prefix,root_dir);
    end

    best_idx = [];
    best_n = -inf;
    best_date = -inf;
    for k = 1:numel(D)
        p = fullfile(D(k).folder,D(k).name);
        f = fullfile(p,required_csv);
        if ~exist(f,'file'), continue; end
        try
            T = readtable(f,'TextType','string');
            if ~ismember('Trial',T.Properties.VariableNames), continue; end
            if ismember('Success',T.Properties.VariableNames)
                trial_ids = unique(T.Trial(T.Success==1));
            else
                trial_ids = unique(T.Trial);
            end
            n = numel(trial_ids);
            % group_var is checked to reject malformed tables.
            if ~ismember(group_var,T.Properties.VariableNames), continue; end
        catch
            continue;
        end

        if require_exact && n ~= expected_runs
            % Still keep as a fallback only if no exact folder is found.
        end

        exact_bonus = double(n==expected_runs)*1e9;
        score = exact_bonus + n;
        if score > best_n || (score==best_n && D(k).datenum>best_date)
            best_idx = k;
            best_n = score;
            best_date = D(k).datenum;
        end
    end

    if isempty(best_idx)
        error('ab_finalize_evidence:NoUsableResultFolder', ...
            'No usable result folder containing %s was found.',required_csv);
    end

    out_dir = fullfile(D(best_idx).folder,D(best_idx).name);
    T = readtable(fullfile(out_dir,required_csv),'TextType','string');
    if ismember('Success',T.Properties.VariableNames)
        n = numel(unique(T.Trial(T.Success==1)));
    else
        n = numel(unique(T.Trial));
    end
    if require_exact && n ~= expected_runs
        error('ab_finalize_evidence:NoExactRunFolder', ...
            'Best available %s folder has %d successful trials; expected %d.', ...
            prefix,n,expected_runs);
    end
end


function ab_fin_validate_required_file(folder,name)
    if ~exist(folder,'dir')
        error('ab_finalize_evidence:MissingFolder','Folder not found: %s',folder);
    end
    if ~exist(fullfile(folder,name),'file')
        error('ab_finalize_evidence:MissingFile', ...
            'Required file not found: %s',fullfile(folder,name));
    end
end


function ab_fin_assert_table_vars(T,names,label)
    missing = names(~ismember(names,T.Properties.VariableNames));
    if ~isempty(missing)
        error('ab_finalize_evidence:MissingVariables', ...
            '%s is missing variables: %s',label,strjoin(missing,', '));
    end
end

%% ========================================================================

function Tout = ab_fin_build_core_mc_table(T)
    variants = ["Full Proposed","No Mismatch Compensation","Simplified Adaptation (B/K)"];
    joints = ["Hip","Knee","Ankle"];
    metric_cols = {'RMSE_deg','JerkRMS_deg_s3','RawResidualRMS_Nm'};

    rows = cell(0,11);
    for v = 1:numel(variants)
        for j = 1:numel(joints)
            idx = T.Success==1 & T.Variant==variants(v) & T.Joint==joints(j);
            n = sum(idx);
            vals = cell(1,6);
            c = 1;
            for m = 1:numel(metric_cols)
                x = T.(metric_cols{m})(idx);
                x = x(isfinite(x));
                vals{c} = mean(x); c=c+1;
                vals{c} = std(x,0); c=c+1;
            end
            rows(end+1,:) = {variants(v),joints(j),n, ...
                vals{1},vals{2},vals{3},vals{4},vals{5},vals{6}, ...
                mean(T.PositionSaturation_pct(idx),'omitnan'), ...
                mean(T.ControllerInputRMS_Nm(idx),'omitnan')}; %#ok<AGROW>
        end
    end

    Tout = cell2table(rows,'VariableNames',{ ...
        'Variant','Joint','N', ...
        'RMSE_Mean_deg','RMSE_Std_deg', ...
        'JerkRMS_Mean_deg_s3','JerkRMS_Std_deg_s3', ...
        'RawResidualRMS_Mean_Nm','RawResidualRMS_Std_Nm', ...
        'PositionSaturation_Mean_pct','ControllerInputRMS_Mean_Nm'});
end

%% ========================================================================

function Tout = ab_fin_build_pairwise_statistics(T)
    full_name = "Full Proposed";
    comparators = ["No Mismatch Compensation","Simplified Adaptation (B/K)"];
    joints = ["Hip","Knee","Ankle"];
    metrics = struct( ...
        'label',{'RMSE','JerkRMS','RawResidualRMS'}, ...
        'column',{'RMSE_deg','JerkRMS_deg_s3','RawResidualRMS_Nm'}, ...
        'unit',{'deg','deg/s^3','N m'});

    rows = cell(0,18);
    for c = 1:numel(comparators)
        for j = 1:numel(joints)
            F = T(T.Success==1 & T.Variant==full_name & T.Joint==joints(j),:);
            C = T(T.Success==1 & T.Variant==comparators(c) & T.Joint==joints(j),:);
            common = intersect(F.Trial,C.Trial);

            for m = 1:numel(metrics)
                xf = nan(numel(common),1);
                xc = nan(numel(common),1);
                for k = 1:numel(common)
                    rf = F.Trial==common(k);
                    rc = C.Trial==common(k);
                    xf(k) = F.(metrics(m).column)(find(rf,1,'first'));
                    xc(k) = C.(metrics(m).column)(find(rc,1,'first'));
                end
                good = isfinite(xf) & isfinite(xc);
                xf = xf(good); xc = xc(good);
                d = xc-xf; % comparator minus Full; positive means Full is lower/better
                n = numel(d);
                if n < 2
                    error('ab_finalize_evidence:TooFewPairs', ...
                        'Too few paired observations for %s / %s / %s.', ...
                        comparators(c),joints(j),metrics(m).label);
                end

                mu = mean(d);
                sd = std(d,0);
                se = sd/sqrt(n);
                [tstat,p] = ab_fin_paired_t_pvalue(d);
                tcrit = ab_fin_tcrit975(n-1);
                ci_lo = mu-tcrit*se;
                ci_hi = mu+tcrit*se;
                if sd>0, dz=mu/sd; else dz=sign(mu)*Inf; end
                full_better_pct = 100*mean(d>0);
                full_worse_pct = 100*mean(d<0);
                tie_pct = 100*mean(d==0);
                denom = max(abs(xc),eps);
                mean_improvement_pct = mean(100*d./denom);
                full_mean = mean(xf);
                comp_mean = mean(xc);
                if mu>0
                    direction = "Full lower";
                elseif mu<0
                    direction = "Comparator lower";
                else
                    direction = "Equal mean";
                end

                rows(end+1,:) = {comparators(c),joints(j),metrics(m).label,metrics(m).unit,n, ...
                    full_mean,comp_mean,mu,ci_lo,ci_hi,tstat,p,dz, ...
                    full_better_pct,full_worse_pct,tie_pct,mean_improvement_pct,direction}; %#ok<AGROW>
            end
        end
    end

    Tout = cell2table(rows,'VariableNames',{ ...
        'Comparator','Joint','Metric','Unit','N', ...
        'FullMean','ComparatorMean','ComparatorMinusFull_Mean', ...
        'Difference95CI_Low','Difference95CI_High','PairedT','P_Raw','CohensDz', ...
        'FullBetterTrialPct','FullWorseTrialPct','TieTrialPct', ...
        'FullImprovementPct_Mean','MeanDirection'});
end


function [tstat,p] = ab_fin_paired_t_pvalue(d)
    d = d(isfinite(d));
    n = numel(d);
    mu = mean(d);
    sd = std(d,0);
    if sd==0
        if mu==0
            tstat=0; p=1;
        else
            tstat=sign(mu)*Inf; p=0;
        end
        return;
    end
    tstat = mu/(sd/sqrt(n));
    v = n-1;
    x = v/(v+tstat^2);
    % Two-sided Student-t p-value. For |t|, p = I_x(v/2,1/2).
    p = betainc(x,v/2,0.5);
    p = min(max(p,0),1);
end


function tc = ab_fin_tcrit975(v)
    % 97.5th percentile of Student-t(v), without Statistics Toolbox.
    % For t>0, two-sided alpha=.05 -> I_x(v/2,1/2)=.05.
    x = betaincinv(0.05,v/2,0.5);
    tc = sqrt(v*(1/x-1));
end


function padj = ab_fin_holm_adjust(p)
    p = p(:);
    m = numel(p);
    [ps,ord] = sort(p,'ascend');
    adj_sorted = zeros(m,1);
    running = 0;
    for k = 1:m
        val = (m-k+1)*ps(k);
        running = max(running,val);
        adj_sorted(k) = min(1,running);
    end
    padj = zeros(m,1);
    padj(ord) = adj_sorted;
end

%% ========================================================================

function Tout = ab_fin_build_nominal_compact(T)
    wanted = {'Variant','Joint','K_mis','Adapt_B','Adapt_K','Adapt_G', ...
        'RMSE_deg','JerkRMS_deg_s3','RawResidualRMS_Nm', ...
        'ControllerInputRMS_Nm','PositionSaturation_pct','MismatchLeakageRMS_Nm'};
    wanted = wanted(ismember(wanted,T.Properties.VariableNames));
    Tout = T(:,wanted);
end


function Tout = ab_fin_build_mismatch_mechanism_table(T)
    wanted = {'K_mis','Joint','N','RMSE_Mean','RMSE_Std', ...
        'JerkRMS_Mean','JerkRMS_Std','RawResidualRMS_Mean','RawResidualRMS_Std', ...
        'ControllerInputRMS_Mean','MismatchRMS_Mean','PathologicalRMS_Mean', ...
        'MismatchLeakageRMS_Mean','CorrMismatchPathological_Mean', ...
        'MeanProductMismatchPathological_Mean','CancellationRatio_Mean', ...
        'EventRMSE_Mean','EventJerkRMS_Mean','EventRawResidualRMS_Mean', ...
        'EventCorrMismatchPathological_Mean','PositionSaturationPct_Mean'};
    wanted = wanted(ismember(wanted,T.Properties.VariableNames));
    Tout = T(:,wanted);
    Tout = sortrows(Tout,{'K_mis','Joint'});
end


function Tout = ab_fin_build_mismatch_endpoint_table(T)
    idx = abs(T.K_mis-0)<1e-12 | abs(T.K_mis-1)<1e-12;
    Tout = ab_fin_build_mismatch_mechanism_table(T(idx,:));
end


function Tout = ab_fin_build_mismatch_nominal_endpoint_table(T)
    idx = abs(T.K_mis-0)<1e-12 | abs(T.K_mis-1)<1e-12;
    wanted = {'K_mis','Joint','RMSE_deg','JerkRMS_deg_s3','RawResidualRMS_Nm', ...
        'ControllerInputRMS_Nm','MismatchRMS_Nm','PathologicalRMS_Nm', ...
        'MismatchLeakageRMS_Nm','Corr_Mismatch_Pathological', ...
        'MeanProduct_Mismatch_Pathological_Nm2','CancellationRatio', ...
        'Event_RMSE_deg','Event_JerkRMS_deg_s3','Event_RawResidualRMS_Nm', ...
        'Event_Corr_Mismatch_Pathological','PositionSaturation_pct'};
    wanted = wanted(ismember(wanted,T.Properties.VariableNames));
    Tout = T(idx,wanted);
    Tout = sortrows(Tout,{'K_mis','Joint'});
end

%% ========================================================================

function ab_fin_write_readme(out_dir,n)
    f = fopen(fullfile(out_dir,'README_FINAL_EVIDENCE.txt'),'w');
    if f<0, return; end
    C = onCleanup(@()fclose(f)); %#ok<NASGU>
    fprintf(f,'FINAL CORE-ABLATION EVIDENCE PACKAGE\n');
    fprintf(f,'====================================\n\n');
    fprintf(f,'This folder aggregates the completed simulations; it does not rerun or retune controllers.\n\n');
    fprintf(f,'Successful strictly paired Monte Carlo trials used: %d\n\n',n);
    fprintf(f,'Primary comparison groups:\n');
    fprintf(f,'  1. Full Proposed: K_mis=1; B_a, K_a, G_a adaptive.\n');
    fprintf(f,'  2. No Mismatch Compensation: K_mis=0; B_a, K_a, G_a adaptive.\n');
    fprintf(f,'  3. Simplified Adaptation (B/K): K_mis=1; B_a and K_a adaptive, G_a fixed.\n\n');
    fprintf(f,'Primary metrics (lower is better): tracking RMSE, jerk RMS, raw interaction-residual RMS.\n');
    fprintf(f,'Pairwise difference convention: Comparator - Full. Positive values mean Full is lower.\n');
    fprintf(f,'Statistics: two-sided paired t test, t-based 95%% CI, Cohen d_z, and Holm adjustment across 18 endpoint tests.\n\n');
    fprintf(f,'Mismatch K_mis={0,0.25,0.5,0.75,1} data are retained only as mechanism/sensitivity evidence.\n');
    fprintf(f,'They are not used to tune K_mis or to redefine the Full Proposed method.\n\n');
    fprintf(f,'Key files:\n');
    fprintf(f,'  TableA_core_ablation_mc_mean_std.csv\n');
    fprintf(f,'  TableB_core_ablation_pairwise_statistics.csv\n');
    fprintf(f,'  TableC_core_ablation_nominal.csv\n');
    fprintf(f,'  TableD_mismatch_sensitivity_mechanism.csv\n');
    fprintf(f,'  TableE_mismatch_endpoint_mechanism.csv\n');
    fprintf(f,'  TableF_mismatch_nominal_endpoint_mechanism.csv\n');
    fprintf(f,'  figures/FigA_core_ablation_mc.*\n');
    fprintf(f,'  figures/FigB_mismatch_tradeoff_and_leakage.*\n');
    fprintf(f,'  figures/FigC_pairwise_effect_forest.*\n');
    fprintf(f,'  figure_data/FigA_core_ablation_mc.csv/.mat\n');
    fprintf(f,'  figure_data/FigB_mismatch_tradeoff_and_leakage.csv/.mat\n');
    fprintf(f,'  figure_data/FigC_pairwise_effect_forest.csv/.mat\n');
    fprintf(f,'\nEvery ablation PDF/PNG/JPEG has an identically named CSV and MAT figure-data snapshot.\n');
end


function ab_final_generate_figures(out_dir,Ttrial,Tstats,Tmis)
%GENERATE_FINAL_ABLATION_FIGURES Paper-facing figures from completed runs.
% Every exported PDF/PNG receives an identically named CSV and MAT in
% figure_data/ containing the values actually plotted.

    fig_dir=fullfile(out_dir,'figures');
    dat_dir=fullfile(out_dir,'figure_data');
    if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
    if ~exist(dat_dir,'dir'), mkdir(dat_dir); end
    variants=["Full Proposed","No Mismatch Compensation","Simplified Adaptation (B/K)"];
    joints=["Hip","Knee","Ankle"];

    %% FigA: final primary core-ablation figure, MC mean +/- standard deviation
    base_name='FigA_core_ablation_mc';
    metrics=struct('col',{'RMSE_deg','JerkRMS_deg_s3','RawResidualRMS_Nm'}, ...
        'title',{'Tracking RMSE','Jerk RMS','Raw interaction-residual RMS'}, ...
        'ylab',{'RMSE (deg)','Jerk RMS (deg/s^3)','Raw residual RMS (N m)'});
    ntrials=numel(unique(Ttrial.Trial(Ttrial.Success==1)));
    [f1,lay1]=paper_figure_utils('create','Final Core Ablation', ...
        [70 45 1180 1180],3,1,sprintf('Core-Module Ablation (%d Paired Monte Carlo Trials)',ntrials));
    axA=gobjects(3,1);
    nr=3*numel(joints)*numel(variants); Metric=strings(nr,1); Joint=strings(nr,1); Variant=strings(nr,1); Mean=nan(nr,1); SD=nan(nr,1); N=zeros(nr,1); row=0;
    for m=1:3
        ax=paper_figure_utils('tile',lay1,m); axA(m)=ax; hold(ax,'on');
        mu=nan(numel(joints),numel(variants)); sd=mu;
        for j=1:numel(joints)
            for v=1:numel(variants)
                idx=Ttrial.Success==1 & Ttrial.Joint==joints(j) & Ttrial.Variant==variants(v);
                x=Ttrial.(metrics(m).col)(idx); x=x(isfinite(x));
                mu(j,v)=mean(x); sd(j,v)=std(x,0);
                row=row+1; Metric(row)=string(metrics(m).title); Joint(row)=joints(j); Variant(row)=variants(v);
                Mean(row)=mu(j,v); SD(row)=sd(j,v); N(row)=numel(x);
            end
        end
        b=bar(ax,1:numel(joints),mu,'grouped');
        for v=1:numel(b)
            try, x=b(v).XEndPoints; catch
                groupwidth=min(0.8,numel(variants)/(numel(variants)+1.5));
                x=(1:numel(joints))-groupwidth/2+(2*v-1)*groupwidth/(2*numel(variants));
            end
            errorbar(ax,x,mu(:,v),sd(:,v),'k.','LineWidth',1.0,'CapSize',5,'HandleVisibility','off');
        end
        xticks(ax,1:numel(joints)); xticklabels(ax,cellstr(joints)); ylabel(ax,metrics(m).ylab);
        if m==3, xlabel(ax,'Joint'); end
        paper_figure_utils('panel_title',ax,sprintf('(%c) %s (Mean +/- SD)',char('a'+m-1),metrics(m).title));
    end
    lgd = gobjects(numel(axA),1);
    for ii = 1:numel(axA)
        lgd(ii)=legend(axA(ii),cellstr(variants),'Location','northeast','Interpreter','none','AutoUpdate','off');
    end
    paper_figure_utils('align_legends',axA,lgd,'NumColumns',2);
    Tfig=table(Metric,Joint,Variant,Mean,SD,N);
    meta=struct('figure',base_name,'description','Final core-ablation paired Monte Carlo mean and standard deviation','successful_trials',ntrials);
    ab_save_figure_data(dat_dir,base_name,Tfig,meta);
    ab_fin_export_pair(f1,fig_dir,base_name);

    %% FigB: compact mismatch mechanism panel
    base_name='FigB_mismatch_tradeoff_and_leakage';
    kmis=sort(unique(Tmis.K_mis)); kmis=kmis(:);
    [f2,lay2]=paper_figure_utils('create','Mismatch Tradeoff and Leakage', ...
        [100 90 1250 560],1,2,'Mismatch-Compensation Mechanism and Sensitivity');

    ax1=paper_figure_utils('tile',lay2,1); hold(ax1,'on');
    Tk=Tmis(Tmis.Joint=="Knee",:); Tk=sortrows(Tk,'K_mis');
    yyaxis(ax1,'left');
    p1=plot(ax1,Tk.K_mis,Tk.RMSE_Mean,'-o','LineWidth',1.5,'MarkerSize',6,'DisplayName','Tracking RMSE');
    ylabel(ax1,'Knee RMSE (deg)');
    yyaxis(ax1,'right');
    p2=plot(ax1,Tk.K_mis,Tk.RawResidualRMS_Mean,'-s','LineWidth',1.5,'MarkerSize',6,'DisplayName','Raw residual RMS');
    ylabel(ax1,'Knee raw residual RMS (N m)');
    xticks(ax1,kmis); xlabel(ax1,'K_{mis}');
    paper_figure_utils('panel_title',ax1,'(a) Knee Tracking-Interaction Trade-Off');
    lgd1=legend(ax1,[p1 p2],{'Tracking RMSE','Raw residual RMS'},'Location','northeast','AutoUpdate','off');

    ax2=paper_figure_utils('tile',lay2,2); hold(ax2,'on');
    markers={'o','s','^'}; leak=nan(numel(kmis),numel(joints));
    for j=1:numel(joints)
        Tj=Tmis(Tmis.Joint==joints(j),:); Tj=sortrows(Tj,'K_mis');
        semilogy(ax2,Tj.K_mis,max(Tj.MismatchLeakageRMS_Mean,eps),['-' markers{j}], ...
            'LineWidth',1.5,'MarkerSize',6,'DisplayName',joints(j));
        [~,loc]=ismember(Tj.K_mis,kmis); leak(loc,j)=Tj.MismatchLeakageRMS_Mean;
    end
    xticks(ax2,kmis); xlabel(ax2,'K_{mis}'); ylabel(ax2,'Mismatch leakage RMS (N m)');
    paper_figure_utils('panel_title',ax2,'(b) Nominal-Mismatch Leakage');
    lgd2=legend(ax2,'Location','northeast','AutoUpdate','off');
    paper_figure_utils('align_legends',[ax1 ax2],[lgd1 lgd2]);

    K_mis=kmis; Knee_RMSE_deg=nan(size(kmis)); Knee_RawResidualRMS_Nm=nan(size(kmis));
    [~,loc]=ismember(Tk.K_mis,kmis); Knee_RMSE_deg(loc)=Tk.RMSE_Mean; Knee_RawResidualRMS_Nm(loc)=Tk.RawResidualRMS_Mean;
    Hip_MismatchLeakageRMS_Nm=leak(:,1); Knee_MismatchLeakageRMS_Nm=leak(:,2); Ankle_MismatchLeakageRMS_Nm=leak(:,3);
    Tfig=table(K_mis,Knee_RMSE_deg,Knee_RawResidualRMS_Nm,Hip_MismatchLeakageRMS_Nm,Knee_MismatchLeakageRMS_Nm,Ankle_MismatchLeakageRMS_Nm);
    meta=struct('figure',base_name,'description','Knee tracking/raw-residual trade-off and three-joint mismatch leakage versus K_mis');
    ab_save_figure_data(dat_dir,base_name,Tfig,meta);
    ab_fin_export_pair(f2,fig_dir,base_name);

    %% FigC: paired effects (Comparator - Full), 95% CI
    base_name='FigC_pairwise_effect_forest';
    metric_names=["RMSE","JerkRMS","RawResidualRMS"];
    comparator_names=["No Mismatch Compensation","Simplified Adaptation (B/K)"];
    [f3,lay3]=paper_figure_utils('create','Paired Ablation Effects', ...
        [130 45 1180 1180],3,1,'Paired Core-Ablation Effects: Comparator - Full Proposed');
    Tfig=Tstats(ismember(Tstats.Metric,metric_names),:);
    for m=1:numel(metric_names)
        ax=paper_figure_utils('tile',lay3,m); hold(ax,'on');
        S=Tstats(Tstats.Metric==metric_names(m),:); y=1:height(S);
        for r=1:height(S)
            x=S.ComparatorMinusFull_Mean(r); lo=S.Difference95CI_Low(r); hi=S.Difference95CI_High(r);
            errorbar(ax,x,y(r),x-lo,hi-x,'horizontal','o','LineWidth',1.2,'CapSize',5,'MarkerSize',5);
        end
        xline(ax,0,'--','LineWidth',1.0);
        labels=strings(height(S),1);
        for r=1:height(S)
            shortcomp=S.Comparator(r);
            if shortcomp==comparator_names(1), shortcomp="No mismatch"; else, shortcomp="B/K only"; end
            labels(r)=S.Joint(r)+" - "+shortcomp;
        end
        yticks(ax,y); yticklabels(ax,cellstr(labels)); set(ax,'YDir','reverse');
        xlabel(ax,['Comparator - Full (' char(S.Unit(1)) ')']);
        paper_figure_utils('panel_title',ax,sprintf('(%c) %s',char('a'+m-1),char(metric_names(m))));
    end
    meta=struct('figure',base_name,'description','Paired comparator-minus-Full effects with 95% confidence intervals; positive means Full is lower');
    ab_save_figure_data(dat_dir,base_name,Tfig,meta);
    ab_fin_export_pair(f3,fig_dir,base_name);

    % Keep figures visible after export, consistent with generate_paper_figures.m.
end

function ab_save_figure_data(dat_dir,base_name,Tfig,metadata)
%AB_SAVE_FIGURE_DATA Save one-to-one plot data for an exported ablation figure.
% CSV contains the values actually plotted. MAT contains the same table plus
% descriptive metadata, preserving a direct PDF/PNG/JPEG <-> CSV/MAT mapping.
    if ~exist(dat_dir,'dir'), mkdir(dat_dir); end
    if nargin < 4 || isempty(metadata), metadata=struct(); end
    metadata.base_name=char(base_name);
    metadata.csv_file=[char(base_name) '.csv'];
    metadata.mat_file=[char(base_name) '.mat'];
    metadata.note='CSV/MAT are figure-level data snapshots generated from the same values used for plotting.';
    csv_path=fullfile(dat_dir,[char(base_name) '.csv']);
    mat_path=fullfile(dat_dir,[char(base_name) '.mat']);
    writetable(Tfig,csv_path);
    figure_data=Tfig; %#ok<NASGU>
    figure_metadata=metadata; %#ok<NASGU>
    save(mat_path,'figure_data','figure_metadata','-v7');
end

function ab_fin_export_pair(fig,fig_dir,base_name)
    paper_figure_utils('export',fig,fig_dir,base_name, ...
        'SavePNG',true,'SaveJPEG',true,'JPEGResolution',1000);
end
