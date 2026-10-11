function monte_carlo_runner(num_runs, base_params, reference_bundle)
%MONTE_CARLO_RUNNER Four-controller strictly paired Monte Carlo workflow.
% Controllers:
%   Fixed | Zhang-2018 literature-derived VA | Cao-2024-derived VA | Proposed IE-AAC
% Primary metrics:
%   Tracking RMSE | Jerk RMS | Raw Interaction Residual RMS
% Supervisory diagnostic:
%   Position-correction saturation at the common Delta-q limit.

    if nargin < 1 || isempty(num_runs), num_runs = 100; end
    fprintf('=== 开始四控制器 Monte-Carlo 仿真（%d 组） ===\n', num_runs);
    addpath(genpath(pwd));

    if nargin < 2 || isempty(base_params), base_params = parameters(); end
    assert(isfield(base_params,'literature_zhang2018') && logical(base_params.literature_zhang2018.frozen), ...
        'monte_carlo_runner:ZhangNotFrozen', 'Zhang-2018 literature comparator 必须在正式 MC 前冻结。');
    assert(isfield(base_params,'literature_cao2024') && logical(base_params.literature_cao2024.frozen), ...
        'monte_carlo_runner:CaoNotFrozen', 'Cao-2024-derived literature comparator 必须在正式 MC 前冻结。');

    %% ================= 统一参考轨迹 =================
    if nargin >= 3 && ~isempty(reference_bundle)
        required_fields = {'q_ref','dq_ref','ddq_ref','t_all','metrics_teaching'};
        for k = 1:numel(required_fields)
            assert(isfield(reference_bundle, required_fields{k}), ...
                'reference_bundle 缺少字段 %s', required_fields{k});
        end
        q_teach = reference_bundle.q_ref;
        dq_teach = reference_bundle.dq_ref;
        ddq_teach = reference_bundle.ddq_ref;
        t_all = reference_bundle.t_all;
        metrics_teaching = reference_bundle.metrics_teaching;
        reference_source = 'shared_nominal_reference';
    else
        [base_params, metrics_teaching, q_teach, dq_teach, ddq_teach, t_all] = ...
            TeachingDataPipeline(base_params, 'Seed', base_params.seed, ...
                'PlotFlag', false, 'SaveFlag', false, 'Verbose', false);
        reference_source = 'standalone_same_teaching_pipeline';
    end

    if isfield(metrics_teaching,'q_controller_ref')
        assert(max(abs(q_teach(:)-metrics_teaching.q_controller_ref(:))) < 1e-14, 'MC q_ref 不一致');
        assert(max(abs(dq_teach(:)-metrics_teaching.dq_controller_ref(:))) < 1e-14, 'MC dq_ref 不一致');
        assert(max(abs(ddq_teach(:)-metrics_teaching.ddq_controller_ref(:))) < 1e-14, 'MC ddq_ref 不一致');
        assert(max(abs(t_all(:)-metrics_teaching.t_controller_ref(:))) < 1e-14, 'MC 时间轴不一致');
    end
    fprintf('✅ MC reference 来源：%s\n', reference_source);

    %% ================= 预分配 =================
    stats = struct();
    stats.fixed = init_metric_block(num_runs);
    stats.zhang2018 = init_metric_block(num_runs);
    stats.cao2024 = init_metric_block(num_runs);
    stats.adaptive = init_metric_block(num_runs);
    stats.success_mask = false(num_runs,1);
    stats.parameter_seed = nan(num_runs,1);
    stats.noise_seed = nan(num_runs,1);
    stats.failure_message = strings(num_runs,1);
    stats.trial_audit = cell(num_runs,1);
    stats.fail_reasons = {};
    stats.reference_source = reference_source;
    stats.reference_seed = metrics_teaching.teaching_seed;
    stats.reference_num_keypoints = metrics_teaching.num_keypoints;
    stats.reference_rmse_deg = metrics_teaching.rmse_global;
    stats.reference_max_key_gap_s = metrics_teaching.max_keypoint_gap_s;

    results_dir = 'results';
    if ~exist(results_dir,'dir'), mkdir(results_dir); end
    mc_timestamp = datestr(now,'yyyymmdd_HHMMSS');
    mc_folder = fullfile(results_dir,['MC_' mc_timestamp]);
    mkdir(mc_folder);
    fprintf('📁 Monte-Carlo 结果保存至：%s\n', mc_folder);

    %% ================= 严格配对主循环 =================
    for r = 1:num_runs
        try
            fprintf('运行 %d/%d ...\n', r, num_runs);
            [params, trial_audit] = mc_trial_utils('sample', base_params, r);
            params.config_date = string(datetime('now'),'yyyyMMdd_HHmmss');
            stats.parameter_seed(r) = trial_audit.parameter_seed;
            stats.noise_seed(r) = trial_audit.noise_seed;
            stats.trial_audit{r} = trial_audit;

            [pre_fb, event_timing] = mc_trial_utils('pre_fb', params, q_teach, t_all, trial_audit.noise_seed);
            params.t_abnormal = event_timing.abnormal_window;

            % Fixed
            [dq_q, dq_v, dq_a, ~, ~, ~, qflag, euse, eraw, diagc] = ...
                run_admittance_controller('fixed', q_teach, dq_teach, params, pre_fb, t_all, ddq_teach);
            qcmd = q_teach + dq_q; dcmd = dq_teach + dq_v; acmd = ddq_teach + dq_a;
            m_fixed = compute_metrics(q_teach,qcmd,dcmd,acmd,eraw,euse,params,diagc);
            validate_trial_outputs(qcmd,dcmd,acmd,eraw,m_fixed,'Fixed');
            stats.fixed = store_metric_block(stats.fixed,r,m_fixed,qflag);

            % Zhang-2018
            [dq_q, dq_v, dq_a, ~, ~, ~, qflag, euse, eraw, diagc] = ...
                run_admittance_controller('zhang2018', q_teach, dq_teach, params, pre_fb, t_all, ddq_teach);
            qcmd = q_teach + dq_q; dcmd = dq_teach + dq_v; acmd = ddq_teach + dq_a;
            m_zhang = compute_metrics(q_teach,qcmd,dcmd,acmd,eraw,euse,params,diagc);
            validate_trial_outputs(qcmd,dcmd,acmd,eraw,m_zhang,'Zhang-2018');
            stats.zhang2018 = store_metric_block(stats.zhang2018,r,m_zhang,qflag);

            % Cao-2024-derived
            [dq_q, dq_v, dq_a, ~, ~, ~, qflag, euse, eraw, diagc] = ...
                run_admittance_controller('cao2024', q_teach, dq_teach, params, pre_fb, t_all, ddq_teach);
            qcmd = q_teach + dq_q; dcmd = dq_teach + dq_v; acmd = ddq_teach + dq_a;
            m_cao = compute_metrics(q_teach,qcmd,dcmd,acmd,eraw,euse,params,diagc);
            validate_trial_outputs(qcmd,dcmd,acmd,eraw,m_cao,'Cao-2024-derived');
            stats.cao2024 = store_metric_block(stats.cao2024,r,m_cao,qflag);

            % Proposed
            [dq_q, dq_v, dq_a, ~, ~, ~, qflag, euse, eraw, diagc] = ...
                run_admittance_controller('proposed', q_teach, dq_teach, params, pre_fb, t_all, ddq_teach);
            qcmd = q_teach + dq_q; dcmd = dq_teach + dq_v; acmd = ddq_teach + dq_a;
            m_adp = compute_metrics(q_teach,qcmd,dcmd,acmd,eraw,euse,params,diagc);
            validate_trial_outputs(qcmd,dcmd,acmd,eraw,m_adp,'Proposed');
            stats.adaptive = store_metric_block(stats.adaptive,r,m_adp,qflag);

            stats.success_mask(r) = true;
            fprintf('✅ 第 %d 次四控制器配对运行成功\n', r);
        catch ME
            stats.failure_message(r) = string(ME.message);
            stats.fail_reasons{end+1} = sprintf('Run %d: %s',r,ME.message); %#ok<AGROW>
            stats.success_mask(r) = false;
            fprintf('❌ 第 %d 次运行出错（该 trial 从四控制器配对统计中排除）：%s\n',r,ME.message);
        end
    end

    %% ================= 有效样本审计 =================
    valid_idx = stats.success_mask(:);
    stats.success_count = sum(valid_idx);
    stats.requested_count = num_runs;
    stats.fail_count = num_runs - stats.success_count;
    stats.success_rate = stats.success_count / num_runs;

    fprintf('\n=== Monte Carlo 有效样本审计 ===\n');
    fprintf('请求 trial: %d | 成功: %d | 失败: %d | 成功率: %.2f%%\n', ...
        num_runs,stats.success_count,stats.fail_count,100*stats.success_rate);
    if stats.success_count == 0
        error('MonteCarlo:NoValidTrials','全部 trial 均失败，禁止生成统计结果。');
    end
    if stats.success_rate < base_params.mc.min_success_rate
        warning('MonteCarlo:LowSuccessRate', ...
            'Monte Carlo 成功率 %.2f%% 低于审计阈值 %.2f%%。', ...
            100*stats.success_rate,100*base_params.mc.min_success_rate);
    end

    %% ================= 参考与 trial 审计保存 =================
    reference_metadata = struct();
    reference_metadata.source = reference_source;
    reference_metadata.teaching_seed = metrics_teaching.teaching_seed;
    reference_metadata.num_keypoints = metrics_teaching.num_keypoints;
    reference_metadata.reconstruction_rmse_deg = metrics_teaching.rmse_global;
    reference_metadata.max_key_gap_s = metrics_teaching.max_keypoint_gap_s;
    reference_metadata.pipeline = metrics_teaching.reference_pipeline;
    reference_metadata.pipeline_version = metrics_teaching.reference_pipeline_version;
    save(fullfile(mc_folder,['MC_reference_' mc_timestamp '.mat']), ...
        'q_teach','dq_teach','ddq_teach','t_all','metrics_teaching','reference_metadata','-v7.3');

    T_reference = table(string(reference_source),metrics_teaching.teaching_seed, ...
        metrics_teaching.num_keypoints,metrics_teaching.rmse_global,metrics_teaching.max_keypoint_gap_s, ...
        'VariableNames',{'ReferenceSource','TeachingSeed','NumKeypoints','ReconstructionRMSE_deg','MaxKeyGap_s'});
    writetable(T_reference,fullfile(mc_folder,['MC_reference_audit_' mc_timestamp '.csv']));

    T_ranges = build_mc_uncertainty_table(base_params);
    writetable(T_ranges,fullfile(mc_folder,['MC_uncertainty_ranges_' mc_timestamp '.csv']));
    T_trials = table((1:num_runs)',stats.success_mask,stats.parameter_seed,stats.noise_seed,stats.failure_message, ...
        'VariableNames',{'Run','Success','ParameterSeed','NoiseSeed','FailureMessage'});
    writetable(T_trials,fullfile(mc_folder,['MC_trial_audit_' mc_timestamp '.csv']));
    save(fullfile(mc_folder,['MC_stats_' mc_timestamp '.mat']),'stats','num_runs','base_params','reference_metadata','T_ranges','T_trials');

    %% ================= 主统计表 =================
    Joints = {'Hip';'Knee';'Ankle'};
    F = extract_valid(stats.fixed,valid_idx);
    Z = extract_valid(stats.zhang2018,valid_idx);
    C = extract_valid(stats.cao2024,valid_idx);
    A = extract_valid(stats.adaptive,valid_idx);
    N_success = repmat(stats.success_count,3,1);
    SuccessRate_pct = repmat(100*stats.success_rate,3,1);

    T = table(Joints,N_success,SuccessRate_pct, ...
        mean(F.RMSE,1)',std(F.RMSE,0,1)', mean(Z.RMSE,1)',std(Z.RMSE,0,1)', mean(C.RMSE,1)',std(C.RMSE,0,1)', mean(A.RMSE,1)',std(A.RMSE,0,1)', ...
        mean(F.JerkRMS,1)'*180/pi,std(F.JerkRMS,0,1)'*180/pi, mean(Z.JerkRMS,1)'*180/pi,std(Z.JerkRMS,0,1)'*180/pi, mean(C.JerkRMS,1)'*180/pi,std(C.JerkRMS,0,1)'*180/pi, mean(A.JerkRMS,1)'*180/pi,std(A.JerkRMS,0,1)'*180/pi, ...
        mean(F.RMS_InteractionTorque,1)',std(F.RMS_InteractionTorque,0,1)', mean(Z.RMS_InteractionTorque,1)',std(Z.RMS_InteractionTorque,0,1)', mean(C.RMS_InteractionTorque,1)',std(C.RMS_InteractionTorque,0,1)', mean(A.RMS_InteractionTorque,1)',std(A.RMS_InteractionTorque,0,1)', ...
        'VariableNames',{'Joint','N_Success','SuccessRate_pct', ...
        'Fixed_RMSE_mean','Fixed_RMSE_std','Zhang2018_RMSE_mean','Zhang2018_RMSE_std','Cao2024_RMSE_mean','Cao2024_RMSE_std','Proposed_RMSE_mean','Proposed_RMSE_std', ...
        'Fixed_Jerk_deg_s3_mean','Fixed_Jerk_deg_s3_std','Zhang2018_Jerk_deg_s3_mean','Zhang2018_Jerk_deg_s3_std','Cao2024_Jerk_deg_s3_mean','Cao2024_Jerk_deg_s3_std','Proposed_Jerk_deg_s3_mean','Proposed_Jerk_deg_s3_std', ...
        'Fixed_RawResidualRMS_mean','Fixed_RawResidualRMS_std','Zhang2018_RawResidualRMS_mean','Zhang2018_RawResidualRMS_std','Cao2024_RawResidualRMS_mean','Cao2024_RawResidualRMS_std','Proposed_RawResidualRMS_mean','Proposed_RawResidualRMS_std'});
    disp(T);
    writetable(T,fullfile(mc_folder,['MC_summary_primary_four_controller_' mc_timestamp '.csv']));

    T_diag = table(Joints, ...
        mean(F.PathologicalRMS,1)',mean(Z.PathologicalRMS,1)',mean(C.PathologicalRMS,1)',mean(A.PathologicalRMS,1)', ...
        mean(F.ControllerInputRMS,1)',mean(Z.ControllerInputRMS,1)',mean(C.ControllerInputRMS,1)',mean(A.ControllerInputRMS,1)', ...
        mean(F.InteractionEnergy,1)',mean(Z.InteractionEnergy,1)',mean(C.InteractionEnergy,1)',mean(A.InteractionEnergy,1)', ...
        mean(F.PeakInteractionTorque,1)',mean(Z.PeakInteractionTorque,1)',mean(C.PeakInteractionTorque,1)',mean(A.PeakInteractionTorque,1)', ...
        mean(F.MaxError,1)',std(F.MaxError,0,1)',mean(Z.MaxError,1)',std(Z.MaxError,0,1)',mean(C.MaxError,1)',std(C.MaxError,0,1)',mean(A.MaxError,1)',std(A.MaxError,0,1)', ...
        mean(F.PositionSaturationPct,1)',std(F.PositionSaturationPct,0,1)',mean(Z.PositionSaturationPct,1)',std(Z.PositionSaturationPct,0,1)',mean(C.PositionSaturationPct,1)',std(C.PositionSaturationPct,0,1)',mean(A.PositionSaturationPct,1)',std(A.PositionSaturationPct,0,1)', ...
        mean(F.PositionMaxRatio,1)',mean(Z.PositionMaxRatio,1)',mean(C.PositionMaxRatio,1)',mean(A.PositionMaxRatio,1)', ...
        'VariableNames',{'Joint','Fixed_PathologicalRMS','Zhang2018_PathologicalRMS','Cao2024_PathologicalRMS','Proposed_PathologicalRMS', ...
        'Fixed_ControllerInputRMS','Zhang2018_ControllerInputRMS','Cao2024_ControllerInputRMS','Proposed_ControllerInputRMS', ...
        'Fixed_RawResidualEnergy','Zhang2018_RawResidualEnergy','Cao2024_RawResidualEnergy','Proposed_RawResidualEnergy', ...
        'Fixed_RawInteractionPeak','Zhang2018_RawInteractionPeak','Cao2024_RawInteractionPeak','Proposed_RawInteractionPeak', ...
        'Fixed_PeakError_mean','Fixed_PeakError_std','Zhang2018_PeakError_mean','Zhang2018_PeakError_std','Cao2024_PeakError_mean','Cao2024_PeakError_std','Proposed_PeakError_mean','Proposed_PeakError_std', ...
        'Fixed_PositionSaturation_pct_mean','Fixed_PositionSaturation_pct_std','Zhang2018_PositionSaturation_pct_mean','Zhang2018_PositionSaturation_pct_std','Cao2024_PositionSaturation_pct_mean','Cao2024_PositionSaturation_pct_std','Proposed_PositionSaturation_pct_mean','Proposed_PositionSaturation_pct_std', ...
        'Fixed_PositionMaxRatio_mean','Zhang2018_PositionMaxRatio_mean','Cao2024_PositionMaxRatio_mean','Proposed_PositionMaxRatio_mean'});
    writetable(T_diag,fullfile(mc_folder,['MC_diagnostic_four_controller_' mc_timestamp '.csv']));

    T_sat = table(Joints,N_success, ...
        mean(F.PositionSaturationPct,1)',std(F.PositionSaturationPct,0,1)', ...
        mean(Z.PositionSaturationPct,1)',std(Z.PositionSaturationPct,0,1)', ...
        mean(C.PositionSaturationPct,1)',std(C.PositionSaturationPct,0,1)', ...
        mean(A.PositionSaturationPct,1)',std(A.PositionSaturationPct,0,1)', ...
        mean(F.PositionMaxRatio,1)',mean(Z.PositionMaxRatio,1)',mean(C.PositionMaxRatio,1)',mean(A.PositionMaxRatio,1)', ...
        'VariableNames',{'Joint','N_Success','Fixed_PositionSaturation_pct_mean','Fixed_PositionSaturation_pct_std', ...
        'Zhang2018_PositionSaturation_pct_mean','Zhang2018_PositionSaturation_pct_std', ...
        'Cao2024_PositionSaturation_pct_mean','Cao2024_PositionSaturation_pct_std', ...
        'Proposed_PositionSaturation_pct_mean','Proposed_PositionSaturation_pct_std', ...
        'Fixed_PositionMaxRatio_mean','Zhang2018_PositionMaxRatio_mean','Cao2024_PositionMaxRatio_mean','Proposed_PositionMaxRatio_mean'});
    writetable(T_sat,fullfile(mc_folder,['MC_position_saturation_summary_' mc_timestamp '.csv']));

    Run = find(valid_idx);
    T_sat_trials = table(Run, ...
        F.PositionSaturationPct(:,1),F.PositionSaturationPct(:,2),F.PositionSaturationPct(:,3), ...
        Z.PositionSaturationPct(:,1),Z.PositionSaturationPct(:,2),Z.PositionSaturationPct(:,3), ...
        C.PositionSaturationPct(:,1),C.PositionSaturationPct(:,2),C.PositionSaturationPct(:,3), ...
        A.PositionSaturationPct(:,1),A.PositionSaturationPct(:,2),A.PositionSaturationPct(:,3), ...
        'VariableNames',{'Run','Fixed_Hip_pct','Fixed_Knee_pct','Fixed_Ankle_pct', ...
        'Zhang2018_Hip_pct','Zhang2018_Knee_pct','Zhang2018_Ankle_pct', ...
        'Cao2024_Hip_pct','Cao2024_Knee_pct','Cao2024_Ankle_pct', ...
        'Proposed_Hip_pct','Proposed_Knee_pct','Proposed_Ankle_pct'});
    writetable(T_sat_trials,fullfile(mc_folder,['MC_position_saturation_trials_' mc_timestamp '.csv']));

    %% ================= Trial-level long data =================
    Run_long = repelem(find(valid_idx),3,1);
    Joint_long = repmat(Joints,stats.success_count,1);
    flatten = @(x) reshape(x.',[],1);
    T_primary_trials = table(Run_long,Joint_long, ...
        flatten(F.RMSE),flatten(Z.RMSE),flatten(C.RMSE),flatten(A.RMSE), ...
        flatten(F.JerkRMS*180/pi),flatten(Z.JerkRMS*180/pi),flatten(C.JerkRMS*180/pi),flatten(A.JerkRMS*180/pi), ...
        flatten(F.RMS_InteractionTorque),flatten(Z.RMS_InteractionTorque),flatten(C.RMS_InteractionTorque),flatten(A.RMS_InteractionTorque), ...
        flatten(F.PositionSaturationPct),flatten(Z.PositionSaturationPct),flatten(C.PositionSaturationPct),flatten(A.PositionSaturationPct), ...
        flatten(F.MaxError),flatten(Z.MaxError),flatten(C.MaxError),flatten(A.MaxError), ...
        'VariableNames',{'Run','Joint','Fixed_RMSE_deg','Zhang2018_RMSE_deg','Cao2024_RMSE_deg','Proposed_RMSE_deg', ...
        'Fixed_JerkRMS_deg_s3','Zhang2018_JerkRMS_deg_s3','Cao2024_JerkRMS_deg_s3','Proposed_JerkRMS_deg_s3', ...
        'Fixed_RawResidualRMS_Nm','Zhang2018_RawResidualRMS_Nm','Cao2024_RawResidualRMS_Nm','Proposed_RawResidualRMS_Nm', ...
        'Fixed_PositionSaturation_pct','Zhang2018_PositionSaturation_pct','Cao2024_PositionSaturation_pct','Proposed_PositionSaturation_pct', ...
        'Fixed_PeakErrorDiagnostic_deg','Zhang2018_PeakErrorDiagnostic_deg','Cao2024_PeakErrorDiagnostic_deg','Proposed_PeakErrorDiagnostic_deg'});
    writetable(T_primary_trials,fullfile(mc_folder,['MC_trial_level_four_controller_' mc_timestamp '.csv']));

    red = @(base,comp) 100*(base-comp)./max(abs(base),eps); % positive = comp lower/better for lower-is-better metrics
    T_change = table(Joints, ...
        mean(red(F.RMSE,Z.RMSE),1)',mean(red(F.RMSE,C.RMSE),1)',mean(red(F.RMSE,A.RMSE),1)',mean(red(Z.RMSE,A.RMSE),1)',mean(red(C.RMSE,A.RMSE),1)', ...
        mean(red(F.JerkRMS,Z.JerkRMS),1)',mean(red(F.JerkRMS,C.JerkRMS),1)',mean(red(F.JerkRMS,A.JerkRMS),1)',mean(red(Z.JerkRMS,A.JerkRMS),1)',mean(red(C.JerkRMS,A.JerkRMS),1)', ...
        mean(red(F.RMS_InteractionTorque,Z.RMS_InteractionTorque),1)',mean(red(F.RMS_InteractionTorque,C.RMS_InteractionTorque),1)',mean(red(F.RMS_InteractionTorque,A.RMS_InteractionTorque),1)',mean(red(Z.RMS_InteractionTorque,A.RMS_InteractionTorque),1)',mean(red(C.RMS_InteractionTorque,A.RMS_InteractionTorque),1)', ...
        'VariableNames',{'Joint', ...
        'Zhang_vs_Fixed_RMSE_reduction_pct','Cao_vs_Fixed_RMSE_reduction_pct','Proposed_vs_Fixed_RMSE_reduction_pct','Proposed_vs_Zhang_RMSE_reduction_pct','Proposed_vs_Cao_RMSE_reduction_pct', ...
        'Zhang_vs_Fixed_Jerk_reduction_pct','Cao_vs_Fixed_Jerk_reduction_pct','Proposed_vs_Fixed_Jerk_reduction_pct','Proposed_vs_Zhang_Jerk_reduction_pct','Proposed_vs_Cao_Jerk_reduction_pct', ...
        'Zhang_vs_Fixed_RawResidual_reduction_pct','Cao_vs_Fixed_RawResidual_reduction_pct','Proposed_vs_Fixed_RawResidual_reduction_pct','Proposed_vs_Zhang_RawResidual_reduction_pct','Proposed_vs_Cao_RawResidual_reduction_pct'});
    writetable(T_change,fullfile(mc_folder,['MC_paired_change_four_controller_' mc_timestamp '.csv']));

    save(fullfile(mc_folder,['MC_summary_tables_' mc_timestamp '.mat']), ...
        'T','T_diag','T_sat','T_sat_trials','T_primary_trials','T_change','-v7.3');

    %% ================= Paper-style MC figures =================
    if stats.success_count >= 2
        labels = {'Fixed','Zhang18','Cao24','Proposed'};
        joints_names = {'Hip','Knee','Ankle'};
        [fig_box, lay_box] = paper_figure_utils('create','Monte Carlo Primary Metrics', ...
            [80 40 1550 1250],3,3,sprintf('Paired Monte Carlo: Four Controllers (n=%d/%d)',stats.success_count,num_runs));
        letters = 'abcdefghi'; pidx=0;
        for j=1:3
            pidx=pidx+1; ax=paper_figure_utils('tile',lay_box,pidx);
            draw_box4(ax,F.RMSE(:,j),Z.RMSE(:,j),C.RMSE(:,j),A.RMSE(:,j),labels);
            ylabel(ax,'RMSE (deg)'); paper_figure_utils('panel_title',ax,sprintf('(%c) %s Tracking RMSE',letters(pidx),joints_names{j}));
        end
        for j=1:3
            pidx=pidx+1; ax=paper_figure_utils('tile',lay_box,pidx);
            draw_box4(ax,F.JerkRMS(:,j)*180/pi,Z.JerkRMS(:,j)*180/pi,C.JerkRMS(:,j)*180/pi,A.JerkRMS(:,j)*180/pi,labels);
            ylabel(ax,'Jerk RMS (deg/s^3)'); paper_figure_utils('panel_title',ax,sprintf('(%c) %s Jerk RMS',letters(pidx),joints_names{j}));
        end
        for j=1:3
            pidx=pidx+1; ax=paper_figure_utils('tile',lay_box,pidx);
            draw_box4(ax,F.RMS_InteractionTorque(:,j),Z.RMS_InteractionTorque(:,j),C.RMS_InteractionTorque(:,j),A.RMS_InteractionTorque(:,j),labels);
            ylabel(ax,'Raw Residual RMS (N m)'); paper_figure_utils('panel_title',ax,sprintf('(%c) %s Raw Residual',letters(pidx),joints_names{j}));
        end
        paper_figure_utils('export',fig_box,mc_folder,['MC_boxplot_four_controller_' mc_timestamp],'SavePNG',true);

        [fig_sat, lay_sat] = paper_figure_utils('create','Position-Correction Saturation', ...
            [120 120 1450 500],1,3,sprintf('Supervisory Diagnostic: Shared Position-Correction Limit (n=%d/%d)',stats.success_count,num_runs));
        for j=1:3
            ax=paper_figure_utils('tile',lay_sat,j);
            draw_box4(ax,F.PositionSaturationPct(:,j),Z.PositionSaturationPct(:,j),C.PositionSaturationPct(:,j),A.PositionSaturationPct(:,j),labels);
            ylabel(ax,'Position saturation (%)'); paper_figure_utils('panel_title',ax,sprintf('(%c) %s Delta-q Limit Activation',char('a'+j-1),joints_names{j}));
        end
        paper_figure_utils('export',fig_sat,mc_folder,['MC_position_saturation_diagnostic_' mc_timestamp],'SavePNG',true);
    else
        warning('MonteCarlo:TooFewValidTrials','有效 trial 少于2个，跳过箱线图。');
    end

    if stats.fail_count > 0
        fprintf('⚠️ %d/%d 个 trial 失败，已从四控制器配对统计中排除。\n',stats.fail_count,num_runs);
        for i=1:numel(stats.fail_reasons), fprintf('  - %s\n',stats.fail_reasons{i}); end
    else
        fprintf('✅ 全部 %d 次四控制器配对运行成功。\n',num_runs);
    end
end

function s = init_metric_block(n)
    names = {'RMSE','JerkRMS','RMS_InteractionTorque','QuietZone','InteractionEnergy', ...
        'PeakInteractionTorque','MaxError','PositionSaturationPct','PositionMaxRatio', ...
        'PathologicalRMS','ControllerInputRMS','ControllerInputEnergy'};
    s = struct();
    for k=1:numel(names), s.(names{k}) = nan(n,3); end
end

function s = store_metric_block(s,r,m,qflag)
    s.RMSE(r,:) = m.RMSE;
    s.JerkRMS(r,:) = m.RMS_jerk;
    s.RMS_InteractionTorque(r,:) = m.RMS_InteractionTorque;
    s.QuietZone(r,:) = mean(qflag,1)*100;
    s.InteractionEnergy(r,:) = m.InteractionEnergy;
    s.PeakInteractionTorque(r,:) = m.PeakInteractionTorque;
    s.MaxError(r,:) = m.MaxError;
    s.PositionSaturationPct(r,:) = m.PositionCorrectionSaturationPct;
    s.PositionMaxRatio(r,:) = m.PositionCorrectionMaxRatio;
    s.PathologicalRMS(r,:) = m.RMS_PathologicalInteraction;
    s.ControllerInputRMS(r,:) = m.RMS_ControllerTorque;
    s.ControllerInputEnergy(r,:) = m.ControllerInputEnergy;
end

function out = extract_valid(s,idx)
    out = struct();
    f = fieldnames(s);
    for k=1:numel(f), out.(f{k}) = s.(f{k})(idx,:); end
end

function draw_box4(ax,f,z,c,a,labels)
    boxplot(ax,[f;z;c;a],[zeros(size(f));ones(size(z));2*ones(size(c));3*ones(size(a))],'Labels',labels);
end

function validate_trial_outputs(q,dq,ddq,e_fb_raw,metrics,label)
    if ~isreal(q) || any(~isfinite(q(:))), error('MonteCarlo:InvalidTrialOutput','%s q 包含 NaN/Inf/复数',label); end
    if ~isreal(dq) || any(~isfinite(dq(:))), error('MonteCarlo:InvalidTrialOutput','%s dq 包含 NaN/Inf/复数',label); end
    if ~isreal(ddq) || any(~isfinite(ddq(:))), error('MonteCarlo:InvalidTrialOutput','%s ddq 包含 NaN/Inf/复数',label); end
    if ~isreal(e_fb_raw) || any(~isfinite(e_fb_raw(:))), error('MonteCarlo:InvalidTrialOutput','%s interaction torque 包含 NaN/Inf/复数',label); end
    required = {'RMSE','RMS_jerk','RMS_InteractionTorque','InteractionEnergy','PeakInteractionTorque','MaxError', ...
        'PositionCorrectionSaturationPct','PositionCorrectionMaxRatio','RMS_PathologicalInteraction','RMS_ControllerTorque','ControllerInputEnergy'};
    for kk=1:numel(required)
        val=metrics.(required{kk});
        if ~isreal(val) || any(~isfinite(val(:)))
            error('MonteCarlo:InvalidTrialMetric','%s 指标 %s 包含 NaN/Inf/复数',label,required{kk});
        end
    end
end

function T = build_mc_uncertainty_table(base_params)
    c = base_params.mc;
    Group = ["Admittance";"Literature baseline";"Literature baseline";"Normal feedback";"Pathology";"Pathology";"Noise";"Event timing";"Quiet zone";"Excluded";"Excluded"];
    Parameter = ["Ma0/Ba0/Ka0/Ga0";"Zhang2018 a1/a2/a3/b1/b2/b3";"Cao2024 zeta_B/zeta_K/eta_B/eta_K";"Kq_fb/Kv_fb/Ke_fb";"spasm thresholds";"spasm gains/Csyn/cfatigue/event_amp/fatigue_slope";"feedback_noise_level";"event_time_offset";"quiet_dq_th/quiet_ddq_th";"model.L1/L2";"B0"];
    Distribution = ["bounded uniform";"fixed after nominal audit";"fixed after structural + nominal audit";"bounded uniform";"bounded uniform";"bounded uniform";"bounded uniform";"bounded uniform";"bounded uniform";"fixed nominal";"fixed nominal"];
    Range = [string(sprintf("nominal ± %.0f%%",100*c.rel_admittance));"0%";"0%";string(sprintf("nominal ± %.0f%%",100*c.rel_normal_feedback));string(sprintf("nominal ± %.0f%%",100*c.rel_spasm_threshold));string(sprintf("nominal ± %.0f%%",100*c.rel_pathology));string(sprintf("nominal ± %.0f%%",100*c.rel_noise));string(sprintf("± %.3f s",c.event_time_offset_abs));string(sprintf("nominal ± %.0f%%",100*c.rel_quiet_threshold));"0%";"0%"];
    Status = ["Perturbed";"Not perturbed";"Not perturbed";"Perturbed";"Perturbed";"Perturbed";"Perturbed";"Perturbed";"Perturbed";"Not perturbed";"Not perturbed"];
    Rationale = ["Shared nominal admittance parameters";"Zhang scheduling coefficients frozen before Monte Carlo";"Cao-derived effective momentum gains frozen before Monte Carlo";"Active interaction-feedback parameters";"Active spasticity gates";"Active pathological disturbance parameters";"Active feedback noise amplitude";"Moves event centers and event metric windows together";"Shared quiet-zone supervisor thresholds";"Link lengths do not enter current joint-space MC metrics";"Constant baseline term cancels in current feedback-error construction"];
    T = table(Group,Parameter,Distribution,Range,Status,Rationale);
end
