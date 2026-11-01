%% ================= 下肢康复机器人柔顺控制仿真主入口 =================
clc; clear; close all;
addpath(genpath(pwd));
paper_figure_utils('defaults');

%% ================= 0. 一键完整工作流配置 =================
% 唯一正式入口：运行 main_simulation 即按顺序执行全部已纳入工程的研究。
% 完整论文级运行默认开启所有实验；各正式实验输出均保留在 results/ 中。
RUN_MONTE_CARLO = true;
MC_NUM_RUNS = 100;
CLEAN_OLD_MC = false;              % 默认不删除历史结果

RUN_ABLATION_WORKFLOW = true;
ABLATION_NUM_RUNS = 100;
MAKE_ABLATION_FIGURES = true;
KMIS_GRID = [0 0.25 0.5 0.75 1];

RUN_EXECUTION_VALIDATION = true;   % 六工况非理想执行层验证
RUN_TRUE_ESTIMATED_MISMATCH = true; % Reviewer 3 C4：真/估参数失配
MODEL_MISMATCH_NUM_RUNS = 100;
MODEL_MISMATCH_SEVERITIES = [0 0.10 0.20];

RUN_JOINT_SPECIFIC = true;         % Reviewer 3 C3：关节特异性标称参数
JOINT_SPECIFIC_NUM_RUNS = 100;

EXPORT_ADDITIONAL_PAPER_EVIDENCE = true; % Fig16、Table6/7、实验目录索引

RUN_DEADZONE_SENSITIVITY = true;   % Reviewer 2 P11: soft dead-zone/noise paired study
RUN_STABILITY_VALIDATION = true;   % R3-C5/R2-P10: read-only numerical evidence, no control-law change
DEADZONE_SENSITIVITY_NUM_TRIALS = 100;

% Fig16/Table6/7 的正式归档以完整的 100-trial 数据为前提。
if EXPORT_ADDITIONAL_PAPER_EVIDENCE
    assert(RUN_MONTE_CARLO && RUN_ABLATION_WORKFLOW && RUN_EXECUTION_VALIDATION && ...
        RUN_TRUE_ESTIMATED_MISMATCH && RUN_JOINT_SPECIFIC, ...
        '论文级统一归档要求所有正式实验开关均开启。');
    assert(MC_NUM_RUNS == 100 && ABLATION_NUM_RUNS == 100 && ...
        MODEL_MISMATCH_NUM_RUNS == 100 && JOINT_SPECIFIC_NUM_RUNS == 100, ...
        '论文级统一归档要求各配对 Monte Carlo 实验均为 100 trials。');
end

%% ================= 1. 参数加载与初始化 =================
params = parameters();
rng(params.seed);

assert(isfield(params,'literature_zhang2018') && isstruct(params.literature_zhang2018), ...
       '缺少 params.literature_zhang2018 配置。');
assert(isfield(params,'literature_cao2024') && isstruct(params.literature_cao2024), ...
       '缺少 params.literature_cao2024 配置。');
assert(logical(params.literature_zhang2018.frozen), ...
       '正式主流程要求 Zhang-2018 comparator 已冻结。');
assert(logical(params.literature_cao2024.frozen), ...
       '正式主流程要求 Cao-2024-derived comparator 已冻结。');

fprintf('=== 下肢康复机器人四控制器柔顺控制仿真开始 ===\n');
fprintf('配置版本: %s | 时间戳: %s | 循环次数: %d\n', ...
        params.version, char(params.config_date), params.loop_num);
fprintf('对照组: Fixed | Zhang-2018 VA | Cao-2024-derived VA | Proposed IE-AAC\n');

%% ================= 2. 统一示教重建入口 =================
[params_updated, metrics_teaching, ...
 q_teach, dq_teach, ddq_teach, t_all] = TeachingDataPipeline(params, ...
    'Seed', params.seed, ...
    'PlotFlag', true, ...
    'SaveFlag', true, ...
    'Verbose', true);
params = params_updated;

%% ================= 3. 唯一权威参考轨迹一致性审计 =================
assert(isequal(size(q_teach), size(metrics_teaching.q_controller_ref)), ...
       '控制器参考轨迹尺寸与示教管线保存轨迹不一致');
assert(max(abs(q_teach(:) - metrics_teaching.q_controller_ref(:))) < 1e-14, ...
       'q_teach 与示教管线评估的控制器参考轨迹不一致');
assert(max(abs(dq_teach(:) - metrics_teaching.dq_controller_ref(:))) < 1e-14, ...
       'dq_teach 与示教管线评估的控制器参考轨迹不一致');
assert(max(abs(ddq_teach(:) - metrics_teaching.ddq_controller_ref(:))) < 1e-14, ...
       'ddq_teach 与示教管线评估的控制器参考轨迹不一致');
assert(max(abs(t_all(:) - metrics_teaching.t_controller_ref(:))) < 1e-14, ...
       't_all 与示教管线评估的控制器时间轴不一致');
fprintf('✅ 示教重建轨迹与四控制器实际参考轨迹严格一致。\n');

%% ================= 4. Monte Carlo：四控制器严格配对 =================
mc_result_dir = '';
primary_mc_statistics = [];


if RUN_MONTE_CARLO
    mc_root = 'results';
    if CLEAN_OLD_MC && exist(mc_root, 'dir')
        old_mc = dir(fullfile(mc_root, 'MC_*'));
        if ~isempty(old_mc)
            fprintf('⚠️ 发现 %d 个旧蒙特卡洛文件夹，将按 CLEAN_OLD_MC=true 清理。\n', length(old_mc));
            for i = 1:length(old_mc)
                rmdir(fullfile(old_mc(i).folder, old_mc(i).name), 's');
            end
            fprintf('✅ 旧蒙特卡洛数据已清除。\n');
        end
    end

    reference_bundle = struct();
    reference_bundle.q_ref = q_teach;
    reference_bundle.dq_ref = dq_teach;
    reference_bundle.ddq_ref = ddq_teach;
    reference_bundle.t_all = t_all;
    reference_bundle.metrics_teaching = metrics_teaching;

    fprintf('=== 自动运行四控制器蒙特卡洛仿真（%d 组，共享 nominal reference）===\n', MC_NUM_RUNS);
    tic;
    mc_result_dir = monte_carlo_runner(MC_NUM_RUNS, params, reference_bundle);
    %% Primary paired inferential statistics (same completed MC dataset)
    trial_files = dir(fullfile(mc_result_dir, 'MC_trial_level_four_controller_*.csv'));
    assert(numel(trial_files)==1, 'Expected exactly one primary MC trial-level CSV.');
    primary_mc_statistics = export_main_mc_pairwise( ...
        fullfile(mc_result_dir, trial_files(1).name), ...
        fullfile(mc_result_dir, 'paired_statistics'));
    fprintf('Primary MC Holm-adjusted statistics: %s\n',primary_mc_statistics.JointCSV);
    elapsed_mc = toc;
    fprintf('=== 蒙特卡洛完成，耗时 %.1f 秒 ===\n\n', elapsed_mc);
end

%% ================= 5. 构造四控制器共享的预生成扰动 =================
[pre_fb, event_timing] = mc_trial_utils('pre_fb', params, q_teach, t_all, params.seed);
params.t_abnormal = event_timing.abnormal_window;

fprintf('   Event centers (s): Hip %.3f | Knee %.3f | Ankle %.3f\n', ...
        event_timing.centers(1), event_timing.centers(2), event_timing.centers(3));
fprintf('   Event-active interval (±%.1fσ): [%.3f, %.3f] s\n', ...
        event_timing.n_sigma, event_timing.abnormal_window(1), event_timing.abnormal_window(2));
fprintf('✅ 四控制器共享 fatigue/event/noise realization 已构造完成。\n');

%% ================= 6. Fixed nominal admittance =================
[delta_q_fixed, delta_dq_fixed, delta_ddq_fixed, ...
 Ba_fixed, Ka_fixed, Ga_fixed, quiet_zone_flag_fixed, e_fb_used_fixed, e_fb_raw_fixed, interaction_diag_fixed] = ...
    run_admittance_controller('fixed', q_teach, dq_teach, params, pre_fb, t_all, ddq_teach);
q_fixed   = q_teach + delta_q_fixed;
dq_fixed  = dq_teach + delta_dq_fixed;
ddq_fixed = ddq_teach + delta_ddq_fixed;

%% ================= 7. Zhang-2018 literature-derived variable admittance =================
[delta_q_zhang, delta_dq_zhang, delta_ddq_zhang, ...
 Ba_zhang, Ka_zhang, Ga_zhang, quiet_zone_flag_zhang, e_fb_used_zhang, e_fb_raw_zhang, interaction_diag_zhang] = ...
    run_admittance_controller('zhang2018', q_teach, dq_teach, params, pre_fb, t_all, ddq_teach);
q_zhang   = q_teach + delta_q_zhang;
dq_zhang  = dq_teach + delta_dq_zhang;
ddq_zhang = ddq_teach + delta_ddq_zhang;

%% ================= 8. Cao-2024-derived momentum variable admittance =================
[delta_q_cao, delta_dq_cao, delta_ddq_cao, ...
 Ba_cao, Ka_cao, Ga_cao, quiet_zone_flag_cao, e_fb_used_cao, e_fb_raw_cao, interaction_diag_cao] = ...
    run_admittance_controller('cao2024', q_teach, dq_teach, params, pre_fb, t_all, ddq_teach);
q_cao   = q_teach + delta_q_cao;
dq_cao  = dq_teach + delta_dq_cao;
ddq_cao = ddq_teach + delta_ddq_cao;

%% ================= 9. Proposed IE-AAC =================
[delta_q_adp, delta_dq_adp, delta_ddq_adp, ...
 Ba_adp, Ka_adp, Ga_adp, quiet_zone_flag_adp, e_fb_used_adp, e_fb_raw_adp, interaction_diag_adp] = ...
    run_admittance_controller('proposed', q_teach, dq_teach, params, pre_fb, t_all, ddq_teach);
q_adaptive   = q_teach + delta_q_adp;
dq_adaptive  = dq_teach + delta_dq_adp;
ddq_adaptive = ddq_teach + delta_ddq_adp;

%% ================= 10. 足部轨迹 =================
[~, ~, ~, heel_teach, toe_teach] = compute_trajectories(q_teach, params);
[~, ~, ~, heel_fixed, toe_fixed] = compute_trajectories(q_fixed, params);
[~, ~, ~, heel_zhang, toe_zhang] = compute_trajectories(q_zhang, params);
[~, ~, ~, heel_cao, toe_cao] = compute_trajectories(q_cao, params);
[~, ~, ~, heel_adaptive, toe_adaptive] = compute_trajectories(q_adaptive, params);

%% ================= 11. 统一性能指标 =================
metrics_fixed = compute_metrics(q_teach, q_fixed, dq_fixed, ddq_fixed, ...
                                e_fb_raw_fixed, e_fb_used_fixed, params, interaction_diag_fixed);
metrics_zhang = compute_metrics(q_teach, q_zhang, dq_zhang, ddq_zhang, ...
                                e_fb_raw_zhang, e_fb_used_zhang, params, interaction_diag_zhang);
metrics_cao = compute_metrics(q_teach, q_cao, dq_cao, ddq_cao, ...
                              e_fb_raw_cao, e_fb_used_cao, params, interaction_diag_cao);
metrics_adaptive = compute_metrics(q_teach, q_adaptive, dq_adaptive, ddq_adaptive, ...
                                   e_fb_raw_adp, e_fb_used_adp, params, interaction_diag_adp);

fprintf('\n=== Nominal primary metrics: Fixed | Zhang-2018 | Cao-2024-derived | Proposed ===\n');
Joint = {'Hip';'Knee';'Ankle'};
T_nominal_primary = table(Joint, ...
    metrics_fixed.RMSE(:), metrics_zhang.RMSE(:), metrics_cao.RMSE(:), metrics_adaptive.RMSE(:), ...
    metrics_fixed.RMS_jerk_deg(:), metrics_zhang.RMS_jerk_deg(:), metrics_cao.RMS_jerk_deg(:), metrics_adaptive.RMS_jerk_deg(:), ...
    metrics_fixed.RMS_InteractionTorque(:), metrics_zhang.RMS_InteractionTorque(:), metrics_cao.RMS_InteractionTorque(:), metrics_adaptive.RMS_InteractionTorque(:), ...
    'VariableNames', {'Joint', ...
    'Fixed_RMSE_deg','Zhang2018_RMSE_deg','Cao2024Derived_RMSE_deg','Proposed_RMSE_deg', ...
    'Fixed_JerkRMS_deg_s3','Zhang2018_JerkRMS_deg_s3','Cao2024Derived_JerkRMS_deg_s3','Proposed_JerkRMS_deg_s3', ...
    'Fixed_RawResidualRMS_Nm','Zhang2018_RawResidualRMS_Nm','Cao2024Derived_RawResidualRMS_Nm','Proposed_RawResidualRMS_Nm'});
disp(T_nominal_primary);

fprintf('\n=== Nominal supervisory diagnostics (not primary performance metrics) ===\n');
T_nominal_diagnostic = table(Joint, ...
    metrics_fixed.MaxError(:), metrics_zhang.MaxError(:), metrics_cao.MaxError(:), metrics_adaptive.MaxError(:), ...
    metrics_fixed.PositionCorrectionSaturationPct(:), metrics_zhang.PositionCorrectionSaturationPct(:), metrics_cao.PositionCorrectionSaturationPct(:), metrics_adaptive.PositionCorrectionSaturationPct(:), ...
    metrics_fixed.PositionCorrectionMaxRatio(:), metrics_zhang.PositionCorrectionMaxRatio(:), metrics_cao.PositionCorrectionMaxRatio(:), metrics_adaptive.PositionCorrectionMaxRatio(:), ...
    'VariableNames', {'Joint', ...
    'Fixed_PeakError_deg','Zhang2018_PeakError_deg','Cao2024Derived_PeakError_deg','Proposed_PeakError_deg', ...
    'Fixed_PositionSaturation_pct','Zhang2018_PositionSaturation_pct','Cao2024Derived_PositionSaturation_pct','Proposed_PositionSaturation_pct', ...
    'Fixed_PositionMaxRatio','Zhang2018_PositionMaxRatio','Cao2024Derived_PositionMaxRatio','Proposed_PositionMaxRatio'});
disp(T_nominal_diagnostic);

T_zhang_audit = nominal_diagnostics('zhang2018', ...
    Ba_zhang, Ka_zhang, Ga_zhang, delta_q_zhang, delta_dq_zhang, delta_ddq_zhang, ...
    quiet_zone_flag_zhang, params, interaction_diag_zhang);
fprintf('\n=== Zhang-2018 literature comparator structural audit ===\n');
disp(T_zhang_audit);

T_cao_audit = nominal_diagnostics('cao2024', ...
    Ba_cao, Ka_cao, Ga_cao, delta_q_cao, delta_dq_cao, delta_ddq_cao, ...
    quiet_zone_flag_cao, params, interaction_diag_cao);
fprintf('\n=== Cao-2024-derived literature comparator structural audit ===\n');
disp(T_cao_audit);

T_four_controller_limits = nominal_diagnostics('motion_limits', ...
    delta_q_fixed, delta_dq_fixed, delta_ddq_fixed, ...
    delta_q_zhang, delta_dq_zhang, delta_ddq_zhang, ...
    delta_q_cao, delta_dq_cao, delta_ddq_cao, ...
    delta_q_adp, delta_dq_adp, delta_ddq_adp, params);
fprintf('\n=== Common motion-limit audit: four formal controllers ===\n');
disp(T_four_controller_limits);

%% ================= 12. 保存完整结果 =================
ctrl_params_fixed = struct('Ba', Ba_fixed, 'Ka', Ka_fixed, 'Ga', Ga_fixed, ...
                           'quiet_zone_flag', quiet_zone_flag_fixed);
ctrl_params_zhang = struct('Ba', Ba_zhang, 'Ka', Ka_zhang, 'Ga', Ga_zhang, ...
                           'quiet_zone_flag', quiet_zone_flag_zhang, ...
                           'source', params.literature_zhang2018.source, ...
                           'doi', params.literature_zhang2018.doi);
ctrl_params_cao = struct('Ba', Ba_cao, 'Ka', Ka_cao, 'Ga', Ga_cao, ...
                         'quiet_zone_flag', quiet_zone_flag_cao, ...
                         'source', params.literature_cao2024.source, ...
                         'doi', params.literature_cao2024.doi);
ctrl_params_adp = struct('Ba', Ba_adp, 'Ka', Ka_adp, 'Ga', Ga_adp, ...
                         'quiet_zone_flag', quiet_zone_flag_adp);

fb_nominal_adp = interaction_diag_adp.fb_nominal_ref;
fb_true_adp = interaction_diag_adp.fb_true;
fb_real_adp = interaction_diag_adp.fb_true;

save_results(params, t_all, ...
    q_teach, q_fixed, q_zhang, q_cao, q_adaptive, ...
    dq_teach, dq_fixed, dq_zhang, dq_cao, dq_adaptive, ...
    ddq_teach, ddq_fixed, ddq_zhang, ddq_cao, ddq_adaptive, ...
    heel_teach, heel_fixed, heel_zhang, heel_cao, heel_adaptive, ...
    toe_teach, toe_fixed, toe_zhang, toe_cao, toe_adaptive, ...
    ctrl_params_fixed, ctrl_params_zhang, ctrl_params_cao, ctrl_params_adp, ...
    metrics_fixed, metrics_zhang, metrics_cao, metrics_adaptive, ...
    e_fb_raw_fixed, e_fb_used_fixed, ...
    e_fb_raw_zhang, e_fb_used_zhang, ...
    e_fb_raw_cao, e_fb_used_cao, ...
    e_fb_raw_adp, e_fb_used_adp, ...
    fb_nominal_adp, fb_true_adp, fb_real_adp, ...
    interaction_diag_fixed, interaction_diag_zhang, interaction_diag_cao, interaction_diag_adp, ...
    T_nominal_primary, T_nominal_diagnostic);

nominal_result_dir = fullfile('results', char(params.config_date));
writetable(T_zhang_audit, fullfile(nominal_result_dir, ...
    ['zhang2018_structural_audit_' char(params.config_date) '.csv']));
writetable(T_cao_audit, fullfile(nominal_result_dir, ...
    ['cao2024_derived_structural_audit_' char(params.config_date) '.csv']));
writetable(T_four_controller_limits, fullfile(nominal_result_dir, ...
    ['four_controller_motion_limit_audit_' char(params.config_date) '.csv']));

%% ================= 12A. R3-C5 / R2-P10: nominal stability evidence =================
% Read-only post-processing of the already computed Proposed trace.
% Fixed-parameter Schur test is NOT an adaptive closed-loop stability proof.
if RUN_STABILITY_VALIDATION
    trace_stability = struct('t',t_all, ...
        'Ba',Ba_adp,'Ka',Ka_adp,'Ga',Ga_adp, ...
        'delta_q',delta_q_adp,'delta_dq',delta_dq_adp, ...
        'delta_ddq',delta_ddq_adp,'e_used',e_fb_used_adp);
    stability_evidence = run_stability_validation(params,trace_stability, ...
        fullfile(nominal_result_dir,'stability_validation'));
    fprintf('R3-C5/R2-P10 stability evidence: %s\n',stability_evidence.OutputDir);
    % Archive stability table sources with the other manuscript figures/tables.
    export_stability_paper_evidence(stability_evidence,nominal_result_dir);
end

%% ================= 13. 论文级图表 =================
if RUN_MONTE_CARLO
    mc_path_for_figures = '';
else
    mc_path_for_figures = '__SKIP_MC__';
end

generate_paper_figures(t_all, q_teach, dq_teach, ddq_teach, ...
    q_fixed, dq_fixed, ddq_fixed, ...
    q_zhang, dq_zhang, ddq_zhang, ...
    q_cao, dq_cao, ddq_cao, ...
    q_adaptive, dq_adaptive, ddq_adaptive, ...
    Ba_zhang, Ka_zhang, Ga_zhang, ...
    Ba_cao, Ka_cao, Ga_cao, ...
    Ba_adp, Ka_adp, Ga_adp, ...
    metrics_fixed, metrics_zhang, metrics_cao, metrics_adaptive, ...
    params, char(params.config_date), mc_path_for_figures);

if RUN_MONTE_CARLO
    fprintf('\n🎉 四控制器 nominal + paired Monte Carlo 流程完成。\n');
else
    fprintf('\n✅ 四控制器 nominal-only 流程完成（未运行 Monte Carlo）。\n');
end

%% ================= 14. 自动执行完整消融与最终证据 =================
R_ablation = [];

if RUN_ABLATION_WORKFLOW
    fprintf('\n============================================================\n');
    fprintf(' Primary workflow completed. Starting ablation automatically.\n');
    fprintf(' Core ablation -> K_mis sensitivity -> finalize\n');
    fprintf(' Ablation paired Monte Carlo runs: %d\n', ABLATION_NUM_RUNS);
    fprintf('============================================================\n\n');

    R_ablation = ablation_study('all', ...
        'NumRuns', ABLATION_NUM_RUNS, ...
        'MakeFigures', MAKE_ABLATION_FIGURES, ...
        'KmisGrid', KMIS_GRID);

    fprintf('\n============================================================\n');
    fprintf(' IE-AAC ONE-CLICK COMPLETE WORKFLOW FINISHED\n');
    fprintf(' Primary result folder: %s\n', fullfile('results', char(params.config_date)));

    if isfield(R_ablation,'core') && isfield(R_ablation.core,'output_dir')
        fprintf(' Core ablation folder: %s\n', R_ablation.core.output_dir);
    end
    if isfield(R_ablation,'mismatch') && isfield(R_ablation.mismatch,'output_dir')
        fprintf(' K_mis sensitivity folder: %s\n', R_ablation.mismatch.output_dir);
    end
    if isfield(R_ablation,'final') && isfield(R_ablation.final,'output_dir')
        fprintf(' Final evidence folder: %s\n', R_ablation.final.output_dir);
    end
    fprintf('============================================================\n');
end

%% ================= 15. 非理想执行层与论文图自动归档 =================
% Uses the SAME nominal reference and common exogenous sequence as this run.
% Execution layer is OFF for original nominal/MC/ablation experiments.
if RUN_EXECUTION_VALIDATION
    fprintf('\n=== Starting six-case nonideal execution-layer validation ===\n');
    execution_evidence = run_execution_dynamics_validation( ...
        params, q_teach, dq_teach, ddq_teach, t_all, char(params.config_date));
    fprintf('Execution-layer figures and data: %s\n', execution_evidence.PaperFigureDir);
end

fprintf('\n=== Primary, ablation and execution experiments completed ===\n');

%% ================= 16. 真/估计标称参数失配验证（Reviewer 3 C4） =================
model_mismatch_result_dir = '';
if RUN_TRUE_ESTIMATED_MISMATCH
    fprintf('\n=== Starting true-vs-estimated interaction model experiment ===\n');
    model_mismatch_result_dir = run_model_mismatch_validation( ...
        MODEL_MISMATCH_NUM_RUNS, MODEL_MISMATCH_SEVERITIES);
    fprintf('Model-mismatch evidence: %s\n', model_mismatch_result_dir);
end

%% ================= 17. 关节特异性标称参数验证（Reviewer 3 C3） =================
joint_specific_result_dir = '';
if RUN_JOINT_SPECIFIC
    fprintf('\n=== Starting joint-specific nominal interaction experiment ===\n');
    joint_specific_result_dir = run_joint_specific_validation(JOINT_SPECIFIC_NUM_RUNS);
    fprintf('Joint-specific evidence: %s\n', joint_specific_result_dir);
end

%% ================= 18. 统一论文图表、表格源数据及实验目录索引 =================
if EXPORT_ADDITIONAL_PAPER_EVIDENCE
    assert(~isempty(model_mismatch_result_dir) && ~isempty(joint_specific_result_dir), ...
        '统一论文证据归档需要先完成真/估失配和关节特异性两项实验。');
    generate_validation_evidence(model_mismatch_result_dir, ...
        joint_specific_result_dir, char(params.config_date), ...
        mc_result_dir, nominal_result_dir, R_ablation);
end

%% ================= 19. Reviewer 2 P11：软死区与模拟反馈噪声敏感性 =================
% Independent full 100-trial experiment; does NOT overwrite primary MC.
% Reuses this run's common reference; appends P11 figure/data to paper archive.
if RUN_DEADZONE_SENSITIVITY
    fprintf('\n=== Starting R2 P11 soft-dead-zone/noise sensitivity ===\n');
    p11_evidence = run_deadzone_sensitivity_validation( ...
        params,q_teach,dq_teach,ddq_teach,t_all, ...
        DEADZONE_SENSITIVITY_NUM_TRIALS,char(params.config_date));
    fprintf('R2 P11 evidence: %s\n',p11_evidence.FigureDataDir);
end

fprintf('\n=== COMPLETE ONE-ENTRY WORKFLOW AND PAPER EVIDENCE FINISHED ===\n');
