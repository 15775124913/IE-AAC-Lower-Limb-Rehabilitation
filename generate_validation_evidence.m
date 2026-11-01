function evidence = generate_validation_evidence(model_dir, joint_dir, sim_timestamp, mc_dir, nominal_dir, ablation_result)
%GENERATE_VALIDATION_EVIDENCE Archive completed validation studies for the paper.
% Reads existing summaries; does not rerun trials, resample data, or modify
% controller/experiment calculations. Outputs follow generate_paper_figures:
% PDF/PNG in paper_figures; CSV/MAT/metadata in figure_data.
    arguments
        model_dir (1,:) char
        joint_dir (1,:) char
        sim_timestamp (1,:) char
        mc_dir (1,:) char
        nominal_dir (1,:) char
        ablation_result
    end
    assert(isfolder(model_dir), 'Missing model-mismatch folder: %s', model_dir);
    assert(isfolder(joint_dir), 'Missing joint-specific folder: %s', joint_dir);
    assert(isfolder(nominal_dir), 'Missing nominal folder: %s', nominal_dir);

    fig_dir = fullfile('results', sim_timestamp, 'paper_figures');
    data_dir = fullfile('results', sim_timestamp, 'figure_data');
    if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
    if ~exist(data_dir,'dir'), mkdir(data_dir); end

    mismatch_summary_path = fullfile(model_dir,'model_mismatch_group_summary.csv');
    mismatch_manifest_path = fullfile(model_dir,'model_mismatch_manifest.csv');
    joint_summary_path = fullfile(joint_dir,'joint_specific_group_summary.csv');
    joint_manifest_path = fullfile(joint_dir,'joint_specific_manifest.csv');
    assert(isfile(mismatch_summary_path) && isfile(mismatch_manifest_path), ...
        'Model-mismatch summary/manifest is incomplete.');
    assert(isfile(joint_summary_path) && isfile(joint_manifest_path), ...
        'Joint-specific summary/manifest is incomplete.');

    Tm = readtable(mismatch_summary_path,'TextType','string');
    Tj = readtable(joint_summary_path,'TextType','string');
    M = readtable(mismatch_manifest_path,'TextType','string');
    J = readtable(joint_manifest_path,'TextType','string');
    assert(height(M)==1 && height(J)==1 && M.NumTrials(1)==100 && J.NumTrials(1)==100, ...
        'Paper evidence export requires completed 100-trial studies.');

    severity = [0;10;20];
    methods = ["Fixed";"Zhang-2018-derived";"Cao-2024-derived"; ...
        "Proposed (Kmis=1)";"Proposed (Kmis=0)"];
    values = nan(3,5); sds = values;
    for a=1:numel(severity)
        for b=1:numel(methods)
            ix = Tm.SeverityPct==severity(a) & Tm.Controller==methods(b) & ...
                Tm.Joint=="ThreeJointMean" & Tm.Metric=="ReferenceRMSEDeg";
            assert(nnz(ix)==1 && Tm.N(ix)==100, ...
                'Missing/duplicate/incomplete Fig16 group at severity=%g, method=%s.',severity(a),methods(b));
            values(a,b)=Tm.Mean(ix); sds(a,b)=Tm.SD(ix);
        end
    end
    assert(all(isfinite(values(:))) && all(isfinite(sds(:))), 'Non-finite Fig16 summary values.');

    paper_figure_utils('defaults');
    [fig,layout] = paper_figure_utils('create','Fig16_model_mismatch_sensitivity', ...
        [110 80 1100 650],1,1, ...
        'True-estimated nominal parameter mismatch: paired simulation','Visible','off');
    ax = paper_figure_utils('tile',layout,1);
    colors = lines(5); hold(ax,'on');
    for b=1:numel(methods)
        plot(ax,severity,values(:,b),'o-','Color',colors(b,:), ...
            'LineWidth',1.8,'MarkerSize',7);
    end
    set(ax,'XTick',severity,'XLim',[-1 21]);
    xlabel(ax,'Nominal parameter uncertainty severity (%)');
    ylabel(ax,'Three-joint mean reference RMSE (deg)');
    legend(ax,cellstr(methods),'Location','best','FontSize',9); grid(ax,'on');
    paper_figure_utils('export',fig,fig_dir,'Fig16_model_mismatch_sensitivity');
    close(fig);

    n=numel(methods)*numel(severity);
    Sev=zeros(n,1); Method=strings(n,1); MeanRMSE=zeros(n,1); SDRMSE=zeros(n,1); N=100*ones(n,1); k=0;
    for a=1:numel(severity)
        for b=1:numel(methods)
            k=k+1; Sev(k)=severity(a); Method(k)=methods(b);
            MeanRMSE(k)=values(a,b); SDRMSE(k)=sds(a,b);
        end
    end
    Fig16=table(Sev,Method,N,MeanRMSE,SDRMSE, ...
        'VariableNames',{'SeverityPct','Controller','N','MeanRMSEDeg','SDRMSEDeg'});
    writetable(Fig16,fullfile(data_dir,'Fig16_model_mismatch_sensitivity_data.csv'));
    save(fullfile(data_dir,'Fig16_model_mismatch_sensitivity_data.mat'), ...
        'Fig16','Tm','M','values','sds','severity','methods','model_dir');
    Figure16Metadata=table("Figure 16",string(model_dir),"100 paired Monte Carlo trials", ...
        "Mean of three joint RMSE values within trial, then mean and SD across trials", ...
        "0/10/20% denote simulated true-model coefficient uncertainty severity", ...
        'VariableNames',{'ManuscriptItem','SourceFolder','Design','Statistic','ScopeLimit'});
    writetable(Figure16Metadata,fullfile(data_dir,'Fig16_model_mismatch_sensitivity_metadata.csv'));

    t6=Tm(Tm.Joint=="ThreeJointMean" & ...
        (Tm.Controller=="Proposed (Kmis=1)" | Tm.Controller=="Proposed (Kmis=0)") & ...
        ismember(Tm.Metric,["ReferenceRMSEDeg","JerkRMSRadPerS3", ...
        "RawResidualRMSNm","CompensatedContaminationRMSNm"]),:);
    assert(height(t6)==24 && all(t6.N==100), ...
        'Table 6 requires 3 severities x 2 Kmis modes x 4 metrics.');
    % Preserve the experiment's predefined order: severity 0/10/20%, then
    % Proposed Kmis=1 (ON) followed by Kmis=0 (OFF), with metrics in source order.
    writetable(t6,fullfile(data_dir,'Table6_model_mismatch_on_off_source.csv'));
    save(fullfile(data_dir,'Table6_model_mismatch_on_off_source.mat'),'t6','Tm');
    Table6Metadata=table("Table 6",string(model_dir),"100 paired trials", ...
        "Three-joint mean within each trial; Jerk RMS is in rad/s^3", ...
        'VariableNames',{'ManuscriptItem','SourceFolder','Design','UnitsNote'});
    writetable(Table6Metadata,fullfile(data_dir,'Table6_model_mismatch_metadata.csv'));

    t7=Tj(Tj.Joint=="ThreeJointMean" & ...
        ismember(Tj.Metric,["ReferenceRMSEDeg","JerkRMSDegPerS3","RawResidualRMSNm"]),:);
    assert(height(t7)==12 && all(t7.N==100), ...
        'Table 7 requires four controllers x three metrics with 100 trials.');
    % Preserve the predefined controller and metric order from the experiment summary.
    writetable(t7,fullfile(data_dir,'Table7_joint_specific_threejoint_source.csv'));
    save(fullfile(data_dir,'Table7_joint_specific_threejoint_source.mat'),'t7','Tj','J');
    Table7Metadata=table("Table 7",string(joint_dir),"100 paired Monte Carlo trials", ...
        "Joint-specific nominal coefficients; simulated truth and estimate are matched", ...
        'VariableNames',{'ManuscriptItem','SourceFolder','Design','ModelAssumption'});
    writetable(Table7Metadata,fullfile(data_dir,'Table7_joint_specific_metadata.csv'));

    %% Section 5.1: C6 reconstruction audit -> unified manuscript evidence.
    % Mirror the already-computed C6 files. Do not rerun reconstruction,
    % recompute a reference, change a controller, or sample random numbers.
    teaching_source_dir = fullfile('results', sim_timestamp);
    c6_quality_source = fullfile(teaching_source_dir, ...
        ['teaching_reconstruction_quality_' sim_timestamp '.csv']);
    c6_checks_source = fullfile(teaching_source_dir, ...
        ['teaching_reconstruction_checks_' sim_timestamp '.csv']);
    c6_metrics_source = fullfile(teaching_source_dir, ...
        ['teaching_metrics_' sim_timestamp '.mat']);
    assert(isfile(c6_quality_source) && isfile(c6_checks_source) && ...
        isfile(c6_metrics_source), ...
        'Missing C6 source quality/checks/metrics for timestamp %s.', sim_timestamp);
    C6Quality = readtable(c6_quality_source,'TextType','string');
    C6Checks = readtable(c6_checks_source,'TextType','string');
    assert(height(C6Quality)==8 && height(C6Checks)==6 && ...
        all(C6Checks.Passed==1), ...
        'C6 manuscript export requires eight quality and six passing check records.');
    assert(nnz(C6Quality.ReferenceSignal=="TeachingPCHIPTarget" & ...
        C6Quality.Joint=="AllJoints")==1 && ...
        nnz(C6Quality.ReferenceSignal=="LatentTruthOffline" & ...
        C6Quality.Joint=="AllJoints")==1, ...
        'C6 audit does not contain the two expected global reference baselines.');
    writetable(C6Quality,fullfile(data_dir, ...
        'Section5_1_reconstruction_quality_source.csv'));
    writetable(C6Checks,fullfile(data_dir, ...
        'Section5_1_reconstruction_checks_source.csv'));
    copyfile(c6_metrics_source,fullfile(data_dir, ...
        'Section5_1_reconstruction_metrics_source.mat'));
    C6Metadata = table("Section 5.1",string(c6_quality_source), ...
        "1001 first-cycle time samples for each of three joints", ...
        "PCHIP teaching proxy for fitting; latent truth is offline only", ...
        "1 degree tolerance applies to pooled fitting RMSE, not maximum error", ...
        'VariableNames',{'ManuscriptItem','SourceCSV', ...
        'SampleScope','ReferenceDefinition','ToleranceDefinition'});
    writetable(C6Metadata,fullfile(data_dir, ...
        'Section5_1_reconstruction_metadata.csv'));

    ablationDirs=strings(3,1);
    ablationFields={"core","mismatch","final"};
    if isstruct(ablation_result)
        for ii=1:3
            f=ablationFields{ii};
            if isfield(ablation_result,f) && isfield(ablation_result.(f),'output_dir')
                ablationDirs(ii)=string(ablation_result.(f).output_dir);
            else
                ablationDirs(ii)="";
            end
        end
    else
        ablationDirs(:)="";
    end
    items=["Reconstruction quality (R3 C6)";"Primary Monte Carlo";"Nominal results";"Core ablation";"K_mis sensitivity"; ...
        "Final ablation evidence";"Execution six-case";"True-estimated model mismatch"; ...
        "Joint-specific nominal study";"Conditional discrete stability audit";"Paper figures";"Paper figure data"];
    dirs=[string(data_dir);string(mc_dir);string(nominal_dir);ablationDirs; ...
        string(fullfile(nominal_dir,'execution_dynamics'));string(model_dir); ...
        string(joint_dir);string(fullfile(nominal_dir,'stability_validation'));string(fig_dir);string(data_dir)];
    keyFiles=["Section5_1_reconstruction_*_source.csv and metadata";"MC_trial_level_four_controller_*.csv";"nominal MAT/CSV + Fig05-Fig11"; ...
        "ablation_mc_trial_metrics_long.csv";"mismatch_mc_sensitivity_compact.csv"; ...
        "TableA-F evidence tables";"execution_dynamics_results.csv"; ...
        "model_mismatch_group_summary.csv";"joint_specific_group_summary.csv"; ...
        "stability_frozen_parameter_boxes.csv; figure_data/Section5_3_1_stability_*"; ...
        "Fig16 PDF/PNG";"Fig16/Table6/Table7 and Section5_1 CSV/MAT/metadata"];
    status=strings(numel(items),1);
    for ii=1:numel(items)
        if strlength(strtrim(dirs(ii)))==0
            status(ii)="NOT_RUN";
        elseif isfolder(char(dirs(ii)))
            status(ii)="FOLDER_PRESENT";
        else
            status(ii)="MISSING";
        end
    end
    Index=table(items,status,dirs,keyFiles, ...
        'VariableNames',{'Experiment','Status','RelativeFolder','KeyEvidence'});
    writetable(Index,fullfile(nominal_dir,'all_experiment_output_index.csv'));
    evidence=struct('PaperFigureDir',fig_dir,'FigureDataDir',data_dir, ...
        'Index',Index,'Figure16Source',Fig16,'Table6Source',t6,'Table7Source',t7, ...
        'C6QualitySource',C6Quality,'C6ChecksSource',C6Checks);
    save(fullfile(data_dir,'unified_paper_evidence.mat'), ...
        'evidence','Tm','Tj','M','J','Fig16','t6','t7','C6Quality','C6Checks','-v7.3');
    UnifiedMetadata=table("Unified evidence archive",string(sim_timestamp), ...
        "Reads completed summaries; no resampling", "main_simulation.m", ...
        'VariableNames',{'ArchiveItem','SimulationTimestamp','DataLineage','EntryPoint'});
    writetable(UnifiedMetadata,fullfile(data_dir,'unified_paper_evidence_metadata.csv'));
    fprintf('Figure 16, Table 6/7 and Section 5.1 C6 evidence saved under %s\n',data_dir);
    fprintf('Unified experiment index saved to %s\n',fullfile(nominal_dir,'all_experiment_output_index.csv'));
end
