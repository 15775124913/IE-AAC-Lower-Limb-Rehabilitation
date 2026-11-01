function evidence = run_deadzone_sensitivity_validation(params, q_ref, dq_ref, ddq_ref, t_all, num_trials, sim_timestamp)
%RUN_DEADZONE_SENSITIVITY_VALIDATION Reviewer 2 P11, noise/dead-zone study.
% Dedicated paired simulation; NO changes to the controller update law.
% Called by main_simulation.m after the previously frozen experiments.
% Can also be called with num_trials=2 for a quick smoke test.
%
% Grid: D = [0 .025 .05 .10] N*m; feedback-noise std = [0 .02 .04] N*m.
% All cells share the SAME reference, sampled non-noise parameters,
% event/fatigue realizations, and underlying unit-Gaussian noise in a trial.
% Exact nominal model, ideal motion execution and K_mis=1 are held fixed.
%
% Returned and saved data are *simulation diagnostics*, not measured forces.
% The matched-state noise transmission metric is an OFFLINE difference:
% soft(e_comp,D) - soft(e_comp - injected_noise,D), with both operands
% evaluated at the noisy-run commanded state. It does NOT establish a
% causal closed-loop noise-to-output transfer or closed-loop stability.

    if nargin < 6 || isempty(num_trials), num_trials = 100; end
    if nargin < 7 || isempty(sim_timestamp)
        sim_timestamp = char(params.config_date);
    end
    validateattributes(num_trials,{'numeric'}, ...
        {'scalar','integer','>=',2,'finite'},mfilename,'num_trials');
    sim_timestamp = char(sim_timestamp);
    assert(size(q_ref,2)==3 && isequal(size(q_ref),size(dq_ref)) && ...
        isequal(size(q_ref),size(ddq_ref)) && numel(t_all)==size(q_ref,1), ...
        'P11:ReferenceSize','Reference q/dq/ddq/time must have compatible sizes.');
    assert(all(isfinite(q_ref(:))) && all(isfinite(dq_ref(:))) && ...
        all(isfinite(ddq_ref(:))) && all(isfinite(t_all(:))), ...
        'P11:ReferenceFinite','Reference contains nonfinite values.');
    assert(~params.execution_dynamics.enabled && ...
        all(params.K_mismatch_comp(:)==1), ...
        'P11:StudyAssumptions','Requires ideal execution and K_mis=1.');
    assert(~isfield(params,'patient_nominal_truth'), ...
        'P11:StudyAssumptions','No true-estimated parameter separation in this study.');

    d_grid = [0, 0.025, 0.05, 0.10];
    noise_grid = [0, 0.02, 0.04];
    nominal_d = 0.05;
    assert(all(abs(params.e_deadzone(:)-nominal_d)<1e-14), ...
        'P11:BaselineDeadzone','The frozen baseline must be D=0.05 Nm at all joints.');
    assert(abs(params.feedback_noise_level-0.02)<1e-14, ...
        'P11:BaselineNoise','The frozen baseline must use sigma=0.02 Nm.');
    joints = ["Hip","Knee","Ankle","ThreeJointMean"];
    metric_keys = ["ReferenceRMSEDeg","JerkRMSRadPerS3", ...
        "RawResidualRMSNm","NoiseTransmissionRMSNm"];

    fig_dir = fullfile('results',sim_timestamp,'paper_figures');
    data_dir = fullfile('results',sim_timestamp,'figure_data');
    if ~isfolder(fig_dir), mkdir(fig_dir); end
    if ~isfolder(data_dir), mkdir(data_dir); end
    fprintf('\n=== Reviewer 2 P11: paired dead-zone/noise sensitivity ===\n');
    fprintf('Trials=%d, D=%s Nm, noise std=%s Nm\n', ...
        num_trials,mat2str(d_grid),mat2str(noise_grid));

    blank = struct('Trial',0,'DeadzoneNm',0,'NoiseStdNm',0, ...
        'Joint',"",'ReferenceRMSEDeg',0,'JerkRMSRadPerS3',0, ...
        'RawResidualRMSNm',0,'CompensatedRMSNm',0, ...
        'ProcessedInputRMSNm',0,'NoiseTransmissionRMSNm',0, ...
        'NoiseTransmissionRatio',nan,'DeadzoneZeroFraction',0, ...
        'CleanEventSignalRetention',0,'ParameterSeed',0,'NoiseSeed',0);
    rows = repmat(blank,num_trials*numel(d_grid)*numel(noise_grid)*4,1);
    at=0;
    for r=1:num_trials
        [sampled,seed_audit] = mc_trial_utils('sample',params,r);
        sampled.execution_dynamics.enabled=false;
        % Force unit variance for a reproducible standard-normal noise draw.
        % MC sampling still applies to ALL other parameters; noise amplitude
        % is an explicitly controlled grid factor rather than sampled +/-20%.
        unit_params = sampled;
        unit_params.feedback_noise_level = 1;
        [unit_fb,timing] = mc_trial_utils('pre_fb', ...
            unit_params,q_ref,t_all,seed_audit.noise_seed);
        unit_noise = unit_fb.noise_all;
        for ns=1:numel(noise_grid)
            noise_std=noise_grid(ns);
            pre_fb=unit_fb;
            pre_fb.noise_all = noise_std .* unit_noise;
            p = sampled;
            p.feedback_noise_level=noise_std;
            p.sigma_N=noise_std; % explanatory alias, not used by controller
            p.t_abnormal=timing.abnormal_window;
            for di=1:numel(d_grid)
                deadzone=d_grid(di);
                p.e_deadzone = repmat(deadzone,1,3);
                [delta_q,delta_dq,delta_ddq,~,~,~,~,used,raw,diag] = ...
                    run_admittance_controller('proposed', ...
                    q_ref,dq_ref,p,pre_fb,t_all,ddq_ref);
                assert(all(isfinite(delta_q(:))) && ...
                    all(isfinite(delta_dq(:))) && all(isfinite(delta_ddq(:))), ...
                    'P11:Nonfinite','Nonfinite correction in trial %d.',r);
                % The controller evaluates interaction at samples 1:N-1;
                % finalize_interaction_diag then copies sample N-1 into N.
                % pre_fb.noise_all(N,:) is an unused independent draw, so
                % audit the effective diagnostic trace with the same
                % terminal-sample convention (without changing inputs).
                expected_noise = pre_fb.noise_all;
                expected_noise(end,:) = expected_noise(end-1,:);
                assert(isequal(size(diag.noise),size(expected_noise)) && ...
                    max(abs(diag.noise(:)-expected_noise(:)))<1e-12, ...
                    'P11:Noise', ...
                    'Noise diagnostics differ from the effective shared realization.');
                assert(max(abs(used(:)-soft_threshold(diag.compensated(:),deadzone)))<1e-10, ...
                    'P11:SoftThreshold','Applied dead-zone differs from exact definition.');
                assert(max(abs(diag.compensated(:)-diag.pathological(:)))<1e-9, ...
                    'P11:ExactModel','Nominal compensation is not exact as assumed.');
                metric = compute_metrics(q_ref,q_ref+delta_q,dq_ref+delta_dq, ...
                    ddq_ref+delta_ddq,raw,used,p,diag);
                clean=diag.compensated-diag.noise;
                matched_clean_used=soft_threshold(clean,deadzone);
                passed_noise=used-matched_clean_used;
                mv=zeros(3,9);
                for j=1:3
                    nrm=@(x)sqrt(mean(x.^2));
                    nj=diag.noise(:,j);
                    cj=clean(:,j);
                    ec=diag.compensated(:,j);
                    mask=timing.event_mask_joint(:,j);
                    assert(nnz(mask)>=2,'P11:EventMask','Event mask empty.');
                    clean_event_rms=nrm(cj(mask));
                    event_retention=nrm(matched_clean_used(mask,j)) / ...
                        max(clean_event_rms,eps);
                    noise_rms=nrm(nj);
                    noise_ratio=nan;
                    if noise_rms>0
                        noise_ratio=nrm(passed_noise(:,j))/noise_rms;
                    end
                    mv(j,:)=[metric.RMSE(j),metric.RMS_jerk(j), ...
                        metric.RMS_InteractionTorque(j),nrm(ec),nrm(used(:,j)), ...
                        nrm(passed_noise(:,j)),noise_ratio, ...
                        mean(abs(ec)<=deadzone),event_retention];
                end
                % A mean of three joint-level metrics, not pooled time-sample RMS.
                mv(4,:)=mean(mv,1,'omitnan');
                for j=1:4
                    at=at+1;
                    rows(at).Trial=r;
                    rows(at).DeadzoneNm=deadzone;
                    rows(at).NoiseStdNm=noise_std;
                    rows(at).Joint=joints(j);
                    rows(at).ReferenceRMSEDeg=mv(j,1);
                    rows(at).JerkRMSRadPerS3=mv(j,2);
                    rows(at).RawResidualRMSNm=mv(j,3);
                    rows(at).CompensatedRMSNm=mv(j,4);
                    rows(at).ProcessedInputRMSNm=mv(j,5);
                    rows(at).NoiseTransmissionRMSNm=mv(j,6);
                    rows(at).NoiseTransmissionRatio=mv(j,7);
                    rows(at).DeadzoneZeroFraction=mv(j,8);
                    rows(at).CleanEventSignalRetention=mv(j,9);
                    rows(at).ParameterSeed=seed_audit.parameter_seed;
                    rows(at).NoiseSeed=seed_audit.noise_seed;
                end
            end
        end
        if mod(r,10)==0 || r==1 || r==num_trials
            fprintf('P11 paired trial %d/%d completed (%d rows)\n',r,num_trials,at);
        end
    end
    assert(at==numel(rows),'P11:RecordCount','Incomplete P11 experiment.');
    trial=struct2table(rows);
    assert(height(trial)==num_trials*4*numel(d_grid)*numel(noise_grid), ...
        'P11:RecordCount','Incorrect P11 trial count.');
    trial_file=fullfile(data_dir,'P11_deadzone_noise_trial_data.csv');
    writetable(trial,trial_file);

    % Condition summary, preserving trial as the independent unit.
    srows=struct('DeadzoneNm',{},'NoiseStdNm',{},'Joint',{}, ...
        'Metric',{},'N',{},'Mean',{},'SD',{});
    summary_keys=[metric_keys, "CompensatedRMSNm","ProcessedInputRMSNm", ...
        "NoiseTransmissionRatio","DeadzoneZeroFraction", ...
        "CleanEventSignalRetention"];
    k=0;
    for ni=1:numel(noise_grid)
        for di=1:numel(d_grid)
            for j=1:4
                idx=trial.NoiseStdNm==noise_grid(ni) & ...
                    trial.DeadzoneNm==d_grid(di) & trial.Joint==joints(j);
                assert(nnz(idx)==num_trials,'P11:Pairing','Incomplete condition.');
                for mi=1:numel(summary_keys)
                    all_vals=trial.(char(summary_keys(mi)));
                    vals=all_vals(idx);
                    vals=vals(isfinite(vals));
                    avg=nan; stdval=nan;
                    if ~isempty(vals)
                        avg=mean(vals); stdval=std(vals,0);
                    end
                    k=k+1;
                    srows(k)=struct('DeadzoneNm',d_grid(di), ...
                        'NoiseStdNm',noise_grid(ni),'Joint',joints(j), ...
                        'Metric',summary_keys(mi),'N',numel(vals), ...
                        'Mean',avg,'SD',stdval); %#ok<AGROW>
                end
            end
        end
    end
    summary=struct2table(srows);
    writetable(summary,fullfile(data_dir,'P11_deadzone_noise_group_summary.csv'));

    % Exploratory paired 95% CIs; no significance claims or unadjusted p-values.
    crows=struct('NoiseStdNm',{},'DeadzoneNm',{},'ReferenceDeadzoneNm',{}, ...
        'Joint',{},'Metric',{},'N',{},'MeanDifference',{}, ...
        'SDDifference',{},'CI95Lower',{},'CI95Upper',{});
    k=0;
    for ni=1:numel(noise_grid)
        for di=1:numel(d_grid)
            if d_grid(di)==nominal_d, continue; end
            for j=1:4
                base_ix=trial.NoiseStdNm==noise_grid(ni) & ...
                    trial.DeadzoneNm==nominal_d & trial.Joint==joints(j);
                alt_ix=trial.NoiseStdNm==noise_grid(ni) & ...
                    trial.DeadzoneNm==d_grid(di) & trial.Joint==joints(j);
                base=trial(base_ix,:); alt=trial(alt_ix,:);
                base=sortrows(base,'Trial'); alt=sortrows(alt,'Trial');
                assert(height(base)==num_trials && height(alt)==num_trials && ...
                    isequal(base.Trial,alt.Trial) && ...
                    isequal(base.ParameterSeed,alt.ParameterSeed) && ...
                    isequal(base.NoiseSeed,alt.NoiseSeed), ...
                    'P11:Pairing','Paired trial / seed mismatch.');
                for mi=1:numel(metric_keys)
                    metric_name=char(metric_keys(mi));
                    d=alt.(metric_name)-base.(metric_name);
                    assert(all(isfinite(d)),'P11:NonfiniteContrast', ...
                        'Nonfinite paired metric.');
                    mu=mean(d); sd=std(d,0);
                    % Fixed design has >=2 trials. Existing project already
                    % uses tinv for paired confidence intervals.
                    half=tinv(0.975,num_trials-1)*sd/sqrt(num_trials);
                    k=k+1;
                    crows(k)=struct('NoiseStdNm',noise_grid(ni), ...
                        'DeadzoneNm',d_grid(di),'ReferenceDeadzoneNm',nominal_d, ...
                        'Joint',joints(j),'Metric',metric_keys(mi), ...
                        'N',num_trials,'MeanDifference',mu,'SDDifference',sd, ...
                        'CI95Lower',mu-half,'CI95Upper',mu+half); %#ok<AGROW>
                end
            end
        end
    end
    contrasts=struct2table(crows);
    writetable(contrasts,fullfile(data_dir,'P11_deadzone_noise_paired_differences.csv'));

    % Analytical noise-only (zero underlying signal) suppression probability.
    % This probability is NOT the retained/suppressed fraction of modeled HRI.
    [dd,ss]=ndgrid(d_grid,noise_grid);
    theoretical=zeros(numel(dd),1);
    for k=1:numel(dd)
        if ss(k)==0
            theoretical(k)=1; % degenerate noise-only signal at zero
        else
            theoretical(k)=erf(dd(k)/(sqrt(2)*ss(k)));
        end
    end
    theoretical_table=table(dd(:),ss(:),theoretical, ...
        'VariableNames',{'DeadzoneNm','NoiseStdNm', ...
        'GaussianNoiseOnlyZeroedProbability'});
    writetable(theoretical_table, ...
        fullfile(data_dir,'P11_deadzone_noise_analytical_reference.csv'));

    % Supplementary figure is tentative until manuscript revision is frozen.
    paper_figure_utils('defaults');
    figure_base='FigS3_deadzone_noise_sensitivity';
    [fig,layout]=paper_figure_utils('create',figure_base, ...
        [110 70 1120 790],2,2, ...
        'Dead-zone threshold and simulated feedback-noise sensitivity', ...
        'Visible','off');
    specs={"ReferenceRMSEDeg","JerkRMSRadPerS3", ...
        "RawResidualRMSNm","NoiseTransmissionRMSNm"};
    ylabels={'Mean joint reference RMSE (deg)', ...
        'Mean joint jerk RMS (rad/s^3)', ...
        'Mean joint raw residual RMS (N*m)', ...
        'Matched-state injected-noise effect (N*m RMS)'};
    colors=lines(numel(noise_grid));
    for ip=1:4
        ax=paper_figure_utils('tile',layout,ip); hold(ax,'on');
        for ni=1:numel(noise_grid)
            vals=nan(size(d_grid));
            for di=1:numel(d_grid)
                ix=summary.NoiseStdNm==noise_grid(ni) & ...
                    summary.DeadzoneNm==d_grid(di) & ...
                    summary.Joint=="ThreeJointMean" & ...
                    summary.Metric==specs{ip};
                assert(nnz(ix)==1,'P11:PlotSource','Nonunique source data.');
                vals(di)=summary.Mean(ix);
            end
            plot(ax,d_grid,vals,'o-','Color',colors(ni,:), ...
                'LineWidth',1.6,'MarkerSize',6, ...
                'DisplayName',sprintf('Noise SD = %.2f Nm',noise_grid(ni)));
        end
        threshold_line=xline(ax,nominal_d,':','Baseline D=0.05');
        if isprop(threshold_line,'HandleVisibility')
            set(threshold_line,'HandleVisibility','off');
        end
        xlabel(ax,'Soft-dead-zone threshold D (N*m)');
        ylabel(ax,ylabels{ip});
        set(ax,'XTick',d_grid,'XLim',[-0.005 0.105]);
        grid(ax,'on');
        if ip==1, legend(ax,'Location','best'); end
    end
    paper_figure_utils('export',fig,fig_dir,figure_base);
    close(fig);

    note1="Noise is simulated additive Gaussian feedback perturbation; not measured interaction forces";
    note2="Noise-only zero probability assumes true input is identically zero";
    note3="Matched-state noise difference is an offline diagnostic, not closed-loop transfer";
    note4="ThreeJointMean is unweighted mean of jointwise metrics within each trial";
    note5="100 trials required for paper claim; smaller counts are smoke tests only";
    metadata=table("Reviewer 2 P11",num_trials,string(sim_timestamp), ...
        string(mat2str(d_grid)),string(mat2str(noise_grid)), ...
        "Exact nominal model; ideal command execution; Kmis=1", ...
        note1,note2,note3,note4,note5, ...
        'VariableNames',{'ReviewerItem','NumPairedTrials','RunTimestamp', ...
        'ThresholdsNm','FeedbackNoiseStdNm','Assumptions', ...
        'SignalScope','NoiseOnlyScope','MatchedStateScope', ...
        'Aggregation','PublicationCondition'});
    writetable(metadata,fullfile(data_dir,'P11_deadzone_noise_metadata.csv'));
    manifest=table(string(sim_timestamp),num_trials, ...
        height(trial),height(summary),height(contrasts), ...
        string(trial_file),string(fig_dir),string(data_dir), ...
        'VariableNames',{'RunTimestamp','Trials','TrialRows', ...
        'SummaryRows','PairedContrastRows','TrialCSV','PaperFigureDir','FigureDataDir'});
    writetable(manifest,fullfile(data_dir,'P11_deadzone_noise_manifest.csv'));
    save(fullfile(data_dir,'P11_deadzone_noise_source.mat'), ...
        'trial','summary','contrasts','theoretical_table','metadata', ...
        'manifest','d_grid','noise_grid','nominal_d','num_trials','-v7.3');

    % Append a new line to the existing one-entry master index only AFTER
    % previous evidence export. A targeted P11 run does not forge that index.
    index_file=fullfile('results',sim_timestamp,'all_experiment_output_index.csv');
    if isfile(index_file)
        idx=readtable(index_file,'TextType','string');
        assert(all(ismember({'Experiment','Status','RelativeFolder','KeyEvidence'}, ...
            idx.Properties.VariableNames)),'P11:IndexSchema', ...
            'Unexpected master-index schema.');
        label="Soft-dead-zone and noise sensitivity (R2 P11)";
        assert(~any(idx.Experiment==label), ...
            'P11:DuplicateIndex','P11 experiment already indexed for this run.');
        next=table(label,"COMPLETE",string(data_dir), ...
            "FigS3 PDF/PNG; P11 trial/summary/paired/metadata CSV and MAT", ...
            'VariableNames',{'Experiment','Status','RelativeFolder','KeyEvidence'});
        idx=[idx; next]; %#ok<AGROW>
        writetable(idx,index_file);
    end
    evidence=struct('TrialCSV',trial_file, ...
        'SummaryCSV',fullfile(data_dir,'P11_deadzone_noise_group_summary.csv'), ...
        'PairedCSV',fullfile(data_dir,'P11_deadzone_noise_paired_differences.csv'), ...
        'FigureDir',fig_dir,'FigureDataDir',data_dir, ...
        'Manifest',manifest,'NumTrials',num_trials);
    fprintf('P11 sensitivity evidence exported to %s\n',data_dir);
    fprintf('P11 %d paired trials, %d records; full-paper status=%d\n', ...
        num_trials,height(trial),num_trials==100);
end

function out=soft_threshold(x,d)
    out=sign(x).*max(abs(x)-d,0);
end
