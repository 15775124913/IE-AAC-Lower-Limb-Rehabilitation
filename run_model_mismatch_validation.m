function output_dir = run_model_mismatch_validation(num_runs, severity_levels)
%RUN_MODEL_MISMATCH_VALIDATION Reviewer 3 Comment 4, isolated experiment.
%  run_model_mismatch_validation(2)     : small workflow/identity smoke check
%  run_model_mismatch_validation(100)   : full 100-trial paired sensitivity
%  run_model_mismatch_validation(100,[0 .1 .2]) : explicit levels
%
% Does NOT run or overwrite the frozen 18-module main_simulation workflow.
% Each trial uses the same pre-generated disturbance stream across all
% severities/controllers. All four controllers run K_mis=1, and Proposed
% additionally runs K_mis=0 under identical patient truth and disturbances.
% Three objectives are kept separate: tracking deviation (lower better),
% jerk RMS (lower better), raw interaction-residual RMS (lower better).
% Also report compensated-signal contamination vs known simulated pathology
% as a model-internal mechanism diagnostic, NOT a measured clinical torque.
%
% IMPORTANT: unrun source code does not establish numerical/statistical results.
    if nargin < 1 || isempty(num_runs), num_runs=100; end
    if nargin < 2 || isempty(severity_levels), severity_levels=[0 0.10 0.20]; end
    validateattributes(num_runs,{'numeric'},{'scalar','integer','>=',2});
    validateattributes(severity_levels,{'numeric'},{'vector','real','finite','>=',0,'<',1});
    severity_levels=unique(double(severity_levels(:)'),'stable');
    assert(any(severity_levels==0),'ModelMismatch:NeedZero','Levels must include the exact-model 0 baseline');

    params=parameters();
    assert(~params.execution_dynamics.enabled,'ModelMismatch:Execution', ...
        'Keep the execution layer off to isolate nominal-model parameter uncertainty');
    params.execution_dynamics.enabled=false;
    [params,teaching,q_ref,dq_ref,ddq_ref,t_all] = TeachingDataPipeline(params, ...
        'Seed',params.seed,'PlotFlag',false,'SaveFlag',false,'Verbose',false);
    assert(size(q_ref,2)==3 && numel(t_all)==size(q_ref,1));
    assert(max(abs(q_ref(:)-teaching.q_controller_ref(:)))<1e-14);
    fprintf('\n=== R3C4 independent true-vs-estimated nominal model study ===\n');
    fprintf('Trials: %d | severity levels: %s | sampling: 100 Hz, 14 cycles\n', ...
        num_runs,mat2str(severity_levels));

    stamp=datestr(now,'yyyymmdd_HHMMSS');
    output_dir=fullfile('results',['MODEL_MISMATCH_' stamp]);
    if ~exist(output_dir,'dir'), mkdir(output_dir); end
    diary_file=fullfile(output_dir,'model_mismatch_run_diary.txt');
    diary(diary_file); diary on;
    diary_cleanup=onCleanup(@() diary('off')); %#ok<NASGU>
    try

    variants = {'fixed','zhang2018','cao2024','proposed','proposed'};
    labels = {'Fixed','Zhang-2018-derived','Cao-2024-derived', ...
        'Proposed (Kmis=1)','Proposed (Kmis=0)'};
    comp_levels=[1 1 1 1 0];
    joints={'Hip','Knee','Ankle'};
    n_total=num_runs*numel(severity_levels)*numel(variants)*3;
    blank=struct('Run',0,'SeverityPct',0,'Controller',string(''), ...
        'Kmis',0,'Joint',string(''), ...
        'ReferenceRMSEDeg',nan,'JerkRMSRadPerS3',nan, ...
        'RawResidualRMSNm',nan,'CompensatedContaminationRMSNm',nan, ...
        'ReferenceBiasRMSNm',nan,'TrueMotionMismatchRMSNm',nan, ...
        'EstimatedCommandMismatchRMSNm',nan,'PredictedErrorRMSNm',nan, ...
        'PathologyRMSNm',nan,'CorrectionSaturationPct',nan, ...
        'ParameterSeed',0,'NoiseSeed',0,'TruthSeed',0, ...
        'EstB0',nan,'TrueB0',nan,'EstKq',nan,'TrueKq',nan, ...
        'EstKv',nan,'TrueKv',nan,'EstKe',nan,'TrueKe',nan);
    records=repmat(blank,n_total,1);
    pos=0;
    for r=1:num_runs
        [sampled,seed_audit]=mc_trial_utils('sample',params,r);
        [pre_fb,event_timing]=mc_trial_utils('pre_fb',sampled,q_ref,t_all,seed_audit.noise_seed);
        sampled.t_abnormal=event_timing.abnormal_window;
        for lev=1:numel(severity_levels)
            severity=severity_levels(lev);
            [truth,truth_audit]=model_mismatch_utils('draw',sampled,r,severity);
            % Same nominal reference, patient truth and exogenous disturbance
            % for all controllers and compensation variants in this pair.
            sampled.patient_nominal_truth=truth;
            for v=1:numel(variants)
                p_run=sampled;
                p_run.K_mismatch_comp=comp_levels(v)*ones(1,3);
                [delta_q,delta_dq,delta_ddq,~,~,~,~,eused,eraw,diag] = ...
                    run_admittance_controller(variants{v}, ...
                    q_ref,dq_ref,p_run,pre_fb,t_all,ddq_ref);
                q=q_ref+delta_q; dq=dq_ref+delta_dq; ddq=ddq_ref+delta_ddq;
                assert(all(isfinite(q(:))) && all(isfinite(dq(:))) && all(isfinite(ddq(:))), ...
                    'ModelMismatch:NonfiniteState','Nonfinite trajectory at run %d level %.2f %s', ...
                    r,severity,labels{v});
                met=compute_metrics(q_ref,q,dq,ddq,eraw,eused,p_run,diag);
                contamination=diag.compensated-diag.pathological;
                comp_equation=diag.mismatch + diag.reference_model_bias - ...
                    bsxfun(@times,diag.predicted_mismatch,p_run.K_mismatch_comp);
                assert(max(abs(contamination(:)-comp_equation(:)))<1e-10, ...
                    'ModelMismatch:Identity','Compensation balance mismatch');
                if severity==0 && comp_levels(v)==1
                    assert(max(abs(contamination(:)))<1e-9, ...
                        'ModelMismatch:NominalRegression', ...
                        'Exact model + ideal execution + Kmis=1 should remove nominal component');
                end
                for j=1:3
                    pos=pos+1;
                    records(pos).Run=r;
                    records(pos).SeverityPct=100*severity;
                    records(pos).Controller=string(labels{v});
                    records(pos).Kmis=comp_levels(v);
                    records(pos).Joint=string(joints{j});
                    records(pos).ReferenceRMSEDeg=met.RMSE(j);
                    records(pos).JerkRMSRadPerS3=met.RMS_jerk(j);
                    records(pos).RawResidualRMSNm=met.RMS_InteractionTorque(j);
                    cr=sqrt(mean(contamination(:,j).^2));
                    if cr<1e-10, cr=0; end % ignore algebraic roundoff
                    records(pos).CompensatedContaminationRMSNm=cr;
                    records(pos).ReferenceBiasRMSNm=...
                        sqrt(mean(diag.reference_model_bias(:,j).^2));
                    records(pos).TrueMotionMismatchRMSNm=...
                        sqrt(mean(diag.mismatch(:,j).^2));
                    records(pos).EstimatedCommandMismatchRMSNm=...
                        sqrt(mean(diag.predicted_mismatch(:,j).^2));
                    records(pos).PredictedErrorRMSNm=...
                        sqrt(mean((diag.mismatch(:,j)-diag.predicted_mismatch(:,j)).^2));
                    records(pos).PathologyRMSNm=met.RMS_PathologicalInteraction(j);
                    records(pos).CorrectionSaturationPct=met.PositionCorrectionSaturationPct(j);
                    records(pos).ParameterSeed=seed_audit.parameter_seed;
                    records(pos).NoiseSeed=seed_audit.noise_seed;
                    records(pos).TruthSeed=truth_audit.seed;
                    records(pos).EstB0=sampled.B0(j);
                    records(pos).TrueB0=truth.B0(j);
                    records(pos).EstKq=sampled.Kq_fb(j);
                    records(pos).TrueKq=truth.Kq_fb(j);
                    records(pos).EstKv=sampled.Kv_fb(j);
                    records(pos).TrueKv=truth.Kv_fb(j);
                    records(pos).EstKe=sampled.Ke_fb(j);
                    records(pos).TrueKe=truth.Ke_fb(j);
                end
            end
        end
        fprintf('R3C4 trial %d/%d complete: %d rows\n',r,num_runs,pos);
    end
    assert(pos==n_total,'ModelMismatch:RowCount','Incomplete paired trials');
    % Add unweighted three-joint mean for each paired trial, without
    % pooling the 100 Hz time samples as independent observations.
    agg=repmat(blank,n_total/3,1);
    metric_names={'ReferenceRMSEDeg','JerkRMSRadPerS3','RawResidualRMSNm', ...
        'CompensatedContaminationRMSNm','ReferenceBiasRMSNm', ...
        'TrueMotionMismatchRMSNm','EstimatedCommandMismatchRMSNm', ...
        'PredictedErrorRMSNm','PathologyRMSNm','CorrectionSaturationPct'};
    par_names={'EstB0','TrueB0','EstKq','TrueKq','EstKv','TrueKv','EstKe','TrueKe'};
    for k=1:numel(agg)
        idx=3*(k-1)+(1:3);
        block=records(idx);
        assert(block(1).Run==block(2).Run && block(2).Run==block(3).Run);
        agg(k)=block(1);
        agg(k).Joint=string('ThreeJointMean');
        for m=1:numel(metric_names)
            f=metric_names{m};
            agg(k).(f)=mean([block.(f)]);
        end
        for m=1:numel(par_names)
            agg(k).(par_names{m})=nan;
        end
    end
    T=struct2table([records;agg]);
    writetable(T,fullfile(output_dir,'model_mismatch_trial_joint_metrics.csv'));
    eval_joints=[joints {'ThreeJointMean'}];
    T_summary=summarize_model_mismatch(T,severity_levels,labels,eval_joints);
    writetable(T_summary,fullfile(output_dir,'model_mismatch_group_summary.csv'));
    T_pairs=paired_model_mismatch(T,severity_levels,labels,eval_joints);
    writetable(T_pairs,fullfile(output_dir,'model_mismatch_paired_tests.csv'));
    manifest=table(num_runs,string(mat2str(severity_levels)),size(q_ref,1),params.dt, ...
        string('true coefficients drawn as estimate*(1+severity*U[-1,1]), 12 independent coefficients per run, common across variants'), ...
        string('paired trial; ideal execution; 4 controllers Kmis=1 and Proposed Kmis=0'), ...
        'VariableNames',{'NumTrials','SeverityLevels','Samples','DtSeconds','Sampling','Design'});
    writetable(manifest,fullfile(output_dir,'model_mismatch_manifest.csv'));
    save(fullfile(output_dir,'model_mismatch_results.mat'), ...
        'T','T_summary','T_pairs','manifest','num_runs','severity_levels','-v7.3');
    fprintf('=== R3C4 finished: %s | rows=%d | comparisons=%d ===\n', ...
        output_dir,height(T),height(T_pairs));
    catch ME
        report=getReport(ME,'extended','hyperlinks','off');
        fid=fopen(fullfile(output_dir,'MODEL_MISMATCH_ERROR_REPORT.txt'),'wt');
        if fid>0
            fprintf(fid,'%s\n',report);
            fclose(fid);
        end
        fprintf(2,'R3C4 validation failed: %s\n',report);
        rethrow(ME);
    end
end

function T_sum = summarize_model_mismatch(T,levels,labels,joints)
    metrics={'ReferenceRMSEDeg','JerkRMSRadPerS3','RawResidualRMSNm', ...
             'CompensatedContaminationRMSNm','ReferenceBiasRMSNm', ...
             'CorrectionSaturationPct'};
    rec=struct('SeverityPct',{},'Controller',{},'Joint',{},'Metric',{}, ...
               'N',{},'Mean',{},'SD',{});
    for lv=1:numel(levels)
        for v=1:numel(labels)
            for j=1:numel(joints)
                idx=T.SeverityPct==100*levels(lv) & ...
                    T.Controller==string(labels{v}) & T.Joint==string(joints{j});
                for m=1:numel(metrics)
                    x=T.(metrics{m})(idx);
                    a=struct('SeverityPct',100*levels(lv),'Controller',string(labels{v}), ...
                        'Joint',string(joints{j}),'Metric',string(metrics{m}), ...
                        'N',numel(x),'Mean',mean(x),'SD',std(x));
                    rec(end+1)=a; %#ok<AGROW>
                end
            end
        end
    end
    T_sum=struct2table(rec);
end

function Tp=paired_model_mismatch(T,levels,labels,joints)
    % Negative mean difference means the first controller/variant is lower.
    metrics={'ReferenceRMSEDeg','JerkRMSRadPerS3','RawResidualRMSNm', ...
             'CompensatedContaminationRMSNm'};
    entry=struct('Family',{},'SeverityPct',{},'Joint',{},'Metric',{}, ...
        'First',{},'Second',{},'N',{},'MeanDifferenceFirstMinusSecond',{}, ...
        'CI95Low',{},'CI95High',{},'PUnadjusted',{},'PHolm',{},'CohenDz',{});
    comparisons=[4 5;4 1;4 2;4 3];
    for lv=1:numel(levels)
        for j=1:numel(joints)
            for m=1:numel(metrics)
                field=metrics{m};
                for c=1:size(comparisons,1)
                    a=comparisons(c,1); b=comparisons(c,2);
                    ia=T.SeverityPct==100*levels(lv) & T.Joint==string(joints{j}) & T.Controller==string(labels{a});
                    ib=T.SeverityPct==100*levels(lv) & T.Joint==string(joints{j}) & T.Controller==string(labels{b});
                    A=sortrows(T(ia,{'Run',field}),'Run');
                    B=sortrows(T(ib,{'Run',field}),'Run');
                    assert(isequal(A.Run,B.Run),'ModelMismatch:Unpaired','Trial indices differ');
                    d=A.(field)-B.(field);
                    if max(abs(d))<1e-10
                        p=nan; ci=[0;0]; % no test for algebraic equality
                    else
                        [~,p,ci]=ttest(d);
                    end
                    dsd=std(d);
                    if dsd==0, dz=nan; else, dz=mean(d)/dsd; end
                    family='Proposed_Kmis_On_vs_Off';
                    if c>1, family='Proposed_vs_Comparators'; end
                    entry(end+1)=struct('Family',string(family), ...
                        'SeverityPct',100*levels(lv),'Joint',string(joints{j}), ...
                        'Metric',string(field),'First',string(labels{a}), ...
                        'Second',string(labels{b}),'N',numel(d), ...
                        'MeanDifferenceFirstMinusSecond',mean(d), ...
                        'CI95Low',ci(1),'CI95High',ci(2), ...
                        'PUnadjusted',p,'PHolm',nan,'CohenDz',dz); %#ok<AGROW>
                end
            end
        end
    end
    % Severity-degradation paired tests for the Proposed Kmis=1, keeping
    % exactly matched trial indices and using the 0%-severity condition.
    for lv=1:numel(levels)
        if levels(lv)==0, continue; end
        for j=1:numel(joints)
            for m=1:numel(metrics)
                field=metrics{m};
                i1=T.SeverityPct==100*levels(lv) & T.Joint==string(joints{j}) & ...
                    T.Controller==string(labels{4});
                i0=T.SeverityPct==0 & T.Joint==string(joints{j}) & ...
                    T.Controller==string(labels{4});
                A=sortrows(T(i1,{'Run',field}),'Run');
                B=sortrows(T(i0,{'Run',field}),'Run');
                assert(isequal(A.Run,B.Run),'ModelMismatch:UnpairedSeverity');
                d=A.(field)-B.(field);
                if max(abs(d))<1e-10
                    p=nan; ci=[0;0]; % numerically identical, not significant
                else
                    [~,p,ci]=ttest(d);
                end
                if std(d)==0, dz=nan; else, dz=mean(d)/std(d); end
                entry(end+1)=struct('Family',string('Severity_vs_Zero_ProposedKmis1'), ...
                    'SeverityPct',100*levels(lv),'Joint',string(joints{j}), ...
                    'Metric',string(field),'First',string(labels{4}), ...
                    'Second',string('Same Proposed at 0%'), ...
                    'N',numel(d),'MeanDifferenceFirstMinusSecond',mean(d), ...
                    'CI95Low',ci(1),'CI95High',ci(2), ...
                    'PUnadjusted',p,'PHolm',nan,'CohenDz',dz); %#ok<AGROW>
            end
        end
    end
    Tp=struct2table(entry);
    families=unique(Tp.Family);
    for f=1:numel(families)
        idx=find(Tp.Family==families(f) & isfinite(Tp.PUnadjusted));
        [p_sorted,order]=sort(Tp.PUnadjusted(idx));
        m=numel(idx); running=0;
        for k=1:m
            running=max(running,(m-k+1)*p_sorted(k));
            Tp.PHolm(idx(order(k)))=min(1,running);
        end
    end
end
