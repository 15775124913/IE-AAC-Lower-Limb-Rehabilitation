function output_dir = run_joint_specific_validation(num_runs)
%RUN_JOINT_SPECIFIC_VALIDATION Reviewer 3 Comment 3, isolated experiment.
%  run_joint_specific_validation(2)   small workflow / identity check
%  run_joint_specific_validation(100) formally paired 100-trial experiment
%
% Tests *deterministically heterogeneous* hip/knee/ankle nominal resistance
% magnitudes, not a signed anatomical torque law, identified patient data,
% or a true-vs-estimated parameter error experiment. Both simulated nominal
% magnitude and controller prediction use the same joint-specific nominal
% coefficients; nonzero true/estimated error is studied separately in Node2.
% The original 20 frozen files are unchanged, including main_simulation.
% All four controllers see identical reference, parameter trial seeds,
% pre-generated noise and events in a given run. 100 Hz, 14 cycles.
%
% The model's nonnegative nominal magnitude is NOT a signed torque and its
% Kv*abs(dq) term is NOT a signed viscous damping torque. The *separate*
% optional execution-layer surrogate uses tau_load=-abs(tau_mag)*tanh(v/eps),
% giving tau_load*v <= 0 for nonnegative load magnitudes. That property of
% the execution surrogate is NOT a passivity proof for the adaptive loop.
%
% No numerical claims should be made until MATLAB executes this function.
    if nargin < 1 || isempty(num_runs), num_runs=100; end
    validateattributes(num_runs,{'numeric'},{'scalar','integer','>=',2});
    base=parameters();
    assert(~base.execution_dynamics.enabled, 'JointSpecific:ExecutionOn', ...
        'Use ideal execution to isolate joint-specific magnitude parameters');
    assert(~isfield(base,'patient_nominal_truth') || isempty(base.patient_nominal_truth), ...
        'JointSpecific:TruthOverride','Do not mix joint specificity with estimated-model bias');

    % Illustrative, positive, predeclared asymmetric triples (Hip/Knee/Ankle).
    % These are sensitivity settings, not anthropometry or clinical estimates.
    % B0 [Nm], Kq [Nm/rad], Kv [Nm*s/rad], Ke [Nm] because phi is unitless.
    base.B0    = [0.60 0.45 0.25];
    base.Kq_fb = [2.50 1.80 1.10];
    base.Kv_fb = [1.50 1.10 0.65];
    base.Ke_fb = [55.0 40.0 25.0];
    names={'B0','Kq_fb','Kv_fb','Ke_fb'};
    for k=1:numel(names)
        v=base.(names{k});
        assert(isnumeric(v) && isequal(size(v),[1 3]) && ...
            all(isfinite(v)) && all(v>0) && numel(unique(v))==3, ...
            'JointSpecific:InvalidParameters','All 3 coefficients must differ in %s',names{k});
    end
    [base,teaching,q_ref,dq_ref,ddq_ref,t_all]=TeachingDataPipeline(base, ...
        'Seed',base.seed,'PlotFlag',false,'SaveFlag',false,'Verbose',false);
    assert(size(q_ref,2)==3 && numel(t_all)==size(q_ref,1));
    assert(max(abs(q_ref(:)-teaching.q_controller_ref(:)))<1e-14);
    assert(max(abs(dq_ref(:)-teaching.dq_controller_ref(:)))<1e-14);
    assert(max(abs(ddq_ref(:)-teaching.ddq_controller_ref(:)))<1e-14);
    modes={'fixed','zhang2018','cao2024','proposed'};
    labels={'Fixed','Zhang-2018-derived','Cao-2024-derived','Proposed'};
    joint_names={'Hip','Knee','Ankle'};
    stamp=datestr(now,'yyyymmdd_HHMMSS');
    output_dir=fullfile('results',['JOINT_SPECIFIC_' stamp]);
    if ~exist(output_dir,'dir'), mkdir(output_dir); end
    diary(fullfile(output_dir,'joint_specific_run_diary.txt')); diary on;
    diary_cleanup=onCleanup(@() diary('off')); %#ok<NASGU>
    try
        fprintf('=== R3C3 joint-specific nominal parameter validation ===\n');
        fprintf('Trials=%d, reference samples=%d, dt=%.5f s\n',num_runs,size(q_ref,1),base.dt);
        for k=1:numel(names)
            fprintf('%s [Hip Knee Ankle]: %s\n',names{k},mat2str(base.(names{k}),5));
        end
        n_joint=num_runs*4*3;
        blank=struct('Run',0,'Controller',string(''),'Joint',string(''), ...
            'ParameterSeed',0,'NoiseSeed',0,'ReferenceRMSEDeg',nan, ...
            'JerkRMSDegPerS3',nan,'RawResidualRMSNm',nan, ...
            'PositionSaturationPct',nan,'CompensatedContaminationRMSNm',nan, ...
            'B0',nan,'Kq',nan,'Kv',nan,'Ke',nan);
        rec=repmat(blank,n_joint,1); pos=0;
        for r=1:num_runs
            [sampled,audit]=mc_trial_utils('sample',base,r);
            assert(~sampled.execution_dynamics.enabled);
            assert(all(sampled.K_mismatch_comp==1));
            [pre_fb,event_timing]=mc_trial_utils('pre_fb',sampled,q_ref,t_all,audit.noise_seed);
            sampled.t_abnormal=event_timing.abnormal_window;
            assert(isequal(sampled.B0,base.B0), 'JointSpecific:B0Modified');
            for v=1:4
                [x,dx,ddx,~,~,~,~,eused,eraw,diag]= ...
                    run_admittance_controller(modes{v},q_ref,dq_ref, ...
                    sampled,pre_fb,t_all,ddq_ref);
                q=q_ref+x; dq=dq_ref+dx; ddq=ddq_ref+ddx;
                assert(all(isfinite(q(:))) && all(isfinite(dq(:))) && ...
                    all(isfinite(ddq(:))),'JointSpecific:NonFinite','Nonfinite state');
                metric=compute_metrics(q_ref,q,dq,ddq,eraw,eused,sampled,diag);
                contamination=diag.compensated-diag.pathological;
                assert(max(abs(contamination(:))) < 1e-8, ...
                    'JointSpecific:NominalIdentity', ...
                    'Exact nominal model plus Kmis=1 must remove magnitude mismatch');
                for j=1:3
                    pos=pos+1;
                    rec(pos).Run=r;
                    rec(pos).Controller=string(labels{v});
                    rec(pos).Joint=string(joint_names{j});
                    rec(pos).ParameterSeed=audit.parameter_seed;
                    rec(pos).NoiseSeed=audit.noise_seed;
                    rec(pos).ReferenceRMSEDeg=metric.RMSE(j);
                    rec(pos).JerkRMSDegPerS3=metric.RMS_jerk(j)*180/pi;
                    rec(pos).RawResidualRMSNm=metric.RMS_InteractionTorque(j);
                    rec(pos).PositionSaturationPct=metric.PositionCorrectionSaturationPct(j);
                    rec(pos).CompensatedContaminationRMSNm= ...
                        sqrt(mean(contamination(:,j).^2));
                    rec(pos).B0=sampled.B0(j);
                    rec(pos).Kq=sampled.Kq_fb(j);
                    rec(pos).Kv=sampled.Kv_fb(j);
                    rec(pos).Ke=sampled.Ke_fb(j);
                end
            end
            fprintf('R3C3 paired trial %d/%d completed, rows=%d\n',r,num_runs,pos);
        end
        assert(pos==n_joint,'JointSpecific:Incomplete','Missing joint records');
        agg=repmat(blank,n_joint/3,1);
        fields={'ReferenceRMSEDeg','JerkRMSDegPerS3', ...
            'RawResidualRMSNm','PositionSaturationPct', ...
            'CompensatedContaminationRMSNm'};
        for k=1:numel(agg)
            idx=3*(k-1)+(1:3); block=rec(idx);
            assert(block(1).Run==block(2).Run && block(2).Run==block(3).Run);
            assert(block(1).Controller==block(2).Controller && ...
                   block(2).Controller==block(3).Controller);
            agg(k)=block(1);
            agg(k).Joint=string('ThreeJointMean');
            for m=1:numel(fields)
                f=fields{m}; agg(k).(f)=mean([block.(f)]);
            end
            agg(k).B0=nan; agg(k).Kq=nan; agg(k).Kv=nan; agg(k).Ke=nan;
        end
        T=struct2table([rec;agg]);
        assert(height(T)==num_runs*16,'JointSpecific:FinalRows');
        writetable(T,fullfile(output_dir,'joint_specific_trial_metrics.csv'));
        Tsum=make_summary(T,labels,[joint_names {'ThreeJointMean'}]);
        writetable(Tsum,fullfile(output_dir,'joint_specific_group_summary.csv'));
        Tp=make_pairwise(T,labels,[joint_names {'ThreeJointMean'}]);
        writetable(Tp,fullfile(output_dir,'joint_specific_pairwise_Holm.csv'));
        M=table(num_runs,size(q_ref,1),base.dt,base.seed, ...
            string('hip-knee-ankle coefficients explicitly different; exact joint-specific estimate and truth'), ...
            string('same MC seeds as Node0; ideal execution; 4 controllers; same exogenous disturbances'), ...
            'VariableNames',{'NumTrials','Samples','DtSeconds','Seed', ...
            'Design','Pairing'});
        writetable(M,fullfile(output_dir,'joint_specific_manifest.csv'));
        P=table(string(joint_names(:)),base.B0(:),base.Kq_fb(:),base.Kv_fb(:),base.Ke_fb(:), ...
            'VariableNames',{'Joint','B0_Nm','Kq_NmPerRad', ...
            'Kv_NmSecondsPerRad','Ke_Nm'});
        writetable(P,fullfile(output_dir,'joint_specific_nominal_parameters.csv'));
        Ta=make_sign_audit(base);
        writetable(Ta,fullfile(output_dir,'joint_specific_sign_convention_audit.csv'));
        assert(all(Ta.PowerW<=1e-12),'JointSpecific:SignedLoadPower');
        save(fullfile(output_dir,'joint_specific_results.mat'), ...
            'T','Tsum','Tp','M','P','Ta','-v7.3');
        fprintf('R3C3 complete; %d trial rows; %d summary rows; %d pairwise rows\n', ...
            height(T),height(Tsum),height(Tp));
        fprintf('All tested execution-proxy load power values <= 0.\n');
    catch ME
        report=getReport(ME,'extended','hyperlinks','off');
        fid=fopen(fullfile(output_dir,'JOINT_SPECIFIC_ERROR_REPORT.txt'),'wt');
        if fid>0, fprintf(fid,'%s\n',report); fclose(fid); end
        fprintf(2,'R3C3 validation failed: %s\n',report);
        rethrow(ME);
    end
end

function Ts=make_summary(T,labels,joints)
    fields={'ReferenceRMSEDeg','JerkRMSDegPerS3','RawResidualRMSNm', ...
        'PositionSaturationPct','CompensatedContaminationRMSNm'};
    rows=struct('Controller',{},'Joint',{},'Metric',{},'N',{},'Mean',{},'SD',{});
    for i=1:numel(labels)
        for j=1:numel(joints)
            ix=T.Controller==string(labels{i}) & T.Joint==string(joints{j});
            for k=1:numel(fields)
                z=T.(fields{k})(ix);
                rows(end+1)=struct('Controller',string(labels{i}), ...
                    'Joint',string(joints{j}),'Metric',string(fields{k}), ...
                    'N',numel(z),'Mean',mean(z),'SD',std(z)); %#ok<AGROW>
            end
        end
    end
    Ts=struct2table(rows);
end

function Tp=make_pairwise(T,labels,joints)
    fields={'ReferenceRMSEDeg','JerkRMSDegPerS3','RawResidualRMSNm'};
    rows=struct('Family',{},'Joint',{},'Metric',{},'First',{}, ...
        'Second',{},'N',{},'MeanDifferenceFirstMinusSecond',{}, ...
        'CI95Low',{},'CI95High',{},'PUnadjusted',{},'PHolm',{},'CohenDz',{});
    for j=1:numel(joints)
        for k=1:numel(fields)
            f=fields{k};
            for c=1:3
                a=sortrows(T(T.Controller==string(labels{4}) & ...
                    T.Joint==string(joints{j}),{'Run',f}),'Run');
                b=sortrows(T(T.Controller==string(labels{c}) & ...
                    T.Joint==string(joints{j}),{'Run',f}),'Run');
                assert(isequal(a.Run,b.Run),'JointSpecific:Unpaired');
                d=a.(f)-b.(f);
                if max(abs(d))<1e-10
                    p=nan; ci=[0;0];
                else
                    [~,p,ci]=ttest(d);
                end
                s=std(d); if s==0, dz=nan; else, dz=mean(d)/s; end
                if j<=3, family='JointSpecific_27'; else, family='ThreeJointMean_9'; end
                rows(end+1)=struct('Family',string(family), ...
                    'Joint',string(joints{j}),'Metric',string(f), ...
                    'First',string(labels{4}),'Second',string(labels{c}), ...
                    'N',numel(d),'MeanDifferenceFirstMinusSecond',mean(d), ...
                    'CI95Low',ci(1),'CI95High',ci(2), ...
                    'PUnadjusted',p,'PHolm',nan,'CohenDz',dz); %#ok<AGROW>
            end
        end
    end
    Tp=struct2table(rows);
    families=unique(Tp.Family);
    for i=1:numel(families)
        idx=find(Tp.Family==families(i) & isfinite(Tp.PUnadjusted));
        [sorted,order]=sort(Tp.PUnadjusted(idx));
        count=numel(idx); current=0;
        for k=1:count
            current=max(current,(count-k+1)*sorted(k));
            Tp.PHolm(idx(order(k)))=min(1,current);
        end
    end
end

function Ta=make_sign_audit(p)
    % The original phenomenological model has no torque sign; this audit
    % applies ONLY to the separate velocity-opposing execution load proxy.
    joint=string({'Hip','Knee','Ankle'});
    velocities=[-2 -1 -.2 0 .2 1 2];
    n=3*numel(velocities);
    Joint=strings(n,1); VelocityRadPerS=zeros(n,1);
    NominalMagnitudeNm=zeros(n,1); SignedExecutionLoadNm=zeros(n,1);
    PowerW=zeros(n,1); row=0;
    epsv=p.execution_dynamics.velocity_sign_eps;
    maxq=[p.q1_max p.q2_max p.q3_max];
    q=0.30; % representative angle, below all three limits
    for j=1:3
        for i=1:numel(velocities)
            v=velocities(i); row=row+1;
            phi=min((abs(q)/maxq(j))^2,1);
            mag=p.B0(j)+p.Kq_fb(j)*abs(q)+p.Kv_fb(j)*abs(v)+p.Ke_fb(j)*phi;
            signed=-abs(mag)*tanh(v/epsv);
            Joint(row)=joint(j);
            VelocityRadPerS(row)=v;
            NominalMagnitudeNm(row)=mag;
            SignedExecutionLoadNm(row)=signed;
            PowerW(row)=signed*v;
        end
    end
    Ta=table(Joint,VelocityRadPerS,NominalMagnitudeNm, ...
        SignedExecutionLoadNm,PowerW);
end
