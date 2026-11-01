function outputs = run_execution_dynamics_validation(params, q_ref, dq_ref, ddq_ref, t_all, sim_timestamp)
%RUN_EXECUTION_DYNAMICS_VALIDATION Six-case execution study and paper evidence.
% Called automatically by main_simulation using its unique nominal reference.
% The model is a reduced-order simulated execution surrogate, not a hardware
% model or stability/passivity proof. Ideal nominal/MC results are unchanged.
% Results use the same folder convention as generate_paper_figures:
%   results/<timestamp>/execution_dynamics/ (metrics, audit, traces)
%   results/<timestamp>/paper_figures/      (PDF and PNG)
%   results/<timestamp>/figure_data/        (MAT, CSV, metadata)
    assert(nargin==6,'This experiment is called by main_simulation with six inputs.');
    assert(all(size(q_ref)==size(dq_ref)) && all(size(q_ref)==size(ddq_ref)));
    assert(numel(t_all)==size(q_ref,1));
    [pre_fb, timing] = mc_trial_utils('pre_fb',params,q_ref,t_all,params.seed);
    params.t_abnormal = timing.abnormal_window;
    save_full_mat = false;  % Paper-specific figure MAT files are always written.
cases = struct('name',{},'enable',{},'delay_samples',{},'coupling_fraction',{});
cases(1) = struct('name','IdealCommandBaseline','enable',false, 'delay_samples',0,'coupling_fraction',0);
cases(2) = struct('name','Servo0msNoCoupling', 'enable',true, 'delay_samples',0,'coupling_fraction',0.00);
cases(3) = struct('name','Servo0msCoupling012','enable',true, 'delay_samples',0,'coupling_fraction',0.12);
cases(4) = struct('name','Servo20msCoupling012','enable',true,'delay_samples',2,'coupling_fraction',0.12);
cases(5) = struct('name','Servo50msCoupling012','enable',true,'delay_samples',5,'coupling_fraction',0.12);
cases(6) = struct('name','Servo50msCoupling020','enable',true,'delay_samples',5,'coupling_fraction',0.20);
methods = {'fixed','zhang2018','cao2024','proposed'};
joints = {'Hip','Knee','Ankle'};
outdir = fullfile('results',sim_timestamp,'execution_dynamics');
fig_dir = fullfile('results',sim_timestamp,'paper_figures');
data_dir = fullfile('results',sim_timestamp,'figure_data');
if ~exist(outdir,'dir'), mkdir(outdir); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
if ~exist(data_dir,'dir'), mkdir(data_dir); end

rows = cell(numel(cases)*numel(methods)*3,14);
audit_rows = cell(numel(cases)*numel(methods),5); audit_idx = 0;
records = struct(); ri=0;
for c = 1:numel(cases)
    p = params;
    p.execution_dynamics.enabled = cases(c).enable;
    p.execution_dynamics.delay_samples = cases(c).delay_samples;
    p.execution_dynamics.coupling_fraction = cases(c).coupling_fraction;
    fprintf('\nScenario %s: enabled=%d, delay=%d samples, inertia coupling=%.2f\n', ...
        cases(c).name, cases(c).enable, cases(c).delay_samples, cases(c).coupling_fraction);
    for m = 1:numel(methods)
        [d_q,d_dq,d_ddq,~,~,~,~,~,raw,diagH] = ...
            run_admittance_controller(methods{m},q_ref,dq_ref,p,pre_fb,t_all,ddq_ref);
        q_cmd = q_ref+d_q;
        if cases(c).enable
            assert(isfield(diagH,'execution'));
            ex = diagH.execution;
            q_real = ex.q;
            dq_real = ex.dq;
            ddq_real = ex.ddq;
            sat = ex.torque_saturation;
        else
            q_real = q_cmd;
            dq_real = dq_ref+d_dq;
            ddq_real = ddq_ref+d_ddq;
            sat = false(size(q_ref));
        end
        deviation = q_real-q_ref;
        servo_error = q_real-q_cmd;
        % Safety metrics here are descriptive finite-horizon observations,
        % NOT formal stability or passivity guarantees.
        for j=1:3
            ri=ri+1;
            rows(ri,:) = {cases(c).name,methods{m},joints{j},cases(c).enable, ...
                cases(c).delay_samples*p.dt*1000,cases(c).coupling_fraction, ...
                sqrt(mean(deviation(:,j).^2))*180/pi, ...
                max(abs(deviation(:,j)))*180/pi, ...
                sqrt(mean(servo_error(:,j).^2))*180/pi, ...
                max(abs(servo_error(:,j)))*180/pi, ...
                100*mean(sat(:,j)), ...
                max(abs(dq_real(:,j))), ...
                max(abs(ddq_real(:,j))), ...
                sqrt(mean(raw(:,j).^2))};
        end
        records(c,m).scenario = cases(c).name;
        records(c,m).method = methods{m};
        records(c,m).q_cmd = q_cmd;
        records(c,m).q_executed = q_real;
        records(c,m).dq_executed = dq_real;
        records(c,m).ddq_executed = ddq_real;
        records(c,m).interaction_diag = diagH;
        audit_idx = audit_idx + 1;
        if cases(c).enable
            assert(max(ex.velocity_solver_residual_Nm)<1e-6, ...
                'Implicit velocity solver residual exceeded tolerance.');
            vel_inconsistent = ex.raw_des_dq - ex.des_dq;
            audit_rows(audit_idx,:) = {cases(c).name,methods{m}, ...
                max(ex.velocity_solver_residual_Nm), ...
                max(ex.velocity_solver_iterations), ...
                max(abs(vel_inconsistent(:)))};
        else
            audit_rows(audit_idx,:) = {cases(c).name,methods{m},NaN,NaN,NaN};
        end
        fprintf('  %-12s joint-mean executed reference RMSE = %.3f deg\n', ...
            methods{m},mean(sqrt(mean(deviation.^2,1))*180/pi));
    end
end
T = cell2table(rows,'VariableNames',{'Scenario','Method','Joint','NonidealEnabled', ...
    'Delay_ms','InertiaCouplingFraction','ExecutedRefRMSE_deg','ExecutedRefMax_deg', ...
    'ServoRMSE_deg','ServoMax_deg','TorqueSaturation_pct', ...
    'MaxAbsJointVelocity_rad_s','MaxAbsJointAcceleration_rad_s2','RawResidualRMS_Nm'});
writetable(T,fullfile(outdir,'execution_dynamics_results.csv'));
Audit = cell2table(audit_rows,'VariableNames',{'Scenario','Method', ...
    'MaxSolverResidual_Nm','MaxSolverIterations','MaxRawVelocityMismatch_rad_s'});
writetable(Audit,fullfile(outdir,'execution_layer_numerical_audit.csv'));
if save_full_mat
    save(fullfile(outdir,'execution_dynamics_complete.mat'), ...
        'T','Audit','cases','records','params','pre_fb','t_all','q_ref','dq_ref','ddq_ref','-v7.3');
end

% Compact 10-Hz stress-case traces: enough for diagnostic inspection;
% users do not need to upload the large optional full MAT.
k_sample = 1:10:numel(t_all);
stress_rows = cell(numel(methods)*3*numel(k_sample),13); sr=0;
for m = 1:numel(methods)
    rr = records(end,m);
    ex = rr.interaction_diag.execution;
    for j = 1:3
        for kk = k_sample
            sr = sr+1;
            stress_rows(sr,:) = {rr.scenario,methods{m},joints{j},t_all(kk), ...
                q_ref(kk,j),rr.q_cmd(kk,j),ex.q(kk,j), ...
                ex.raw_des_dq(kk,j),ex.des_dq(kk,j),ex.dq(kk,j), ...
                ex.ddq(kk,j),ex.tau_act_hist(kk,j),ex.tau_hri_hist(kk,j)};
        end
    end
end
Trace = cell2table(stress_rows,'VariableNames',{'Scenario','Method','Joint', ...
    'Time_s','ReferencePosition_rad','CommandPosition_rad', ...
    'ExecutedPosition_rad','RawOuterVelocity_rad_s','ServoTargetVelocity_rad_s', ...
    'ExecutedVelocity_rad_s','ExecutedAcceleration_rad_s2', ...
    'ActuatorTorque_Nm','SignedLoad_Nm'});
writetable(Trace,fullfile(outdir,'execution_stress_signals_10Hz.csv'));


% Preserve exact scenario definitions, parameters, and the input-reference source.
ScenarioConfig = struct2table(cases);
writetable(ScenarioConfig,fullfile(outdir,'execution_scenarios.csv'));
paramNames = {'dt','inertia_diag','viscous_damping','kp','kd', ...
    'torque_limit','actuator_time_constant','velocity_sign_eps'};
fid = fopen(fullfile(outdir,'execution_parameters.txt'),'w');
assert(fid>=0,'Could not create execution parameter manifest.');
fprintf(fid,'Reference: main_simulation TeachingDataPipeline output; 14 cycles, common disturbances.\n');
fprintf(fid,'Simulation only; not identified patient/robot dynamics or a stability proof.\n');
for k=1:numel(paramNames)
    key=paramNames{k};
    if strcmp(key,'dt'), value=params.dt; else, value=params.execution_dynamics.(key); end
    fprintf(fid,'%s: %s\n',key,mat2str(value));
end
fclose(fid);

% --- Fig12: time-domain response for all three joints in the worst case.
methodNames = {'Fixed','Zhang-2018-derived','Cao-2024-derived','Proposed'};
methodColors = lines(4);
worst_command = zeros(numel(t_all),3,4);
worst_executed = zeros(numel(t_all),3,4);
for m=1:4
    worst_command(:,:,m) = records(end,m).q_cmd;
    worst_executed(:,:,m) = records(end,m).q_executed;
end
[fig,layout] = paper_figure_utils('create','Fig12_executed_joint_tracking', ...
    [120 80 1150 940],3,1,'Executed joint trajectories: 50 ms delay and coupling 0.20', ...
    'Visible','off');
for j=1:3
    ax=paper_figure_utils('tile',layout,j);
    plot(ax,t_all,q_ref(:,j)*180/pi,'k--','LineWidth',1.3); hold(ax,'on');
    for m=1:4
        plot(ax,t_all,worst_executed(:,j,m)*180/pi,'Color',methodColors(m,:),'LineWidth',1.0);
    end
    ylabel(ax,[joints{j} ' (deg)']); xlabel(ax,'Time (s)'); grid(ax,'on');
    legend(ax,[{'Reference'},methodNames],'Location','best','FontSize',8);
end
paper_figure_utils('export',fig,fig_dir,'Fig12_executed_joint_tracking');
close(fig);

% --- Fig13: six scenarios x four methods: unweighted three-joint RMSE.
caseNames={cases.name}; meanRMSE=zeros(numel(cases),4);
for c=1:numel(cases)
    for m=1:4
        ix=strcmp(string(T.Scenario),cases(c).name) & strcmp(string(T.Method),methods{m});
        assert(sum(ix)==3,'Each scenario-method must provide exactly three joints.');
        meanRMSE(c,m)=mean(T.ExecutedRefRMSE_deg(ix));
    end
end
[fig,layout]=paper_figure_utils('create','Fig13_execution_scenario_rmse', ...
    [100 100 1200 630],1,1,'Execution-layer reference-deviation RMSE (simulation only)', ...
    'Visible','off');
ax=paper_figure_utils('tile',layout,1);
bar(ax,meanRMSE,'grouped'); grid(ax,'on');
set(ax,'XTick',1:numel(cases),'XTickLabel', ...
    {'Ideal','0 ms / 0','0 ms / .12','20 ms / .12','50 ms / .12','50 ms / .20'});
ylabel(ax,'Three-joint mean RMSE (deg)'); xlabel(ax,'Delay / inertia coupling');
legend(ax,methodNames,'Location','best');
paper_figure_utils('export',fig,fig_dir,'Fig13_execution_scenario_rmse');
close(fig);

% --- Supplementary: ankle/knee high-stress execution traces.
[fig,layout]=paper_figure_utils('create','FigS2_execution_stress_knee_ankle', ...
    [120 80 1150 710],2,1,'Nonideal execution: knee and ankle stress trajectories', ...
    'Visible','off');
for jj=1:2
    j=jj+1; ax=paper_figure_utils('tile',layout,jj);
    plot(ax,t_all,q_ref(:,j)*180/pi,'k--','LineWidth',1.3); hold(ax,'on');
    for m=1:4
        plot(ax,t_all,worst_executed(:,j,m)*180/pi,'Color',methodColors(m,:),'LineWidth',1.0);
    end
    ylabel(ax,[joints{j} ' (deg)']); xlabel(ax,'Time (s)'); grid(ax,'on');
    legend(ax,[{'Reference'},methodNames],'Location','best','FontSize',8);
end
paper_figure_utils('export',fig,fig_dir,'FigS2_execution_stress_knee_ankle');
close(fig);

% Figures' underlying numerical values are archived analogously to
% generate_paper_figures.m, not just as raster images.
trackData=table(t_all(:),'VariableNames',{'Time_s'});
for j=1:3
    trackData.(['Reference_' joints{j} '_deg'])=q_ref(:,j)*180/pi;
    for m=1:4
        trackData.([methods{m} '_' joints{j} '_executed_deg'])= ...
            worst_executed(:,j,m)*180/pi;
    end
end
writetable(trackData,fullfile(data_dir,'Fig12_executed_joint_tracking_data.csv'));
save(fullfile(data_dir,'Fig12_executed_joint_tracking_data.mat'), ...
    't_all','q_ref','worst_command','worst_executed','methodNames','joints','cases');
fig13Table=array2table(meanRMSE,'VariableNames',methods);
fig13Table=addvars(fig13Table,string(caseNames(:)),'Before',1,'NewVariableNames','Scenario');
writetable(fig13Table,fullfile(data_dir,'Fig13_execution_scenario_rmse_data.csv'));
save(fullfile(data_dir,'Fig13_execution_scenario_rmse_data.mat'), ...
    'T','meanRMSE','cases','methods','joints');
% Use the same complete trace for the supplementary figure, maintaining
% an unambiguous mapping between reference, method and joint.
stressFigData=trackData(:,[1,7:16]);
writetable(stressFigData,fullfile(data_dir,'FigS2_execution_stress_knee_ankle_data.csv'));
save(fullfile(data_dir,'FigS2_execution_stress_knee_ankle_data.mat'), ...
    't_all','q_ref','worst_executed','methods','cases');
Meta = table(string('50 ms; inertia coupling 0.20'), string('Simulation only'), ...
    string('Joint angular position is the realized execution state'), ...
    'VariableNames',{'WorstCase','ValidationScope','SignalDefinition'});
writetable(Meta,fullfile(data_dir,'Fig12_execution_metadata.csv'));
writetable(Meta,fullfile(data_dir,'FigS2_execution_metadata.csv'));
writetable(ScenarioConfig,fullfile(data_dir,'Fig13_execution_scenario_metadata.csv'));

assert(all(isfinite(T.ExecutedRefRMSE_deg)), 'Nonfinite executed-state metric');
assert(all(isfinite(T.ServoRMSE_deg)), 'Nonfinite inner-loop tracking metric');
outputs=struct('ResultDir',outdir,'PaperFigureDir',fig_dir,'FigureDataDir',data_dir, ...
    'ScenarioResults',T,'NumericalAudit',Audit,'ScenarioConfig',ScenarioConfig);
fprintf('Execution evidence written to: %s\n',outdir);
fprintf('NOTE: simulation only; no formal stability or passivity claim.\n');
end
