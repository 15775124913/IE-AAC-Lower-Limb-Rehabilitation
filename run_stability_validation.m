function evidence = run_stability_validation(params, trace, output_dir)
%RUN_STABILITY_VALIDATION R3-C5 / R2-P10 non-invasive stability evidence audit.
%
% Modes:
%   run_stability_validation
%       Read parameters(), export only the frozen-(B,K) box Jury audit.
%   run_stability_validation(params, trace, output_dir)
%       Also check an EXISTING Proposed controller trajectory, with fields
%       t, Ba, Ka, Ga, delta_q, delta_dq, delta_ddq, e_used.
%
% The audit DOES NOT change controller parameters, commands, interactions,
% random streams, or the motion supervisor. The Schur/Jury test applies to
% each unforced, unsupervised semi-implicit Euler system with FROZEN B/K.
% Pointwise Schur checks for changing B/K do NOT establish switched-system
% stability, ISS, recursive feasibility, passivity, or physical safety.
%
% Read-only evidence: three-joint parameter-box CSV, optionally nominal
% trajectory summary and every-step CSV, and limitation statements.

    if nargin < 1 || isempty(params)
        params = parameters();
    end
    if nargin < 2
        trace = [];
    end
    if nargin < 3 || isempty(output_dir)
        stamp = datestr(now,'yyyymmdd_HHMMSS');
        output_dir = fullfile('results',['stability_audit_' stamp]);
    end
    if exist(output_dir,'dir') ~= 7
        mkdir(output_dir);
    end

    h = params.dt;
    M = to3(params.Ma0,'Ma0');
    B0 = to3(params.Ba0,'Ba0');
    K0 = to3(params.Ka0,'Ka0');
    G0 = to3(params.Ga0,'Ga0');
    validateattributes(h,{'numeric'},{'scalar','finite','positive'},mfilename,'dt');
    assert(all(M>0 & B0>0 & K0>0 & G0>0), ...
        'StabilityAudit:InvalidParameters','M, B0, K0 and G0 must be positive.');

    % Identical numerical projection bounds to proposed_controller_core in
    % run_admittance_controller.m (two-stage projection and quiet scaling).
    Blo = 0.10*B0; Bhi = 4.0*B0;
    Klo = 0.01*K0; Khi = 5.0*K0;
    Glo = 0.08*G0; Ghi = 2.5*G0;
    jointNames = {'Hip';'Knee';'Ankle'};

    r_lo = h*Blo./M; r_hi = h*Bhi./M;
    s_lo = (h*h)*Klo./M; s_hi = (h*h)*Khi./M;
    % Jury/Schur conditions: r>0, s>0 and 4-2r-s>0.
    % Worst r/s for positive rectangular B/K projection ranges:
    j1 = r_lo;
    j2 = s_lo;
    j3 = 4-2*r_hi-s_hi;
    boxPassed = (j1>0) & (j2>0) & (j3>0);
    rhoCorner = zeros(1,3);
    for j=1:3
        rr = [r_lo(j), r_hi(j)]; ss = [s_lo(j), s_hi(j)];
        rhomax = 0;
        for ir=1:2
            for is=1:2
                rhomax = max(rhomax,spectral_radius(rr(ir),ss(is)));
            end
        end
        rhoCorner(j) = rhomax;
    end
    Tbox = table(jointNames,M(:),Blo(:),Bhi(:),Klo(:),Khi(:),Glo(:),Ghi(:), ...
        r_lo(:),s_lo(:),j3(:),rhoCorner(:),boxPassed(:), ...
        'VariableNames',{'Joint','M','B_lower','B_upper','K_lower','K_upper', ...
        'G_lower','G_upper','Minimum_r','Minimum_s','Minimum_4minus2rminuss', ...
        'MaxSpectralRadius_4Corners','FrozenParameterBoxSchurPASS'});
    writetable(Tbox,fullfile(output_dir,'stability_frozen_parameter_boxes.csv'));

    evidence = struct();
    evidence.OutputDir = output_dir;
    evidence.FrozenParameterBox = Tbox;
    evidence.FrozenParameterBoxPASS = all(boxPassed);
    evidence.HasTrajectory = ~isempty(trace);
    if ~evidence.FrozenParameterBoxPASS
        error('StabilityAudit:BoxJuryFailed', ...
            'The frozen-parameter Jury test failed for the selected parameter box.');
    end

    nsteps = 0;
    intervention_count = 0;
    if ~isempty(trace)
        required = {'t','Ba','Ka','Ga','delta_q','delta_dq','delta_ddq','e_used'};
        for k=1:numel(required)
            assert(isfield(trace,required{k}), ...
                'StabilityAudit:MissingTraceField', ...
                'Missing trace field: %s',required{k});
        end
        B = trace.Ba; K = trace.Ka; G = trace.Ga;
        q = trace.delta_q; v = trace.delta_dq; a = trace.delta_ddq;
        u = trace.e_used; t = trace.t(:);
        N = size(q,1);
        assert(N>=2 && size(t,1)==N, 'StabilityAudit:BadLength', ...
            'Expected N>=2 and N time samples.');
        for A = {B,K,G,q,v,a,u}
            assert(isequal(size(A{1}),[N,3]) && all(isfinite(A{1}(:))), ...
                'StabilityAudit:BadArray','All sampled arrays must be finite Nx3.');
        end
        assert(all(isfinite(t)) && max(abs(diff(t)-h))<1e-8, ...
            'StabilityAudit:BadTimeGrid','Time samples must use params.dt.');
        nsteps=N-1;
        % Last row is a terminal copy in the controller; exclude it from
        % parameter-step checks and keep only actual k -> k+1 transitions.
        B=B(1:nsteps,:); K=K(1:nsteps,:); G=G(1:nsteps,:);
        r=h*bsxfun(@rdivide,B,M);
        s=h*h*bsxfun(@rdivide,K,M);
        jury=4-2*r-s;
        pbox = bsxfun(@ge,B,Blo-1e-9) & bsxfun(@le,B,Bhi+1e-9) & ...
               bsxfun(@ge,K,Klo-1e-9) & bsxfun(@le,K,Khi+1e-9) & ...
               bsxfun(@ge,G,Glo-1e-9) & bsxfun(@le,G,Ghi+1e-9);
        pointSchur = (r>0)&(s>0)&(jury>0);
        % Independently check measured state recurrence. These residuals
        % audit discrete integration consistency; they are not a proof.
        predicted_v = v(1:nsteps,:)+h*a(2:N,:);
        predicted_q = q(1:nsteps,:)+h*v(2:N,:);
        err_v = v(2:N,:)-predicted_v;
        err_q = q(2:N,:)-predicted_q;
        desired_a = bsxfun(@rdivide, ...
            G.*u(1:nsteps,:)-B.*v(1:nsteps,:)-K.*q(1:nsteps,:),M);
        applied_a = a(2:N,:);
        correction = applied_a-desired_a;
        intervention = abs(correction)>1e-8;
        intervention_count = sum(intervention(:));

        % Match the controller's accepted limit convention: a shared
        % positive scalar is broadcast to Hip/Knee/Ankle, or an explicit
        % three-joint vector is used unchanged.  Do NOT modify params.
        qlim=motion_limit3(params.delta_q_limit,'delta_q_limit');
        vlim=motion_limit3(params.delta_dq_limit,'delta_dq_limit');
        alim=motion_limit3(params.delta_ddq_limit,'delta_ddq_limit');
        jlim=motion_limit3(params.jerk_limit,'jerk_limit');
        q_ok = bsxfun(@le,abs(q(2:N,:)),qlim+2e-8);
        v_ok = bsxfun(@le,abs(v(2:N,:)),vlim+2e-8);
        a_ok = bsxfun(@le,abs(a(2:N,:)),alim+2e-8);
        j_ok = bsxfun(@le,abs(diff(a,1,1))/h,jlim+2e-8);
        motion_ok = q_ok&v_ok&a_ok&j_ok;

        rho = zeros(nsteps,3);
        for j=1:3
            trA = 2-r(:,j)-s(:,j);
            detA=1-r(:,j);
            disc=complex(trA.^2-4*detA,0);
            sqrtDisc=sqrt(disc);
            rho(:,j)=max(abs((trA+sqrtDisc)/2),abs((trA-sqrtDisc)/2));
        end
        Tsummary=table(jointNames, min(B,[],1)',max(B,[],1)', ...
            min(K,[],1)',max(K,[],1)',min(G,[],1)',max(G,[],1)', ...
            min(jury,[],1)',max(rho,[],1)', ...
            sum(intervention,1)', max(abs(correction),[],1)', ...
            max(abs(err_q),[],1)', max(abs(err_v),[],1)', ...
            all(pbox,1)',all(pointSchur,1)',all(motion_ok,1)', ...
            'VariableNames',{'Joint','MinB','MaxB','MinK','MaxK','MinG','MaxG', ...
            'MinPointwiseJuryMargin','MaxPointwiseFrozenRho', ...
            'SupervisorChangedSteps','MaxSupervisorCorrection_rad_s2', ...
            'MaxQRecurrenceResidual_rad','MaxVRecurrenceResidual_rad_s', ...
            'ParametersInProjectionBox','PointwiseFrozenSchurPASS', ...
            'ObservedMotionLimitsPASS'});
        writetable(Tsummary,fullfile(output_dir,'stability_nominal_trajectory_summary.csv'));

        % Long-form CSV is limited to recorded nominal states, not random
        % Monte Carlo realizations or a complete viability-domain sweep.
        total=nsteps*3;
        SampleIndex=zeros(total,1); Time_s=zeros(total,1);
        Joint=cell(total,1); B_used=zeros(total,1); K_used=zeros(total,1);
        G_used=zeros(total,1); JuryMargin=zeros(total,1);
        PointwiseRho=zeros(total,1); DesiredA=zeros(total,1);
        AppliedA=zeros(total,1); SupervisorCorrection=zeros(total,1);
        DeltaQ=zeros(total,1); DeltaDQ=zeros(total,1);
        DeltaDDQ_prev=zeros(total,1);
        ProjectionOK=false(total,1); PointwiseSchurOK=false(total,1);
        MotionOK=false(total,1);
        for j=1:3
            ix=(j-1)*nsteps+(1:nsteps);
            SampleIndex(ix)=(1:nsteps)'; Time_s(ix)=t(1:nsteps);
            Joint(ix)=repmat(jointNames(j),nsteps,1);
            B_used(ix)=B(:,j); K_used(ix)=K(:,j); G_used(ix)=G(:,j);
            JuryMargin(ix)=jury(:,j); PointwiseRho(ix)=rho(:,j);
            DesiredA(ix)=desired_a(:,j); AppliedA(ix)=applied_a(:,j);
            SupervisorCorrection(ix)=correction(:,j);
            DeltaQ(ix)=q(1:nsteps,j); DeltaDQ(ix)=v(1:nsteps,j);
            DeltaDDQ_prev(ix)=a(1:nsteps,j);
            ProjectionOK(ix)=pbox(:,j); PointwiseSchurOK(ix)=pointSchur(:,j);
            MotionOK(ix)=motion_ok(:,j);
        end
        Tsteps=table(SampleIndex,Time_s,Joint,B_used,K_used,G_used, ...
            JuryMargin,PointwiseRho,DesiredA,AppliedA,SupervisorCorrection, ...
            DeltaQ,DeltaDQ,DeltaDDQ_prev,ProjectionOK,PointwiseSchurOK,MotionOK);
        writetable(Tsteps,fullfile(output_dir,'stability_nominal_per_step_audit.csv'));

        evidence.NominalSummary=Tsummary;
        evidence.RecordedSteps=nsteps;
        evidence.ObservedProjectionPASS=all(pbox(:));
        evidence.PointwiseFrozenSchurPASS=all(pointSchur(:));
        evidence.ObservedMotionLimitsPASS=all(motion_ok(:));
        evidence.MaxQRecurrenceResidual=max(abs(err_q(:)));
        evidence.MaxVRecurrenceResidual=max(abs(err_v(:)));
        evidence.SupervisorChangedSteps=intervention_count;
        assert(evidence.ObservedProjectionPASS && evidence.PointwiseFrozenSchurPASS && ...
            evidence.ObservedMotionLimitsPASS && ...
            evidence.MaxQRecurrenceResidual<=1e-8 && ...
            evidence.MaxVRecurrenceResidual<=1e-8, ...
            'StabilityAudit:ObservedConsistencyFailed', ...
            'At least one parameter-box, Jury, kinematic or motion-limit check failed.');
    end

    fid=fopen(fullfile(output_dir,'stability_scope_and_checks.txt'),'w');
    assert(fid~=-1,'StabilityAudit:OutputFile','Cannot open stability scope file.');
    cleaner=onCleanup(@() fclose(fid));
    fprintf(fid,'R3-C5 / R2-P10 STABILITY EVIDENCE (MATLAB R2024b expected)\n');
    fprintf(fid,'Frozen parameter box Jury checks passed: %d\n',evidence.FrozenParameterBoxPASS);
    fprintf(fid,'Observed nominal trajectory present: %d\n',evidence.HasTrajectory);
    fprintf(fid,'Observed nominal time steps: %d\n',nsteps);
    fprintf(fid,'Observed supervisor correction count (threshold 1e-8 rad/s^2): %d\n',intervention_count);
    fprintf(fid,'Projection bounds: 0.1*Ba0..4*Ba0, 0.01*Ka0..5*Ka0, 0.08*Ga0..2.5*Ga0.\n');
    fprintf(fid,'Semi-implicit Euler: v_next=v+h*a_next; q_next=q+h*v_next.\n');
    fprintf(fid,'Jury: r=h*B/M>0; s=h^2*K/M>0; 4-2*r-s>0.\n');
    fprintf(fid,'PROVEN: each fixed positive B/K pair within the listed projection box is Schur stable for zero input and an INACTIVE supervisor.\n');
    fprintf(fid,'OBSERVED ONLY: parameter-box membership, sampled Jury margins, bounded sampled motion, state recurrence and supervisor intervention counts.\n');
    fprintf(fid,'NOT PROVEN: arbitrary switched/adaptive closed-loop stability, complete ISS, recursive feasibility, virtual/physical port passivity or patient safety.\n');
    fprintf(fid,'An imposed input and the nonlinear supervisor invalidate identification of the running trajectory with the unforced frozen homogeneous system.\n');
    clear cleaner;
    fprintf('R3-C5 / R2-P10 stability evidence saved: %s\n',output_dir);
    fprintf('Frozen parameter box minimum 4-2r-s: [%.4f %.4f %.4f]\n',j3);
    if evidence.HasTrajectory
        fprintf('Nominal observed supervisor corrections (>1e-8): %d\n',intervention_count);
    end
end

function x=to3(x,label)
    assert(isnumeric(x)&&all(isfinite(x(:)))&&numel(x)==3, ...
        'StabilityAudit:VectorSize','%s must contain 3 finite values.',label);
    x=reshape(x,1,3);
end

function x=motion_limit3(x,label)
%MOTION_LIMIT3 Accept the SAME scalar-or-3-vector limits as the controller.
% parameters.m: delta_q_limit=0.08, delta_dq_limit=1.5,
% delta_ddq_limit=5.0 (scalars); jerk_limit=[20 20 15] (vector).
    assert(isnumeric(x) && isreal(x) && isvector(x) && ...
        (isscalar(x) || numel(x)==3) && all(isfinite(x(:))) && ...
        all(x(:)>0), ...
        'StabilityAudit:MotionLimitSize', ...
        '%s must be one positive finite scalar or three positive finite values.',label);
    if isscalar(x)
        x=repmat(x,1,3);
    else
        x=reshape(x,1,3);
    end
end

function rho=spectral_radius(r,s)
    trA=2-r-s;
    detA=1-r;
    d=sqrt(complex(trA^2-4*detA,0));
    rho=max(abs((trA+d)/2),abs((trA-d)/2));
end
