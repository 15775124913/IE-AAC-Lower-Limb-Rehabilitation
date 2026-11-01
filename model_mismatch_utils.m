function varargout = model_mismatch_utils(mode, varargin)
%MODEL_MISMATCH_UTILS Independent simulated-patient nominal model support.
%  truth = model_mismatch_utils('resolve',params) reads the optional
%  params.patient_nominal_truth struct (B0,Kq_fb,Kv_fb,Ke_fb). If absent,
%  the plant uses the original same-as-controller nominal parameters.
%
%  [truth,audit] = model_mismatch_utils('draw',params,run_index,severity)
%  deterministically perturbs all four true-model parameter groups using
%  separate seeded uniform signs, independently for Hip/Knee/Ankle.
%
%  This phenomenological positive-magnitude model is NOT a signed rigid-
%  body torque model. These tests assess parametric error, not stability.
    switch lower(char(mode))
        case 'resolve'
            varargout{1} = resolve_truth(varargin{1});
        case 'draw'
            [varargout{1:nargout}] = draw_truth(varargin{:});
        otherwise
            error('ModelMismatch:Mode','Unsupported mode %s',char(mode));
    end
end

function truth = resolve_truth(params)
    names = {'B0','Kq_fb','Kv_fb','Ke_fb'};
    if isfield(params,'patient_nominal_truth') && ~isempty(params.patient_nominal_truth)
        truth = params.patient_nominal_truth;
    else
        truth = struct();
        for k=1:numel(names)
            truth.(names{k}) = params.(names{k});
        end
    end
    for k=1:numel(names)
        f=names{k};
        assert(isfield(truth,f), 'ModelMismatch:MissingField', ...
            'Patient true nominal model lacks %s',f);
        v=truth.(f);
        assert(isnumeric(v) && isreal(v) && numel(v)==3 && ...
               all(isfinite(v(:))) && all(v(:)>0), ...
            'ModelMismatch:InvalidTruth', 'Invalid positive 3-vector %s',f);
        truth.(f)=reshape(double(v),1,3);
    end
end

function [truth,audit] = draw_truth(params,run_index,severity)
    assert(isscalar(run_index) && isnumeric(run_index) && isfinite(run_index) && ...
        run_index>=1 && run_index==round(run_index),'ModelMismatch:BadTrial','Invalid trial');
    assert(isscalar(severity) && isnumeric(severity) && isfinite(severity) && ...
        severity>=0 && severity<1,'ModelMismatch:BadSeverity','Severity must be in [0,1)');
    names={'B0','Kq_fb','Kv_fb','Ke_fb'};
    % 30000 is distinct from the main MC sampling/noise stream offsets.
    seed=double(params.seed+30000+run_index);
    stream=RandStream('mt19937ar','Seed',seed);
    truth=struct(); audit=struct('seed',seed,'severity',severity);
    for k=1:numel(names)
        f=names{k};
        signed_draw=2*rand(stream,1,3)-1; % same signs per trial at every severity
        v=reshape(params.(f),1,3);
        truth.(f)=v.*(1+severity*signed_draw);
        audit.([f '_signed_unit_draw'])=signed_draw;
        audit.([f '_relative_error'])=severity*signed_draw;
    end
    audit.true_model=truth;
    % Run the full type/positivity check before the plant receives the model.
    probe=params;
    probe.patient_nominal_truth=truth;
    truth=resolve_truth(probe);
end
