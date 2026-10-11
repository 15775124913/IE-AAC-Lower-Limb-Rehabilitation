function varargout = mc_trial_utils(mode, varargin)
%MC_TRIAL_UTILS Shared Monte-Carlo trial utilities for primary and ablation studies.
% This file centralizes the bounded parameter sampler and pregenerated
% exogenous disturbance construction so the primary Monte Carlo experiment,
% core ablation, and mismatch-sensitivity study use the same trial policy.
%
% Usage:
%   [params_trial, audit] = mc_trial_utils('sample', base_params, run_index);
%   pre_fb = mc_trial_utils('pre_fb', params_trial, q_ref, t_all, noise_seed);

    switch lower(char(mode))
        case 'sample'
            [varargout{1:nargout}] = sample_trial_params(varargin{:});
        case 'pre_fb'
            [varargout{1:nargout}] = build_pre_fb(varargin{:});
        otherwise
            error('mc_trial_utils:UnknownMode','Unknown mode: %s',char(mode));
    end
end

function [params, audit] = sample_trial_params(base_params, run_index)
%SAMPLE_TRIAL_PARAMS Bounded sampling around the nominal configuration.
    assert(isfield(base_params,'mc'), 'base_params.mc does not exist');
    cfg = base_params.mc;

    parameter_seed = double(base_params.seed + cfg.parameter_seed_offset + run_index);
    noise_seed = double(base_params.seed + cfg.noise_seed_offset + run_index);
    stream = RandStream('mt19937ar','Seed',parameter_seed);
    params = base_params;

    params.Ma0 = rel_uniform(stream,base_params.Ma0,cfg.rel_admittance);
    params.Ba0 = rel_uniform(stream,base_params.Ba0,cfg.rel_admittance);
    params.Ka0 = rel_uniform(stream,base_params.Ka0,cfg.rel_admittance);
    params.Ga0 = rel_uniform(stream,base_params.Ga0,cfg.rel_admittance);

    params.Kq_fb = rel_uniform(stream,base_params.Kq_fb,cfg.rel_normal_feedback);
    params.Kv_fb = rel_uniform(stream,base_params.Kv_fb,cfg.rel_normal_feedback);
    params.Ke_fb = rel_uniform(stream,base_params.Ke_fb,cfg.rel_normal_feedback);

    params.spasm_q_hip_th = rel_uniform(stream,base_params.spasm_q_hip_th,cfg.rel_spasm_threshold);
    params.spasm_q_knee_th = rel_uniform(stream,base_params.spasm_q_knee_th,cfg.rel_spasm_threshold);

    params.spasm_c1 = rel_uniform(stream,base_params.spasm_c1,cfg.rel_pathology);
    params.spasm_c2 = rel_uniform(stream,base_params.spasm_c2,cfg.rel_pathology);
    params.spasm_c3 = rel_uniform(stream,base_params.spasm_c3,cfg.rel_pathology);
    params.Csyn = rel_uniform(stream,base_params.Csyn,cfg.rel_pathology);
    params.cfatigue = rel_uniform(stream,base_params.cfatigue,cfg.rel_pathology);
    params.event_amp = rel_uniform(stream,base_params.event_amp,cfg.rel_pathology);
    params.fatigue_slope = rel_uniform(stream,base_params.fatigue_slope,cfg.rel_pathology);

    params.feedback_noise_level = rel_uniform(stream,base_params.feedback_noise_level,cfg.rel_noise);
    params.event_time_offset = cfg.event_time_offset_abs .* (2*rand(stream,1,3)-1);
    params.quiet_dq_th = rel_uniform(stream,base_params.quiet_dq_th,cfg.rel_quiet_threshold);
    params.quiet_ddq_th = rel_uniform(stream,base_params.quiet_ddq_th,cfg.rel_quiet_threshold);

    % Explicitly fixed in the current joint-space Monte Carlo design.
    params.model.L1 = base_params.model.L1;
    params.model.L2 = base_params.model.L2;
    params.B0 = base_params.B0;

    audit = struct();
    audit.parameter_seed = parameter_seed;
    audit.noise_seed = noise_seed;
    audit.distribution = cfg.distribution;
    audit.Ma0 = params.Ma0;
    audit.Ba0 = params.Ba0;
    audit.Ka0 = params.Ka0;
    audit.Ga0 = params.Ga0;
    audit.Kq_fb = params.Kq_fb;
    audit.Kv_fb = params.Kv_fb;
    audit.Ke_fb = params.Ke_fb;
    audit.spasm_q = [params.spasm_q_hip_th, params.spasm_q_knee_th];
    audit.spasm_c = [params.spasm_c1, params.spasm_c2, params.spasm_c3];
    audit.Csyn = params.Csyn;
    audit.cfatigue = params.cfatigue;
    audit.event_amp = params.event_amp;
    audit.event_time_offset = params.event_time_offset;
    audit.fatigue_slope = params.fatigue_slope;
    audit.feedback_noise_level = params.feedback_noise_level;
    audit.quiet_dq_th = params.quiet_dq_th;
    audit.quiet_ddq_th = params.quiet_ddq_th;
end

function [pre_fb, event_timing] = build_pre_fb(params, q_ref, t_all, noise_seed)
%BUILD_PRE_FB Shared pregenerated fatigue/event/noise realization.
    N = size(q_ref,1);
    pre_fb = struct();
    pre_fb.fatigue_trend = linspace(0,1,N)';
    pre_fb.fatigue_slope = params.fatigue_slope;
    [pre_fb.event_all, event_timing] = event_model(t_all, params);
    noise_stream = RandStream('mt19937ar','Seed',double(noise_seed));
    pre_fb.noise_all = params.feedback_noise_level .* randn(noise_stream,N,3);
end

function y = rel_uniform(stream,x,rel_range)
    scale = 1 + rel_range .* (2*rand(stream,size(x))-1);
    y = x .* scale;
end
