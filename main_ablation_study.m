function R = main_ablation_study(action, varargin)
%MAIN_ABLATION_STUDY Legacy/selective ablation wrapper.
%
% Normal workflow:
%   Run main_simulation.m ONCE. It automatically runs primary + core +
%   mismatch sensitivity + finalize, so this file is not required.
%
% This wrapper is retained only when an ablation stage must be run separately:
%   main_ablation_study('core','NumRuns',100)
%   main_ablation_study('mismatch','NumRuns',100)
%   main_ablation_study('finalize','ExpectedRuns',100,'RequireExactRuns',true)
%   main_ablation_study('all','NumRuns',100)

    if nargin < 1 || isempty(action)
        action = 'core';
    end

    action = lower(char(action));
    if strcmp(action,'ablation_all') || strcmp(action,'ablation')
        action = 'all';
    end

    warning('IEAAC:LegacyEntry', ...
        ['main_ablation_study.m is optional. The standard workflow is to run ' ...
         'main_simulation.m once; it automatically includes the full ablation workflow.']);

    R = ablation_study(action, varargin{:});
end
