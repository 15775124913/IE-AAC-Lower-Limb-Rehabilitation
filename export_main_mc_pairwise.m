function outputs = export_main_mc_pairwise(inputCsv, outputDir)
%EXPORT_MAIN_MC_PAIRWISE 100-run paired analysis, no controller changes.
%   outputs = export_main_mc_pairwise()
%   outputs = export_main_mc_pairwise('path/to/MC_trial_level_four_controller_*.csv')
%   outputs = export_main_mc_pairwise(inputCsv, outputDir)
%
%  Analysis unit: the Monte Carlo trial (N=100), NOT the 14001 time samples.
%  Positive comparator-minus-Proposed difference means Proposed is LOWER.
%  All 27 per-joint comparisons form one Holm family; the nine three-joint
%  means form a SEPARATE Holm family. The family choices are prespecified.
%  Uses MATLAB betainc/betaincinv; Statistics Toolbox is not required.
%
    if nargin < 1 || isempty(inputCsv)
        [fn, pn] = uigetfile('MC_trial_level_four_controller_*.csv', ...
            'Select the 100-run primary Monte Carlo trial CSV');
        if isequal(fn,0)
            error('Input selection canceled.');
        end
        inputCsv = fullfile(pn,fn);
    end
    inputCsv = char(inputCsv);
    assert(exist(inputCsv,'file')==2,'CSV does not exist: %s',inputCsv);
    if nargin < 2 || isempty(outputDir)
        outputDir = fullfile(fileparts(inputCsv),'MainMC_PairedStatistics');
    end
    outputDir = char(outputDir);
    if ~exist(outputDir,'dir'), mkdir(outputDir); end

    data = readtable(inputCsv);
    joints = {'Hip','Knee','Ankle'};
    comps = {'Fixed','Zhang2018','Cao2024'};
    metrics = {'RMSE_deg','JerkRMS_deg_s3','RawResidualRMS_Nm'};
    units = {'deg','deg/s^3','Nm'};
    keys = {'Run','Joint'};
    required = keys;
    for m = 1:numel(metrics)
        required{end+1} = ['Proposed_' metrics{m}]; %#ok<AGROW>
        for c = 1:numel(comps)
            required{end+1} = [comps{c} '_' metrics{m}]; %#ok<AGROW>
        end
    end
    absent = setdiff(required,data.Properties.VariableNames);
    assert(isempty(absent),'Missing required column(s): %s',strjoin(absent,', '));
    runs = sort(unique(double(data.Run)));
    assert(isequal(runs(:),(1:100)'), ...
        'Expected exactly 100 unique Monte Carlo runs indexed 1:100.');
    assert(height(data)==300, ...
        'Expected 300 rows (100 trial x 3 joints), got %d.',height(data));

    jRows = struct([]); overallRows = struct([]);
    for m = 1:numel(metrics)
        proposedByJoint = zeros(100,3);
        compByJoint = zeros(100,3,3);
        for j = 1:3
            block = data(strcmp(string(data.Joint),joints{j}),:);
            block = sortrows(block,'Run');
            assert(height(block)==100 && isequal(double(block.Run(:)), runs(:)), ...
                'Missing/duplicate trials in joint %s.',joints{j});
            prop = double(block.(['Proposed_' metrics{m}]));
            assert(all(isfinite(prop)),'Non-finite proposed metrics: %s',joints{j});
            proposedByJoint(:,j) = prop;
            for c = 1:3
                comp = double(block.([comps{c} '_' metrics{m}]));
                assert(all(isfinite(comp)),'Non-finite comparator metrics: %s',comps{c});
                compByJoint(:,j,c) = comp;
                s = onePairedSummary(comp,prop,joints{j},comps{c},metrics{m},units{m});
                if isempty(jRows), jRows = s; else, jRows(end+1) = s; end %#ok<AGROW>
            end
        end
        propMean = mean(proposedByJoint,2);
        for c = 1:3
            compMean = mean(compByJoint(:,:,c),2);
            s = onePairedSummary(compMean,propMean, ...
                'ThreeJointArithmeticMean',comps{c},metrics{m},units{m});
            if isempty(overallRows), overallRows = s; else, overallRows(end+1)=s; end %#ok<AGROW>
        end
    end
    jointResults = holmAdjust(struct2table(jRows));
    overallResults = holmAdjust(struct2table(overallRows));
    assert(height(jointResults)==27 && height(overallResults)==9);

    output27 = fullfile(outputDir,'MainMC_27_joint_pairwise_Holm.csv');
    output9  = fullfile(outputDir,'MainMC_9_threejoint_pairwise_Holm.csv');
    writetable(jointResults,output27);
    writetable(overallResults,output9);
    readme = fullfile(outputDir,'METHODS_AND_PROVENANCE.txt');
    fid = fopen(readme,'w');
    assert(fid >= 0,'Cannot open provenance text file.');
    fprintf(fid,'Dataset: primary MC, expected 100 paired runs.\n');
    fprintf(fid,'Input file: %s\n',inputCsv);
    fprintf(fid,'Unit of analysis: matched Monte Carlo trial.\n');
    fprintf(fid,'Direction: comparator minus Proposed, positive = Proposed lower.\n');
    fprintf(fid,'Inference: two-sided paired t test, t-based 95%% CI, Cohens dz.\n');
    fprintf(fid,'Holm multiplicity: 27 comparisons per-joint; 9 three-joint comparisons in a separate family.\n');
    fprintf(fid,'Joint mean: first arithmetic mean across Hip/Knee/Ankle within each run, then paired inference across trials.\n');
    fprintf(fid,'ImprovementPct_RatioOfMeans: 100 * (mean comparator - mean proposed)/mean comparator.\n');
    fprintf(fid,'No inference on execution-layer stability or physical patient validity is implied.\n');
    fclose(fid);
    outputs = struct('JointCSV',output27,'ThreeJointCSV',output9,'ProvenanceTXT',readme, ...
        'JointTable',jointResults,'ThreeJointTable',overallResults);
    fprintf('Primary MC paired statistics exported.\n');
    fprintf('   27 per-joint tests: %d significant after Holm (directions may differ).\n', ...
        sum(jointResults.Significant_Holm_0p05));
    fprintf('    9 three-joint tests: %d significant after Holm.\n', ...
        sum(overallResults.Significant_Holm_0p05));
    fprintf('CSV folder: %s\n',outputDir);
end

function s = onePairedSummary(comp,proposed,group,comparison,metric,unit)
    n = numel(comp);
    assert(n==numel(proposed) && n>1,'Paired sample size invalid.');
    delta = comp(:) - proposed(:);
    mu = mean(delta);
    sd = std(delta,0);  % sample standard deviation, N-1 denominator
    df = n-1;
    % Two-sided Student-t CI and p-values via base-MATLAB incomplete beta.
    bcrit = betaincinv(0.05,df/2,0.5);
    criticalT = sqrt(df*(1/bcrit-1));
    if sd==0
        tstat = 0; p = double(mu==0);
        if mu~=0, tstat = sign(mu)*Inf; end
        cohendz = 0;
        if mu~=0, cohendz = sign(mu)*Inf; end
        halfWidth = 0;
    else
        tstat = mu/(sd/sqrt(n));
        p = betainc(df/(df+tstat^2),df/2,0.5);
        cohendz = mu/sd;
        halfWidth = criticalT*sd/sqrt(n);
    end
    s = struct( ...
        'Group',string(group),'Comparator',string(comparison), ...
        'Metric',string(metric),'Unit',string(unit),'N',n, ...
        'ProposedMean',mean(proposed), 'ComparatorMean',mean(comp), ...
        'ComparatorMinusProposed_Mean',mu, ...
        'Difference95CI_Low',mu-halfWidth,'Difference95CI_High',mu+halfWidth, ...
        'PairedT',tstat,'P_Raw',p,'CohensDz',cohendz, ...
        'ProposedBetterTrialPct',100*mean(delta>0), ...
        'ProposedWorseTrialPct',100*mean(delta<0), ...
        'ProposedImprovementPct_RatioOfMeans',100*mu/mean(comp), ...
        'MeanDirection',string(directionText(mu)));
end

function label = directionText(d)
    if d>0
        label = 'Proposed lower';
    elseif d<0
        label = 'Comparator lower';
    else
        label = 'Tie';
    end
end

function tableOut = holmAdjust(tableIn)
    tableOut = tableIn;
    m = height(tableIn);
    [pSort, idx] = sort(tableIn.P_Raw);
    adjSort = min(1,cummax(((m:-1:1)').*pSort));
    adj = nan(m,1);
    adj(idx) = adjSort;
    tableOut.P_Holm = adj;
    tableOut.Significant_Holm_0p05 = adj<0.05;
end
