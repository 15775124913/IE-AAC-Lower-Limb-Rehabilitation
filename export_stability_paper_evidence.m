function outputs = export_stability_paper_evidence(evidence, nominal_result_dir)
%EXPORT_STABILITY_PAPER_EVIDENCE Archive validated stability tables for paper.
% Called only by main_simulation.m after run_stability_validation.
% Stable naming and common storage with generate_paper_figures.m:
%   results/<config_date>/stability_validation/ : complete raw audit
%   results/<config_date>/figure_data/           : paper table CSV + metadata
% No new simulation, graph, control-law change or additional entry point.
% The fixed-parameter Jury certificate is NOT a proof for the switched
% adaptive closed loop; sampled quantities are observed diagnostics only.
    assert(isstruct(evidence) && isfield(evidence,'OutputDir') && ...
        isfield(evidence,'FrozenParameterBoxPASS') && ...
        logical(evidence.FrozenParameterBoxPASS), ...
        'StabilityPaper:Incomplete','Stability audit is not complete.');
    assert(ischar(nominal_result_dir) || isstring(nominal_result_dir), ...
        'StabilityPaper:Folder','Invalid nominal folder.');
    nominal_result_dir=char(nominal_result_dir);
    audit_dir=char(evidence.OutputDir);
    assert(isfolder(audit_dir) && isfolder(nominal_result_dir), ...
        'StabilityPaper:Folder','Stability source folder is missing.');
    fig_data=fullfile(nominal_result_dir,'figure_data');
    if ~exist(fig_data,'dir'), mkdir(fig_data); end
    source_box=fullfile(audit_dir,'stability_frozen_parameter_boxes.csv');
    source_nom=fullfile(audit_dir,'stability_nominal_trajectory_summary.csv');
    assert(isfile(source_box) && isfile(source_nom), ...
        'StabilityPaper:MissingCSV','Stability audit CSV is missing.');
    box=readtable(source_box);
    nom=readtable(source_nom);
    assert(height(box)==3 && height(nom)==3 && ...
        all(box.FrozenParameterBoxSchurPASS) && ...
        all(nom.ParametersInProjectionBox) && ...
        all(nom.PointwiseFrozenSchurPASS) && ...
        all(nom.ObservedMotionLimitsPASS), ...
        'StabilityPaper:Checks','Cannot export failed/incomplete evidence.');
    out_box=fullfile(fig_data,'Section5_3_1_stability_frozen_box_source.csv');
    out_nom=fullfile(fig_data,'Section5_3_1_stability_nominal_source.csv');
    writetable(box,out_box);
    writetable(nom,out_nom);
    % Data lineage and logical scope accompany the table source data.
    ManuscriptSection=["Section 5.3.1";"Section 5.3.1"];
    SourceFile=[string(source_box);string(source_nom)];
    PaperSource=[string(out_box);string(out_nom)];
    DataType=["Analytical frozen-parameter box";"Observed nominal Proposed trajectory"];
    SampleScope=["Hip/Knee/Ankle parameter box, no input, inactive supervisor"; ...
        "14,000 transitions; sampled switched B/K/G; not a stability proof"];
    EvidenceScope=["Each fixed parameter pair within the box is Schur stable"; ...
        "Parameter bounds, pointwise frozen Jury, motion and recurrence audit only"];
    Tmeta=table(ManuscriptSection,SourceFile,PaperSource,DataType,SampleScope,EvidenceScope);
    metadata=fullfile(fig_data,'Section5_3_1_stability_metadata.csv');
    writetable(Tmeta,metadata);
    outputs=struct('FigureDataDir',fig_data,'BoxTableCSV',out_box, ...
        'NominalTableCSV',out_nom,'MetadataCSV',metadata);
    fprintf('Stability Section 5.3.1 source CSV and metadata archived: %s\n',fig_data);
end
