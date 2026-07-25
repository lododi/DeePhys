function classify(ctc)
% CLASSIFY  Top-level classification entry point.
%
% Routes to the appropriate classification method based on
% Parameters.Ensemble.Enabled:
%   true  -> classifyUnitsEnsemble()  (majority vote across multiple seeds)
%   false -> classifyUnits()          (single-seed classification)
%
% Call identifyGroundTruthUnits() and generateTrainLabels() before this.
%
% Sets: ctc.UnitLabels, ctc.UnitConfidence, ctc.UnitGraphConnectivity

arguments
    ctc CellTypeClassifier
end

if ctc.Parameters.Ensemble.Enabled
    fprintf('classify: ensemble mode (seeds: %s)\n', ...
        num2str(ctc.Parameters.Ensemble.Seeds));
    ctc.classifyUnitsEnsemble();
else
    fprintf('classify: single-seed mode (seed: %d)\n', ...
        ctc.Parameters.RNGSeed);
    ctc.classifyUnits();
end

end
