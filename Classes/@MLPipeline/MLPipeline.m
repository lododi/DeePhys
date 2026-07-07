classdef MLPipeline
    % MLPIPELINE  Static ML helpers: classifier/regressor creation, CV splits, label pooling.
    %
    %   All methods are static — no instantiation needed.
    %   Access default hyperparameters via MLPipeline.returnDefaultParams().

    methods (Static)

        function params = returnDefaultParams()
            % RETURNDEFAULTPARAMS  Default ML hyperparameters.
            %   Replaces hardcoded values scattered across RecordingGroup methods.
            params.RF.NumCycles = 500;
            params.RF.MinLeafSize = 5;  % ≥5 prevents leaf-level memorisation on typical neuro datasets
            params.RF.NumVariablesToSample = 'auto';  % auto = sqrt(F) for classif, F/3 for regress (Breiman, 2001)
            params.RF.Surrogate = 'on';
            params.RF.Reproducible = true;
            params.RF.Prior = 'empirical';  % 'uniform' gives equal weight to each class regardless of frequency
            params.RF.HyperKFold = [];  % [] = fitcensemble/fitrensemble default (5) for the Bayesian search's internal CV

            params.UMAP.NNeighbors = 100;
            params.UMAP.NNeighborsCulture = 10;
            params.UMAP.MinDist = 1;
            params.UMAP.Spread = 5;
            params.UMAP.SGDTasks = 20;
            params.UMAP.ClusterDetail = 'adaptive';

            params.PCA = struct();
            params.tSNE = struct();
        end

        function [clf, train_acc] = createClassifier(X_train, Y_train, alg, N_hyper, params)
            % CREATECLASSIFIER  Train a classifier with optional hyperparameter optimization.
            %
            % INPUTS:
            %   X_train - (N x F) training feature matrix
            %   Y_train - (N x 1) training labels
            %   alg     - "rf", "svm", "cnb", "knn"
            %   N_hyper - number of hyperparameter optimization evaluations (0 = none)
            %   params  - (optional) struct with RF hyperparameters
            arguments
                X_train {isnumeric}
                Y_train
                alg string = "rf"
                N_hyper (1,1) double = 0
                params struct = MLPipeline.returnDefaultParams()
            end

            if N_hyper > 0
                opt_opts = struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
                    'MaxObjectiveEvaluations', N_hyper, 'ShowPlots', false, 'Verbose', 0);
                % The Bayesian search scores each candidate via its own internal
                % K-fold CV on X_train — a plain random split, NOT grouped by
                % CVGroups/recording like the outer CV in Classifier.classify. That
                % can let a unit's siblings straddle the inner train/validation
                % split. Defaults to fitcensemble/fitcsvm's own default (5-fold)
                % unless params.RF.HyperKFold overrides it (e.g. to reduce cost).
                if isfield(params, 'RF') && isfield(params.RF, 'HyperKFold') && ~isempty(params.RF.HyperKFold)
                    opt_opts.Kfold = params.RF.HyperKFold;
                end
                switch alg
                    case 'svm'
                        clf = fitcsvm(X_train, Y_train, 'Prior', params.RF.Prior, ...
                            'OptimizeHyperparameters', 'all', ...
                            'HyperparameterOptimizationOptions', opt_opts);
                    case 'cnb'
                        clf = fitcnb(X_train, Y_train, 'Prior', params.RF.Prior, ...
                            'OptimizeHyperparameters', 'all', ...
                            'HyperparameterOptimizationOptions', opt_opts);
                    case 'knn'
                        clf = fitcknn(X_train, Y_train, 'Prior', params.RF.Prior, ...
                            'OptimizeHyperparameters', 'all', ...
                            'HyperparameterOptimizationOptions', opt_opts);
                    case 'rf'
                        % NumLearningCycles is fixed (not searched): tree count should
                        % be an explicit memory/time choice (see Classifier opts.NumTrees),
                        % not something bayesopt spends its budget exploring. Surrogate is
                        % likewise fixed so opts.Surrogate keeps working when NHyper > 0 —
                        % previously both were silently ignored here, defeating the RF
                        % memory knobs added for large Unit-level classification.
                        hyperparams = {'MinLeafSize', 'MaxNumSplits', 'SplitCriterion', 'NumVariablesToSample'};
                        t = templateTree('Reproducible', true, 'Surrogate', params.RF.Surrogate);
                        clf = fitcensemble(X_train, Y_train, 'Method', 'Bag', 'Prior', params.RF.Prior, ...
                            'NumLearningCycles', params.RF.NumCycles, ...
                            'OptimizeHyperparameters', hyperparams, 'Learners', t, ...
                            'HyperparameterOptimizationOptions', opt_opts);
                end
                % MinObjective is error rate → convert to accuracy. Guard against a
                % degenerate/aborted search (empty results table) rather than crashing —
                % this has been observed to happen on some folds while others succeed.
                if isprop(clf, 'HyperparameterOptimizationResults') && ...
                        ~isempty(clf.HyperparameterOptimizationResults) && ...
                        ~isempty(clf.HyperparameterOptimizationResults.MinObjective)
                    train_acc = 1 - clf.HyperparameterOptimizationResults.MinObjective;
                else
                    warning('MLPipeline:optimizationResultsMissing', ...
                        'Hyperparameter optimization produced no usable results table for this fold — falling back to OOB/resubstitution accuracy.');
                    if alg == "rf"
                        train_acc = 1 - oobLoss(clf, 'LossFun', 'classiferror');
                    else
                        train_acc = 1 - resubLoss(clf, 'LossFun', 'classiferror');
                    end
                end
            else
                switch alg
                    case 'svm'
                        clf = fitcsvm(X_train, Y_train, 'Prior', params.RF.Prior);
                    case 'cnb'
                        clf = fitcnb(X_train, Y_train, 'Prior', params.RF.Prior);
                    case 'knn'
                        clf = fitcknn(X_train, Y_train, 'Prior', params.RF.Prior);
                    case 'rf'
                        if isnumeric(params.RF.NumVariablesToSample)
                            nvts = params.RF.NumVariablesToSample;
                        elseif params.RF.NumVariablesToSample == "all"
                            nvts = 'all';
                        else  % 'auto' or unrecognized — apply Breiman default
                            F = size(X_train, 2);
                            nvts = max(1, round(sqrt(F)));  % classification default; use F/3 for regression variant
                        end
                        t = templateTree('Surrogate', params.RF.Surrogate, ...
                            'MinLeafSize', params.RF.MinLeafSize, ...
                            'NumVariablesToSample', nvts, ...
                            'Reproducible', params.RF.Reproducible);
                        clf = fitcensemble(X_train, Y_train, 'Method', 'Bag', 'Prior', params.RF.Prior, ...
                            'NumLearningCycles', params.RF.NumCycles, ...
                            'Learners', t, 'Options', statset("UseParallel", true));
                end
                if alg == "rf"
                    train_acc = 1 - oobLoss(clf, 'LossFun', 'classiferror');  % OOB ≈ held-out performance, not true training accuracy
                else
                    train_acc = 1 - resubLoss(clf, 'LossFun', 'classiferror');  % resubstitution; optimistic (~1.0 for SVM/KNN)
                end
            end
        end

        function [mdl, train_r2] = createRegressor(X_train, Y_train, alg, N_hyper)
            % CREATEREGRESSOR  Train a regression model with optional hyperparameter optimization.
            %
            % INPUTS:
            %   X_train - (N x F) training feature matrix
            %   Y_train - (N x 1) numeric training targets
            %   alg     - "rf", "svm", "knn"
            %   N_hyper - Bayesian optimization evaluations (0 = none)
            %
            % OUTPUTS:
            %   mdl      - trained MATLAB regression model
            %   train_r2 - training-set R² (in-bag for RF, resubstitution for SVM/KNN)
            arguments
                X_train {isnumeric}
                Y_train {isnumeric}
                alg     string = "rf"
                N_hyper (1,1) double = 0
            end

            params = MLPipeline.returnDefaultParams();
            opt_opts = struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
                'MaxObjectiveEvaluations', N_hyper, 'ShowPlots', false, 'Verbose', 0);

            if N_hyper > 0
                switch alg
                    case 'svm'
                        mdl = fitrsvm(X_train, Y_train, ...
                            'OptimizeHyperparameters', 'all', ...
                            'HyperparameterOptimizationOptions', opt_opts);
                    case 'knn'
                        mdl = fitrknn(X_train, Y_train, ...
                            'OptimizeHyperparameters', 'all', ...
                            'HyperparameterOptimizationOptions', opt_opts);
                    case 'rf'
                        % NumLearningCycles/Surrogate fixed, not searched — see the
                        % matching fix in createClassifier's 'rf' branch for why.
                        hyperparams = {'MinLeafSize', 'MaxNumSplits', 'NumVariablesToSample'};
                        t = templateTree('Reproducible', true, 'Surrogate', params.RF.Surrogate);
                        mdl = fitrensemble(X_train, Y_train, 'Method', 'Bag', ...
                            'Learners', t, ...
                            'NumLearningCycles', params.RF.NumCycles, ...
                            'OptimizeHyperparameters', hyperparams, ...
                            'HyperparameterOptimizationOptions', opt_opts);
                    otherwise
                        error('MLPipeline:createRegressor', 'Unknown algorithm "%s". Use rf, svm, or knn.', alg);
                end
            else
                switch alg
                    case 'svm'
                        mdl = fitrsvm(X_train, Y_train);
                    case 'knn'
                        mdl = fitrknn(X_train, Y_train);
                    case 'rf'
                        t = templateTree('Surrogate', params.RF.Surrogate, ...
                            'MinLeafSize', params.RF.MinLeafSize, ...
                            'NumVariablesToSample', params.RF.NumVariablesToSample, ...
                            'Reproducible', params.RF.Reproducible);
                        mdl = fitrensemble(X_train, Y_train, 'Method', 'Bag', ...
                            'NumLearningCycles', params.RF.NumCycles, ...
                            'Learners', t, 'Options', statset('UseParallel', true));
                    otherwise
                        error('MLPipeline:createRegressor', 'Unknown algorithm "%s". Use rf, svm, or knn.', alg);
                end
            end

            % Training R² (OOB estimate for RF; resubstitution for SVM/KNN)
            if alg == "rf"
                Y_train_pred = oobPredict(mdl);
            else
                Y_train_pred = predict(mdl, X_train);
            end
            ss_res = sum((Y_train - Y_train_pred).^2);
            ss_tot = sum((Y_train - mean(Y_train)).^2);
            if ss_tot == 0
                train_r2 = NaN;
            else
                train_r2 = 1 - ss_res / ss_tot;
            end
        end

        function [Y_train, Y_test, train_idx, test_idx] = cvSplit(Y, cv, k)
            % CVSPLIT  Extract train/test split for fold k.
            arguments
                Y
                cv cvpartition
                k double
            end
            if iscell(Y)
                Y_train = [Y{cv.training(k)}];
                Y_test  = [Y{cv.test(k)}];
                k_train = cv.training(k);
                train_idx = arrayfun(@(x) ones(size(Y{x})) * k_train(x), 1:length(Y), 'un', 0);
                train_idx = [train_idx{:}];
                test_idx = ~train_idx;
            else
                Y_train = Y(cv.training(k));
                Y_test  = Y(cv.test(k));
                train_idx = cv.training(k);
                test_idx  = cv.test(k);
            end
        end

        function [new_group_idx, new_group_labels] = poolMetadataValues(group_idx, group_labels, classification_val)
            % POOLMETADATAVALUES  Pool metadata values into binary groups.
            arguments
                group_idx double
                group_labels string
                classification_val string
            end
            clf_group_idx = find(contains(group_labels, classification_val));
            new_group_idx = ismember(group_idx, clf_group_idx) * 1;
            new_group_labels(1) = join(group_labels(clf_group_idx), '/');
            clf_group_idx = find(~ismember(group_labels, classification_val));
            new_group_idx(new_group_idx == 0) = 2;
            new_group_labels(2) = join(group_labels(clf_group_idx), '/');
        end

    end
end
