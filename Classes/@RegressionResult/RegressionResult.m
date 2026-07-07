classdef RegressionResult < MLResult
% REGRESSIONRESULT  Result container for one cross-validation fold of regression.
%
%   Stores regression model, predictions, feature importance, and all
%   parameters used — decoupled from RecordingGroup state.
%
%   Inherits shared properties (Mdl, Y_pred, Y_test, objects, train_acc,
%   predImp, Features, Parameters) and summarizeImportance from MLResult.

    properties
        mse_train   % Training mean squared error (resubstitution loss)
    end

    methods

        function r = RegressionResult(varargin)
        % Constructor. Accepts name-value pairs or no arguments.
            r = r@MLResult(varargin{:});   % parent handles shared properties
            if nargin > 0
                for i = 1:2:length(varargin)
                    if strcmp(varargin{i}, 'mse_train')
                        r.mse_train = varargin{i+1};
                    end
                end
            end
        end

        function metrics = computeMetrics(r)
        % COMPUTEMETRICS  Compute MSE, MAE, R², correlation from predictions.
            Y_pred    = r.Y_pred(:);
            Y_test    = r.Y_test(:);
            residuals = Y_pred - Y_test;

            metrics.MSE  = mean(residuals.^2);
            metrics.RMSE = sqrt(metrics.MSE);
            metrics.MAE  = mean(abs(residuals));

            ss_tot = sum((Y_test - mean(Y_test)).^2);
            if ss_tot == 0
                metrics.R2 = NaN;
            else
                metrics.R2 = 1 - sum(residuals.^2) / ss_tot;
            end
            metrics.Correlation = corr(Y_pred, Y_test);
            metrics.N           = length(Y_test);
            metrics.mse_train   = r.mse_train;
        end

    end

    methods (Static)

        function T = results2table(result_array)
        % RESULTS2TABLE  Convert a RegressionResult array to a flat table.
        %
        % Returns a table with one row per test object across all folds:
        %   Fold, Y_pred, Y_test, Residual
            rows = {};
            for k = 1:length(result_array)
                r        = result_array(k);
                n        = length(r.Y_test);
                fold_col = repmat(k, n, 1);
                residual = r.Y_pred(:) - r.Y_test(:);

                rows{k} = table(fold_col, r.Y_pred(:), r.Y_test(:), residual, ...
                    'VariableNames', {'Fold', 'Y_pred', 'Y_test', 'Residual'}); %#ok<AGROW>
            end
            T = vertcat(rows{:});
        end

        function summary = summarizeFolds(result_array)
        % SUMMARIZEFOLDS  Compute per-fold and aggregate regression metrics.
        %
        % Returns struct with: per_fold (table), mean_MSE, mean_R2, mean_Correlation.
        %
        % A fold whose computeMetrics() errors (e.g. a degenerate model from a
        % failed hyperparameter search, or corr() on non-finite predictions)
        % contributes NaN rather than aborting the whole summary.
            n_folds = length(result_array);
            mses = nan(n_folds, 1);
            r2s  = nan(n_folds, 1);
            cors = nan(n_folds, 1);
            for k = 1:n_folds
                try
                    m      = result_array(k).computeMetrics();
                    mses(k) = m.MSE;
                    r2s(k)  = m.R2;
                    cors(k) = m.Correlation;
                catch ME
                    warning('RegressionResult:foldMetricsFailed', ...
                        'Fold %d: computeMetrics() failed (%s) — reporting NaN for this fold.', ...
                        k, ME.message);
                end
            end
            summary.per_fold         = table((1:n_folds)', mses, r2s, cors, ...
                'VariableNames', {'Fold', 'MSE', 'R2', 'Correlation'});
            summary.mean_MSE         = mean(mses, 'omitnan');
            summary.mean_R2          = mean(r2s, 'omitnan');
            summary.mean_Correlation = mean(cors, 'omitnan');
        end

    end
end
