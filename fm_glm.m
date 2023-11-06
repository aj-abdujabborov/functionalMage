classdef fm_glm < matlab.mixin.Copyable
%FM_GLM Fit a General Linear Model on fMRI data
% Fits fMRI data with a model encapsulated inside an fm_eventList vector
%
% Input properties
%   <fmriData> is a vector of fm_data objects
%   <eventList> is an fm_eventList vector representing the design matrix
%
%   <framework> is the approach taken to analyzing the data. It can be one
%   of three values
%     'findBestHrfs': Perform a Least Squares All analysis on the data
%       for the purpose of identifying the best-fitting HRF. When 'hrf' is
%       set to 'derivative-1' or 'derivative-2', this means collapsing all
%       events into one regressor. Otherwise, this means identifying which
%       HRF from the provided or pre-existing set best fits each voxel by
%       running a separate model for every HRF on all the voxels.
%     'OLS': Perform an Ordinary Least Squares analysis on the data. This
%       is the default and should be used unless you're doing a Least
%       Squares Separate analysis. 
%     'LSS': Least Squares Separate analysis, where each event's magnitude
%       of activity is estimated by running a model where that event gets a
%       dedicated regressor and all other events are collapsed into one or
%       a small number of regressors. "help fm_designMatrix" for more info.
%
%       This paper gives an overview:
%       Arco, J. E., González-García, C., Díaz-Gutiérrez, P., Ramírez, J.,
%       & Ruz, M. (2018). Influence of activation pattern estimates and
%       statistical significance tests in fMRI decoding analysis. Journal
%       of Neuroscience Methods, 308, 248–260
%       
%   <lssEventsToDo> LSS analyses can be time-consuming. With this property,
%     you can specify which events the estimates should be computed for.
%     'lssEventsToDo' is a cell array of column vectors, each cell
%     corresponding to a run. Each column vector should contain the
%     *indices* of the events in the corresponding eventList property for
%     which you need to get an estimate. 
% 
%     For example, if you only want estimates for IDs 4 and 5, each cell of
%     this property could contain: find(any(eventList(runNo).ID == [4 5]))
%
%   <hrf> sets the analytical HRF model. Possible values:
%       'canonical': double gamma canonical HRF from SPM
%       'nsd': the 20 HRFs extracted from the Natural Scenes Dataset
%       'nsd+canonical: NSD HRFs + canonical (21 in total)
%       'derivative-1': model the HRF with the canonical HRF and its
%         temporal derivative
%       'derivative-2': model the HRF with the canonical HRF and its
%         temporal and dispersion derivatives
%       'custom': provide your own HRFs inside 'hrfsMatrix' property
%   <hrfsMatrix> If you set 'hrf' to 'custom', provide a 2D matrix of your
%     own analytical HRFs where each column is an HRF. Do ensure the
%     temporal resolution of the HRFs matches the TR of the design matrix.
%     Note that the design matrix may be downsampled during the analysis if
%     it's TR is higher than the data's. 'hrfsMatrix' should match this
%     downsampled design matrix.
%   <hrfsIdx> When you request that more than one HRF is supplied to
%     analyze the data, this property should be a vector of size [1 x
%     numVoxels] and each value should indicate the column index of
%     'hrfsMatrix' that should be used to analyze that voxel.
%     If you set 'framework' to 'findBestHrfs', this property needs to be
%     empty. Once you perform the analysis with that framework, 'hrfsIdx'
%     will be automatically filled.
%
%   <nuissances> can be one of two values:
%     'baselines': For each run, have a baseline regressor. This is the
%       default.
%     'custom': Custom nuissances. See next property.
%   <nuissancesData> A vector of fm_data objects containing your own custom
%     nuissances. In case you want to regress out slow drifts in the data,
%     with polynomials for example, you'll need to use this.
%
%   <ttestContrasts> A matrix indicating ttest contrasts. Each row is a
%     contrast and columns correspond to the IDs in the supplied eventList.
%     Empty by default.
%   
%   <saveFits> Whether to save the fit time series. False by default.
%
% Output properties
%   <results> is a structure array containing several possible fields:
%     'betas': fm_data object containing the parameter estimates. The ID
%       in eventList corresponds to the row number here.
%     'biasCorrectedBetas': if the derivative approach was used, this
%       contains bias-corrected betas as seen in Lindquist et al. (2008),
%       "Modeling the hemodynamic...".
%     'fits': fm_data object containing the fit time series if 'saveFits'
%       is true. Columns correspond to voxels.
%     'RSqs': fm_data object containig R-squared value for each voxel. 
%     'tstats': fm_data object with t-stat values. Each row
%       corresponds to a contrast specified in 'ttestContrasts'.
%
% Constructors
%   > obj = fm_glm()
%   > obj = fm_glm(fmriData)
%   > obj = fm_glm(fmriData, eventList)
%
% Methods
%   > obj.go() runs the analysis
%   > gt = obj.getGroundTruthInBetasForm(patternPerEvent) is a special
%     method for a ground-truth-neural-activity to betas correlation done
%     with fm_mvpa. 'patternPerEvent' can be supplied with the
%     'totalNeuralActivityPerEvent' property on the fm_simulation object.
%     The output 'gt' is a vector of fm_data objects containing the
%     ground-truth neural patterns in the same format as the betas
%
% Example
%   % Supposed 'dm' is an fm_designMatrix object and 'simulation' is an
%   % fm_simulation object.
%   % Fist identify the best-fitting HRFs
%   glm = fm_glm(simulation.boldTimeSeries, dm.getGlmLSA());
%   glm.hrf = "nsd+canonical";
%   glm.framework = "findbesthrfs";
%   glm.go();
%   disp(glm.results.hrfsIdx) % indices of best-fitting HRFs
%   % Now do an LSA analysis
%   glm.framework = 'ols';
%   glm.saveFits = true;
%   glm.go();
%   plot(glm.results(1).fits.data(:,1)) % plot fits from voxel 1 of the first run
%   % Now do an LSS analysis
%   glm.framework = 'lss';
%   glm.eventList = dm.getGlmLSS2();
%   glm.go();
%   disp(glm.results(1).betas) % betas from the 1st run
%
% Part of package funkyMage. November 2023.
% https://github.com/aj-abdujabborov/funkyMage

% Not implemented
%   <saveDesignMatrices> Whether to save the design matrices. False by
%   default.
% Code notes
% * Derivative works with frameworks 'OLS' and 'findBestHrfs' but not with
% 'LSS'.

    properties
        fmriData (1,:) fm_data;
        eventList (1,:) fm_eventList;

        framework string {matlab.system.mustBeMember(framework, {'findBestHrfs', 'ols', 'lss'})} = "OLS";
        lssEventsToDo;

        hrf (1,1) string {matlab.system.mustBeMember(hrf,...
                   {'nsd', 'nsd+canonical', 'canonical',...
                    'derivative-1', 'derivative-2', 'custom'})} = "canonical";
        hrfsMatrix (:,:) {mustBeFinite};
        hrfsIdx (1,:) {mustBeInteger, mustBePositive};

        nuissances (1,1) string {matlab.system.mustBeMember(nuissances, {'baselines', 'custom'})} = "baselines";
        nuissancesData (1,:) fm_data;

        ttestContrasts {mustBeFinite} = [];

        saveFits (1,1) logical = false;
        saveDesignMatrices (1,1) logical = false;
    end

    properties (SetAccess = protected)
        results; 
    end

    properties (Access = protected)
        numRuns;
        numRegressors;

        dataTR;
        desMatTR;
        bDownsampleDesMat;

        privateLssEventsToDo;

        saveTstats;
    end

    methods
        function obj = fm_glm(fmriData, eventList)
            if nargin > 0
                obj.fmriData = fmriData;
            end
            if nargin > 1
                obj.eventList = eventList;
            end
        end
    end

    methods
        function go(obj)
            obj.checkProperties();
            EL = obj.eventList;

            if lower(obj.framework) == "findbesthrfs"
                if any(contains(obj.hrf, "derivative"))
                    desMat = EL.cat().collapseIntoSuperEvent().computeDesignMatrix(obj.dataTR);
                    obj.results = obj.analyzeWithDerivatives(desMat, diagCat(obj.nuissancesData), cat(obj.fmriData));
                else
                    desMat = EL.cat().openIntoUniqueEvents().computeDesignMatrix(obj.dataTR);
                    obj.results = obj.analyzeWithLibrary(desMat, diagCat(obj.nuissancesData), cat(obj.fmriData));
                    obj.hrfsIdx = obj.results.hrfsIdx;
                end

            elseif lower(obj.framework) == "ols"
                if any(contains(obj.hrf, 'derivative'))
                    obj.results = struct('betas', fm_data.empty,...
                                         'biasCorrectedBetas', fm_data.empty,...
                                         'fits', fm_data.empty,...
                                         'tstats', fm_data.empty,...
                                         'RSqs', []);

                    for i = obj.numRuns:-1:1
                        desMat = EL(i).computeDesignMatrix(obj.dataTR);
                        obj.results(i) = obj.analyzeWithDerivatives(desMat, obj.nuissancesData(i), obj.fmriData(i));
                    end

                else
                    obj.results = struct('betas', fm_data.empty,...
                                         'fits', fm_data.empty,...
                                         'tstats', fm_data.empty,...
                                         'RSqs', []);

                    for i = obj.numRuns:-1:1
                        desMat = EL(i).computeDesignMatrix(obj.dataTR);
                        obj.results(i) = obj.analyzeWithGivenHrfs(desMat, obj.nuissancesData(i), obj.fmriData(i));
                    end
                end

            elseif lower(obj.framework) == "lss"
                obj.results = struct('betas', fm_data.empty);

                for i = obj.numRuns:-1:1
                    obj.results(i) = obj.performLSS(EL(i), obj.nuissancesData(i), obj.fmriData(i), obj.lssEventsToDo{i});
                end
            end
        end

        function betas = getGroundTruthInBetasForm(obj, patternPerEvent)
            arguments
                obj;
                patternPerEvent (1,:) fm_data;
            end

            obj.checkProperties();
            assert(lower(obj.framework) ~= "findbesthrfs",...
                  "getGroundTruthInBetasForm only works with 'LSS' and 'OLS' frameworks");

            if lower(obj.framework) == "lss"
                betas = patternPerEvent;
                return;
            end
            
            EL = obj.eventList;
            for i = obj.numRuns:-1:1
                numEvs = height(EL(i));
                numRegIDs = max(EL(i).ID);

                idx = sub2ind([numEvs, numRegIDs], 1:numEvs, EL(i).ID(:).');
                X = zeros(numEvs, numRegIDs);
                X(idx) = 1;

                betas(i) = fm_data(obj.computeOLS(X, patternPerEvent(i).data));
            end
        end
    end

    methods (Access = private)
        function out = performLSS(obj, EL, nuissancesData, fmriData, lssEventsToDo)
            obj.saveFits = false;
            obj.saveTstats = false;

            out = obj.preallocateOutputStructure(height(EL), fmriData.numTRs, ["betas"]);

            numEventsPerID = sum(EL.ID == unique(EL.ID)', 1);
            for i = lssEventsToDo(:)'
                ELCopy = EL;

                if numEventsPerID(EL.ID(i)) > 1
                    targetEventID = max(EL.ID) + 1;
                    ELCopy.ID(i) = targetEventID;
                else
                    targetEventID = EL.ID(i);
                end

                desMat = ELCopy.computeDesignMatrix(obj.desMatTR);
                tmp = obj.analyzeWithGivenHrfs(desMat, nuissancesData, fmriData);
                out.betas.data(i, :) = tmp.betas.data(targetEventID, :);
            end
        end

        function [out] = analyzeWithLibrary(obj, desMat, nuissancesData, fmriData)
            numReg = width(desMat);
            out = obj.preallocateOutputStructure(numReg, fmriData.numTRs, ["Rsqs", "betas", "fits"]);
            out.hrfsIdx = nan(1, fmriData.numVoxels);
            
            for h = 1:width(obj.hrfsMatrix)
                boldDesMat = fm_glm.conv(desMat, obj.hrfsMatrix(:,h));
                if obj.bDownsampleDesMat
                    boldDesMat = obj.bDownsampleDesMatToFmriResolution(boldDesMat);
                end
                [currentBetas, SSRegression] ...
                    = obj.computeOLS([boldDesMat, nuissancesData.data], fmriData.data);
                [~,            SSNullModel ] ... 
                    = obj.computeOLS(nuissancesData.data, fmriData.data);
                currentRSqs = 1 - (SSRegression ./ SSNullModel);

                betterFitVoxels = currentRSqs > out.RSqs.data;
                
                out.hrfsIdx(betterFitVoxels)       = h;
                out.RSqs.data(betterFitVoxels)     = currentRSqs(betterFitVoxels);
                out.betas.data(:, betterFitVoxels) = currentBetas(1:numReg, betterFitVoxels);
                if obj.saveFits
                    out.fits.data(:, betterFitVoxels) = [boldDesMat, nuissancesData.data] * currentBetas(:, betterFitVoxels);
                end
            end
        end

        function out = analyzeWithGivenHrfs(obj, desMat, nuissancesData, fmriData)
            numReg = width(desMat);
            out = obj.preallocateOutputStructure(numReg, fmriData.numTRs, ["fits", "tstats"]);

            betas = nan(width(desMat)+width(nuissancesData), fmriData.numVoxels);
            SSRegression = nan(1, fmriData.numVoxels);
            SSNullModel = nan(1, fmriData.numVoxels);
            degFreedom = fmriData.numTRs - rank(desMat);

            for h = unique(obj.hrfsIdx)
                currVox = obj.hrfsIdx == h;

                boldDesMat = fm_glm.conv(desMat, obj.hrfsMatrix(:,h));
                if obj.bDownsampleDesMat
                    boldDesMat = obj.bDownsampleDesMatToFmriResolution(boldDesMat);
                end

                X = [boldDesMat, nuissancesData.data];
                [betas(:,currVox), SSRegression(:,currVox)] = obj.computeOLS(...
                    X,...
                    fmriData.data(:,currVox));
                [~,                SSNullModel(:,currVox)]  = obj.computeOLS(...
                    nuissancesData.data,...
                    fmriData.data(:,currVox));

                if obj.saveFits
                    out.fits.data(:,currVox) = X * betas(:,currVox);
                end

                if obj.saveTstats
                    out.tstats(:,currVox) = obj.ttestCalculator(...
                        X, betas(:,currVox), SSRegression(:,currVox), degFreedom);
                    % TODO: you're inverting X twice, one during regression
                    % and one during t-test
                end
            end

            out.RSqs = 1 - (SSRegression ./ SSNullModel);
            out.betas = fm_data(betas(1:width(desMat), :));
        end

        function out = analyzeWithDerivatives(obj, desMat, nuissancesData, fmriData)
            numReg = width(desMat);
            out = obj.preallocateOutputStructure(numReg, fmriData.numTRs, ["tstats"]);

            boldDesMat = {};
            for h = width(obj.hrfsMatrix):-1:1
                boldDesMat{h} = fm_glm.conv(desMat, obj.hrfsMatrix(:,h));
                if obj.bDownsampleDesMat
                    boldDesMat{h} = obj.bDownsampleDesMatToFmriResolution(boldDesMat{h});
                end
            end
            boldDesMat = cat(2, boldDesMat{:});

            [out.betas, SSRegression] = obj.computeOLS(...
                    [boldDesMat, nuissancesData.data],...
                    fmriData.data);
            [~,      SSNullModel] = obj.computeOLS(...
                    nuissancesData.data,...
                    fmriData.data);

            out.RSqs = 1 - (SSRegression ./ SSNullModel);
            if obj.saveFits
                out.fits = fm_data([boldDesMat, nuissancesData.data] * out.betas);
            end

            
            out.biasCorrectedBetas = fm_data(nan(numReg, fmriData.numVoxels));
            for b = 1:numReg
                out.biasCorrectedBetas.data(b,:) = ...
                    sign(out.betas(b,:)) .* ...
                    sqrt(sum(out.betas(b:numReg:width(boldDesMat), :).^2, 1));
                    % Lindquist et al. (2008), "Modeling the hemodynamic..."
            end

            if obj.saveTstats
                out.tstats = fm_data( ...
                    obj.ttestCalculatorForDerivative(biasCorrectedBetas,...
                                                     SSRegression,...
                                                     fmriData.numTRs) ...
                                     );
            end

            out.betas = fm_data(out.betas(1:width(boldDesMat), :));
        end


        function checkProperties(obj)
            assert(~isempty(obj.fmriData), "fmriData property is empty");
            assert(~isempty(obj.eventList), "eventList (design matrix) propety is empty");

            checkFrameworkAndHrfAreCompatible();

            obj.dataTR = obj.fmriData.TR;
            obj.numRuns = length(obj.fmriData);
            obj.numRegressors = max(obj.eventList.cat().ID);

            obj.hrfsMatrix                               = getHrfs(obj.dataTR);
            obj.hrfsIdx                                  = getHrfsIdx(obj.hrfsMatrix);
            obj.lssEventsToDo                            = getLssEventsToDo(obj.lssEventsToDo);
            [obj.desMatTR, obj.bDownsampleDesMat]        = calculateDesMatResolution();
            obj.nuissancesData                           = getNuissances();

            validateTtestContrasts(obj.numRegressors);
            
                        
            %%%
            function checkFrameworkAndHrfAreCompatible()
                if lower(obj.framework) == "lss"
                    assert(~contains(obj.hrf, 'derivative'),...
                           "LSS cannot be implemented with HRF derivative");
                end
            end

            function hrfsMatrix = getHrfs(TR)
                switch lower(obj.hrf)
                    case 'nsd+canonical'
                        hrfsMatrix = [fm_hrf.getAllNsdHrfs(TR), fm_hrf.getCanonicalHrf(TR)];
                    case 'nsd'
                        hrfsMatrix = fm_hrf.getAllNsdHrfs(TR);
                    case 'canonical'
                        hrfsMatrix = fm_hrf.getCanonicalHrf(TR);
                    case 'derivative-1'
                        hrfsMatrix = fm_hrf.getDerivativeHrfs(TR, 1);
                    case 'derivative-2'
                        hrfsMatrix = fm_hrf.getDerivativeHrfs(TR, 2);
                    case 'custom'
                    otherwise
                        error("'%s' for property hrf is invalid", obj.hrf);
                end
            end

            function hrfsIdx = getHrfsIdx(hrfsMatrix)
                hrfsIdx = obj.hrfsIdx;
                if lower(obj.framework) == "findbesthrfs"
                    assert(isempty(obj.hrfsIdx), "hrfsIdx must be empty when framework is 'findBestHrfs'");
                elseif contains(obj.hrf, 'derivative')
                    assert(isempty(obj.hrfsIdx), "hrfsIdx must be empty when the hrf method is 'derivative'");
                elseif width(hrfsMatrix) == 1
                    hrfsIdx = ones(1, obj.fmriData(1).numVoxels);
                else
                    assert(~isempty(obj.hrfsIdx), "When the basis set has more than one HRF, " + ...
                                                  "hrfsIdx must be provided");
                end
            end

            function lssEventsToDo = getLssEventsToDo(lssEventsToDo)
                if lower(obj.framework) ~= "lss"
                    return;
                end

                if isempty(lssEventsToDo)
                    for i = length(obj.eventList):-1:1
                        lssEventsToDo{i}(:,1) = 1:height(obj.eventList(i));
                    end
                    return;
                end
    
                if length(obj.eventList) > 1
                    assert(iscell(lssEventsToDo), ...
                           "For multiple runs, eventList must be a cell array");
                    assert(length(obj.eventList) == length(lssEventsToDo),...
                           "Number of cells in 'lssEventsToDo' must match the number of runs");
                elseif ~iscell(lssEventsToDo)
                    lssEventsToDo = {lssEventsToDo};
                end
    
                for i = length(lssEventsToDo):-1:1
                    assert(width(lssEventsToDo{i}) == 1, "Input should be a column vector");
                    assert(max(lssEventsToDo{i}) <= height(obj.eventList(i)) & ...
                           min(lssEventsToDo{i}) >= 1 & ...
                           isRound(lssEventsToDo{i}), ...
                           "Input should contain indices of events in eventList " + ...
                           "that are of interest");
                end
            end

            function [desMatTR, bDownsampleDesMat] = calculateDesMatResolution()
                timeInfo = cat(1, obj.eventList.Duration, obj.eventList.Onset, obj.dataTR);
                timeInfo = unique(timeInfo);
    
                multFactor = 4;
                decimallessTimeInfo = timeInfo * 10^(multFactor);
                assert(isRound(decimallessTimeInfo),...
                    "The temporal resolution of the design matrix is too fine.");
                
                desMatTR = gcdOfMany(decimallessTimeInfo) / 10^(multFactor);
                bDownsampleDesMat = desMatTR ~= obj.dataTR;
            end

            function validateTtestContrasts(numRegressors)
                if isempty(obj.ttestContrasts), return; end
                assert(width(obj.ttestContrasts) == numRegressors,...
                       "Number of columns in contrasts matrix should be equal " + ...
                       "to the number of event types. Do not include nuissances.");
            end

            function nuissancesData = getNuissances()
                if strcmpi(obj.nuissances, 'baselines')
                    nuissancesData = obj.getBaselineNuissances([obj.fmriData.numTRs]);
                else
                    assert(~isempty(obj.nuissancesData) && ~isempty(obj.nuissancesData.data),...
                           "'nuissancesData' cannot be empty unless 'nuissances' is 'baseline'");
                    for i = 1:obj.numRuns
                        assert(obj.nuissancesData(i).numTRs == obj.fmriData(i).numTRs,...
                               "Number of TRs in nuissancesData does not match the fMRI data");
                    end
                    nuissancesData = obj.nuissancesData;
                end
            end

        end

        function out = preallocateOutputStructure(obj, numReg, numTRs, neededFields)
            numVoxels = obj.fmriData.getNumVoxels();

            if any(strcmpi(neededFields, "RSqs"))
                out.RSqs = fm_data(zeros(1, numVoxels));
            end

            if any(strcmpi(neededFields, "betas"))
                out.betas = fm_data(nan(numReg, numVoxels));
            end

            if any(strcmpi(neededFields, "fits"))
                out.fits = fm_data.empty(1,0);
                if obj.saveFits
                    out.fits = fm_data(nan(numTRs, numVoxels));
                end
            end
            
            if any(strcmpi(neededFields, "tstats"))
                out.tstats = fm_data.empty(1,0);
                if obj.saveTstats
                    out.tstats = fm_data(nan(height(obj.ttestContrasts), numVoxels));
                end
            end
        end

        function data = downsampleDesMatToFmriResolution(obj, data)
            idx = 1 : (obj.dataTR/obj.desMatTR) : height(data); 
            data = data(idx, :);
        end

        function tstats = ttestCalculator(obj, X, betas, SSRegression, degFreedom)
            variance = SSRegression/degFreedom;
            for c = height(obj.ttestContrasts):-1:1
                currContrast = obj.ttestContrasts(c,:) * betas;
                standardError = sqrt(currContrast*inv(X'*X)*currContrast'*variance); %#ok<MINV>
                tstats(c, :) = currContrast ./ standardError;
                % See Appendix A from Handbook of Functional MRI by Poldrack et al.
                
                % In evoked response (1-condition) univariate analyses, this could inflate tstat
            end            
        end

        function tstats = ttestCalculatorForDerivative(obj, betas, SSRegression, numTRs)
            errorTerms = sqrt(SSRegression ./ (numTRs * (numTRs - 1)));        
            for c = height(obj.ttestContrasts):-1:1
                contrastTerms = obj.ttestContrasts(c,:) * betas;
                tstats(c,:) = contrastTerms ./ errorTerms;
            end
            % Calhoun et al. (2003), "fMRI analysis with the general linear model..."
        end
    end

    methods (Static = true)
        function [betas, SS] = computeOLS(X, data)
            valid = ~all(X == 0, 1);
            betas(valid, :) = inv(X(:,valid)' * X(:,valid)) * X(:,valid)' * data; %#ok<MINV>
            betas(~valid, :) = nan;
            SS = sum((data-X*betas).^2, 1);
        end

        function Y = conv(data, HRF)
            arguments
                data {mustBeFinite};
                HRF (:,1) {mustBeFinite};
            end
            Y = conv2(data, HRF);
            Y = Y(1 : end-length(HRF)+1, :);
        end

        function baselines = getBaselineNuissances(numTRsPerRun)
            for i = length(numTRsPerRun):-1:1
                tmp = ones(numTRsPerRun(i), 1);
                baselines(i) = fm_data(tmp);
            end
        end
    end

    methods % Set and get methods
        function set.fmriData(obj, fmriData)
            assert(length([fmriData.TR]) == length(fmriData),...
                   "TR of all the input data should be set to the same value");
            assert(length([fmriData.numVoxels]) == length(fmriData),...
                   "Number of voxels of all the input data should be the same value");
            obj.fmriData = fmriData;
        end
    end
end