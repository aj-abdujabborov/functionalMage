% Guide notes

% Code notes
% * Derivative works with frameworks 'OLS' and 'findBestHrfs' but not with
% 'LSS'.

classdef fm_glm < matlab.mixin.Copyable
    properties
        fmriData (1,:) fm_data;
        eventList (1,:) fm_eventList;

        framework string {matlab.system.mustBeMember(framework, {'findBestHrfs', 'ols', 'lss'})} = "OLS";
        lssEventsToDo;

        hrf (1,1) string {matlab.system.mustBeMember(hrf,...
                   {'nsd', 'nsd+canonical', 'canonical',...
                    'derivative-1', 'derivative-2', 'custom'})} = "canonical";
        hrfsMatrix (:,:) uint32 {mustBeFinite};
        hrfsIdx (1,:) {mustBeInteger, mustBePositive};

        nuissances (1,1) string {matlab.system.mustBeMember(nuissances, {'baselines', 'custom'})} = "baselines";
        nuissancesData (1,:) fm_data;

        ttestContrasts {mustBeFinite};

        saveFits (1,1) logical = true;
        saveDesignMatrices (1,1) logical = true;
            % NOT IMPLEMENTED YET
    end

    properties (SetAccess = protected)
        results; 
    end

    properties (Access = private)
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
        function execute(obj)
            obj.checkProperties();
            EL = obj.eventList;

            if lower(obj.framework) == "findbesthrfs"
                if any(contains(obj.hrf, "derivative"))
                    desMat = EL.cat().collapseIntoSuperEvent().computeDesignMatrix(obj.dataTR);
                    obj.results = obj.analyzeWithDerivatives(desMat, diagCat(obj.nuissancesData), cat(obj.fmriData));
                else
                    desMat = EL.cat().openIntoUniqueEvents().computeDesignMatrix(obj.dataTR);
                    obj.results = obj.analyzeWithLibrary(desMat, diagCat(obj.nuissancesData), cat(obj.fmriData));
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

        function out = performLSS(obj, EL, nuissancesData, fmriData, lssEventsToDo)
            obj.saveFits = false;
            obj.saveTstats = false;

            out = obj.preallocateOutputStructure(height(EL), fmriData.numTRs, ["betas"]);

            numEventsPerID = sum(EL.ID == unique(EL.ID)', 1);
            for i = lssEventsToDo
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

    methods (Access = private)
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
                bDownsampleDesMat = obj.desMatTR ~= obj.dataTR;
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
end