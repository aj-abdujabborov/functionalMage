
classdef fm_glm < matlab.mixin.Copyable
    properties
        fmriData (1,:) fm_data;
        eventList (1,:) fm_eventList;

        framework {matlab.system.mustBeMember(framework, {'ols', 'lss'})} = 'OLS';
        lssEventsToDo;

        hrf {matlab.system.mustBeMember(hrf,...
            {'nsd', 'nsd+canonical', 'canonical',...
            'derivative-1', 'derivative-2', 'custom'})} = 'canonical';
        hrfsMatrix;
        hrfsIdx (1,:) {mustBeInteger, mustBePositive};

        nuissances {matlab.system.mustBeMember(nuissances, {'baselines', 'custom'})} = 'baselines';
        nuissancesData;

        ttestContrasts {mustBeFinite};

        findBestHrfsWithSuperEvent = true;

        saveFits (1,1) logical;
        saveRsqs (1,1) logical;
        saveDesignMatrices (1,1) logical;
    end

    properties (SetAccess = protected)
        betas;
        fits;
        tstats;
    end

    properties (Access = private)
        numVoxels;
        numRuns;
        numRegressors;

        dataTR;
        desMatTR;
        downsampleDesMat;

        privateLssEventsToDo;

        bTtestContrasts;
    end

    methods
        function obj = fm_glm(fmriData, eventList, varargin)
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
            obj.runPreparations();
            EL = obj.eventList;

            if obj.findBestHrfsWithSuperEvent
                switch lower(obj.hrf)
                    case {'nsd+canonical', 'nsd', 'custom'}
                        desMat = EL.cat().openIntoUniqueEvents().computeDesignMatrix(obj.dataTR);
                        [obj.hrfsIdx] = obj.performLibrary(desMat, cat(obj.nuissancesData), cat(obj.fmriData));
                        % TODO: make computeDesignMatrix return an fm_data
                    case {'derivative-1', 'derivative-2'}
                        desMat = EL.cat().collapseIntoSuperEvent().computeDesignMatrix(obj.dataTR);
                        obj.doHrfDerivative(desMat, cat(obj.nuissancesData), cat(obj.fmriData));
                end
            else
                switch lower(obj.framework)
                    % TODO: you should be able to do derivative without collapsing
                    case 'ols'
                        for i = obj.numRuns:-1:1
                            desMat = EL(i).computeDesignMatrix(obj.dataTR);
                            obj.performGlmWithGivenHRFs(desMat, obj.nuissancesData(i), obj.fmriData(i));
                        end
                    case 'lss'
                        for i = obj.numRuns:-1:1
                            obj.doLSS(EL(i), obj.nuissancesData(i), obj.fmriData(i));
                        end
                end
            end
        end

        function [betas, fits] = doLSS(obj, EL, nuissancesData, fmriData)
            betas = nan(height(EL), fmriData.numVoxels);

            numEventsPerID = sum(EL.ID == unique(EL.ID)', 1);
            for i = obj.lssEventsToDo
                ELCopy = EL;

                if numEventsPerID(EL.ID(i)) > 1
                    targetEventID = max(EL.ID) + 1;
                    ELCopy.ID(i) = targetEventID;
                else
                    targetEventID = EL.ID(i);
                end

                desMat = ELCopy.computeDesignMatrix(obj.desMatTR);
                [betasTmp, fitTmp, rsqsTmp] = obj.performGlmWithGivenHRFs(desMat, nuissancesData, fmriData);
                betas(i, :) = betasTmp(targetEventID, :);
            end


            if length(obj.lssEventsToDo) == height(EL.ID)
                fits = EL.openIntoUniqueEvents().computeDesignMatrix(obj.desMatTR) * betas;
                % TODO: you forgot to convolve before you got the fits
            end
        end

        function [bestHrfsIdx, highestRsqs] = performLibrary(obj, desMat, nuissancesData, fmriData)
            % TODO: what if we have no nuissance data? would be
            % incompatible with a size check

            highestRsqs = zeros(1, fmriData.numVoxels);
            bestHrfsIdx = nan(1, fmriData.numVoxels);
            for h = 1:width(obj.hrfsMatrix)
                boldDesMat = fm_glm.conv(desMat, obj.hrfsMatrix(:,h));
                if obj.downsampleDesMat
                    boldDesMat = obj.downsampleDesMatToFmriResolution(boldDesMat);
                end
                [~, SSRegression] = obj.computeOLS([boldDesMat, nuissancesData.data], fmriData.data);
                [~, SSNullModel]  = obj.computeOLS(nuissancesData.data, fmriData.data);
                RSqs = 1 - (SSRegression ./ SSNullModel);

                betterFitVoxels = RSqs > highestRsqs;
                highestRsqs(betterFitVoxels) = RSqs(betterFitVoxels);
                bestHrfsIdx(betterFitVoxels) = h;
            end
        end

        function [betas, fit, RSqs] = performGlmWithGivenHRFs(obj, desMat, nuissancesData, fmriData)
            % For univariate analysis, we'll probably concatenate before we
            % call this function.
            betas = nan(width(desMat)+width(nuissancesData), fmriData.numVoxels);
            SSRegression = nan(1, fmriData.numVoxels);
            SSNullModel = nan(1, fmriData.numVoxels);
            fit = nan(size(fmriData.data));
            degFreedom = fmriData.numTRs - rank(desMat);

            for h = unique(obj.hrfsIdx)
                currVox = obj.hrfsIdx == h;

                boldDesMat = fm_glm.conv(desMat, obj.hrfsMatrix(:,h));
                if obj.downsampleDesMat
                    boldDesMat = obj.downsampleDesMatToFmriResolution(boldDesMat);
                end

                X = [boldDesMat, nuissancesData.data];
                [betas(:,currVox), SSRegression(:,currVox)] = obj.computeOLS(...
                    X,...
                    fmriData.data(:,currVox));
                [~,                SSNullModel(:,currVox)]  = obj.computeOLS(...
                    nuissancesData.data,...
                    fmriData.data(:,currVox));

                fit(:,currVox) = [boldDesMat, nuissancesData.data] * betas(:,currVox);


                if obj.bTtestContrasts
                    tstats(:,currVox) = obj.ttestCalculator(...
                        X, betas(:,currVox), SSRegression(:,currVox), degFreedom);
                    % TODO: you're inverting X twice, one during regression
                    % and one during analysis
                end
            end
            
            RSqs = 1 - (SSRegression ./ SSNullModel);
        end

        function doHrfDerivative(obj, desMat, nuissancesData, fmriData)
            boldDesMat = {};
            for h = width(obj.hrfsMatrix):-1:1
                boldDesMat{h} = fm_glm.conv(desMat, obj.hrfsMatrix(:,h));
                if obj.downsampleDesMat
                    boldDesMat{h} = obj.downsampleDesMatToFmriResolution(boldDesMat{h});
                end
            end
            boldDesMat = cat(2, boldDesMat{:});

            [betas, SSRegression] = obj.computeOLS(...
                    [boldDesMat, nuissancesData.data],...
                    fmriData.data);
            [~,      SSNullModel] = obj.computeOLS(...
                    nuissancesData.data,...
                    fmriData.data);

            RSqs = 1 - (SSRegression ./ SSNullModel);
            fit = [boldDesMat, nuissancesData.data] * betas;


            numRegressors = width(desMat);
            biasCorrectedBetas = nan(numRegressors, fmriData.numVoxels);
            for b = 1:numRegressors
                biasCorrectedBetas(b,:) = ...
                    sign(betas(b,:)) .* ...
                    sqrt(sum(betas(b:numRegressors:width(boldDesMat), :).^2, 1));
                    % Lindquist et al. (2008), "Modeling the hemodynamic..."
            end

            if obj.bTtestContrasts
                tstats = obj.ttestCalculatorForDerivative(...
                    biasCorrectedBetas, SSRegression, fmriData.numTRs);
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

    methods (Access = private)
        function runPreparations(obj)
            assert(~isempty(obj.fmriData), "fmriData property is empty");
            assert(~isempty(obj.eventList), "eventList (design matrix) propety is empty");

            obj.dataTR = obj.fmriData.TR;
            obj.numRuns = length(obj.fmriData);
            obj.numRegressors = max(obj.eventList.cat().ID);


            timeInfo = cat(1, obj.eventList.Duration, obj.eventList.Onset, obj.dataTR);
            timeInfo = unique(timeInfo);

            multFactor = 4;
            decimallessTimeInfo = timeInfo * 10^(multFactor);
            assert(isRound(decimallessTimeInfo),...
                "The temporal resolution of the design matrix is too fine.");
            
            obj.desMatTR = gcdOfMany(decimallessTimeInfo) / 10^(multFactor);
            obj.downsampleDesMat = obj.desMatTR ~= obj.dataTR;


            TR = obj.dataTR;
            switch lower(obj.hrf)
                case 'nsd+canonical'
                    obj.hrfsMatrix = [fm_hrf.getAllNsdHrfs(TR), fm_hrf.getCanonicalHrf(TR)];
                case 'nsd'
                    obj.hrfsMatrix = fm_hrf.getAllNsdHrfs(TR);
                case 'canonical'
                    obj.hrfsMatrix = fm_hrf.getCanonicalHrf(TR);
                    obj.hrfsIdx = ones(1, obj.fmriData(1).numVoxels);
                        % TODO: do any other HRF methods need this?
                case 'derivative-1'
                    obj.hrfsMatrix = fm_hrf.getDerivativeHrfs(TR, 1);
                case 'derivative-2'
                    obj.hrfsMatrix = fm_hrf.getDerivativeHrfs(TR, 2);
                case 'custom'
                otherwise
                    error("Provided hrf parameter %s is invalid", obj.hrf);
            end


            if ~isempty(obj.ttestContrasts)
                assert(width(obj.ttestContrasts) == obj.numRegressors,...
                       "Number of columns in contrasts matrix should be equal " + ...
                       "to the number of event types. Do not include nuissances.");
            end
            


            if isempty(obj.hrfsIdx)
                obj.hrfsIdx = 1:width(obj.hrfsMatrix);
            end



            if strcmpi(obj.nuissances, 'baselines')
                obj.nuissancesData = obj.getBaselineNuissances([obj.fmriData.numTRs]);
            end


            %{
            if isempty(obj.lssEventsToDo)
                for i = length(obj.eventList):-1:1
                    obj.lssEventsToDo{i} = 1:height(obj.eventList);
                end
                % return; % return when this is in a function
            end

            if length(obj.eventList) > 1
                assert(iscell(obj.lssEventsToDo), ...
                       "For multiple runs, eventList must be a cell array");
            elseif ~iscell(obj.lssEventsToDo)
                obj.lssEventsToDo = {obj.lssEventsToDo};
            end

            for i = length(obj.lssEventsToDo):-1:1
                assert(width(obj.lssEventsToDo) == 1, "Input should be a column vector");
                assert(max(obj.lssEventsToDo) <= height(obj.eventList{i}) & ...
                       min(obj.lssEventsToDo) >= 1 & ...
                       isRound(obj.lssEventsToDo), ...
                       "Input should contain indices of events in eventList " + ...
                       "that are of interest");
            end
            %}





            % TODO: check that all parameter inputs are compatible
                % ideally, setup properties such that incompatible
                % parameters are not possible to specify
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