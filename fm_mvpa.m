classdef fm_mvpa < matlab.mixin.Copyable
    properties
        betas (1,:) fm_data;
        labels;
        groups;

        doNotClassifyGroupIDs = [0];
    end

    properties (Access = private)
        numRuns;
        numGroups;
        groupIDs;

        runSelector;
        groupSelector;
        labelSelector;
    end

    properties (Access = private, Constant = true)
        classifCfg = struct('classifier', 'lda', ...
                            'metric', 'accuracy', ...
                            'feedback', 0, ...
                            'cv', 'none');
    end

    methods
        function obj = fm_mvpa(betas, varargin)
            if nargin == 0
                return;
            elseif nargin >= 1
                obj.betas = betas;
            end

            if isempty(varargin)
                return;
            end

            if length(varargin) == 1
                obj.labels = {varargin{1}.Labels};
                obj.groups = {varargin{1}.Groups};
            elseif length(varargin) == 2
                obj.labels = varargin{1};
                obj.groups = varargin{2};
            else
                error("Invalid fm_mvpa constructor call");
            end
        end

        function ca = getClassificationAccuracy(obj)
            obj.checkProperties();

            allBetas = cat(obj.betas).data;
            
            for g = obj.numGroups:-1:1
                caIterations = nan(obj.numRuns, 1);

                currLabelSelector = obj.labelSelector;
                currLabelSelector(obj.groupSelector ~= g) = nan;
                currLabelSelector = obj.mapValuesInto1ToN(currLabelSelector);
    
                for r = 1:obj.numRuns
                    trainSelector = obj.groupSelector == g & obj.runSelector ~= r;
                    testSelector  = obj.groupSelector == g & obj.runSelector == r;
    
                    trainData   = allBetas(trainSelector, :);
                    testData    = allBetas(testSelector, :);

                    trainLabels = currLabelSelector(trainSelector, :);
                    testLabels  = currLabelSelector(testSelector, :);
    
                    perf = mv_classify(obj.classifCfg, trainData, trainLabels, testData, testLabels);
                    caIterations(r) = perf;
                end

                ca(g) = mean(caIterations);
            end
        end

        function pattSim = getPatternSimilarity(obj, groundTruthInBetasForm)
            arguments
                obj;
                groundTruthInBetasForm (1,:) fm_data;
            end

            obj.checkProperties();
            validateInputs();

            corrList(:,1) = correlateCorrespondingRows(...
                                cat(obj.betas).data,...
                                cat(groundTruthInBetasForm).data);

            for g = obj.numGroups:-1:1
                currLabelSelector = obj.labelSelector;
                currLabelSelector(obj.groupSelector ~= g) = nan;
                currLabelSelector = obj.mapValuesInto1ToN(currLabelSelector);
                numLabels = max(currLabelSelector);

                pattSimRuns = nan(obj.numRuns, numLabels);
                for r = 1:obj.numRuns
                    for l = 1:numLabels
                        currEstsIdx = obj.groupSelector == g & obj.runSelector == r & currLabelSelector == l;
                        pattSimRuns(r, l) = mean(corrList(currEstsIdx));
                    end
                end
                pattSim(g) = mean(pattSimRuns, "all");
            end

            function validateInputs()
                assert(isequal([groundTruthInBetasForm.numRows], [obj.betas.numRows]), ...
                  "Number of parameters in groundTruth input and betas do not match");
                assert(isequal([groundTruthInBetasForm.numVoxels], [obj.betas.numVoxels]), ...
                      "Number of voxels in groundTruth input and betas do not match");
            end
            
            function correlation = correlateCorrespondingRows(A, B)
                correlation = nan(1, height(A));
                for i = 1:height(A)
                    correlation(i) = corr(A(i,:)', B(i,:)');
                end
            end
        end
    end

    methods (Access = private)
        function checkProperties(obj)
            obj.numRuns = length(obj.betas);
            for r = 1:obj.numRuns
                assert(obj.betas(r).numRows == height(obj.groups{r}) ...
                     & obj.betas(r).numRows == height(obj.labels{r}),...
                       "Number of rows in betas, groups and labels must match");
            end

            obj.groupIDs = unique(cat(1, obj.groups{:}));
            obj.groupIDs(ismember(obj.groupIDs, obj.doNotClassifyGroupIDs)) = [];
            obj.numGroups = length(obj.groupIDs);

            obj.runSelector = getRunSelector();
            obj.groupSelector = getGroupSelector();
            obj.labelSelector = cat(1, obj.labels{:});

            function runSelector = getRunSelector()
                runSelector = nan(sum([obj.betas.numRows]), 1);
                startInd = 1;
                for i = 1:obj.numRuns
                    endInd = startInd + obj.betas(i).numRows - 1;
                    runSelector(startInd : endInd) = i;
                    startInd = endInd + 1;
                end
            end

            function groupSelector = getGroupSelector()
                groupSelector = cat(1, obj.groups{:});
                groupSelector = obj.mapValuesInto1ToN(groupSelector, obj.groupIDs(:)');
            end
        end
    end

    methods % Set methods
        function set.labels(obj, labels)
            if isnumeric(labels)
                obj.labels = {labels(:)};
            elseif iscell(labels)
                obj.labels = cellfun(@(x) x(:), labels, 'UniformOutput', false);
            else
                error("Invalid input. Labels should be a [1 x numRuns] " + ...
                       "cell array of column vectors");
            end
        end

        function set.groups(obj, groups)
            if isnumeric(groups)
                obj.groups = {groups(:)};
            elseif iscell(groups)
                obj.groups = cellfun(@(x) x(:), groups, 'UniformOutput', false);
            else
                error("Invalid input. Groups should be a [1 x numRuns] " + ...
                       "cell array of column vectors");
            end
        end
    end


    methods (Access = private, Static = true)
        function columnVec = mapValuesInto1ToN(columnVec, valsToMap)
            nanIdx = isnan(columnVec);
            cleaned = columnVec(~nanIdx, :);
            
            if ~exist('valsToMap', 'var')
                valsToMap = unique(cleaned);
            end
            numValsToMap = length(valsToMap);
            cleaned = sum((cleaned == valsToMap(:)') .* (1:numValsToMap), 2);

            columnVec(~nanIdx, :) = cleaned;
        end
    end
end