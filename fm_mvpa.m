classdef fm_mvpa < matlab.mixin.Copyable
%FM_MVPA Perform multivariate pattern analysis
% Take neural activity data as a vector of fm_data objects along with
% their labels and classification groups and output classification accuracy
% or "pattern similarity" score
% 
% Input properties
%   <betas> A vector of fm_data objects containing estimates of neural
%     activity. This is generally taken from an fm_glm object.
%   <labels> A cell array of column vectors indicating the classification
%     label of each estimate. Vector values should be integers.
%   <groups> A cell array of column vectors consisting of integers
%     indicating the classification group of each estimate. Only estimates
%     with the same 'group' values are classified together.
%   <labelsAndGroups> You can fill both labels and groups directly if you
%     supply this property with the output from a getMvpa*() method from
%     an fm_designMatrix object.
%   <doNotClassifyGroupIDs> If you want to omit classifying certain events,
%     add their 'group' values to this vector. By default, this property
%     consists of only a zero. Therefore, if you used fm_task and set
%     'ClassificationGroup' to 0, you do not need to modify this property.
%
% Methods
%   > obj = fm_mvpa(betas, labelsAndGroupsStruct) will return an fm_mvpa
%     object. 'labelsAndGroupsStruct' is a [1 x numRuns] structure array with
%     fields 'Labels' and 'Groups' (case-sensitive).
%   > obj = fm_mvpa(betas, labels, groups) is another way to create an
%     fm_mvpa object. 'labels' and 'groups' are as they are defined above
%     in the 'Input properties' section.
%   > obj = fm_mvpa() will return an fm_mvpa object without input
%     properties set.
%   > ca = obj.getClassificationAccuracy(), where ca is a [1 x numGroups]
%     vector of classifcation accuracy values. MVPA-Light toolbox is used
%     for a Linear Discriminant Analysis. Data from all but one run is used
%     for training and the last run is used for testing.
%     https://github.com/treder/MVPA-Light
%   > pattSim = obj.getPatternSimilarity(gtInBetasForm), where
%     pattSim is a [1 x numGroups] vector of Pearson's correlation values
%     between ground-truth neural patterns and the betas. gtInBetasForm is
%     an fm_data vector the same size as the 'betas' properties containing
%     the ground-truth neural activity patterns. You can get it by running
%     the getGroundTruthInBetasForm() method on the fm_glm object.
%
% Protected methods
%   > perf = runClassification(obj, trainData, trainLabels, testData,
%   testLabels) runs one iteration of the classification analysis. You can
%   make a new class that inherits from fm_mvpa and override this method to
%   implement a custom classification algorithm.
%
% Example
%   % add MVPA-light toolbox
%   addpath('~/Downloads/MVPA-Light-master/startup')
%   startup_MVPA_Light()
%   % classify
%   mvpa = fm_mvpa([glm.results.betas], dm.getMvpaLSA());
%       % glm is an fm_glm object and dm is an fm_designMatrix object
%   ca = mvpa.getClassificationAccuracy();
%   disp(ca) shows a classification accuracy value for each group
%
% Part of package funkyMage. November 2023.
% https://github.com/aj-abdujabborov/funkyMage

    properties
        betas (1,:) fm_data;
    end

    properties (Dependent = true)
        labelsAndGroups;
    end

    properties
        labels;
        groups;

        doNotClassifyGroupIDs = [0];
    end

    properties (Access = protected)
        numRuns;
        numGroups;
        groupIDs;

        runSelector;
        groupSelector;
        labelSelector;
    end

    properties (Access = protected, Constant = true)
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
                obj.labelsAndGroups = varargin{1};
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
    
                    perf = obj.runClassification(trainData, trainLabels, testData, testLabels);
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

    methods (Access = protected)
        function perf = runClassification(obj, trainData, trainLabels, testData, testLabels)
            % varargin{1} -> trainData
            % varargin{2} -> trainLabels
            % varargin{3} -> testData
            % varargin{4} -> testLabels
            perf = mv_classify(obj.classifCfg, trainData, trainLabels, testData, testLabels);
        end
    end

    methods (Access = protected)
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

            obj.runSelector   = getRunSelector();
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

        function set.labelsAndGroups(obj, labelsAndGroups)
            obj.labels = {labelsAndGroups.Labels};
            obj.groups = {labelsAndGroups.Groups};
        end
    end


    methods (Access = protected, Static = true)
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