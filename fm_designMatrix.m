classdef fm_designMatrix < matlab.mixin.Copyable
%FM_DESIGNMATRIX Do all design-matrix related maneuvers
% Combine an 'fm_task' object and an associated 'fm_eventList' object
% (acquired from 'fm_sequence') and produce design matrices for simulation,
% General Linear Modeling, and classification analysis.
%
% Input properties
%   <task> A filled fm_task object
%   <unqEventList> A vector of 'fm_eventList' objects acquired from
%   'fm_sequence'
%   <runwiseAnalysisIDs> The AnalysisIDs of events (from 'fm_task') that
%   you'd like to analyze run-wise no matter what (even with a Least
%   Squares All or trial-wise analysis). Empty by default.
%   <doNotClassifyGroupIDs> ClassificationGroup values of events (from
%   'fm_task') that you do not want to classify. 0 by default. In other
%   words, set a 0 to all events you do not want classified in
%   'ClassificationGroups'.
% 
% Methods
%   > dm = fm_designMatrix(task, unqEventList) returns an fm_designMatrix
%   object.
%   > getNeuralPatternIDEventList() returns a vector of fm_eventList
%   objects that are supplied to fm_simulation.
%   > getGlmLSA() returns a vector of fm_eventList objects to be used by
%   fm_glm for a Least Squares All analysis.
%   > getMvpaLSA() returns a structure with fields 'label' and 'group'
%   indicating the classification label and ClassificationGroup of each
%   beta estimate produced by fm_glm. It's to be used with fm_mvpa. In
%   other words, once you use getGlmLSA() on fm_glm, you'd use getMvpaLSA()
%   on fm_mvpa.
%   > getGlmLSS1() for a Least Squares Separate analysis, which takes a
%   "pull one out" approach. For each event, a separate OLS General Linear
%   Model is run where that event gets its own regressor and the other
%   events are all combined into another regressor. 
%   > getGlmLSS2() Also a Least Squares Separate, but the "other" events
%   are combined according to their AnalysisID rather than into one
%   regressor.
%   > getMvpaLSS() for Least Squares Separate classifications.
%   > getGlmLSU() for a Least Squares Unitary analysis, where rather than
%   get an estimate for each event, an estimate is produced for each event
%   type (here indicated by AnalysisID).
%   > getMvpaLSU() for a Least Squares Unitary classification analysis.
%
% Examples
%   dm = fm_designMatrix(task, seq.eventList) % seq is an 'fm_sequence' object
%   dm.runwiseAnalysisIDs = [1 2]; 
%        % collapse multiple events into two regressors
%   dm.getGlmLSA() % for LSA GLM
%   dm.getMvpaLSA() % for MVPA analysis using fm_glm output
%
% Part of package funkyMage. November 2023.
% https://github.com/aj-abdujabborov/funkyMage

    properties (Dependent = true)
        task {mustBeA(task, 'fm_task')};
    end

    properties
        unqIDEventList (1,:) {mustBeA(unqIDEventList, 'fm_eventList')} = fm_eventList.empty;
        runwiseAnalysisIDs (1,:) {mustBeInteger} = [];
        doNotClassifyGroupIDs = [0];
    end

    properties (Access = protected)
        privateTaskTable;

        % 1D identifier vectors
        NeuralIntensity;
        NeuralPatternIDs; numNeuralPatterns;
        AnalysisIDs;
        ClassificationGroups;
        EventIDs;

        idLessEventList = fm_eventList.empty;
        trialSequence;

        eventNoIntoAnalysisID;
        eventNoIntoRegID;
        regIDIntoAnalysisID;
        regIDIntoClassifLabel;
        regIDIntoClassifGroup;

        lastMethod string {matlab.system.mustBeMember(lastMethod, {'lsa', 'lss2', 'lss1', 'lsu', ''})} = "";
        numRuns;

        prepSuccess = false;
    end

    methods
        %%%
        function obj = fm_designMatrix(task, eventList)
            if nargin == 0
                return;
            end
            if nargin >= 1
                obj.task = task;
            end
            if nargin >= 2
                obj.unqIDEventList = eventList;
            end
        end

        %%%
        function neuralPatternIDEventList = getNeuralPatternIDEventList(obj)
            obj.prepare();

            for i = obj.numRuns:-1:1
                neuralPatternIDEventList(i)          = obj.unqIDEventList(i);
                neuralPatternIDEventList(i).Activity = obj.NeuralIntensity(obj.unqIDEventList(i).ID);
                neuralPatternIDEventList(i).ID       = obj.NeuralPatternIDs(obj.unqIDEventList(i).ID);
            end
        end

        %%% LSA
        function eventList = getGlmLSA(obj)
            obj.prepare();

            for i = obj.numRuns:-1:1
                obj.eventNoIntoRegID{i} = obj.transformToRunwiseAndEventwise(...
                    obj.eventNoIntoAnalysisID{i},...
                    obj.runwiseAnalysisIDs);
            end

            eventList = obj.replaceIDsOfEventList(obj.idLessEventList, obj.eventNoIntoRegID);

            obj.lastMethod = "lsa";
        end

        function labelsAndGroups = getMvpaLSA(obj)
            if ~(obj.lastMethod == "lsa")
                obj.getGlmLSA();
            end
            
            obj.computeLabelsAndGroupsAdjustedToRegIDs();
            
            labelsAndGroups = struct('Labels', obj.regIDIntoClassifLabel, ...
                                     'Groups', obj.regIDIntoClassifGroup);
        end

        %%% LSS
        function eventList = getGlmLSS2(obj)
            obj.prepare();

            obj.computeRegIDWithAllEventsRunwise();
            eventList = obj.replaceIDsOfEventList(obj.idLessEventList, obj.eventNoIntoRegID);

            obj.lastMethod = "lss2";
        end

        function eventList = getGlmLSS1(obj)
            obj.prepare();
            
            for i = obj.numRuns:-1:1
                eventList(i) = obj.idLessEventList(i);
                eventList(i).ID = 1;
            end

            obj.lastMethod = "lss1";
        end

        function labelsAndGroups = getMvpaLSS(obj)
            obj.prepare();
            
            for i = obj.numRuns:-1:1
                obj.regIDIntoClassifGroup{i} = obj.replaceValues(...
                    obj.unqIDEventList(i).ID, ...
                    obj.EventIDs, ...
                    obj.ClassificationGroups);
                nonClaIdx = ismember(obj.regIDIntoClassifGroup{i},...
                                     obj.doNotClassifyGroupIDs);

                obj.regIDIntoClassifLabel{i} = obj.replaceValues(...
                    obj.unqIDEventList(i).ID, ...
                    obj.EventIDs, ...
                    obj.NeuralPatternIDs);
                obj.regIDIntoClassifLabel{i}(nonClaIdx) = nan;
            end

            labelsAndGroups = struct('Labels', obj.regIDIntoClassifLabel, ...
                                     'Groups', obj.regIDIntoClassifGroup);

        end

        %%% LSU
        function eventList = getGlmLSU(obj)
            obj.prepare();

            obj.computeRegIDWithAllEventsRunwise();
            eventList = obj.replaceIDsOfEventList(obj.idLessEventList, obj.eventNoIntoRegID);

            obj.lastMethod = "lsu";
        end


        function labelsAndGroups = getMvpaLSU(obj)
            if ~(obj.lastMethod == "lsu")
                obj.getGlmLSU();
            end
            
            obj.computeLabelsAndGroupsAdjustedToRegIDs();
            
            labelsAndGroups = struct('Labels', obj.regIDIntoClassifLabel, ...
                                     'Groups', obj.regIDIntoClassifGroup);
        end
    end

    methods (Access = private)
        function prepare(obj)
            if obj.prepSuccess
                return;
            end
            
            assert(~isempty(obj.task), "task property is not set");
            assert(~isempty(obj.unqIDEventList), "unqIDEventList should be set");
            obj.numRuns = length(obj.unqIDEventList);

            for i = obj.numRuns:-1:1
                obj.eventNoIntoAnalysisID{i} = obj.AnalysisIDs(obj.unqIDEventList(i).ID);

                obj.idLessEventList(i) = obj.unqIDEventList(i);
                obj.idLessEventList(i).ID = 999;
            end

            obj.prepSuccess = true;
        end

        function computeRegIDWithAllEventsRunwise(obj)
            for i = obj.numRuns:-1:1
                obj.eventNoIntoRegID{i} = obj.transformToRunwiseAndEventwise(...
                    obj.eventNoIntoAnalysisID{i},...
                    unique(obj.AnalysisIDs(:)'));
            end
        end

        function computeLabelsAndGroupsAdjustedToRegIDs(obj)
            for i = obj.numRuns:-1:1 
                obj.regIDIntoAnalysisID{i} = obj.transform2AIntoB(...
                    obj.eventNoIntoRegID{i},...
                    obj.eventNoIntoAnalysisID{i});

                obj.regIDIntoClassifGroup{i} = obj.replaceValues(...
                    obj.regIDIntoAnalysisID{i}, ...
                    obj.AnalysisIDs, ...
                    obj.ClassificationGroups);
                nonClaIdx = ismember(obj.regIDIntoClassifGroup{i},...
                                     obj.doNotClassifyGroupIDs);

                obj.regIDIntoClassifLabel{i} = obj.replaceValues(...
                    obj.regIDIntoAnalysisID{i}, ...
                    obj.AnalysisIDs, ...
                    obj.NeuralPatternIDs);
                obj.regIDIntoClassifLabel{i}(nonClaIdx) = nan;
            end
        end
    end

    methods (Access = protected, Static = true)
        function transformed = transformToRunwiseAndEventwise(IDs, IDsToMakeRunwise)
            transformed = nan(size(IDs));

            assignID = 1;
            runwiseEventIdx = IDs == IDsToMakeRunwise;
            for j = 1:length(IDsToMakeRunwise)
                transformed(runwiseEventIdx(:,j)) = assignID;
                assignID = assignID + 1;
            end

            nonRunwiseEventIdx = ~any(runwiseEventIdx, 2);
            transformed(nonRunwiseEventIdx) = assignID : (assignID + sum(nonRunwiseEventIdx) - 1);
        end

        function AIntoB = transform2AIntoB(A, B)
            unqRegIDs = unique(A(:)');
            for unqID = unqRegIDs(end:-1:1)
                tmp = B(A == unqID);
                AIntoB(unqID) = unique(tmp); %#ok<AGROW>
            end
            AIntoB = AIntoB(:);
        end

        function eventList = replaceIDsOfEventList(eventList, IDs)
            for i = length(eventList):-1:1
                eventList(i).ID = IDs{i};
            end
        end

        function replaced = replaceValues(data, A, B)
            assert(isequal(size(A), size(B)), "The two vectors must be equal sizes");
            
            replaced = nan(size(data));
            unqA = unique(A);
            for currUnqA = unqA(:).'
                correspondingB = unique(B(A == currUnqA));
                assert(numel(correspondingB) <= 1, "All elements of the same value in vecA must map onto elments of the same value in vecB");
                replaced(data == currUnqA) = correspondingB;
            end
        end
    end

    methods % Set and get methods
        function set.runwiseAnalysisIDs(obj, runwiseAnalysisIDs)
            obj.runwiseAnalysisIDs = obj.runwiseAnalysisIDs(:)';
            obj.runwiseAnalysisIDs = unique(runwiseAnalysisIDs);
        end

        function set.task(obj, task)
            obj.prepSuccess = false;

            obj.NeuralIntensity = collapseCellArray(task.contentNumerical.NeuralIntensity);
            obj.NeuralPatternIDs = collapseCellArray(task.contentNumerical.NeuralPatternIDs);
            obj.AnalysisIDs = collapseCellArray(task.contentNumerical.AnalysisIDs);
            obj.ClassificationGroups = collapseCellArray(task.contentNumerical.ClassificationGroups);
            obj.EventIDs = collapseCellArray(task.contentNumerical.EventIDs);

            obj.numNeuralPatterns = length(unique(obj.NeuralPatternIDs));
            % obj.timings = task.content(:, ["Durations", "Onsets"]);

            obj.privateTaskTable = task;

            function collapsed = collapseCellArray(cellArray)
                collapsed = [cellArray{:}];
                collapsed = collapsed(:);
            end
        end

        function task = get.task(obj)
            task = obj.privateTaskTable;
        end
    end
end