% TODO: add storage of timing and add capabiilty to update the timings

% TODO:
% * check mapability from one vector to another e.g., the
% same AnalysisID values cannot map to
% different NeuralPatternIDs (unless they are chosen to not
% be classified). same from AnalysisID to ClassificationGroups

classdef fm_designMatrix < matlab.mixin.Copyable
    properties (Dependent = true)
        taskTable {mustBeA(taskTable, 'fm_taskTable')};
        simProperties {mustBeA(simProperties, 'fm_simulationProperties')};
    end

    properties
        runwiseAnalysisIDs (1,:) {mustBeInteger} = [];
        timings;
    end

    properties (Hidden = true)
        doNotClassifyGroupIDs = [0];
    end

    properties (Access = private)
        privateSimProperties;
        privateTaskTable;

        % 1D identifier vectors
        NeuralIntensity;
        NeuralPatternIDs; numNeuralPatterns;
        AnalysisIDs;
        ClassificationGroups;
        EventIDs;

        unqIDEventList = fm_eventList.empty;
        idLessEventList = fm_eventList.empty;
        trialSequence;

        eventNoIntoAnalysisID;
        eventNoIntoRegID;
        regIDIntoAnalysisID;
        regIDIntoClassifLabel;
        regIDIntoClassifGroup;

        lastMethod string {matlab.system.mustBeMember(lastMethod, {'lsa', 'lss2', 'lss1', 'lsu', ''})} = "";
        numRuns;

        bEventListGenerated = false;
    end

    methods
        %%%
        function obj = fm_designMatrix(taskTable, simProperties)
            if nargin == 0
                return;
            end
            if nargin >= 1
                obj.taskTable = taskTable;
            end
            if nargin >= 2
                obj.simProperties = simProperties;
            end
        end

        %%%
        function regenerateEventList(obj)
            obj.bEventListGenerated = false;
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
            if obj.bEventListGenerated
                return;
            end
            
            assert(~isempty(obj.simProperties), "simProperties property is not set");
            assert(~isempty(obj.taskTable), "taskTable property is not set");

            obj.generateEventList();

            for i = obj.numRuns:-1:1
                obj.eventNoIntoAnalysisID{i} = obj.AnalysisIDs(obj.unqIDEventList(i).ID);

                obj.idLessEventList(i) = obj.unqIDEventList(i);
                obj.idLessEventList(i).ID = 999;
            end

            obj.bEventListGenerated = true;
        end

        function generateEventList(obj)
            for i = obj.numRuns:-1:1
                [eventList, trialSeq] = makefmriseq(...
                    obj.taskTable.contentNumerical.Durations(:)',...
                    obj.taskTable.contentNumerical.EventIDs(:)',...
                    obj.taskTable.contentNumerical.Probability(:)',...
                    obj.simProperties.runDuration,...
                    1,...
                    obj.simProperties.itiModel,...
                    obj.simProperties.itiParams,...
                    obj.simProperties.TR,...
                    'addExtraTrials', 0);
                obj.unqIDEventList(i) = fm_eventList(eventList{1}, obj.simProperties.runDuration);
                obj.trialSequence{i} = trialSeq{1};
            end
            % TODO: onset times are not used at all
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

    methods (Access = private, Static = true)
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

        %%% destined to trash:
        function bEmpty = isFieldClear(structIn, fieldname)
            bEmpty = isempty(structIn) || ~isfield(structIn, fieldname) || isempty(structIn(1).(fieldname));
        end
    end

    methods % Set and get methods
        function set.runwiseAnalysisIDs(obj, runwiseAnalysisIDs)
            obj.runwiseAnalysisIDs = obj.runwiseAnalysisIDs(:)';
            obj.runwiseAnalysisIDs = unique(runwiseAnalysisIDs);
        end

        function set.taskTable(obj, taskTable)
            if obj.bEventListGenerated
                error("Cannot set new taskTable unless you run the regenerateEventList() function");
            end

            obj.NeuralIntensity = collapseCellArray(taskTable.contentNumerical.NeuralIntensity);
            obj.NeuralPatternIDs = collapseCellArray(taskTable.contentNumerical.NeuralPatternIDs);
            obj.AnalysisIDs = collapseCellArray(taskTable.contentNumerical.AnalysisIDs);
            obj.ClassificationGroups = collapseCellArray(taskTable.contentNumerical.ClassificationGroups);
            obj.EventIDs = collapseCellArray(taskTable.contentNumerical.EventIDs);

            obj.numNeuralPatterns = length(unique(obj.NeuralPatternIDs));
            obj.timings = taskTable.content(:, ["Durations", "Onsets"]);

            obj.privateTaskTable = taskTable;

            function collapsed = collapseCellArray(cellArray)
                collapsed = [cellArray{:}];
                collapsed = collapsed(:);
            end
        end

        function taskTable = get.taskTable(obj)
            taskTable = obj.privateTaskTable;
        end

        function set.simProperties(obj, simProperties)
            if obj.bEventListGenerated
                error("Cannot set new simProperties unless you run the regenerateEventList() function");
            end

            obj.numRuns = simProperties.numRuns;
            obj.privateSimProperties = simProperties;
        end

        function simProperties = get.simProperties(obj)
            simProperties = obj.privateSimProperties;
        end
    end
end