% TODO: add storage of timing and add capabiilty to update the timings

classdef fm_designMatrix < matlab.mixin.Copyable
    properties
        runwiseAnalysisIDs;

        neuralPatternIDEventList = fm_eventList.empty;

        idVectors;
        timings;

        doNotClassifyGroupIDs = [0];
    end

    properties (Access = private)
        unqIDEventList = fm_eventList.empty;
        trialSequence;

        NeuralIntensity;
        NeuralPatternIDs; numNeuralPatterns;
        AnalysisIDs;
        ClassificationGroups;
        EventIDs;

        numRuns;

        lsa (1,:) struct = fm_designMatrix.getGlmInfoStruct();
    end

    methods
        %%%
        function obj = fm_designMatrix(taskTable, simProperties)
            if nargin == 0
                return;
            end

            % TODO: check that inputs are of correct classes.
            extractMappingVectors();
            obj.NeuralIntensity = collapseCellArray(taskTable.contentNumerical.NeuralIntensity);
            obj.numNeuralPatterns = numUniqueElements(obj.NeuralPatternIDs);
            obj.timings = taskTable.content(:, ["Durations", "Onsets"]);
            obj.numRuns = simProperties.numRuns;

            for i = obj.numRuns:-1:1
                [obj.unqIDEventList(i), obj.trialSequence{i}] ...
                    = generateEventList(taskTable, simProperties);

                obj.neuralPatternIDEventList(i) = obj.unqIDEventList(i);
                obj.neuralPatternIDEventList(i).Activity = obj.NeuralIntensity(obj.neuralPatternIDEventList(i).ID);
                obj.neuralPatternIDEventList(i).ID = obj.NeuralPatternIDs(obj.neuralPatternIDEventList(i).ID);
            end

            function extractMappingVectors()
                obj.NeuralPatternIDs = collapseCellArray(taskTable.contentNumerical.NeuralPatternIDs);
                obj.AnalysisIDs = collapseCellArray(taskTable.contentNumerical.AnalysisIDs);
                obj.ClassificationGroups = collapseCellArray(taskTable.contentNumerical.ClassificationGroups);
                obj.EventIDs = collapseCellArray(taskTable.contentNumerical.EventIDs);
            end
            
            function [eventList, trialSequence] = generateEventList(taskTable, simProperties)
                [eventList, trialSequence] = makefmriseq(...
                    taskTable.contentNumerical.Durations(:)',...
                    taskTable.contentNumerical.EventIDs(:)',...
                    taskTable.contentNumerical.Probability(:)',...
                    simProperties.runDuration,...
                    1,...
                    simProperties.itiModel,...
                    simProperties.itiParams,...
                    simProperties.TR,...
                    'addExtraTrials', 0);
                eventList = fm_eventList(eventList{1}, simProperties.runDuration);
                trialSequence = trialSequence{1};
            end

            function collapsed = collapseCellArray(cellArray)
                collapsed = [cellArray{:}];
                collapsed = collapsed(:);
            end

            function numElements = numUniqueElements(vector)
                numElements = length(unique(vector));
            end
        end        
    end

    methods
        function glmLSA = getGlmLSA(obj)
            if ~obj.isFieldClear(obj.lsa, 'regEventList')
                glmLSA = [obj.lsa.regEventList];
                return;
            end

            for i = obj.numRuns:-1:1
                obj.lsa(i).eventNoIntoAnalysisID = obj.AnalysisIDs(obj.unqIDEventList(i).ID);
                
                obj.lsa(i).eventNoIntoRegID = nan(size(obj.unqIDEventList(i).ID));
                assignID = 1;
                runwiseEventIdx = obj.lsa(i).eventNoIntoAnalysisID == obj.runwiseAnalysisIDs;
                for j = 1:length(obj.runwiseAnalysisIDs)
                    obj.lsa(i).eventNoIntoRegID(runwiseEventIdx(:,j)) = assignID;
                    assignID = assignID + 1;
                end

                nonRunwiseEventIdx = ~any(runwiseEventIdx, 2);
                obj.lsa(i).eventNoIntoRegID(nonRunwiseEventIdx) ...
                    = assignID : (assignID + sum(nonRunwiseEventIdx) - 1);

                obj.lsa(i).regEventList = ...
                    obj.unqIDEventList(i);
                obj.lsa(i).regEventList.ID = ...
                    obj.lsa(i).eventNoIntoRegID;
            end

            glmLSA = [obj.lsa.regEventList];
        end

        function mvpaLSA = getMvpaLSA(obj)
            if obj.isFieldClear(obj.lsa, 'regEventList')
                obj.getGlmLSA();
            end

            if obj.isFieldClear(obj.lsa, 'regIDIntoAnalysisID')
                for i = obj.numRuns:-1:1 
                    obj.lsa(i).regIDIntoAnalysisID = getRegIDIntoAnalysisID(...
                        obj.lsa(i).eventNoIntoRegID,...
                        obj.lsa(i).eventNoIntoAnalysisID);

                    obj.lsa(i).regIDIntoClassifGroup = getRegIDIntoClassifGroup(...
                        obj.lsa(i).regIDIntoAnalysisID);
                    nonClaIdx = ismember(obj.lsa(i).regIDIntoClassifGroup,...
                                         obj.doNotClassifyGroupIDs);

                    obj.lsa(i).regIDIntoClassifLabel = getRegIDIntoClassifLabel(...
                        obj.lsa(i).regIDIntoAnalysisID);
                    obj.lsa(i).regIDIntoClassifLabel(nonClaIdx) = nan;
                end
            end

            mvpaLSA = struct('Labels', {obj.lsa.regIDIntoClassifLabel}, ...
                             'Groups', {obj.lsa.regIDIntoClassifGroup});

            function reg2analysis = getRegIDIntoAnalysisID(regressionID, analysisID)
                unqRegIDs = unique(regressionID(:)');
                for unqID = unqRegIDs(end:-1:1)
                    tmp = analysisID(regressionID == unqID);
                    reg2analysis(unqID) = unique(tmp); %#ok<AGROW>
                end
                reg2analysis = reg2analysis(:);
            end

            function classifLabels = getRegIDIntoClassifLabel(analysisID)
                classifLabels = obj.replaceValues(analysisID, obj.AnalysisIDs, obj.NeuralPatternIDs);
            end

            function classifGroups = getRegIDIntoClassifGroup(analysisID)
                classifGroups = obj.replaceValues(analysisID, obj.AnalysisIDs, obj.ClassificationGroups);
            end
        end
    end

    methods % Set methods
        function set.runwiseAnalysisIDs(obj, runwiseAnalysisIDs)
            obj.runwiseAnalysisIDs = obj.runwiseAnalysisIDs(:)';
            obj.runwiseAnalysisIDs = unique(runwiseAnalysisIDs);
        end
    end

    methods (Access = private, Static = true)
        function bEmpty = isFieldClear(structIn, fieldname)
            bEmpty = isempty(structIn) || ~isfield(structIn, fieldname) || isempty(structIn(1).(fieldname));
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

        function outStruct = getGlmInfoStruct()
            outStruct = struct('eventNoIntoAnalysisID', [], ...
                               'eventNoIntoRegID',      [], ...
                               'regEventList',          fm_eventList.empty, ...
                               'regIDIntoAnalysisID',   [], ...
                               'regIDIntoClassifLabel', [], ...
                               'regIDIntoClassifGroup', []);
        end

        %%% destined to trash
        function replaced = replaceValues2(data, old, new)
            [found, idx] = ismember(data, old);
            replaced = data;
            replaced(found) = new(idx(found));
        end
    end

end