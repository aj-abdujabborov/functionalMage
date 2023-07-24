% TODO: add storage of timing and add capabiilty to update the timings

classdef fm_designMatrix < matlab.mixin.Copyable
    properties
        runwiseAnalysisIDs;

        neuralPatternIDEventList;

        idVectors;
        timings;

        doNotClassifyGroupIDs = [0]; 
    end

    properties (Dependent = true, SetAccess = private)
        glmLSA;
        mvpaLSA;
    end

    properties (Access = private)
        unqIDEventList;
        trialSequence;

        NeuralIntensity;
        NeuralPatternIDs; numNeuralPatterns;
        AnalysisIDs;
        ClassificationGroups;
        EventIDs;

        numRuns;

        privateLsa = [];
    end

    methods
        %%%
        function obj = fm_designMatrix(taskTable, simProperties)
            % TODO: check that inputs are of correct classes.
            extractMappingVectors();
            obj.NeuralIntensity = collapseCellArray(taskTable.contentNumerical.NeuralIntensity);
            obj.numNeuralPatterns = numUniqueElements(obj.NeuralPatternIDs);
            obj.timings = taskTable.content(:, ["Durations", "Onsets"]);
            obj.numRuns = simProperties.numRuns;

            for i = obj.numRuns:-1:1
                [obj.unqIDEventList{i}, obj.trialSequence{i}] ...
                    = generateEventList(taskTable, simProperties);

                obj.neuralPatternIDEventList{i} = obj.unqIDEventList{i};
                obj.neuralPatternIDEventList{i}.Activity = obj.NeuralIntensity(obj.neuralPatternIDEventList{i}.ID);
                obj.neuralPatternIDEventList{i}.ID = obj.NeuralPatternIDs(obj.neuralPatternIDEventList{i}.ID);
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
                eventList = eventList{1};
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

    methods % Get methods
        function glmLSA = get.glmLSA(obj)
            if isempty(obj.privateLsa)
                for i = obj.numRuns:-1:1
                    obj.privateLsa.unqID2AnalysisID{i} = obj.AnalysisIDs(obj.unqIDEventList{i}.ID);
                    
                    
                    obj.privateLsa.unqID2RegressionID{i} = nan(size(obj.unqIDEventList{i}.ID));
                    assignID = 1;
                    runwiseEventIdx = obj.privateLsa.unqID2AnalysisID{i} == obj.runwiseAnalysisIDs;
                    for j = 1:length(obj.runwiseAnalysisIDs)
                        obj.privateLsa.unqID2RegressionID{i}(runwiseEventIdx(:,j)) = assignID;
                        assignID = assignID + 1;
                    end

                    nonRunwiseEventIdx = ~any(runwiseEventIdx, 2);
                    obj.privateLsa.unqID2RegressionID{i}(nonRunwiseEventIdx) ...
                        = assignID : (assignID + sum(nonRunwiseEventIdx) - 1);


                    obj.privateLsa.regressionDesignMatrix{i} = ...
                        obj.neuralPatternIDEventList{i};
                    obj.privateLsa.regressionDesignMatrix{i}.ID = ...
                        obj.privateLsa.unqID2RegressionID{i};
                end
            end

            glmLSA = obj.privateLsa.regressionDesignMatrix;
        end

        function mvpaLSA = get.mvpaLSA(obj)
            if obj.isFieldClear(obj.privateLsa, 'regressionDesignMatrix')
                obj.glmLSA();
            end
            PL = obj.privateLsa;

            if obj.isFieldClear(PL, 'regressionID2AnalysisID')
                for i = obj.numRuns:-1:1 
                    PL.regressionID2AnalysisID{i} = getRegressionID2AnalysisID(...
                        PL.unqID2RegressionID{i},...
                        PL.unqID2AnalysisID{i});

                    PL.regressionID2ClassifLabel{i} = getRegressionID2ClassifLabel(...
                        PL.regressionID2AnalysisID{i});

                    PL.regressionID2ClassifGroup{i} = getRegressionID2ClassifGroup(...
                        PL.regressionID2AnalysisID{i});
                end
            end

            obj.privateLsa = PL;

            mvpaLSA = [];
            mvpaLSA.Labels = PL.regressionID2ClassifLabel;
            mvpaLSA.Groups = PL.regressionID2ClassifGroup;

            function reg2analysis = getRegressionID2AnalysisID(regressionID, analysisID)
                unqRegIDs = unique(regressionID(:)');
                for unqID = unqRegIDs(end:-1:1)
                    tmp = analysisID(regressionID == unqID);
                    reg2analysis(unqID) = unique(tmp); %#ok<AGROW>
                end
                reg2analysis = reg2analysis(:);
            end

            function classifLabels = getRegressionID2ClassifLabel(analysisID)
                nonClassif = ismember(obj.ClassificationGroups, obj.doNotClassifyGroupIDs);
                old = obj.AnalysisIDs;      old(nonClassif) = 0;
                new = obj.NeuralPatternIDs; new(nonClassif) = 0;
                classifLabels = obj.replaceValues(analysisID, old, new);
            end

            function classifGroups = getRegressionID2ClassifGroup(analysisID)
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
        function bEmpty = isFieldClear(structVar, fieldname)
            bEmpty = isempty(structVar) || ~isfield(structVar, fieldname) || isempty(structVar.(fieldname));
        end

        function data = replaceValues(data, A, B)
            assert(isequal(size(A), size(B)), "The two vectors must be equal sizes");
            
            unqVecA = unique(A);
            for vecAValue = unqVecA(:).'
                valInVecB = unique(B(A == vecAValue));
                assert(numel(valInVecB) <= 1, "All elements of the same value in vecA must map onto elments of the same value in vecB");
                data(data == vecAValue) = valInVecB;
            end
        end

        %%% destined to trash
        function replaced = replaceValues2(data, old, new)
            [found, idx] = ismember(data, old);
            replaced = data;
            replaced(found) = new(idx(found));
        end
    end

end