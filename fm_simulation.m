
classdef fm_simulation < matlab.mixin.Copyable
    properties
        taskTable;
        simProperties;

        groundTruthPatterns;
        groundTruthHRFs;
        eventList;
        trialSequence;

        eventBasePattern;
        eventNoisePattern;
        eventGroundTruth;
        neuralTimeSeries;
    end
 
    methods
        %%%
        function obj = fm_simulation(taskTable, simulationProperties)
            if (nargin > 0)
                obj.taskTable = taskTable;
                obj.simProperties = simulationProperties;
            end
        end

        %%%
        function generate(obj)
            % TODO: verify taskTable and simProperties are set

            obj.makeGroundTruthPatterns();
            obj.makeEventList();
            obj.scaleEventListByActivity();
            obj.makeNeuralTimeSeries();
            keyboard;

            obj.makeGroundTruthHRFs();
           
            % createNTS
            % convolve with HRF
            % add noise

            
        end

        %%%
        function makeGroundTruthPatterns(obj)
            obj.groundTruthPatterns = rand(obj.taskTable.numNeuralPatterns, obj.simProperties.numVoxels);
        end

        %%%
        function makeEventList(obj)
            taskTableContent = obj.taskTable.contentNumerical;
            obj.eventList = cell(1, obj.simProperties.numRuns);
            obj.trialSequence = cell(1, obj.simProperties.numRuns);
            for i = 1:obj.simProperties.numRuns
                [obj.eventList(i), obj.trialSequence(i)] = makefmriseq(...
                    taskTableContent.Durations(:)',...
                    taskTableContent.EventIDs(:)',...
                    taskTableContent.Probability(:)',...
                    obj.simProperties.runDuration,...
                    1,...
                    obj.simProperties.itiModel,...
                    obj.simProperties.itiParams,...
                    obj.simProperties.TR,...
                    'addExtraTrials', 0);
            end
        end

        %%%
        function scaleEventListByActivity(obj)
            for i = 1:obj.simProperties.numRuns
                obj.eventList{i}.Activity = obj.taskTable.NeuralActivity(obj.eventList{i}.ID)';
            end
        end

        %%%
        function makeNeuralTimeSeries(obj)
            obj.eventBasePattern = cell(1, obj.simProperties.numRuns);
            obj.eventNoisePattern = cell(1, obj.simProperties.numRuns);
            obj.eventGroundTruth = cell(1, obj.simProperties.numRuns);
            obj.neuralTimeSeries = cell(1, obj.simProperties.numRuns);

            for i = 1:obj.simProperties.numRuns
                obj.eventBasePattern{i} = generateBasePattern(obj.eventList{i});
                obj.eventNoisePattern{i} = generateEventNoisePattern(obj.eventList{i});
                obj.eventGroundTruth{i} = obj.eventBasePattern{i} + obj.eventNoisePattern{i};
                obj.neuralTimeSeries{i} = computeGroundTruthTimeSeries(obj.eventList{i}, obj.eventGroundTruth{i});
            end

            function eventBasePatternRun = generateBasePattern(eventListRun)
                neuralPatternPerEvent = obj.taskTable.NeuralPatternIDs(eventListRun.ID);
                eventBasePatternRun = obj.groundTruthPatterns(neuralPatternPerEvent, :);
                eventBasePatternRun = eventBasePatternRun .* eventListRun.Activity;
            end

            function eventNoisePatternRun = generateEventNoisePattern(eventListRun)                
                noiseScale = eventListRun.Activity * obj.simProperties.eventToEventNoiseAmount;
                if obj.simProperties.eventToEventNoiseInOnlyClassifiedEvents
                    eventsToClassify = obj.taskTable.ClassificationGroups((eventListRun.ID)) ~= obj.taskTable.NON_CLASSIFIED_EVENT;
                    noiseScale = noiseScale .* eventsToClassify(:);
                    eventsToCorrelate = eventsToClassify(:);
                else
                    eventsToCorrelate = ones(size(noiseScale));
                end
                eventsToCorrelate = logical(eventsToCorrelate);

                eventNoisePatternRun = rand(height(eventListRun), obj.simProperties.numVoxels) - 0.5;
                eventNoisePatternRun = eventNoisePatternRun .* noiseScale;
                eventNoisePatternRun(eventsToCorrelate, :) = makecorrelated(...
                    eventNoisePatternRun(eventsToCorrelate, :),...
                    obj.simProperties.eventToEventNoiseCoherence);
            end

            function neuralTimeSeriesRun = computeGroundTruthTimeSeries(eventListRun, eventPatterns)
                eventListRun.ID = (1:height(eventListRun))';
                timePointsOfEachEvent = computeDesignMatrix(...
                    eventListRun,...
                    obj.simProperties.runDuration,...
                    obj.simProperties.TR);

                neuralTimeSeriesRun = zeros(obj.simProperties.numTRs, obj.simProperties.numVoxels);
                for ev = 1:height(eventListRun)
                    thisEventTimePoints = logical(timePointsOfEachEvent(:,ev));
                    neuralTimeSeriesRun(thisEventTimePoints, :) = ...
                        repmat(eventPatterns(ev,:), sum(thisEventTimePoints), 1);
                end
            end
        end

        %%%
        function makeGroundTruthHRFs(obj)
            cacheLocation = fullfile(functionalMage.getCacheDirectory(), 'hrfDatabase.mat');
            hrfLen = 40;
            obj.groundTruthHRFs = getHrfDb(...
                obj.simProperties.hrfsCorrelationWithDoubleGammaCanonical,...
                obj.simProperties.numVoxels,...
                obj.simProperties.TR,...
                obj.simProperties.hrfLibrary,...
                cacheLocation, ...
                hrfLen);
        end
    end

end