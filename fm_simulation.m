
classdef fm_simulation < matlab.mixin.Copyable
    properties
        taskTable;
        simProperties;

        eventList;
        trialSequence;

        groundTruthPatterns;
        eventBasePattern;
        eventNoisePattern;
        eventGroundTruth;
        neuralTimeSeries;

        groundTruthHRFs;
        boldTimeSeries;
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

            keyboard;
            obj.groundTruthPatterns = obj.makeGroundTruthPatterns();
            obj.groundTruthHRFs     = obj.generateGroundTruthHRFs();

            for i = obj.simProperties.numRuns:-1:1
                [obj.eventList{i}, obj.trialSequence{i}] = obj.makeEventList();
                obj.eventList{i}                         = obj.scaleEventListByActivity(obj.eventList{i});

                obj.eventBasePatterns{i}  = obj.computeEventBasePatterns(obj.eventList{i}, obj.groundTruthPatterns);
                obj.eventNoisePatterns{i} = obj.generateEventNoisePatterns(obj.eventList{i});
                obj.eventGroundTruth{i}   = obj.eventBasePatterns{i} + obj.eventNoisePatterns{i};
                obj.neuralTimeSeries{i}   = obj.computeNeuralTimeSeries(obj.eventList{i}, obj.eventGroundTruth{i});

                obj.boldTimeSeries{i} = obj.convolveWithHRFs(obj.neuralTimeSeries{i}, obj.groundTruthHRFs);

                % obj.makeNoise();
                % obj.addNoise();
            end
        end

        %%%
        function [eventList, trialSequence] = makeEventList(obj)
            [eventList, trialSequence] = makefmriseq(...
                obj.taskTable.contentNumerical.Durations(:)',...
                obj.taskTable.contentNumerical.EventIDs(:)',...
                obj.taskTable.contentNumerical.Probability(:)',...
                obj.simProperties.runDuration,...
                1,...
                obj.simProperties.itiModel,...
                obj.simProperties.itiParams,...
                obj.simProperties.TR,...
                'addExtraTrials', 0);
        end

        %%%
        function eventList = scaleEventListByActivity(obj, eventList)
            eventList.Activity = obj.taskTable.NeuralActivity(eventList.ID)';
        end

        %%%
        function groundTruthPatterns = makeGroundTruthPatterns(obj)
            groundTruthPatterns = rand(obj.taskTable.numNeuralPatterns, obj.simProperties.numVoxels);
        end

        %%%
        function eventBasePatterns = computeEventBasePatterns(obj, eventList, groundTruthPatterns)
            neuralPatternPerEvent = obj.taskTable.NeuralPatternIDs(eventList.ID);
            eventBasePatterns = groundTruthPatterns(neuralPatternPerEvent, :);
            eventBasePatterns = eventBasePatterns .* eventListRun.Activity;
        end

        %%%
        function eventNoisePatterns = generateEventNoisePatterns(obj, eventList)
            noiseScale = eventList.Activity * obj.simProperties.eventToEventNoiseAmount;
            if obj.simProperties.eventToEventNoiseInOnlyClassifiedEvents
                eventsToClassify = obj.taskTable.ClassificationGroups(eventList.ID) ~= obj.taskTable.NON_CLASSIFIED_EVENT;
                noiseScale = noiseScale .* eventsToClassify(:);
                eventsToCorrelate = eventsToClassify(:);
            else
                eventsToCorrelate = ones(size(noiseScale));
            end
            eventsToCorrelate = logical(eventsToCorrelate);

            eventNoisePatterns = rand(height(eventList), obj.simProperties.numVoxels) - 0.5;
            eventNoisePatterns = eventNoisePatterns .* noiseScale;
            eventNoisePatterns(eventsToCorrelate, :) = makecorrelated(...
                eventNoisePatterns(eventsToCorrelate, :),...
                obj.simProperties.eventToEventNoiseCoherence);
        end

        %%%
        function neuralTimeSeries = computeNeuralTimeSeries(obj, eventList, eventPatterns)
            eventList.ID = (1:height(eventList))';
            timePointsOfEachEvent = computeDesignMatrix(...
                eventList,...
                obj.simProperties.runDuration,...
                obj.simProperties.TR);

            neuralTimeSeries = zeros(obj.simProperties.numTRs, obj.simProperties.numVoxels);
            for ev = 1:height(eventList)
                thisEventTimePoints = logical(timePointsOfEachEvent(:,ev));
                neuralTimeSeries(thisEventTimePoints, :) = ...
                    repmat(eventPatterns(ev,:), sum(thisEventTimePoints), 1);
            end
        end

        %%%
        function HRFs = generateGroundTruthHRFs(obj)
            cacheLocation = fullfile(functionalMage.getCacheDirectory(), 'hrfDatabase.mat');
            hrfLen = 40;
            HRFs = getHrfDb(...
                obj.simProperties.hrfsCorrelationWithDoubleGammaCanonical,...
                obj.simProperties.numVoxels,...
                obj.simProperties.TR,...
                obj.simProperties.hrfLibrary,...
                cacheLocation, ...
                hrfLen);
        end

        %%%
        function boldTimeSeries = convolveWithHRFs(~, timeSeries, HRFs)
            boldTimeSeries = convolveByColumn(timeSeries, HRFs);
        end
    end

end