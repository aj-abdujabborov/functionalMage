
classdef fm_simulation < matlab.mixin.Copyable
    properties
        taskTable;
        simProperties;

        eventList;
        trialSequence;

        groundTruthPatterns;

        basePatternPerEvent;
        noisePatternPerEvent;
        groundTruthPatternPerEvent;
        groundTruthTimeSeries;

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
            obj.groundTruthPatterns = obj.generateGroundTruthPatterns();
            obj.groundTruthHRFs     = obj.generateGroundTruthHRFs();

            for i = obj.simProperties.numRuns:-1:1
                [obj.eventList{i}, obj.trialSequence{i}] = obj.generateEventList();
                obj.eventList{i}                         = obj.scaleEventListByActivity(obj.eventList{i});

                obj.basePatternPerEvent{i}        = obj.computeBasePatternPerEvent(obj.eventList{i}, obj.groundTruthPatterns);
                obj.noisePatternPerEvent{i}       = obj.generateNoisePatternPerEvent(obj.eventList{i});
                obj.groundTruthPatternPerEvent{i} = obj.basePatternPerEvent{i} + obj.noisePatternPerEvent{i};
                obj.groundTruthTimeSeries{i}      = obj.computeGroundTruthTimeSeries(obj.eventList{i}, obj.groundTruthPatternPerEvent{i});

                obj.boldTimeSeries{i} = obj.convolveWithHRFs(obj.groundTruthTimeSeries{i}, obj.groundTruthHRFs);

                % obj.makeNoise();
                % obj.addNoise();
            end
        end

        %%%
        function [eventList, trialSequence] = generateEventList(obj)
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
            eventList = eventList{1};
            trialSequence = trialSequence{1};
        end

        %%%
        function eventList = scaleEventListByActivity(obj, eventList)
            eventList.Activity = obj.taskTable.NeuralActivity(eventList.ID)';
        end

        %%%
        function groundTruthPatterns = generateGroundTruthPatterns(obj)
            groundTruthPatterns = rand(obj.taskTable.numNeuralPatterns, obj.simProperties.numVoxels);
        end

        %%%
        function basePatternPerEvent = computeBasePatternPerEvent(obj, eventList, groundTruthPatterns)
            neuralPatternPerEvent = obj.taskTable.NeuralPatternIDs(eventList.ID);
            basePatternPerEvent = groundTruthPatterns(neuralPatternPerEvent, :);
            basePatternPerEvent = basePatternPerEvent .* eventList.Activity;
        end

        %%%
        function noisePatternPerEvent = generateNoisePatternPerEvent(obj, eventList)
            noiseScale = eventList.Activity * obj.simProperties.eventToEventNoiseAmount;
            if obj.simProperties.eventToEventNoiseInOnlyClassifiedEvents
                eventsToClassify = obj.taskTable.ClassificationGroups(eventList.ID) ~= obj.taskTable.NON_CLASSIFIED_EVENT;
                noiseScale = noiseScale .* eventsToClassify(:);
                eventsToCorrelate = eventsToClassify(:);
            else
                eventsToCorrelate = ones(size(noiseScale));
            end
            eventsToCorrelate = logical(eventsToCorrelate);

            noisePatternPerEvent = rand(height(eventList), obj.simProperties.numVoxels) - 0.5;
            noisePatternPerEvent = noisePatternPerEvent .* noiseScale;
            noisePatternPerEvent(eventsToCorrelate, :) = makecorrelated(...
                noisePatternPerEvent(eventsToCorrelate, :),...
                obj.simProperties.eventToEventNoiseCoherence);
        end

        %%%
        function groundTruthTimeSeries = computeGroundTruthTimeSeries(obj, eventList, eventPatterns)
            eventList.ID = (1:height(eventList))';
            timePointsOfEachEvent = computeDesignMatrix(...
                eventList,...
                obj.simProperties.runDuration,...
                obj.simProperties.TR);

            groundTruthTimeSeries = zeros(obj.simProperties.numTRs, obj.simProperties.numVoxels);
            for ev = 1:height(eventList)
                thisEventTimePoints = logical(timePointsOfEachEvent(:,ev));
                groundTruthTimeSeries(thisEventTimePoints, :) = ...
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