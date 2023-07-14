
classdef fm_simulation < matlab.mixin.Copyable
    properties
        taskTable;
        simProperties;

        eventList;
        trialSequence;

        neuralPatterns;
        HRFs;

        neuralPatternPerEvent;
        neuralFluctuationPerEvent;
        totalNeuralActivityPerEvent;

        neuralTimeSeries;
        noNoiseBoldTimeSeries;
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
            obj.neuralPatterns = obj.generateNeuralPatterns();
            obj.HRFs     = obj.generateHRFs();

            for i = obj.simProperties.numRuns:-1:1
                [obj.eventList{i}, obj.trialSequence{i}] = obj.generateEventList();
                obj.eventList{i}                         = obj.scaleEventListByActivity(obj.eventList{i});

                obj.neuralPatternPerEvent{i}       = obj.computeNeuralPatternPerEvent(obj.eventList{i}, obj.neuralPatterns);
                obj.neuralFluctuationPerEvent{i}   = obj.generateNeuralFluctuationPerEvent(obj.eventList{i});
                obj.totalNeuralActivityPerEvent{i} = obj.neuralPatternPerEvent{i} + obj.neuralFluctuationPerEvent{i};
                obj.neuralTimeSeries{i}            = obj.computeNeuralTimeSeries(obj.eventList{i}, obj.totalNeuralActivityPerEvent{i});

                obj.noNoiseBoldTimeSeries{i} = obj.convolveWithHRFs(obj.neuralTimeSeries{i}, obj.HRFs);
                obj.boldTimeSeries{i} = obj.noNoiseBoldTimeSeries{i} + obj.generateNoise();
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
            eventList.Activity = obj.taskTable.NeuralIntensity(eventList.ID)';
        end

        %%%
        function neuralPatterns = generateNeuralPatterns(obj)
            neuralPatterns = rand(obj.taskTable.numNeuralPatterns, obj.simProperties.numVoxels);
        end

        %%%
        function neuralPatternPerEvent = computeNeuralPatternPerEvent(obj, eventList, neuralPatterns)
            neuralPatternPerEvent = obj.taskTable.NeuralPatternIDs(eventList.ID);
            neuralPatternPerEvent = neuralPatterns(neuralPatternPerEvent, :);
            neuralPatternPerEvent = neuralPatternPerEvent .* eventList.Activity;
        end

        %%%
        function neuralFluctuationPerEvent = generateNeuralFluctuationPerEvent(obj, eventList)
            noiseScale = eventList.Activity * obj.simProperties.neuralFluctuationAmount;
            if obj.simProperties.neuralFluctuationInOnlyClassifiedEvents
                eventsToClassify = obj.taskTable.ClassificationGroups(eventList.ID) ~= obj.taskTable.NON_CLASSIFIED_EVENT;
                noiseScale = noiseScale .* eventsToClassify(:);
                eventsToCorrelate = eventsToClassify(:);
            else
                eventsToCorrelate = ones(size(noiseScale));
            end
            eventsToCorrelate = logical(eventsToCorrelate);

            neuralFluctuationPerEvent = rand(height(eventList), obj.simProperties.numVoxels) - 0.5;
            neuralFluctuationPerEvent = neuralFluctuationPerEvent .* noiseScale;
            neuralFluctuationPerEvent(eventsToCorrelate, :) = makecorrelated(...
                neuralFluctuationPerEvent(eventsToCorrelate, :),...
                obj.simProperties.neuralFluctuationCoherence);
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
        function HRFs = generateHRFs(obj)
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
        function noNoiseBoldTimeSeries = convolveWithHRFs(~, timeSeries, HRFs)
            noNoiseBoldTimeSeries = convolveByColumn(timeSeries, HRFs);
        end

        %%%
        function noise = generateNoise(obj)
            noise = fm_noiseMachine( ...
                obj.simProperties.numTRs,...
                obj.simProperties.numVoxels,...
                obj.simProperties.TR,...
                obj.simProperties.noiseSD,...
                obj.simProperties.noiseSources{:} ...
                ).generateNoise();
        end
    end

end