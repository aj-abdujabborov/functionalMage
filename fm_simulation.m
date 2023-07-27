
classdef fm_simulation < matlab.mixin.Copyable
    properties
        eventList;
        simProperties;

        neuralPatterns;
        HRFs;

        neuralPatternPerEvent;
        neuralFluctuationPerEvent;
        totalNeuralActivityPerEvent;

        neuralTimeSeries;
        noNoiseBoldTimeSeries;
        boldTimeSeries = fm_data.empty();
    end
 
    methods
        %%%
        function obj = fm_simulation(eventList, simulationProperties)
            if (nargin > 0)
                obj.eventList = eventList;
                obj.simProperties = simulationProperties;
            end
        end

        %%%
        function generate(obj)
            % TODO: verify eventList and simProperties are set
            obj.neuralPatterns = obj.generateNeuralPatterns();
            if isempty(obj.simProperties.hrfSet)
                obj.HRFs = obj.generateHRFs();
            else
                obj.HRFs = obj.simProperties.hrfSet;
            end

            for i = obj.simProperties.numRuns:-1:1
                obj.neuralPatternPerEvent{i}       = obj.computeNeuralPatternPerEvent(obj.eventList(i), obj.neuralPatterns);
                obj.neuralFluctuationPerEvent{i}   = obj.generateNeuralFluctuationPerEvent(obj.eventList(i));
                obj.totalNeuralActivityPerEvent{i} = obj.neuralPatternPerEvent{i} + obj.neuralFluctuationPerEvent{i};
                obj.neuralTimeSeries{i}            = obj.computeNeuralTimeSeries(obj.eventList(i), obj.totalNeuralActivityPerEvent{i});

                obj.noNoiseBoldTimeSeries{i} = obj.convolveWithHRFs(obj.neuralTimeSeries{i}, obj.HRFs);
                
                obj.boldTimeSeries(i) = fm_data(...
                    obj.noNoiseBoldTimeSeries{i} + obj.generateNoise(),...
                    obj.simProperties.TR);
            end
        end

        %%%
        function neuralPatterns = generateNeuralPatterns(obj)
            numNeuralPatterns = max(cat(1, obj.eventList.ID));
            neuralPatterns = rand(numNeuralPatterns, obj.simProperties.numVoxels);
        end

        %%%
        function neuralPatternPerEvent = computeNeuralPatternPerEvent(~, eventList, neuralPatterns)
            neuralPatternPerEvent = neuralPatterns(eventList.ID, :);
            neuralPatternPerEvent = neuralPatternPerEvent .* eventList.Activity;
        end

        %%%
        function neuralFluctuationPerEvent = generateNeuralFluctuationPerEvent(obj, eventList)
            noiseScale = eventList.Activity * obj.simProperties.neuralFluctuationAmount;
            neuralFluctuationPerEvent = rand(height(eventList), obj.simProperties.numVoxels) - 0.5;
            neuralFluctuationPerEvent = neuralFluctuationPerEvent .* noiseScale;
            neuralFluctuationPerEvent = makecorrelated(...
                neuralFluctuationPerEvent,...
                obj.simProperties.neuralFluctuationCoherence);
        end

        %%%
        function neuralTimeSeries = computeNeuralTimeSeries(obj, eventList, eventPatterns)
            eventList.ID = (1:height(eventList))';
            timePointsOfEachEvent = eventList.computeDesignMatrix(...
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
            opts = [];
            opts.correlationWithCanonical = obj.simProperties.hrfsCorrelationWithDoubleGammaCanonical;
            opts.library = obj.simProperties.hrfLibrary;

            cacheLocation = fullfile(functionalMage.getCacheDirectory(), 'hrfDatabase.mat');

            HRFs = fm_hrf.getHrfsCache(obj.simProperties.TR,...
                                       obj.simProperties.numVoxels,...
                                       cacheLocation,...
                                       opts);
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