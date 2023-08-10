
classdef fm_simulation < matlab.mixin.Copyable
    properties
        eventList fm_eventList;
        simProperties fm_simulationProperties;

        neuralPatterns (:,:) {mustBeFinite};
        HRFs (:,:) {mustBeFinite};

        neuralPatternPerEvent (1,:) fm_data;
        neuralFluctuationPerEvent (1,:) fm_data;
        totalNeuralActivityPerEvent (1,:) fm_data;

        neuralTimeSeries (1,:) fm_data;
        noNoiseBoldTimeSeries (1,:) fm_data;
        boldTimeSeries (1,:) fm_data;
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
            assert(~isempty(obj.eventList), "Set eventList property");
            assert(~isempty(obj.simProperties), "Set simProperties property");

            obj.neuralPatterns = obj.generateNeuralPatterns();
            if isempty(obj.simProperties.hrfSet)
                obj.HRFs = obj.generateHRFs();
            else
                obj.HRFs = obj.simProperties.hrfSet;
            end

            for i = obj.simProperties.numRuns:-1:1
                obj.neuralPatternPerEvent(i)       = obj.computeNeuralPatternPerEvent(obj.eventList(i), obj.neuralPatterns);
                obj.neuralFluctuationPerEvent(i)   = obj.generateNeuralFluctuationPerEvent(obj.eventList(i));
                obj.totalNeuralActivityPerEvent(i) = obj.neuralPatternPerEvent(i) + obj.neuralFluctuationPerEvent(i);
                obj.neuralTimeSeries(i)            = obj.computeNeuralTimeSeries(obj.eventList(i), obj.totalNeuralActivityPerEvent(i));

                obj.noNoiseBoldTimeSeries(i) = obj.convolveWithHRFs(obj.neuralTimeSeries(i), obj.HRFs);
                
                obj.boldTimeSeries(i) = obj.noNoiseBoldTimeSeries(i) + obj.generateNoise();
                obj.boldTimeSeries(i).TR = obj.simProperties.TR;
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
            neuralPatternPerEvent = fm_data(neuralPatternPerEvent .* eventList.Activity);
        end

        %%%
        function neuralFluctuationPerEvent = generateNeuralFluctuationPerEvent(obj, eventList)
            if obj.simProperties.neuralFluctuationAmount == 0
                neuralFluctuationPerEvent = fm_data(zeros(height(eventList), obj.simProperties.numVoxels));
                return;
            end
            
            noiseScale = eventList.Activity * obj.simProperties.neuralFluctuationAmount;
            neuralFluctuationPerEvent = rand(height(eventList), obj.simProperties.numVoxels) - 0.5;
            neuralFluctuationPerEvent = neuralFluctuationPerEvent .* noiseScale;
            neuralFluctuationPerEvent = fm_data( ...
                    makecorrelated(...
                        neuralFluctuationPerEvent,...
                        obj.simProperties.neuralFluctuationCoherence) ...
                        );
        end

        %%%
        function neuralTimeSeries = computeNeuralTimeSeries(obj, eventList, eventPatterns)
            timePointsOfEachEvent = eventList.openIntoUniqueEvents().computeDesignMatrix(obj.simProperties.TR);
            
            neuralTimeSeries = fm_data(...
                zeros(obj.simProperties.numTRs, obj.simProperties.numVoxels),...
                obj.simProperties.TR);

            for ev = 1:height(eventList)
                thisEventTimePoints = logical(timePointsOfEachEvent(:,ev));
                neuralTimeSeries.data(thisEventTimePoints, :) = ...
                    repmat(eventPatterns.data(ev,:), sum(thisEventTimePoints), 1);
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
        function convTimeSeries = convolveWithHRFs(~, timeSeries, HRFs)
            convTimeSeries = fm_data(convolveByColumn(timeSeries.data, HRFs));
            convTimeSeries.TR = timeSeries.TR;
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