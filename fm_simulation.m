classdef fm_simulation < matlab.mixin.Copyable
%FM_SIMULATION Simulate fMRI data
% Generate fMRI data with variable HRFs and various noise sources
%
% Input properties
%   <neurIDEventList> is a vector of fm_eventList objects where each object
%     represents a run. Each ID value in neurIDEventList is a random
%     pattern of neural activity and each 'Activity' value sets the
%     intensity of neural activity. 'neurIDEventList' is extracted from an
%     'fm_designMatrix' object using the getNeuralPatternIDEventList()
%     method.
%   
%   <TR> of simulated data. Default is 1.
%   <numVoxels> Default is 30.
%   <hrfLibrary> can be 'SimTB', 'NSD', or 'custom' (HRFs from Natural
%     Scenes Dataset paper). Default is 'SimTB'.
%   <hrfSet> is a matrix where each column is an HRF. If 'hrfLibrary' is
%     custom, you need to set this yourself. Default is empty.
%   <hrfsCorrelationWithDoubleGammaCanonical> is a 2-element vector with
%     values between 0 and 1, indicating the range of correlations that
%     'hrfSet' HRFs sould have with the double gamma canonical HRF. The
%     correlations are distributed roughly uniformly. Default: [0.6 1.0]
%
%   <noiseSD> The standard deviation of noise added to the simulated data.
%     A value of 2-4 is *approximately* realistic for an event-related
%     working memory study I based it off. Ideally you'd simulate under
%     various noise levels and verify that your conclusion is qualitatively
%     the same. Default is 3.0
%   <noiseSources> A cell array of noise sources to include. The entries
%     come in pairs. The first is the name of the source and second the
%     relative quantity. Noise sources are 'AR1', 'Gaussian', 'Rician',
%     'Physiological'. Do "help fm_noise" for more info.
%
%     Example: {'AR1', 1, 'Gaussian', 2}
%     Default is {'AR1', 1}.
%
%   <neuralFluctuationAmount> is the amount of random event-to-event
%     fluctuation (or "neural noise") to generate. A value of 2 means that
%     every event will have a pattern of noise drawn from a uniform
%     distribution whose minimum is [its neural intensity * -2] and maximum
%     is [its neural intensity * 2]. Default is 2.0
%   <neuralFluctuationCoherence> is how inter-correlated the neural
%     fluctuation is among voxels. This ranges from 0 to 1. Default is 0.5
%
% Output properties
%   <boldTimeSeries> A vector of fm_data objects containing the final
%     simulated data, including noise
%
% Other output properties
%   <neuralPatterns> is an N x numVoxels matrix, where each row contains
%     ground-truth neural patterns. N is equal to the number of unique IDs
%     in the neurIDEventList property.
%   <neuralPatternPerEvent> is an fm_data vector, where each object
%     contains the neural pattern for each event within a run.
%   <neuralFluctuationPerEvent> is an fm_data vector, where each object
%     contains the random fluctuation patterns for each event within a run
%   <totalNeuralActivityPerEvent> is an fm_data vector, where each object
%     contains the final neural activity patterns for each event within a
%     run
%   <neuralTimeSeries> is an fm_data vector containing "neural" time
%     series data (stage prior to convolution)
%   <noNoiseBoldTimeSeries> A vector of fm_data objects containing
%     convolved BOLD time series before noise is added
%
% Constructors
%   > obj = fm_simulation()
%   > obj = fm_simulation(neurIDEventList)
%
% Methods
%   > obj.go() launches the simulation
%
% Example
%   neurIDEventList = dm.getNeuralPatternIDEventList(); 
%       % 'dm' is a fm_designMatrix object
%   sim = fm_simulation();
%   sim.neurIDEventList = neurIDEventList;
%   sim.TR = 0.5;
%   sim.noiseSD = 4;
%   sim.go();
%   disp(sim.boldTimeSeries()) % outputs the generated fMRI time series
%
% Part of package funkyMage. November 2023.
% https://github.com/aj-abdujabborov/funkyMage

    properties
        neurIDEventList fm_eventList;
        TR (1,1) {mustBePositive} = 1;
        numVoxels (1,1) {mustBePositive} = 30;
        hrfLibrary (1,1) {matlab.system.mustBeMember(hrfLibrary, {'simtb', 'nsd'})} = "simtb";
        hrfSet (:,:) = [];
        hrfsCorrelationWithDoubleGammaCanonical (1,2) = [0.6 1];
        noiseSD (1,1) {mustBeNonnegative} = 3;
        noiseSources (1,:) {mustBeA(noiseSources, 'cell')} = {'AR1', 1};
        neuralFluctuationAmount (1,1) {mustBeNonnegative} = 2;
        neuralFluctuationCoherence (1,1) {mustBeNonnegative} = 0.5;

        neuralPatterns (:,:) {mustBeFinite};

        neuralPatternPerEvent (1,:) fm_data;
        neuralFluctuationPerEvent (1,:) fm_data;
        totalNeuralActivityPerEvent (1,:) fm_data;

        neuralTimeSeries (1,:) fm_data;
        noNoiseBoldTimeSeries (1,:) fm_data;
        boldTimeSeries (1,:) fm_data;
    end

    properties (Access = protected)
        numRuns;
        runDuration;
    end

    properties (Access = protected, Dependent = true)
        numTRs;
    end
 
    methods
        %%%
        function obj = fm_simulation(eventList)
            if (nargin > 0)
                obj.neurIDEventList = eventList;
            end
        end

        %%%
        function go(obj)
            assert(~isempty(obj.neurIDEventList), "Set neurIDEventList property");

            obj.neuralPatterns = obj.generateNeuralPatterns();
            if isempty(obj.hrfSet)
                obj.hrfSet = obj.generateHRFs();
            end

            obj.runDuration = obj.neurIDEventList(1).runDuration;
            obj.numRuns = length(obj.neurIDEventList);
            for i = obj.numRuns:-1:1
                obj.neuralPatternPerEvent(i)       = obj.computeNeuralPatternPerEvent(obj.neurIDEventList(i), obj.neuralPatterns);
                obj.neuralFluctuationPerEvent(i)   = obj.generateNeuralFluctuationPerEvent(obj.neurIDEventList(i));
                obj.totalNeuralActivityPerEvent(i) = obj.neuralPatternPerEvent(i) + obj.neuralFluctuationPerEvent(i);
                obj.neuralTimeSeries(i)            = obj.computeNeuralTimeSeries(obj.neurIDEventList(i), obj.totalNeuralActivityPerEvent(i));

                obj.noNoiseBoldTimeSeries(i) = obj.convolveWithHRFs(obj.neuralTimeSeries(i), obj.hrfSet);
                
                obj.boldTimeSeries(i) = obj.noNoiseBoldTimeSeries(i) + obj.generateNoise();
                obj.boldTimeSeries(i).TR = obj.TR;
            end
        end

        %%%
        function neuralPatterns = generateNeuralPatterns(obj)
            numNeuralPatterns = max(cat(1, obj.neurIDEventList.ID));
            neuralPatterns = rand(numNeuralPatterns, obj.numVoxels);
        end

        %%%
        function neuralPatternPerEvent = computeNeuralPatternPerEvent(~, eventList, neuralPatterns)
            neuralPatternPerEvent = neuralPatterns(eventList.ID, :);
            neuralPatternPerEvent = fm_data(neuralPatternPerEvent .* eventList.Activity);
        end

        %%%
        function neuralFluctuationPerEvent = generateNeuralFluctuationPerEvent(obj, eventList)
            if obj.neuralFluctuationAmount == 0
                neuralFluctuationPerEvent = fm_data(zeros(height(eventList), obj.numVoxels));
                return;
            end
            
            noiseScale = eventList.Activity * obj.neuralFluctuationAmount;
            neuralFluctuationPerEvent = rand(height(eventList), obj.numVoxels) - 0.5;
            neuralFluctuationPerEvent = neuralFluctuationPerEvent .* noiseScale;
            neuralFluctuationPerEvent = fm_data( ...
                    makecorrelated(...
                        neuralFluctuationPerEvent,...
                        obj.neuralFluctuationCoherence) ...
                        );
        end

        %%%
        function neuralTimeSeries = computeNeuralTimeSeries(obj, eventList, eventPatterns)
            timePointsOfEachEvent = eventList.openIntoUniqueEvents().computeDesignMatrix(obj.TR);
            
            neuralTimeSeries = fm_data(...
                zeros(obj.numTRs, obj.numVoxels),...
                obj.TR);

            for ev = 1:height(eventList)
                thisEventTimePoints = logical(timePointsOfEachEvent(:,ev));
                neuralTimeSeries.data(thisEventTimePoints, :) = ...
                    repmat(eventPatterns.data(ev,:), sum(thisEventTimePoints), 1);
            end
        end

        %%%
        function HRFs = generateHRFs(obj)
            opts = [];
            opts.correlationWithCanonical = obj.hrfsCorrelationWithDoubleGammaCanonical;
            opts.library = obj.hrfLibrary;

            cacheLocation = fullfile(funkyMage.getCacheDirectory(), 'hrfDatabase.mat');

            HRFs = fm_hrf.getHrfsCache(obj.TR,...
                                       obj.numVoxels,...
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
            noise = fm_noise( ...
                obj.numTRs,...
                obj.numVoxels,...
                obj.TR,...
                obj.noiseSD,...
                obj.noiseSources{:} ...
                ).generateNoise();
        end

        %%%
        function numTRs = get.numTRs(obj)
            numTRs = obj.runDuration ./ obj.TR;
        end
    end

end