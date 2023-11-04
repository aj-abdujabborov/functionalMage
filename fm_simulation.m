classdef fm_simulation < matlab.mixin.Copyable
%FM_SIMULATION Simulate fMRI data
% Generate an fMRI dataset of with TR, variable ground-truth HRFs and noise
%
% Input properties
%   <eventList> is a vector of fm_eventList objects where each object
%   represents a run. Each 'ID' value in eventList is a random pattern of
%   neural activity and each 'Activity' value is used to set the intensity of
%   of neural activity. Generally, 'eventList' should be extracted from an
%   'fm_designMatrix' object using the getNeuralPatternIDEventList()
%   method.
%   
%   <TR> of simulated data
%   <numRuns>
%   <numVoxels>
%   <runDuration> in seconds
%   <hrfLibrary> can be 'SimTB', 'NSD', or 'custom' (HRFs from Natural Scenes Dataset paper)
%   <hrfSet> is a matrix where each column is an HRF. If 'hrfLibrary' is
%   custom, you need to set this yourself.
%   <hrfsCorrelationWithDoubleGammaCanonical> is a 1x2 vector of values
%   between 0 and 1, indicating the range of correlations that 'hrfSet'
%   HRFs sould have with the double gamma canonical HRF. These correlations
%   will be distributed roughly uniformly.
%
%   <noiseSD> The standard deviation of noise added to the simulated data.
%   A value of 2-4 is *approximately* realistic for an event-related
%   working memory study I based it off, but what's realistic will vary.
%   Ideally you'd simulate under various noise levels and verify that your
%   conclusion is qualitatively the same. Note that the default noiseSD
%   assumes neural intensity in 'fm_task' will be a maximum of 1.
%   <noiseSources> is a cell array of noiseSources to include, with its
%   contents coming in pairs. The first value is the name of the source of
%   noise and the second the relative quantity.
%       Example: {'AR1', 1, 'Gaussian', 2}
%       This will create twice as much Gaussian noise as it will AR1 noise
%       and in combination they'll have a standard deviation of 'noiseSD'.
%       
%       Noise sources: 'AR1', 'Gaussian', 'Rician', 'Physiological'. See
%       fm_noiseMachine for more info.
%
%   <neuralFluctuationAmount> is the amount of random event-to-event
%   fluctuation (or "neural noise") to generate. A value of 2 means that
%   every event will have a pattern of noise drawn from a uniform
%   distribution whose minimum is (-1 * its neural intensity * 2) and
%   maximum is (its neural intensity * 2).
%   <neuralFluctuationCoherence> is how inter-correlated the neural
%   fluctuation is among voxels. This ranges from 0 to 1.
%
% Output properties
%   <boldTimeSeries> A vector of fm_data objects containing the final
%   simulated data, including noise.
%
% Other output properties
%   <neuralPatterns> is an N x numVoxels matrix, where each row contains
%   ground-truth neural patterns. N is equal to the number of unique IDs in
%   the input eventList property.
%   <neuralPatternPerEvent> is an fm_data vector, each object containing
%   the neural patterns of each event. Thus, the number of rows in this
%   data are the same as in eventList.
%   <neuralFluctuationPerEvent> contains the random fluctuation patterns
%   per event.
%   <totalNeuralActivityPerEvent> contains the final neural activity per
%   event (including neural fluctuations).
%   <neuralTimeSeries> A vector of fm_data objects containing "neural" time
%   series data (generation before convolution and noise-addition).
%   <noNoiseBoldTimeSeries> A vector of fm_data objects containing
%   convolved data but before noise is added.
%
% Methods
%   > sim = fm_simulation(eventList) returns an fm_simulation object.
%   > generate() will simulate the fMRI data with all the specified
%   properties. It should be launched once all the properties are set.
%
% Examples
%   eventList = dm.getNeuralPatternIDEventList(); 'dm' is a fm_designMatrix object
%   sim = fm_simulation();
%   sim.eventList = eventList;
%   sim.TR = 0.5;
%   sim.noiseSD = 4;
%   sim.generate();
%   disp(sim.boldTimeSeries()) % outputs the generated fMRI time series
%
% Part of package funkyMage. November 2023.
% https://github.com/aj-abdujabborov/funkyMage

    properties
        eventList fm_eventList;
        TR (1,1) {mustBePositive} = 1;
        numVoxels (1,1) {mustBePositive} = 30;
        hrfLibrary (1,1) {matlab.system.mustBeMember(hrfLibrary, {'simtb', 'nsd'})} = "simtb";
        hrfSet (:,:) = [];
        hrfsCorrelationWithDoubleGammaCanonical (1,2) = [0.6 1];
        noiseSD (1,1) {mustBeNonnegative} = 2;
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

    properties (Access = private)
        numRuns;
        runDuration;
    end

    properties (Access = private, Dependent = true)
        numTRs;
    end
 
    methods
        %%%
        function obj = fm_simulation(eventList)
            if (nargin > 0)
                obj.eventList = eventList;
            end
        end

        %%%
        function generate(obj)
            assert(~isempty(obj.eventList), "Set eventList property");

            obj.neuralPatterns = obj.generateNeuralPatterns();
            if isempty(obj.hrfSet)
                obj.hrfSet = obj.generateHRFs();
            end

            obj.runDuration = obj.eventList(1).runDuration;
            obj.numRuns = length(obj.eventList);
            for i = obj.numRuns:-1:1
                obj.neuralPatternPerEvent(i)       = obj.computeNeuralPatternPerEvent(obj.eventList(i), obj.neuralPatterns);
                obj.neuralFluctuationPerEvent(i)   = obj.generateNeuralFluctuationPerEvent(obj.eventList(i));
                obj.totalNeuralActivityPerEvent(i) = obj.neuralPatternPerEvent(i) + obj.neuralFluctuationPerEvent(i);
                obj.neuralTimeSeries(i)            = obj.computeNeuralTimeSeries(obj.eventList(i), obj.totalNeuralActivityPerEvent(i));

                obj.noNoiseBoldTimeSeries(i) = obj.convolveWithHRFs(obj.neuralTimeSeries(i), obj.hrfSet);
                
                obj.boldTimeSeries(i) = obj.noNoiseBoldTimeSeries(i) + obj.generateNoise();
                obj.boldTimeSeries(i).TR = obj.TR;
            end
        end

        %%%
        function neuralPatterns = generateNeuralPatterns(obj)
            numNeuralPatterns = max(cat(1, obj.eventList.ID));
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

            cacheLocation = fullfile(functionalMage.getCacheDirectory(), 'hrfDatabase.mat');

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