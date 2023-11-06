classdef fm_noise < handle
%FM_NOISE Generate fMRI noise
% Generate a combination of AR1, Gaussian, Rician and Physiological noise
%
% Input properties
%   <numTRs>
%   <numVoxels>
%   <TR>
%   <goalNoiseSD> Desired standard deviation of the simulated noise
%   
%   <AR1> Relative proportion of first-order autoregressive Guassian noise.
%     Default is 1. The autoregressive parameter RHO is set to 0.12.
%   <Gaussian> Relative proportion of Gaussian noise. Default is 0
%   <Rician> Relative proportion of Rician noise. Default is 0
%   <Physiological> Relative proportion of physiological noise. Default is 0
%
% Constructors
%   > obj = fm_noise()
%   > obj = fm_noise(numTRs, numVoxels, TR, goalNoiseSD)
%   > obj = fm_noise(numTRs, numVoxels, TR, goalNoiseSD, noise1Name,
%     noise1Amount, noise2Name, noise2Amount...)
%
% Methods
%   > noise = obj.generateNoise() returns a [numTRs by numVoxels] matrix
%     with randomly generated noise
%
% Static methods
%     These allow you to generate the noise sources directly
%   > AR1 = generateAR1Noise(numTRs, numVoxels, ~, rho) returns AR1
%     noise with SD of 1. Optionally you can supply the parameter RHO
%   > Gauss = generateGaussianNoise(numTRs, numVoxels) returns Gaussian
%     noise with SD of 1
%   > Rician = generateRicianNoise(numTRs, numVoxels) returns Rician noise
%     with SD of 1
%   > Phys = generatePhysiologicalNoise(numTRs, numVoxels, TR) returns
%     physiological noise with SD of 1
%
%   References
%   + AR1 process is based on: Mumford, J. A., Turner, B. O., Ashby, F. G., &
%     Poldrack, R. A. (2012). Deconvolving BOLD activation in event-related
%     designs for multivoxel pattern classification analyses. NeuroImage,
%     59(3), 2636–2643
%   + Rician and physiological are based on: Welvaert, M., Durnez, J.,
%     Moerkerke, B., Verdoolaege, G., & Rosseel, Y. (2011). neuRosim: An R
%     Package for Generating fMRI Data. Journal of statistical software,
%     44, 1-18.
%
% Part of package funkyMage. November 2023.
% https://github.com/aj-abdujabborov/funkyMage

    properties
        numTRs {mustBeNonnegative};
        numVoxels {mustBePositive};
        TR;
        goalNoiseSD {mustBeNonnegative};

        AR1 = 0;
        Gaussian = 0;
        Rician = 0;
        Physiological = 0;
    end


    properties (Access = protected, Constant = true)
        noiseNames = {'Rician', 'Gaussian', 'AR1', 'Physiological'};
        generatorHandles = {@fm_noise.generateRicianNoise,...
                            @fm_noise.generateGaussianNoise,...
                            @fm_noise.generateAR1Noise, ...
                            @fm_noise.generatePhysiologicalNoise};
    end


    methods 
        function obj = fm_noise(numTRs, numVoxels, TR, goalNoiseSD, varargin)
            if nargin == 0
                return;
            end

            obj.numTRs = numTRs;
            obj.numVoxels = numVoxels;
            obj.goalNoiseSD = goalNoiseSD;
            obj.TR = TR;
            
            p = inputParser();
            for i = 1:length(obj.noiseNames)
                p.addParameter(obj.noiseNames{i}, 0);
            end
            p.parse(varargin{:});

            for i = 1:length(obj.noiseNames)
                obj.(obj.noiseNames{i}) = p.Results.(obj.noiseNames{i});
            end
        end

        function outNoise = generateNoise(obj)
            for i = length(obj.noiseNames):-1:1
                noiseProportions(i) = obj.(obj.noiseNames{i});
            end
            noiseProportions = noiseProportions ./ sum(noiseProportions);

            noiseSourceIdx = find(noiseProportions > 0);
            numNoiseSources = length(noiseSourceIdx);
            
            outNoise = nan(obj.numTRs, obj.numVoxels, numNoiseSources);
            for i = 1:numNoiseSources
                outNoise(:,:,i) = obj.generatorHandles{noiseSourceIdx(i)}(obj.numTRs, obj.numVoxels, obj.TR);
            end
            outNoise = outNoise .* moveTo3rdDimension(noiseProportions(noiseSourceIdx));
            outNoise = fm_data(sum(outNoise, 3) * obj.goalNoiseSD);

            function out = moveTo3rdDimension(in)
                out = permute(in(:), [2 3 1]);
            end
        end
    end


    methods % Set methods
        function set.numTRs(obj, numTRs)
            assert(numel(numTRs) == 1 & isRound(numTRs) & numTRs > 0,...
                "Number of TRs should be a whole number larger than 0.")
            obj.numTRs = numTRs;
        end

        function set.numVoxels(obj, numVoxels)
            assert(numel(numVoxels) == 1 & isRound(numVoxels) & numVoxels > 0,...
                "Number of voxels should be a whole number larger than 0.")
            obj.numVoxels = numVoxels;
        end

        function set.goalNoiseSD(obj, noiseSD)
            assert(numel(noiseSD) == 1 & noiseSD >= 0, "goalNoiseSD should be a single non-negative value.")
            if (noiseSD == 0)
                warning("goalNoiseSD being exactly 0 can cause analysis issues later on. Setting it to 0.001");
                noiseSD = 0.001;
            end
            obj.goalNoiseSD = noiseSD;
        end

        function set.TR(obj, TR)
            assert(numel(TR) == 1 & TR >= 0, "TR should be a single non-negative value.")
            obj.TR = TR;
        end

        function set.Rician(obj, noiseAmount)
            obj.Rician = fm_noise.validateNoiseAmount('Rician', noiseAmount);
        end
        function set.Gaussian(obj, noiseAmount)
            obj.Gaussian = fm_noise.validateNoiseAmount('Gaussian', noiseAmount);
        end
        function set.AR1(obj, noiseAmount)
            obj.AR1 = fm_noise.validateNoiseAmount('AR1', noiseAmount);
        end
        function set.Physiological(obj, noiseAmount)
            obj.Physiological = fm_noise.validateNoiseAmount('Physiological', noiseAmount);
        end
    end


    methods (Static = true) % all outputs below will have SD = 1
        function outNoise = generateAR1Noise(numTRs, numVoxels, ~, rho)
            if ~exist('rho', 'var')
                rho = 0.12;
            end
            correlationMatrix = rho.^(toeplitz(1:numTRs) - 1);
            diagonalMatrix = diag(ones(1, numTRs));
            covarianceMatrix = diagonalMatrix * correlationMatrix * diagonalMatrix;

            means = zeros(numVoxels, numTRs);
            outNoise = mvnrnd(means, covarianceMatrix)';
            outNoise = 1/mean(std(outNoise)) * outNoise;
        end

        function outNoise = generateGaussianNoise(numTRs, numVoxels, ~)
            outNoise = randn([numTRs, numVoxels]);
            outNoise = 1/mean(std(outNoise)) * outNoise;
        end

        function outNoise = generateRicianNoise(numTRs, numVoxels, ~)
            outSize = [numTRs, numVoxels];
            theta = 1; vee = 1; % non-centrality parameter; 1 is the default in neuRosim
            noise1 = randn(outSize) + vee*cos(theta);
            noise2 = randn(outSize) + vee*sin(theta);
            outNoise = sqrt(noise1.^2 + noise2.^2);
            outNoise = 1/mean(std(outNoise)) * outNoise;
        end

        function outNoise = generatePhysiologicalNoise(numTRs, numVoxels, TR)
            heartHz = 1.17;
            respirationHz = 0.2;
            timePts = 1:numTRs;
            outNoise = sin(2 * pi * heartHz * TR * timePts) + cos(2 * pi * respirationHz * TR * timePts);
            outNoise = 1/std(outNoise) * outNoise;
            outNoise = repmat(outNoise', [1 numVoxels]);
        end
    end


    methods (Static = true, Access = protected)
        function noiseAmount = validateNoiseAmount(noiseName, noiseAmount)
            assert(numel(noiseAmount) == 1 & noiseAmount >= 0,...
                noiseName + " should be a single non-negative value.");
        end
    end

end