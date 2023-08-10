
classdef fm_noiseMachine < handle
    properties
        numTRs;
        numVoxels;
        TR;
        goalNoiseSD;

        Rician = 0;
        Gaussian = 0;
        AR1 = 0;
        Physiological = 0;
    end


    properties (Access = protected, Constant = true)
        noiseNames = {'Rician', 'Gaussian', 'AR1', 'Physiological'};
        generatorHandles = {@fm_noiseMachine.generateRicianNoise,...
                           @fm_noiseMachine.generateGaussianNoise,...
                           @fm_noiseMachine.generateAR1Noise, ...
                           @fm_noiseMachine.generatePhysiologicalNoise};
    end


    methods 
        function obj = fm_noiseMachine(numTRs, numVoxels, TR, goalNoiseSD, varargin)
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
            obj.Rician = fm_noiseMachine.validateNoiseAmount('Rician', noiseAmount);
        end
        function set.Gaussian(obj, noiseAmount)
            obj.Gaussian = fm_noiseMachine.validateNoiseAmount('Gaussian', noiseAmount);
        end
        function set.AR1(obj, noiseAmount)
            obj.AR1 = fm_noiseMachine.validateNoiseAmount('AR1', noiseAmount);
        end
        function set.Physiological(obj, noiseAmount)
            obj.Physiological = fm_noiseMachine.validateNoiseAmount('Physiological', noiseAmount);
        end
    end


    methods (Static = true) % all outputs below will have SD = 1
        function outNoise = generateRicianNoise(numTRs, numVoxels, ~)
            outSize = [numTRs, numVoxels];
            theta = 1; vee = 1; % non-centrality parameter; 1 is the default in neuRosim
            noise1 = randn(outSize) + vee*cos(theta);
            noise2 = randn(outSize) + vee*sin(theta);
            outNoise = sqrt(noise1.^2 + noise2.^2);
            outNoise = 1/mean(std(outNoise)) * outNoise;
        end

        function outNoise = generateGaussianNoise(numTRs, numVoxels, ~)
            outNoise = randn([numTRs, numVoxels]);
            outNoise = 1/mean(std(outNoise)) * outNoise;
        end

        function outNoise = generateAR1Noise(numTRs, numVoxels, ~)
            rho = 0.3; % Set to 0.12 in Mumford et al. ()
            correlationMatrix = rho.^(toeplitz(1:numTRs) - 1);
            diagonalMatrix = diag(ones(1, numTRs));
            covarianceMatrix = diagonalMatrix * correlationMatrix * diagonalMatrix;

            means = zeros(numVoxels, numTRs);
            outNoise = mvnrnd(means, covarianceMatrix)';
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


    methods (Static = true, Access = private)
        function noiseAmount = validateNoiseAmount(noiseName, noiseAmount)
            assert(numel(noiseAmount) == 1 & noiseAmount >= 0,...
                noiseName + " should be a single non-negative value.");
        end
    end

end