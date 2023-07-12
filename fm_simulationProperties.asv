
classdef fm_simulationProperties < matlab.mixin.Copyable
    properties (GetAccess = public, SetAccess = public)
        TR = 1;
        numRuns = 5;
        numVoxels = 30;
        runDuration = 300;

        hrfLibrary = 'SimTB';
        hrfsCorrelationWithDoubleGammaCanonical = [0.6 1];

        itiModel = 'fixed';
        itiParams = 5;

        noiseSd = 2;
        noiseSources = {'tempcorr'};

        eventToEventNoiseAmount = 2;
        eventToEventNoiseCoherence = 0.5;
        eventToEventNoiseInOnlyClassifiedEvents = false;
    end

    properties (GetAccess = public, SetAccess = private, Dependent)
        numTRs;
    end
 
    methods % Set methods
        function set.TR(obj, TR)
            assert(TR > 0, "TR has to be more than 0 seconds.");
            obj.TR = TR;
        end
        function set.numRuns(obj, nRuns)
            assert(nRuns > 0, "Number of runs should be more than 0.");
            obj.numRuns = nRuns;
        end
    end

    methods % Get methods
        function numTRs = get.numTRs(obj)
            numTRs = obj.runDuration ./ obj.TR;
        end
    end
end