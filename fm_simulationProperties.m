
classdef fm_simulationProperties < matlab.mixin.Copyable
    properties (GetAccess = public, SetAccess = public)
        TR = 1;
        numRuns = 5;
        numVoxels = 30;
        runDuration = 300;

        hrfLibrary = 'SimTB';
        hrfSet = [];
        hrfsCorrelationWithDoubleGammaCanonical = [0.6 1];

        itiModel = 'fixed';
        itiParams = 5;

        noiseSd = 2;
        noiseSources = {'tempcorr'};

        eventToEventNoiseAmount = 2; % rename eventToEventNoise to eventNoise EVERYWHERE
        eventToEventNoiseCoherence = 0.5;
        eventToEventNoiseInOnlyClassifiedEvents = false;
    end

    properties (GetAccess = public, SetAccess = private, Dependent)
        numTRs; % TODO: could set complementary set.runDuration, set.numTRs and set.TR methods that update the others
    end

    methods
        function obj = fm_simulationProperties()
            % TODO: implement construction through feeding a structure: obj.("name") works
            %       implement construction through feeding parameters
            %       
        end
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
        function set.hrfLibrary(obj, hrfLibrary)
            if ~any(strcmpi(hrfLibrary, {'simtb', 'kendrickkay', 'custom'}))
                error("hrfLibrary options are 'SimTB', 'KendrickKay' or 'custom");
            end
            if strcmpi(hrfLibrary, 'custom') && isempty(obj.hrfSet)
                error('stuff');
            end
            
            obj.hrfLibrary = hrfLibrary;
            % TODO: rename KendrickKay to maybe GLMSingle or NSD or
            % something.
        end
    end

    methods % Get methods
        function numTRs = get.numTRs(obj)
            numTRs = obj.runDuration ./ obj.TR;
        end
    end
end