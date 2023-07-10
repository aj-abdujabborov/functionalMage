
classdef fm_simulationProperties < matlab.mixin.Copyable
    properties (GetAccess = public, SetAccess = public)
        TR = 1;
        nRuns = 5;
        nVoxels = 30;
        runDur = 300;

        hrfLibrary = 'SimTB';

        itiModel = 'fixed';
        itiParameters = 5;

        noiseSd = 2;
        noiseSources = {'tempcorr'};
    end
 
    methods % Set methods
        function set.TR(obj, TR)
            assert(TR > 0, "TR has to be more than 0 seconds.");
            obj.TR = TR;
        end
        function set.nRuns(obj, nRuns)
            assert(nRuns > 0, "Number of runs should be more than 0.");
            obj.nRuns = nRuns;
        end
    end
end