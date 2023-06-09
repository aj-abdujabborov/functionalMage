

classdef functionalMage < matlab.mixin.Copyable % inherits from handle
    properties (Access = public)
        simulationProperties = struct('TR', 1, ...
            'numVoxels', 30, ...
            'numRuns', 5);
        hrfLibrary = 'SimTB';

        % analysis = fm_analysis();
        simulation = fm_simulation();
    end
    
    properties (GetAccess = public, SetAccess = private)
        analysisType;
        heavy;
    end

    properties (GetAccess = public, SetAccess = private, NonCopyable)
        results;
    end

    methods
        function obj = functionalMage()
        end

        function props = get.simulationProperties(obj)
            props = obj.simulationProperties;
        end

        function set.hrfLibrary(obj, hrfLib)
            obj.validateStringInput(hrfLib, {'SimTB', 'KendrickKay'});
        end

        function set.simulationProperties(obj, props)
            if (~obj.validateSimulationProperties(props))
                error('Invalid input.');
            end
            obj.simulationProperties = props;
        end
    end

    methods (Access = public)
        function run(obj)
            obj.results = repelem(obj.simulationProperties.numRuns, 1, 5);
        end

        function loadUp(obj)
            obj.heavy = randn(10000);
        end

        function loadDown(obj)
            obj.heavy = [];
        end
    end

    methods (Access = protected, Static = true)
        function choice = validateStringInput(string, options)
            match = strcmpi(string, options);
            if ~(any(match))
                error('Input is invalid.')
            end
            choice = options{match};
        end

        function output = validateSimulationProperties(props)
            output = numel(props.numVoxels) == 1;
        end
    end
end