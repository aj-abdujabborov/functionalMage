
classdef functionalMage < matlab.mixin.Copyable
    properties (Access = public, Dependent)
        taskTable;
    end
    
    properties (Access = public)
        simProperties;
        fmriData;

        nSubj;
    end

    properties (Access = private)
        taskTableReal;
        simulation;
    end
    
    properties (GetAccess = public, SetAccess = private, NonCopyable)
    end

    methods
        function obj = functionalMage()
            obj.taskTableReal = fm_taskTable();
            obj.simProperties = fm_simulationProperties();
            obj.nSubj = 10;
        end

        function simulate(obj)
            % Check that taskTable and simProperties are ready
            obj.simulation = fm_simulation(obj.taskTableReal, obj.simProperties);
            obj.simulation.generate();
        end
    end

    %%% Get and Set Methods
    methods
        function output = get.taskTable(obj)
            output = obj.taskTableReal.content;
        end

        function set.taskTable(obj, newTaskTable)
            obj.taskTableReal.content = newTaskTable;
        end
    end

    %%%
    methods (Access = protected)
        
    end

    %%% Static methods
    methods (Static = true)
        function output = getCacheDirectory()
            output = fullfile(fileparts(which('functionalMage.m')), 'cache/');
            if ~exist(output, 'dir')
                mkdir(output);
            end
        end
    end


end