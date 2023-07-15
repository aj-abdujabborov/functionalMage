
classdef functionalMage < matlab.mixin.Copyable
    properties (Access = public, Dependent)
        taskTable;
    end
    
    properties (Access = public)
        simProperties;
        simulation;

        nSubj;
    end

    properties (Access = private)
        privateTaskTable;
    end

    methods
        function obj = functionalMage()
            obj.privateTaskTable = fm_taskTable();
            obj.simProperties = fm_simulationProperties();
            obj.nSubj = 10;
        end

        function simulate(obj)
            % Check that taskTable and simProperties are ready
            obj.simulation = fm_simulation(obj.privateTaskTable, obj.simProperties);
            obj.simulation.generate();
        end
    end

    %%% Get and Set Methods
    methods
        function output = get.taskTable(obj)
            output = obj.privateTaskTable.content;
        end

        function set.taskTable(obj, newTaskTable)
            obj.privateTaskTable.content = newTaskTable;
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