
classdef functionalMage < matlab.mixin.Copyable % inherits from handle
    properties (Access = public, Dependent)
        taskTable;
    end
    
    properties (Access = public)
        simProperties = fm_simulationProperties();
        fmriData;

        nSubj = 10;
    end

    properties (Access = private)
        taskTableReal = fm_taskTable();
    end
    
    properties (GetAccess = public, SetAccess = private, NonCopyable)
    end

    methods
        function simulate(obj)
            
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
    methods (Access = protected, Static = true)

    end


end