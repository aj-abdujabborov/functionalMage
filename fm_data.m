classdef fm_data < matlab.mixin.Copyable
    properties
        data;
        IDs;
        TR;

        rowName;
        dataName;
    end

    properties (SetAccess = private, Dependent = true)
        numVoxels;
        numTRs;
        duration;
        numRows;
    end

    methods
        function obj = fm_data(data, TR, IDs)
            if nargin >= 1, obj.data = data; end
            if nargin >= 2, obj.TR = TR;     end
            if nargin >= 3, obj.IDs = IDs;   end
        end

        function save2File(obj)

        end
    end

    methods % Get methods
        function numVoxels = get.numVoxels(obj)
            numVoxels = size(obj.data, 2);
        end

        function numRows = get.numRows(obj)
            numRows = size(obj.data, 1);
        end
        
        function numTRs = get.numTRs(obj)
            numTRs = obj.numRows;
        end

        function runDuration = get.duration(obj)
            runDuration = obj.numTRs * obj.TR;
        end
    end

    methods % Operator overloading
        % TODO: finish if overloading becomes very beneficial
        %{
        function outObj = plus(inObj1, inObj2)
            outObj = inObj1;

            if isempty(inObj1.TR)
                outObj.TR = inObj2.TR;
            elseif isempty(inObj2.TR)
                outObj.TR = inObj1.TR;
            else
                assert(inObj1.TR == inObj1.TR,...
                    "TRs of two fm_data objects are different.");
                outObj.TR = inObj1.TR;
            end

            if isempty(inObj1.IDs)
                outObj.IDs = inObj2.IDs;
            elseif isempty(inObj2.TR)
                outObj.IDs = inObj1.IDs;
            else
                assert(inObj1.IDs == inObj1.IDs,...
                    "IDs of two fm_data objects are different.");
                outObj.IDs = inObj1.IDs;
            end

        end
        %}
    end

    methods (Static = true)
        function loadFromFile(obj)
            
        end
    end
end
