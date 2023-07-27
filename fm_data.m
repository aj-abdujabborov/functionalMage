classdef fm_data < matlab.mixin.Copyable
    properties
        data;
        TR;

        rowName;
        dataName;
    end

    properties (Dependent = true)
        IDs;
    end

    properties (Dependent = true, SetAccess = private)
        numVoxels;
        numTRs;
        duration;
        numRows;
    end

    properties (Access = private)
        privateIDs;
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

    methods % Get and Set methods
        function set.IDs(obj, IDs)
            if isempty(IDs), return; end
            
            assert(height(IDs) == obj.numRows,...
                  "The number of rows in IDs should be the same as in data");
            assert(width(IDs) == 1,...
                   "IDs should be a column vector");
            obj.privateIDs = IDs;
        end

        function IDs = get.IDs(obj)
            IDs = obj.privateIDs;
        end

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

    methods % Overload MATLAB functions
        function dataCombined = cat(dataVector)
            assert(isa(dataVector, 'fm_data'), "Input is not of fm_data type.");
            
            uniqueTR = unique([dataVector.TR]);
            assert(length(uniqueTR) == 1, "All TRs must be the same.");
            assert(has1UniqueValue([dataVector.numVoxels]), "Number of voxels must be the same for all elements.");

            dataCombined = fm_data(cat(1, dataVector.data), ...
                                   uniqueTR,...
                                   cat(1, dataVector.IDs));
            if has1UniqueString({dataVector.rowName})
                dataCombined.rowName = dataVector(1).rowName;
            end
            if has1UniqueString({dataVector.dataName})
                dataCombined.dataName = dataVector(1).dataName;
            end
            
            function bool = has1UniqueValue(vec)
                bool = length(unique(vec)) == 1;
            end

            function bool = has1UniqueString(vec)
                if isempty(vec{1})
                    bool = false;
                else
                    vec = string(vec);
                    bool = has1UniqueValue(vec);
                end
            end
        end

        function dataCombined = plus(inObj1, inObj2)
            dataCombined = cat([inObj1, inObj2]);
        end
    end

    methods (Static = true)
        function loadFromFile(obj)
            
        end
    end
end
