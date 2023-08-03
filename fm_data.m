classdef fm_data < matlab.mixin.Copyable
    properties
        data {mustBeNumeric};
        TR {mustBePositive};

        rowName {mustBeText};
        dataName {mustBeText};
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
            if nargin >= 2, obj.TR   = TR;   end
            if nargin >= 3, obj.IDs  = IDs;  end
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

        function set.data(obj, data)
            assert(width(data) > 0, "Input data must have more than 0 columns.");
            assert(height(data) > 0, "Input data must have more than 0 rows.");
            obj.data = data;
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
        function TR = getTR(objVector)
            TR = objVector(1).TR;
        end

        function dataCombined = cat(objVector)
            uniqueTR = unique([objVector.TR]);
            assert(length(uniqueTR) <= 1, "All TRs must be the same.");
            assert(has1UniqueValue([objVector.numVoxels]), "Number of voxels must be the same for all elements.");

            dataCombined = fm_data(cat(1, objVector.data), ...
                                   uniqueTR,...
                                   cat(1, objVector.IDs));

            if has1UniqueString({objVector.rowName})
                dataCombined.rowName = objVector(1).rowName;
            end
            if has1UniqueString({objVector.dataName})
                dataCombined.dataName = objVector(1).dataName;
            end
            
            %________________
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
