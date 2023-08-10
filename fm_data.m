classdef fm_data < matlab.mixin.Copyable
    properties
        data {mustBeNumeric};
        TR {mustBePositive};
    end

    properties (Dependent = true)
        IDs;
    end

    properties
        rowName {mustBeText} = 'Undefined';
        dataName {mustBeText} = 'Undefined';
    end

    properties (Dependent = true, SetAccess = private)
        numTRs;
        numRows;
        numCols;
        numVoxels;
        duration;
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
            assert(width(data) >= 0, "Input data must have 0 or more columns.");
            assert(height(data) >= 0, "Input data must have 0 or more rows.");
            obj.data = data;
        end

        function IDs = get.IDs(obj)
            IDs = obj.privateIDs;
        end
        
        function numTRs = get.numTRs(obj)
            numTRs = obj.numRows;
        end

        function numRows = get.numRows(obj)
            numRows = size(obj.data, 1);
        end
        
        function numVoxels = get.numVoxels(obj)
            numVoxels = obj.numCols;
        end

        function numCols = get.numCols(obj)
            numCols = size(obj.data, 2);
        end

        function runDuration = get.duration(obj)
            runDuration = obj.numTRs * obj.TR;
        end
    end

    methods
        function TR = getTR(objVector)
            assert(fm_data.areAllSame([objVector.TR], length(objVector)),...
                   "fm_data objects do not all have the same TR");
            TR = objVector(1).TR;
        end

        function numVoxels = getNumVoxels(objVector)
            assert(fm_data.areAllSame([objVector.numVoxels], length(objVector)),...
                   "fm_data objects do not all have the same number of voxels");
            numVoxels = objVector(1).numVoxels;
        end

        function dataCombined = cat(objVector)
            dataCombined = combineCompatibleParameters(objVector);
            dataCombined.data = cat(1, objVector.data);
            dataCombined.IDs  = cat(1, objVector.IDs);
        end

        function dataCombined = diagCat(objVector)
            dataCombined = combineCompatibleParameters(objVector);
            dataCombined.data = zeros(sum([objVector.numTRs]), sum([objVector.numVoxels]));
            
            rowStart = 1;
            colStart = 1;
            for r = 1:length(objVector)
                rowEnd = rowStart + objVector(r).numRows - 1;
                colEnd = colStart + objVector(r).numVoxels - 1;
                dataCombined.data(rowStart:rowEnd, colStart:colEnd) = objVector(r).data;
                colStart = colEnd + 1;
                rowStart = rowEnd + 1;
            end

            dataCombined.IDs  = cat(1, objVector.IDs);
        end

        function dataCombined = plus(inObj1, inObj2)
            dataCombined = combineCompatibleParameters([inObj1, inObj2]);
            assert(inObj1.numRows == inObj2.numRows,...
                   "Number of rows / TRs in fm_data objects do not match");
            dataCombined.data = inObj1.data + inObj2.data;
            if isequal(inObj1.IDs, inObj2.IDs)
                dataCombined.IDs = inObj1.IDs;
            end
        end
    end

    methods (Access = private)
        function dataCombined = combineCompatibleParameters(objVector)
            if isempty(objVector)
                error("fm_data vector is empty");
            end

            uniqueTR = unique([objVector.TR]);
            assert(length(uniqueTR) <= 1, "All TRs must be the same.");
            assert(has1UniqueValue([objVector.numVoxels]), "Number of voxels must be the same for all elements.");

            dataCombined = fm_data([], uniqueTR, []);

            if has1UniqueString({objVector.rowName})
                dataCombined.rowName = objVector(1).rowName;
            end
            if has1UniqueString({objVector.dataName})
                dataCombined.dataName = objVector(1).dataName;
            end
            
            %%%
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
    end

    methods (Static = true)
        function loadFromFile(obj)
            
        end
    end

    methods (Static = true, Access = private)
        function bool = areAllSame(valueVec, numValues)
            bool = false;
            if isempty(valueVec) || (length(unique(valueVec)) == 1 && length(valueVec) == numValues)
                bool = true;
            end
        end
    end
end
