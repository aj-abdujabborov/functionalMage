classdef fm_data < matlab.mixin.Copyable
%FM_DATA Store any 2D data matrix
%
% Properties
%   <data> The 2D data matrix
%   <TR> Optional repetition time value for time series data
%   <ID> Optional ID values for every row of data
%   <rowName> Optional string to describe the row axis
%   <colName> Optional string to describe the column axis
%   <dataName> Optional string to describe the data
%
% Read-only properties
%   <numRows>
%   <numTRs> returns 'numRows'
%   <numCols>
%   <numVoxels> returns numCols
%   <duration> duration of the data if TR is provided, calculated as
%     numTRs x TR
%
% Constructors
%   > obj = fm_data()
%   > obj = fm_data(data)
%   > obj = fm_data(data, TR)
%   > obj = fm_data(data, TR, ID)
%
% Methods
%   > TR = objVec.getTR(), where objVec is a vector of fm_data objects,
%     will return the TR if all of them have the same TR. Otherwise it will
%     throw an error
%   > numVoxels = objVec.getNumVoxels() will return the number of voxels
%     if all objects have the same number. Otherwise it will throw an
%     error
%   > newObj = objVec.cat() will vertically concatenate all the fm_data
%     objects
%   > newObj = objVec.diagCat() will concatenate all the fm_data objects
%     "diagonally", meaning that
%       newObj.numCols = the sum of numCols of 'objVec' objects
%       newObj.rowCols = the sum of rowCols of 'objVec' objects
%   > newObj = obj1 + obj2 will sum the 'data' properties of both objects
%     and return it in 'newObj'
%
% Part of package funkyMage. November 2023.
% https://github.com/aj-abdujabborov/funkyMage

    properties
        data {mustBeNumeric};
        TR {mustBePositive};
    end

    properties (Dependent = true)
        ID;
    end

    properties
        rowName {mustBeText} = "";
        colName {mustBeText} = "";
        dataName {mustBeText} = "";
    end

    properties (Dependent = true, SetAccess = protected)
        numTRs;
        numRows;
        numCols;
        numVoxels;
        duration;
    end

    properties (Access = protected)
        privateID;
    end

    methods
        function obj = fm_data(data, TR, ID)
            if nargin > 0, obj.data = data; end
            if nargin > 1, obj.TR   = TR;   end
            if nargin > 2, obj.ID  = ID;  end
        end

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
            dataCombined.ID  = cat(1, objVector.ID);
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

            dataCombined.ID  = cat(1, objVector.ID);
        end

        function dataCombined = plus(inObj1, inObj2)
            dataCombined = combineCompatibleParameters([inObj1, inObj2]);
            assert(inObj1.numRows == inObj2.numRows,...
                   "Number of rows / TRs in fm_data objects do not match");
            dataCombined.data = inObj1.data + inObj2.data;
            if isequal(inObj1.ID, inObj2.ID)
                dataCombined.ID = inObj1.ID;
            end
        end
    end

    methods % Get and Set methods
        function set.ID(obj, ID)
            if isempty(ID), return; end
            
            assert(height(ID) == obj.numRows,...
                  "The number of rows in ID should be the same as in data");
            assert(width(ID) == 1,...
                   "ID should be a column vector");
            obj.privateID = ID;
        end

        function set.data(obj, data)
            assert(width(data) >= 0, "Input data must have 0 or more columns.");
            assert(height(data) >= 0, "Input data must have 0 or more rows.");
            obj.data = data;
        end

        function ID = get.ID(obj)
            ID = obj.privateID;
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

    methods (Access = protected)
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
            if has1UniqueString({objVector.colName})
                dataCombined.colName = objVector(1).colName;
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

    methods (Static = true, Access = protected)
        function bool = areAllSame(valueVec, numValues)
            bool = false;
            if isempty(valueVec) || (length(unique(valueVec)) == 1 && length(valueVec) == numValues)
                bool = true;
            end
        end
    end
end
