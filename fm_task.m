
classdef fm_task < matlab.mixin.Copyable
    properties (Dependent)
        content;
    end

    properties
        expand (1,1) logical = false;
    end

    properties (SetAccess = private)
        contentNumerical;
    end

    properties (Access = private)
        privateContent = fm_task.getEmptyTaskTable();
    end

    properties (Access = private, Constant = true)
        colNames = ["Probability", "Durations", "Onsets",...
                    "NeuralIntensity", "NeuralPatternIDs",...
                    "AnalysisIDs", "ClassificationGroups"]';
        colTypes = ["double", "string", "string",...
                    "string", "string",...
                    "string", "string"]';
        colForms = ["", "floating", "floating", ...
                    "floating", "integerPrefixG", ...
                    "integerPrefixG", "whole"];
        necessaryColumns = [true, true, false,...
                            false, true,...
                            true, true];
        numberColumn = 1;
        stringColumns = [2:7];
        numColumns = 7;
        
        floatingExpr = "^(\d+|\d+\.\d+|\.\d+)$";
        integerExpr = "^(\d+)$";
        integerPrefixGExpr = "^(\d+|[gG]\d+)$";

        numExpansionConds = 2;
    end

    methods
        function obj = fm_task(content, expand)
            if nargin == 0
                return;
            end
            if nargin >= 1
                obj.content = content;
            end
            if nargin == 2
                obj.expand = expand;
            end
        end

        function set.content(obj, content)
            numRows = height(content);

            if numRows == 0
                obj.privateContent = obj.getEmptyTaskTable();
                obj.contentNumerical = table.empty;
                return;
            end

            content = handleInputBeingCell(content);

            checkColumnNamesAndOrder(content);

            content{:,obj.stringColumns} = obj.cleanUpString(content{:,obj.stringColumns});
            obj.checkStringEntriesAreParsable(content);
            content = replacePrefixGWithNegativeSign(content);

            [contentNumer, entryLengths] = obj.convertStringTableToNumerical(content);
            
            checkFieldSizes(entryLengths);

            contentNumer = fillOptionalFields(contentNumer, entryLengths);

            contentNumer = obj.expandTable(contentNumer);
            obj.expand = false;
            
            obj.privateContent = obj.convertNumericalTableToString(contentNumer);
            obj.contentNumerical = [contentNumer, obj.calculateEventIDs(contentNumer)];
            
            function content = handleInputBeingCell(content)
                if iscell(content)
                    assert(width(content) == obj.numColumns, ...
                          "Cell array must have a width of %d", obj.numColumns);
                    content = [obj.getEmptyTaskTable(); content];
                end
            end

            function checkColumnNamesAndOrder(content)
                if ~isequal(string(content.Properties.VariableNames(:)), obj.colNames)
                    error("Inputted table columns do not match what is expected in name or order");
                end
            end
            
            function content = replacePrefixGWithNegativeSign(content)
                integerPrefixGCols = strcmpi(fm_task.colForms, "integerPrefixG");
                content{:,integerPrefixGCols} = strrep(content{:,integerPrefixGCols}, "g", "-");
                content{:,integerPrefixGCols} = strrep(content{:,integerPrefixGCols}, "G", "-");
            end

            function checkFieldSizes(entryLengths)
                for r = 1:numRows
                    emptyFields = entryLengths(r, :) == 0;
                    unfilledColumns = emptyFields & obj.necessaryColumns;
                    if any(unfilledColumns)
                        error("%s on row %d cannot be empty", obj.getColumnNamesAsString(unfilledColumns), r);
                    end

                    fieldsMustBeSameSize = entryLengths(r, :) > 0;
                    fieldsMustBeSameSize(obj.numberColumn) = false;
                    requiredLength = entryLengths(:, find(fieldsMustBeSameSize, 1));
                    if any(requiredLength(r) ~= entryLengths(r, fieldsMustBeSameSize))
                        error("All fields of length more than 0 must be the same length\n" ...
                            + "within each row of the table (except for Probability).\n" ...
                            + "This is not true on row %d", r);
                    end    
                end
            end

            function contentNumer = fillOptionalFields(contentNumer, entryLengths)
                for r = 1:numRows
                    if entryLengths(r, strcmpi(obj.colNames, 'Onsets')) == 0
                        contentNumer.Onsets{r} = [0 cumsum(contentNumer.Durations{r}(1:end-1))];
                    end

                    if entryLengths(r, strcmpi(obj.colNames, 'NeuralIntensity')) == 0
                        contentNumer.NeuralIntensity{r} = ones(1, length(contentNumer.Durations{r}));
                    end
                end
            end
        end

        function content = get.content(obj)
            content = obj.privateContent;
        end

        function contentNumer = expandTable(obj, contentNumer)
            NeuralPatternIDs = contentNumer.NeuralPatternIDs;
            AnalysisIDs      = contentNumer.AnalysisIDs;
            neurPatVals      = cat(2, NeuralPatternIDs{:});
            anlysVals        = cat(2, AnalysisIDs{:});

            neurPatHasNegative = any(neurPatVals < 0);
            anysHasNegative    = any(anlysVals < 0);
            
            if ~anysHasNegative && ~neurPatHasNegative && ~obj.expand
                return;
            end

            [newNeuralPatternIDs, numCombos] = makeNewNeuralPatternIDs();
            newAnalysisIDs                   = makeNewAnalysisIDs();

            newNeuralPatternIDs = cellfun(@(x) abs(x), newNeuralPatternIDs, 'UniformOutput', false);
            newAnalysisIDs      = cellfun(@(x) abs(x), newAnalysisIDs, 'UniformOutput', false);

            contentNumer = repmat(contentNumer, numCombos, 1);
            contentNumer.NeuralPatternIDs = newNeuralPatternIDs;
            contentNumer.AnalysisIDs      = newAnalysisIDs;

            function [newNeuralPatternIDs, numCombos] = makeNewNeuralPatternIDs()
                unqPosNeurPatVals = unique(neurPatVals(neurPatVals >= 0), 'stable');
                numVaryingEvents = length(unqPosNeurPatVals);

                newIDs = obj.findSmallestAvailableValues(...
                    abs(neurPatVals(neurPatVals < 0)), ...
                    fm_task.numExpansionConds);
                
                indepCombos = getIndependentCombinations(numVaryingEvents);
                indepCombos = fm_task.replaceValues( ...
                    indepCombos,...
                    1:fm_task.numExpansionConds, ...
                    newIDs);
                numCombos = height(indepCombos);
                
                newNeuralPatternIDs = cell(length(NeuralPatternIDs)*numCombos, 1);
                cntr = 1;
                for i_cmb = 1:height(indepCombos)
                    for j_entr = 1:length(NeuralPatternIDs)
                        newNeuralPatternIDs{cntr} = fm_task.replaceValues(...
                            NeuralPatternIDs{j_entr}, ...
                            unqPosNeurPatVals, ...
                            indepCombos(i_cmb,:));

                        cntr = cntr + 1;
                    end
                end
            end

            function newAnalysisIDs = makeNewAnalysisIDs()
                unqPosAnlysVals = unique(anlysVals(anlysVals >= 0), 'stable');
                numVaryingEvents = length(unqPosAnlysVals);
                newIDs = obj.findSmallestAvailableValues(abs(unique(anlysVals)), numVaryingEvents*(numCombos-1));
                newIDs = reshape(newIDs, [], numCombos-1)';
                newIDs = [unqPosAnlysVals; newIDs];
                
                newAnalysisIDs = cell(length(AnalysisIDs)*numCombos, 1);
                cntr = 1;
                for i_cmb = 1:height(newIDs)
                    for j_entr = 1:length(AnalysisIDs)
                        newAnalysisIDs{cntr} = fm_task.replaceValues(...
                            AnalysisIDs{j_entr},...
                            unqPosAnlysVals,...
                            newIDs(i_cmb, :));

                        cntr = cntr + 1;
                    end
                end
            end

            function combinations = getIndependentCombinations(numIndependentEvents)
                numTotalCombinations = obj.numExpansionConds ^ numIndependentEvents;

                combinations = cell(1, numIndependentEvents);
                [combinations{:}] = ind2sub(obj.numExpansionConds * ones(1, numIndependentEvents), 1:numTotalCombinations);
                combinations = cat(1, combinations{end:-1:1})';
            end
        end
    end

    methods (Access = private, Static = true)
        function content = getEmptyTaskTable()
            content = table('Size', [0, fm_task.numColumns],... 
	            'VariableNames', fm_task.colNames,...
	            'VariableTypes', fm_task.colTypes);
        end

        function [contentNumerical, lengths] = convertStringTableToNumerical(content)
            contentNumerical = table2cell(content);
            lengths = zeros(size(contentNumerical));
            lengths(:,fm_task.numberColumn) = 1;

            for i = 1:height(contentNumerical)
                for j = fm_task.stringColumns
                    tmp = str2doubleMe(split(contentNumerical{i,j}, ", ")');
                    lengths(i,j) = width(tmp);
                    contentNumerical{i, j} = {tmp};
                        % concatenatable contents are always concatenated by
                        % cell2table, so double enclosure makes sure we always get
                        % a table with cell arrays
                end
            end

            contentNumerical = cell2table(contentNumerical, ...
                'VariableNames', fm_task.colNames);

            function output = str2doubleMe(str)
                if isempty(str{1})
                    output = [];
                else
                    output = str2double(str);
                end
            end
        end

        function content = convertNumericalTableToString(contentNumerical)
            colNames = fm_task.colNames;

            content = table.empty();
            content.Probability = contentNumerical.Probability;
            for c = fm_task.stringColumns
                tmp = cellfun(@(x) join(string(x), ", "),...
                              contentNumerical.(colNames(c)), ...
                              'UniformOutput', false);
                content.(colNames(c)) = cat(1, tmp{:});
            end
        end

        function numbers = findSmallestAvailableValues(A, numNeeded)
            numbers = nan(1, numNeeded);
            
            A = [sort([0 A]) Inf];
            df = diff(A);
            dfGt1 = df > 1;

            numSaved = 0;
            for i = find(dfGt1)
                startAt     = A(i) + 1;
                numToSave   = min(A(i+1) - startAt, numNeeded-numSaved);
                tmp         = startAt : startAt+numToSave-1;
                
                numbers(numSaved+1 : numSaved+numToSave) = tmp;
                numSaved = numSaved + numToSave;

                if numSaved == numNeeded
                    return;
                end
            end
        end


        function EventIDs = calculateEventIDs(contentNumerical)
            numEventsPerTask = cellfun(@length, contentNumerical.AnalysisIDs);
            idRangesPerTask = [0 cumsum(numEventsPerTask(:)')];
            eventIDs = cell(length(numEventsPerTask), 1);
            for i = 1:length(numEventsPerTask)
                eventIDs{i} = idRangesPerTask(i)+1 : idRangesPerTask(i+1);
            end
            EventIDs = table(eventIDs, 'VariableNames', {'EventIDs'});
        end

        function strVec = cleanUpString(strVec)
            arguments
                strVec {mustBeA(strVec, 'string')};
            end

            strVec = strrep(strVec, ",", " ");
            strVec = strip(strVec);
            strVec = regexprep(strVec, " +", ", ");
        end

        function checkStringEntriesAreParsable(content)
            assert(all(content.Probability >= 0), "Probabilities must be non-negative.");
            for c = fm_task.stringColumns
                switch fm_task.colForms(c)
                    case 'floating'
                        expr = fm_task.floatingExpr;
                    case 'integer'
                        expr = fm_task.integerExpr;
                    case 'integerPrefixG'
                        expr = fm_task.integerPrefixGExpr;
                end
                
                if ~fm_task.isConvertibleToRegexp(content{:,c}, expr)
                    error("%s entries are invalid", fm_task.colNames(c));
                end
            end
        end

        function bool = isConvertibleToRegexp(strVec, expression)
            bool = true;
            for i = 1:length(strVec)
                if isempty(strVec{i}), continue; end

                match = regexp(split(strVec(i), ", "), expression);
                    % string is either ddd, ddd.ddd or .ddd

                if ~iscell(match), match = {match}; end
                if any(cellfun(@isempty, match))
                    bool = false;
                    return;
                end
            end
        end

        function strOut = getColumnNamesAsString(idx)
            strOut = strjoin(fm_task.colNames(idx), ", ");
        end

        function replaced = replaceValues(data, old, new)
            [found, idx] = ismember(data, old);
            replaced = data;
            replaced(found) = new(idx(found));
        end

    end
end