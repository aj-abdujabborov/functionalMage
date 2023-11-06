classdef fm_task < matlab.mixin.Copyable
%FM_TASK Contains information about task design
% Builds a table containing crucial information about your experimental
% task, including its timings and core simulation and analysis parameters.
%
% Properties:
%   <content> is a table in which each row specifies a condition of your
%   experimental design. It has the following columns:
%       <Probability> The relative proportion of the experiment the
%       condition should constitute. Should be a whole, nonnegative
%       number.
%       <Durations> A string of numbers, separated by spaces, indicating
%       the duration of each event within that condition. For example, to
%       make trials consisting of a 3s and 6s epoch: "2, 6".
%       <Onsets>: A string of numbers (separated by spaces) indicating the
%       onset time of each event with 0 indicating that the event should
%       begin at the start of the trial. This property can be left as an
%       empty string (i.e., ""), in which case it's assumed that you want
%       the first event to begin at 0s and for the events to be temporally
%       contiguous.
%       <NeuralIntensity> A string of numbers indicating the intensity of
%       neural activity to assign to the events during simulation. It is
%       designed to range between 0 and 1, but any value (including
%       negatives) will work. If left as an empty string (i.e., "") I'll
%       fill them with 1s.
%       <NeuralPatternIDs> A string of positive integers indicating the
%       "ID" of each event. All events of the same ID are assigned the same
%       pattern of neural activity during simulation. This column also
%       accepts the prefix "G" in front of any number when using the
%       expansion feature (see expand property below).
%       <AnalysisIDs> A string of positive integers indicating the "ID" of
%       each event during analysis. In a univariate regression, events of
%       the same ID are combined into the same column of the design matrix.
%       This also supports the expansion feature.
%       <ClassificationGroups> A string of positive integers. Events of the
%       same group (same number) will be classified against each other
%       e.g., if there are 3 events having the value '1', classification
%       will be later done across three conditions within this group.
%
%       Across each row, every column (with the exception of Probability or
%       empty string columns) must have the same number of elements.
%       Note that strings are encapsulated by double quotes. In MATLAB,
%       single quotes create character arrays, not strings.
%
%   <expand> In certain cases, you will want to set up a task table where
%   NeuralPatternIDs of events will vary independently between two
%   possibilities. However, manually writing this out is tedious. If you
%   set 'expand' to true, fm_task will automatically independently vary
%   every event within NeuralPatternIDs that does *not* have the prefix
%   'G'. Meanwhile, AnalysisIDs will follow the same ID structure as what
%   you wrote except make it unique for every new generated set. An example
%   will demonstrate it best:
%       Assume you write a table with one condition:
%       (1) NeuralPatternIDs: "g1 2 3". AnalysisIDs: "g1 2 3"
%       It will unfold into 4 conditions:
%       (1) NeuralPatternIDs: "1 2 2". AnalysisIDs: "1 2 3"
%       (2) NeuralPatternIDs: "1 2 3". AnalysisIDs: "1 4 5"
%       (3) NeuralPatternIDs: "1 3 2". AnalysisIDs: "1 6 7"
%       (4) NeuralPatternIDs: "1 3 3". AnalysisIDs: "1 8 9"
%       Non-g-prefixed values are varying between two states in
%       NeuralPatternIDs.
%       Non-g-prefix values are unique in every condition in AnalysisIDs.
%
%       Expansion can work when you specify multiple conditions too (for
%       example if you have multiple types of trials and you want to
%       independently vary them as a whole set).
%
%       'expand' is false by default, but it will automatically go true if
%       the "G" prefix is used in NeuralPatternIDs or AnalysisIDs. Note
%       that if you want expansion, you need to set 'expand' to true before
%       you assign rows to the 'content' property.
%
%       (If you'd really like to expand into 3 or more states, you can set
%       the 'numExpansionConds' property within fm_task.m itself.)
%
%   <contentNumerical> is spun out by using the 'content' property. You
%   should not need to change this property (though you can), but
%   potentially supply it to other funkyMage objects.
%
% Constructor:
%   task = fm_task(taskTable, expand);
%
% Examples:
%   task = fm_task(); % get an empty fm_task object
%   task.expand = true; % enable expansion if wanted
%   task.content = [task.content;
%       {2, "3 2", "", "1 0.5", "1 2", "1 2", "1 0"};...
%       {1, "3",   "", "1",     "1",   "1",   "1"}];
%       % append each row as a cell array of strings
%   % 'task' is now ready to be used
%
%   % Make table first and initialize fm_task object with it
%   table = [{2, "3 2", "", "1 0.5", "1 2", "1 2", "1 0"};...
%            {1, "3",   "", "1",     "1",   "1",   "1"}];
%   task = fm_task(table, true); % second input is 'expand' and optional
%   task.content % 'task' is now ready to be used
%
% Part of package funkyMage. November 2023.
% https://github.com/aj-abdujabborov/funkyMage

    properties (Dependent)
        content;
    end

    properties
        expand (1,1) logical = false;
        contentNumerical;
    end

    properties (Access = protected)
        privateContent = fm_task.getEmptyTaskTable();
    end

    properties (Access = protected, Constant = true)
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
            if nargin == 2
                obj.expand = expand; % expand has to be set first
            end
            obj.content = content;
        end

        function set.content(obj, content)
            numRowsOfContentInitial = height(content);

            if numRowsOfContentInitial == 0
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
            obj.contentNumerical = [...
                contentNumer,...
                obj.calculateEventIDs(contentNumer),...
                computeTotalDurationPerCondition(contentNumer)];
            
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
                for r = 1:numRowsOfContentInitial
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
                for r = 1:numRowsOfContentInitial
                    if entryLengths(r, strcmpi(obj.colNames, 'Onsets')) == 0
                        contentNumer.Onsets{r} = [0 cumsum(contentNumer.Durations{r}(1:end-1))];
                    end

                    if entryLengths(r, strcmpi(obj.colNames, 'NeuralIntensity')) == 0
                        contentNumer.NeuralIntensity{r} = ones(1, length(contentNumer.Durations{r}));
                    end
                end
            end

            function totalDuration = computeTotalDurationPerCondition(contentNumer)
                for r = height(contentNumer):-1:1
                    totalDuration(r) = max(contentNumer.Durations{r} + contentNumer.Onsets{r});
                end
                totalDuration = table(totalDuration(:), 'VariableNames', {'TotalDuration'});
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