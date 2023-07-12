
classdef fm_taskTable < matlab.mixin.Copyable
    properties
        content;
        contentNumerical;

        NeuralActivity;
        NeuralPatternIDs;
        RegressionIDs;
        ClassificationGroups;
        EventIDs;

        numNeuralPatterns;
    end


    properties (Access = private, Constant = true)
        columnNamesTypes = [
                        ["Probability", "double"]; ...
			            ["Durations", "string"]; ...
			            ["Onsets", "string"];
                        ["NeuralActivity", "string"]; ...
                        ["NeuralPatternIDs", "string"]; ...
			            ["RegressionIDs", "string"]; ...
                        ["ClassificationGroups", "string"]; ...
                        ];
        stringColumns = [2:4, 6:7];
        letterColumn = 5;
    end


    methods
        function obj = fm_taskTable()
            obj.content = obj.getEmptyTaskTable();
        end

        function set.content(obj, newContent)
            if isempty(newContent)
                obj.content = newContent;
                return;
            end

            for i = 1:size(newContent, 1)
                % TODO:
                % * Certain entries must be non-empty
                % * Number of elements in all entries of the same row (except
                % Probability) must match, except for empty fields that are
                % allowed to be empty.
                % * Onsets and NeuralActivity may be empty, in which case the
                % numbers must be filled in.
                % * Trailing spaces must be removed. Trailing commas should
                % throw an error.
                % * Entrie with any non-space & non-comma & non-number
                % characters must throw an error (except for
                % NeuralPatternIDs, which should allow letters (can't have
                % multiple letters in sequence without a break).
            end

            obj.content = newContent;
            obj.setContentNumerical();
            obj.setEventIDs();
            obj.setIDVectors();
            obj.setNumNeuralPatterns();
        end
    end


    methods (Access = private)
        function setContentNumerical(obj)
            obj.contentNumerical = obj.convertToNumericalTable(obj.content);
            obj.contentNumerical.NeuralPatternIDs = transformLettersToNumbers(obj.contentNumerical.NeuralPatternIDs);

            function cellArrayNumbers = transformLettersToNumbers(cellArrayLetters)
                letters = unique([cellArrayLetters{:}]);
                numbers = 1:length(letters);
                for i = length(cellArrayLetters):-1:1
                    [~, idx] = ismember(cellArrayLetters{i}, letters);
                    cellArrayNumbers{i} = numbers(idx);
                end
                cellArrayNumbers = cellArrayNumbers';
            end
        end

        function setEventIDs(obj)
            numEventsPerTask = cellfun(@length, obj.contentNumerical.RegressionIDs);
            idRangesPerTask = [0 cumsum(numEventsPerTask(:)')];
            eventIDs = cell(length(numEventsPerTask), 1);
            for i = 1:length(numEventsPerTask)
                eventIDs{i} = idRangesPerTask(i)+1 : idRangesPerTask(i+1);
            end
            obj.contentNumerical = [obj.contentNumerical...
                                    table(eventIDs, 'VariableNames', {'EventIDs'})];
        end

        function setIDVectors(obj)
            obj.NeuralActivity = collapse(obj.contentNumerical.NeuralActivity);
            obj.NeuralPatternIDs = collapse(obj.contentNumerical.NeuralPatternIDs);
            obj.RegressionIDs = collapse(obj.contentNumerical.RegressionIDs);
            obj.ClassificationGroups = collapse(obj.contentNumerical.ClassificationGroups);
            obj.EventIDs = collapse(obj.contentNumerical.EventIDs);

            function output = collapse(cellArray)
                output = [cellArray{:}];
            end
        end

        function setNumNeuralPatterns(obj)
            obj.numNeuralPatterns = length(unique(obj.NeuralPatternIDs));
        end
    end


    methods (Access = protected, Static = true)
        function content = getEmptyTaskTable()
            content = table('Size', [0, size(fm_taskTable.columnNamesTypes,1)],... 
	            'VariableNames', fm_taskTable.columnNamesTypes(:,1),...
	            'VariableTypes', fm_taskTable.columnNamesTypes(:,2));
        end

        function output = convertToNumericalTable(content)
            output = table2cell(content);
            output(:,fm_taskTable.stringColumns) = cellfun(@(x) str2double(split(x))', output(:,fm_taskTable.stringColumns), 'UniformOutput', 0);
            output(:,fm_taskTable.letterColumn) = cellfun(@(x) split(x)', output(:,fm_taskTable.letterColumn), 'UniformOutput', 0);
            output = cell2table(output, 'VariableNames', fm_taskTable.columnNamesTypes(:,1));
        end
    end
end