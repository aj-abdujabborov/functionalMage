
classdef fm_taskTable < matlab.mixin.Copyable
    properties (Dependent)
        content;
    end

    properties
        contentNumerical;
        NeuralIntensity;
        NeuralPatternIDs;
        RegressionIDs;
        ClassificationGroups;
        EventIDs;

        numNeuralPatterns;
    end

    properties (Access = private)
        privateContent;
    end

    properties (Access = private, Constant = true)
        columnNamesTypes = [
                        ["Probability", "double"]; ...
			            ["Durations", "string"]; ...
			            ["Onsets", "string"];
                        ["NeuralIntensity", "string"]; ...
                        ["NeuralPatternIDs", "string"]; ...
			            ["RegressionIDs", "string"]; ...
                        ["ClassificationGroups", "string"]; ...
                        ];
        stringColumns = [2:7];
    end

    properties (Access = public, Constant = true, Hidden = true)
        NON_CLASSIFIED_EVENT = 0; % 0 in ClassificationGroups = do not classify
    end


    methods
        %%%
        function obj = fm_taskTable()
            obj.content = [];
        end
    end

    methods % Get and Set methods
        %%%
        function set.content(obj, content)
            if isempty(content)
                obj.privateContent = obj.getEmptyTaskTable();
                return;
                % PROBLEM: what happens if content is empty but other
                % fields are already set? they have to be cleared too.
            end

            %{
            assert(all(content.Probability > 0), "Probabilities must be larger than 0.");
            content.Durations

            for i = 1:size(content, 1)
                keyboard;
                tokens = split(strrep(content.Durations(i), ",", " "));
                    % do this and save it back into content.
                x = regexp(tokens, "^(\d+|\d+\.\d+|\.\d+)$")
                if all(struct2array(x)==1), inputs are correct; end
                        % checks that the number is either xxx, xxx.xxx or .xxx

                % TODO:
                % * Certain entries must be non-empty
                % * Number of elements in all entries of the same row (except
                % Probability) must match, except for empty fields that are
                % allowed to be empty.
                % * Onsets and NeuralIntensity may be empty, in which case the
                % numbers must be filled in.
                % * Trailing spaces must be removed. Trailing commas should
                % throw an error.
                % * Entrie with any non-space & non-comma & non-number
                % characters must throw an error (except for
                % NeuralPatternIDs, which should allow letters (can't have
                % multiple letters in sequence without a break).
            end
            %}

            obj.privateContent = content;
            obj.contentNumerical = obj.convertContentToNumerical(obj.privateContent);
            obj.contentNumerical = [obj.contentNumerical, obj.calculateEventIDs(obj.contentNumerical)];
            

            obj.NeuralIntensity = obj.collapseCellArray(obj.contentNumerical.NeuralIntensity);
            obj.NeuralPatternIDs = obj.collapseCellArray(obj.contentNumerical.NeuralPatternIDs);
            obj.RegressionIDs = obj.collapseCellArray(obj.contentNumerical.RegressionIDs);
            obj.ClassificationGroups = obj.collapseCellArray(obj.contentNumerical.ClassificationGroups);
            obj.EventIDs = obj.collapseCellArray(obj.contentNumerical.EventIDs);

            obj.numNeuralPatterns = obj.numUniqueElements(obj.NeuralPatternIDs);
        end

        %%%
        function content = get.content(obj)
            content = obj.privateContent;
        end
    end

    methods (Access = private, Static = true)
        %%%
        function content = getEmptyTaskTable()
            content = table('Size', [0, size(fm_taskTable.columnNamesTypes,1)],... 
	            'VariableNames', fm_taskTable.columnNamesTypes(:,1),...
	            'VariableTypes', fm_taskTable.columnNamesTypes(:,2));
        end

        %%%
        function contentNumerical = convertContentToNumerical(content)
            contentNumerical = table2cell(content);

            contentNumerical(:,fm_taskTable.stringColumns) = cellfun(...
                @(x) str2double(split(strrep(x, ",", " "))'), ...
                contentNumerical(:,fm_taskTable.stringColumns),...
                'UniformOutput', 0);

            contentNumerical(:,fm_taskTable.stringColumns) = cellfun(...
                @(x) {x},...
                contentNumerical(:,fm_taskTable.stringColumns),...
                'UniformOutput', false);
                % concatenatable contents are always concatenated by
                % cell2table, so double enclosure makes sure we always get
                % a table with cell arrays

            contentNumerical = cell2table(...
                contentNumerical,...
                'VariableNames', fm_taskTable.columnNamesTypes(:,1));
        end

        %%%
        function EventIDs = calculateEventIDs(contentNumerical)
            numEventsPerTask = cellfun(@length, contentNumerical.RegressionIDs);
            idRangesPerTask = [0 cumsum(numEventsPerTask(:)')];
            eventIDs = cell(length(numEventsPerTask), 1);
            for i = 1:length(numEventsPerTask)
                eventIDs{i} = idRangesPerTask(i)+1 : idRangesPerTask(i+1);
            end
            EventIDs = table(eventIDs, 'VariableNames', {'EventIDs'});
        end

        %%%
        function numElements = numUniqueElements(vector)
            numElements = length(unique(vector));
        end

        %%%
        function collapsed = collapseCellArray(cellArray)
            collapsed = [cellArray{:}];
        end
    end
end