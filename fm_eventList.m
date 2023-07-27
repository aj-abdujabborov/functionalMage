classdef fm_eventList
    properties (Dependent = true)
        ID;
        Activity;
        Duration;
        Onset;
    end

    properties
        content;
        runDuration;
    end    

    methods % Set and Get methods
        function obj = set.content(obj, content)
            assert(istable(content), "Content needs to be a table");
            tableFields = summary(content);
            if ~all(isfield(tableFields, {'ID', 'Activity', 'Duration', 'Onset'}))
                error("A necessary field is missing from the table");
            end
            obj.content = content;
        end

        function obj = set.ID(obj, ID)
            assert(height(obj.content) == height(ID),...
                   "Number of rows does not match.");
            obj.content.ID = ID;
        end

        function obj = set.Activity(obj, Activity)
            assert(height(obj.content) == height(Activity),...
                   "Number of rows does not match.");
            obj.content.Activity = Activity;
        end

        function obj = set.Duration(obj, Duration)
            assert(height(obj.content) == height(Duration),...
                   "Number of rows does not match.");
            obj.content.Duration = Duration;
        end

        function obj = set.Onset(obj, Onset)
            assert(height(obj.content) == height(Onset),...
                   "Number of rows does not match.");
            obj.content.Onset = Onset;
        end

        function ID = get.ID(obj)
            ID = obj.content.ID;
        end

        function Activity = get.Activity(obj)
            Activity = obj.content.Activity;
        end

        function Duration = get.Duration(obj)
            Duration = obj.content.Duration;
        end

        function Onset = get.Onset(obj)
            Onset = obj.content.Onset;
        end
    end

    methods
        function obj = fm_eventList(content, duration)
            if nargin > 0
                obj.content = content;
                obj.runDuration = duration;
            end
        end

        function elCombo = cat(elVec)
            if length(elVec) == 1
                elCombo = elVec;
                return;
            end

            for i = 2:length(elVec)
                elVec(i).Onset = elVec(i).Onset + sum([elVec(1:i-1).runDuration]);
            end
            elCombo = fm_eventList(cat(1, elVec.content), sum([elVec.runDuration]));
        end

        function elCombo = vertcat(elVec)
            elCombo = cat(elVec);
        end

        function h = height(eventList)
            h = height(eventList.ID);
        end

        function eventList = collapseIntoSuperEvent(eventList)
            eventList.ID = 1;
        end

        function eventList = openIntoUniqueEvents(eventList)
            eventList.ID = 1:height(eventList.ID);
        end

        function desMat = computeDesignMatrix(eventList, TR, nColumns)
            if ~(exist('nColumns', 'var') && ~isempty(nColumns))
                nColumns = max(eventList.ID);
            end

            allTimings = [eventList.Onset; eventList.Duration];
            if ~isRound(allTimings ./ TR)
                error('Timings in eventList should be divisible by the TR');
            end

            desMatResolution = 1/TR;
            desMat = zeros(eventList.runDuration * desMatResolution, nColumns);
            
            for i = 1:height(eventList)
                if eventList.ID(i) < 1
                    continue;
                end
            
                fromRow = eventList.Onset(i) * desMatResolution + 1;
                numRows = eventList.Duration(i) * desMatResolution;
                toRow = fromRow + numRows - 1;
            
                activityColumn = eventList.Activity(i) * ones(numRows, 1);
            
                desMat(fromRow:toRow, eventList.ID(i)) = desMat(fromRow:toRow, eventList.ID(i)) + activityColumn;
            end
        end
    end
end