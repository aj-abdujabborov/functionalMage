classdef fm_eventList
    properties (Dependent = true)
        ID {mustBePositive, mustBeInteger};
        Activity {mustBeFinite};
        Duration {mustBeNonnegative};
        Onset {mustBeNonnegative};
        content (:,4) table;
        runDuration (1,1) {mustBeNonnegative};
    end

    properties (Access = private)
        privateContent;
        privateRunDuration;
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
            arguments
                eventList;
                TR (1,1) {mustBePositive};
                nColumns (1,1) {mustBePositive} = max(eventList.ID);
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

    methods % Set and Get methods
        %%%
        function obj = set.ID(obj, ID)
            obj.privateContent.ID(:) = ID(:);
        end

        function ID = get.ID(obj)
            ID = obj.content.ID;
        end

        %%%
        function obj = set.Activity(obj, Activity)
            obj.privateContent.Activity(:) = Activity(:);
        end

        function Activity = get.Activity(obj)
            Activity = obj.content.Activity;
        end

        %%%
        function obj = set.Duration(obj, Duration)
            obj.privateContent.Duration(:) = Duration(:);
        end

        function Duration = get.Duration(obj)
            Duration = obj.content.Duration;
        end

        %%%
        function obj = set.Onset(obj, Onset)
            obj.privateContent.Onset(:) = Onset(:);
        end

        function Onset = get.Onset(obj)
            Onset = obj.content.Onset;
        end

        %%%
        function obj = set.content(obj, content)
            tableFields = summary(content);
            if ~all(isfield(tableFields, {'ID', 'Activity', 'Duration', 'Onset'}))
                error("A necessary field is missing from the table");
            end
            obj.privateContent = table();
            obj.ID = content.ID;
            obj.Activity = content.Activity;
            obj.Duration = content.Duration;
            obj.Onset = content.Onset;
        end

        function content = get.content(obj)
            content = obj.privateContent;
        end

        %%%
        function obj = set.runDuration(obj, runDuration)
            assert(runDuration >= max(obj.Onset + obj.Duration),...
                   "Event timings exceed specified run duration.");
            obj.privateRunDuration = runDuration;
        end

        function runDuration = get.runDuration(obj)
            runDuration = obj.privateRunDuration;
        end
    end
end