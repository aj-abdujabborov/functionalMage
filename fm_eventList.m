classdef fm_eventList
    properties (Dependent = true)
        ID;
        Activity;
        Duration;
        Onset;
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
            end
            if nargin > 1
                obj.runDuration = duration;
            end
        end

        function elCombo = cat(elVec)
            if length(elVec) == 1
                elCombo = elVec;
                return;
            end

            elVec(1).validate();
            for i = 2:length(elVec)
                elVec(i).validate();
                elVec(i).Onset = elVec(i).Onset + sum([elVec(1:i-1).runDuration]);
            end
            elCombo = fm_eventList(cat(1, elVec.content), sum([elVec.runDuration]));
        end

        function elCombo = vertcat(elVec)
            elCombo = cat(elVec);
        end

        function h = height(obj)
            h = height(obj.ID);
        end

        function obj = collapseIntoSuperEvent(obj)
            obj.ID = 1;
        end

        function obj = openIntoUniqueEvents(obj)
            obj.ID = 1:height(obj.ID);
        end

        function desMat = computeDesignMatrix(obj, TR, nColumns)
            arguments
                obj;
                TR (1,1) {mustBePositive};
                nColumns (1,1) {mustBePositive} = max(obj.ID);
            end

            obj.validate();

            allTimings = [obj.Onset; obj.Duration];
            if ~isRound(allTimings ./ TR)
                error('Timings in eventList should be divisible by the TR');
            end

            desMatResolution = 1/TR;
            desMat = zeros(obj.runDuration * desMatResolution, nColumns);
            
            for i = 1:height(obj)
                if obj.ID(i) < 1
                    continue;
                end
            
                fromRow = obj.Onset(i) * desMatResolution + 1;
                numRows = obj.Duration(i) * desMatResolution;
                toRow = fromRow + numRows - 1;
            
                activityColumn = obj.Activity(i) * ones(numRows, 1);
            
                desMat(fromRow:toRow, obj.ID(i)) = desMat(fromRow:toRow, obj.ID(i)) + activityColumn;
            end
        end

        function validate(obj)
            assert(obj.runDuration >= max(obj.Onset + obj.Duration),...
                   "Event timings exceed specified run duration.");
            mustBePositive(obj.ID);
            mustBeInteger(obj.ID);
            mustBeNonnegative(obj.Duration);
            mustBeNonnegative(obj.Onset);
        end
    end

    methods % Set and Get methods
        function obj = set.ID(obj, ID)
            obj.privateContent.ID(:) = ID(:);
        end

        function ID = get.ID(obj)
            ID = obj.content.ID;
        end

        function obj = set.Activity(obj, Activity)
            obj.privateContent.Activity(:) = Activity(:);
        end

        function Activity = get.Activity(obj)
            Activity = obj.content.Activity;
        end

        function obj = set.Duration(obj, Duration)
            obj.privateContent.Duration(:) = Duration(:);
        end

        function Duration = get.Duration(obj)
            Duration = obj.content.Duration;
        end

        function obj = set.Onset(obj, Onset)
            obj.privateContent.Onset(:) = Onset(:);
        end

        function Onset = get.Onset(obj)
            Onset = obj.content.Onset;
        end

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

        function obj = set.runDuration(obj, runDuration)
            obj.privateRunDuration = runDuration;
        end

        function runDuration = get.runDuration(obj)
            runDuration = obj.privateRunDuration;
        end
    end

    methods (Static = true)
        function eventList = preallocate(height, runDuration)
            T = table('Size', [height, 4],...
                'VariableNames', {'ID', 'Activity', 'Duration', 'Onset'},...
                'VariableTypes', {'double', 'double', 'double', 'double'});
            T{:,:} = nan;

            if exist('runDuration', 'var')
                eventList = fm_eventList(T, runDuration);
            else
                eventList = fm_eventList(T);
            end
        end
    end
end