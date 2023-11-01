classdef fm_eventList
    properties (Dependent = true)
        Trial;
        ID;
        Activity;
        Duration;
        Onset;
        content (:,:) table;
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
            obj.privateContent = table();
        end

        function elCombo = cat(elVec)
            if length(elVec) == 1
                elCombo = elVec;
                return;
            end

            trialColExists = false(1, length(elVec));
            for i = 1:length(elVec)
                trialColExists(i) = any("Trial" == string(elVec(i).content.Properties.VariableNames));
            end
            bCatTrial = all(trialColExists);

            elVec(1).validate();
            for i = 2:length(elVec)
                elVec(i).validate();
                elVec(i).Onset = elVec(i).Onset + sum([elVec(1:i-1).runDuration]);
                if bCatTrial
                    elVec(i).Trial = elVec(i).Trial + elVec(i-1).Trial(end);
                end
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
        function obj = set.Trial(obj, Trial)
            obj.privateContent.Trial(:) = Trial(:);
        end

        function Trial = get.Trial(obj)
            Trial = obj.privateContent.Trial;
        end

        function obj = set.ID(obj, ID)
            obj.privateContent.ID(:) = ID(:);
        end

        function ID = get.ID(obj)
            ID = obj.privateContent.ID;
        end

        function obj = set.Activity(obj, Activity)
            obj.privateContent.Activity(:) = Activity(:);
        end

        function Activity = get.Activity(obj)
            Activity = obj.privateContent.Activity;
        end

        function obj = set.Duration(obj, Duration)
            obj.privateContent.Duration(:) = Duration(:);
        end

        function Duration = get.Duration(obj)
            Duration = obj.privateContent.Duration;
        end

        function obj = set.Onset(obj, Onset)
            obj.privateContent.Onset(:) = Onset(:);
        end

        function Onset = get.Onset(obj)
            Onset = obj.privateContent.Onset;
        end

        function obj = set.content(obj, content)
            tableFields = string(content.Properties.VariableNames);
            obj.privateContent = table();

            if any(tableFields == "Trial")
                obj.Trial = nan(height(content), 1);
            end
            if any(tableFields == "ID")
                obj.ID = content.ID;
            end
            if any(tableFields == "Activity")
                obj.Activity = content.Activity;
            end
            if any(tableFields == "Duration")
                obj.Duration = content.Duration;
            end
            if any(tableFields == "Onset")
                obj.Onset = content.Onset;
            end
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
            T = table('Size', [height, 5],...
                'VariableNames', {'Trial', 'ID', 'Activity', 'Duration', 'Onset'},...
                'VariableTypes', {'double', 'double', 'double', 'double', 'double'});
            T{:,:} = nan;

            if exist('runDuration', 'var')
                eventList = fm_eventList(T, runDuration);
            else
                eventList = fm_eventList(T);
            end
        end
    end
end