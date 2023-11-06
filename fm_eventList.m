classdef fm_eventList
%FM_EVENTLIST Tabular representation of a design matrix
% Contains ID, timing and activity information for an fMRI sequence and
% provides methods to quickly build a design matrix
%
% Properties
%   <content> stores a table where each row contains information about an
%   event. It has some or all of these columns:
%     <Trial> Trial number in this sequence
%     <ID> Can be AnalysisID, NeuralPatternID, EventID or any other ID
%     <Activity> Neural intensity for this event
%     <Duration> Event length in seconds
%     <Onset> Event start time in seconds relative to beginning of sequence
%       
%     fm_eventList contains a wrapper so that the above subproperties can
%     be accessed more quickly. Thus, instead of obj.content.ID, you can do
%     obj.ID.
%
%   <runDuration> in seconds
%
% Constructors
%   > obj = fm_eventList()
%   > obj = fm_eventList(content), where 'content' is an existing table with
%     the requisite columns
%   > obj = fm_eventList(content, runDuration)
%
% Methods
%   > newObj = objVec.cat(), where objVec is a vector of fm_eventList
%     objects, concatenates their info together while updating onset data
%     and trial numbers so that they are relative to the beginning of the
%     new, concatenated sequence
%   > h = obj.height(), where h is the number of events in the table
%   > newObj = obj.collapseIntoSuperEvent() will change all ID values to 1
%     and return the new object
%   > newObj = obj.openIntoUniqueEvents() will give each event a unique ID
%     so that there are number-of-events IDs in total
%   > desMat = obj.computeDesignMatrix(TR, numColumns) will transform the
%     eventList into the 2D matrix representation of the design matrix
%     i.e., [number-of-TRs by number-of-IDs]. The 'numColumns' argument is
%     optional
%   > obj.validate() checks that the object's data are internally
%     consistent e.g., the event timings do not exceed the run duration
%   > obj = fm_eventList.preallocate(height, runDuration) is a static
%     method that returns an NaN-filled fm_eventList object with 'height'
%     number of rows. 'runDuration' is optional
% 
% Part of package funkyMage. November 2023.
% https://github.com/aj-abdujabborov/funkyMage


    properties (Dependent = true)
        Trial;
        ID;
        Activity;
        Duration;
        Onset;
        content (:,:) table;
        runDuration (1,1) {mustBeNonnegative};
    end

    properties (Access = protected)
        privateContent;
        privateRunDuration;
    end

    methods
        function obj = fm_eventList(content, duration)
            obj.privateContent = table();
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
                obj.Trial = content.Trial;
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