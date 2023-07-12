function [designMatrix] = computeDesignMatrix(eventList, runDuration, TR, varargin)
%COMPUTEDESIGNMATRIX Compute a design matrix from an eventList
% designMatrix = COMPUTEDESIGNMATRIX(eventList, runDuration, TR) creates a
% separate column for each unique event ID in eventList and assigns it to
% column index event ID in designMatrix.
% 
% Parameters
% 'nColumns' sets the number of columns in the design matrix. Default is
% the maximum ID in eventList.

%% Run
% Check inputs
p = inputParser;
addRequired(p, 'eventList', @(x) istable(x) & all(x.ID > 0) & isequal(round(x.ID), x.ID));
addRequired(p, 'runDuration', @(x) numel(x) == 1);
addRequired(p, 'TR', @(x) numel(x) == 1);
addParameter(p, 'nColumns', [], @(x) numel(x) == 1 || isempty(x));
parse(p, eventList, runDuration, TR, varargin{:});

nColumns = p.Results.nColumns;
if isempty(nColumns)
    nColumns = max(eventList.ID);
end
if max(eventList.Onset + eventList.Duration) > runDuration
    error('Timings in eventList exceed provided run duration');
end
allTimings = [eventList.Onset; eventList.Duration];
if ~isequal(allTimings ./ TR, round(allTimings ./ TR))
    error('Timings in eventList should be divisible by the temporal resolution');
end

% Compute
designMatrixResolution = 1/TR;
designMatrix = zeros(runDuration * designMatrixResolution, nColumns);

for i = 1:height(eventList)
    if eventList.ID(i) < 1
        continue;
    end

    fromRow = eventList.Onset(i) * designMatrixResolution + 1;
    numRows = eventList.Duration(i) * designMatrixResolution;
    toRow = fromRow + numRows - 1;

    activityColumn = eventList.Activity(i) * ones(numRows, 1);

    designMatrix(fromRow:toRow, eventList.ID(i)) = designMatrix(fromRow:toRow, eventList.ID(i)) + activityColumn;
end
    
