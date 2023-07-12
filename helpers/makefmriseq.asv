function [runOutlinesAxRuns, trialSeqAxRuns, fullSetPercent] = makefmriseq(taskDurs, epochIDs, condProp, runDur, nRuns, itiModel, itiParams, res, varargin)
%MAKEFMRISEQ Generate a random task sequence for fMRI study
% [runOutlines, trialSequences] = MAKEFMRISEQ(taskDurs, epochIDs, condProp,
% runDur, nRuns, itiModel, itiParams, res, varargin).
%
% <runOutline> is a [1 nRuns] cell array in which every cell contains
%   information about the sequence. The first column is the event ID.
%   Second is neural activity (useful mostly for simulations). Third is
%   event duration (seconds). Fourth is event onset (s).
% <trialSequences> is a [1 nRuns] cell array in which every cell is a [1
%   nTrials] matrix that indicates the condition ID of every trial.
%
% <taskDurs> is a [1 nTrialTypes] cell array where each cell contains a
%   vector specifying the durations of epochs e.g., an experiment where one
%   trial type has two contiguous epochs of 2 and 4 s and another trial has
%   one epoch of 3 s would be: {[2 4], [3]}.
% <epochIDs> is a cell array of the same dimensions as <taskDurs>. It
%   specifies the identity of each epoch (to which condition each epoch
%   belongs, since epochs from different trial types can be in the same
%   condition e.g., if in the above task the 2 s and 3 s event were the
%   same event, epochIDs might be {[1 2], [1]}). These values will carry to
%   the first column of <runOutlines>.
% <condProp> is a [1 nTrialTypes] vector that specifies the proportion of
%   each trial type e.g, [2 1] if we wanted the first trial type from above
%   occuring 2/3 of the time. (See <nRuns> for important note.)
% <runDur> is the desired run duration in seconds.
% <nRuns> is the number of desired runs. The program finds the minimum set
%   of trials per condition needed to achieve the desired condition
%   proportions, then sequentially fits this set as many times as possible,
%   disregarding run boundaries. For example, if nRuns = 10, and 5 runs are
%   required to achieve a full set, then the first 5 and the last 5 runs
%   will both achieve the requested proportions. Conversely, if 0.4 runs
%   are required to create a set, then the first run will contain two sets
%   as well as some trials from the next set (but these will be shuffled).
%   Thus, if you want to ensure that each run contains the specified and
%   same condition proportions, call makefmriseq with nRuns = 1 separately
%   for each run. nRuns can be a range ([minRuns maxRuns]. If so, the runs
%   are fit with as many trials as possible and empty runs are discarded.
% <itiModel> is a character array that takes one of three values:
%   'fixed', 'uniform', and 'exponential', specifying the model of the
%   inter-trial intervals.
% <itiParams> are the parameters of the ITI. If itiModel is 'fixed', this
%   is just one number. If 'uniform' or 'exponential', this is a vector
%   [minimum ITI, maximum ITI].
% <res> is the resolution of the design matrix. Durations in all your
%   inputs should be divisible by res. The TR should be divisible by and no
%   smaller than res.
%
% Parameters:
% Set 'addExtraTrials' to 1 to fill the time that remains after achieving
%   the requested condition proportions with random trials.
% Set 'lambda' to a particular positive value when you want to modulate the
%   "extremeness" of exponentially distributed ITIs (lower lambda means more
%   short-duration ITIs). Default = 3.
%
% Version 3
% 
% Written by AJ (aa5313@nyu.edu).
%
% Some parts of this function were based on the Python package Neurodesign
% (https://github.com/neuropower/neurodesign).

%% Parse
itiModelNames = {'fixed', 'uniform', 'exponential'};
itiParamsPerModel = [1 2 2];

p = inputParser;
addRequired(p, 'taskDurs', @iscell);
addRequired(p, 'epochIDs', @iscell);
addRequired(p, 'condProp', @ismatrix);
addRequired(p, 'runDur', @(x) isvector(x) & length(x) == 1);
addRequired(p, 'nRuns', @(x) isvector(x) && any(size(x, 2) == [1 2]));
if length(nRuns) == 2
    assert(nRuns(2) >= nRuns(1), 'Second value in nRuns range must be equal to or greater than the first.');
end
addRequired(p, 'itiModel', @(x) any(strcmpi(x, itiModelNames)));
addRequired(p, 'itiParams', @(x) isvector(x) && any(length(x) == [1 2 3]) && all(x >= 0));
addRequired(p, 'res', @(x) isvector(x) & length(x) == 1);
addParameter(p, 'lambda', 3, @(x) isvector(x) & length(x) == 1);
addParameter(p, 'addExtraTrials', 1, @(x) any(x == [1 0]));
parse(p, taskDurs, epochIDs, condProp, runDur, nRuns, itiModel, itiParams, res, varargin{:});

% retrieve optional parameters
lambda = p.Results.lambda;
addExtraTrials = p.Results.addExtraTrials;

% additional checks
assert(all(length(taskDurs) == [length(epochIDs), length(condProp)]),...
    'taskDurs and epochIDs (cell array) and condProp (vector) all need to have the same length.');
itiModelIdx = strcmpi(itiModel, itiModelNames);
assert(itiParamsPerModel(itiModelIdx) == length(itiParams),...
    'Make sure that itiModel and the number of parameters itiParams match.');
    % check works b/c 'fixed', 'uniform', and 'exponential' have 1, 2, 2 parameters respectively

%% Variables
condIDs = 1:length(taskDurs);
epochsPCond = cell2mat(cellfun(@length, taskDurs, 'UniformOutput', 0));

% Normalize condProp so that sum is 1
condProp = condProp ./ sum(condProp);
condDurs = cell2mat(cellfun(@sum, taskDurs, 'UniformOutput', 0));

%% Generate sequence
nConds = length(taskDurs);

% Find itiMean
switch lower(itiModel) % for uniform find itiMean, otherwise it should be given
    case 'fixed'
        itiMean = itiParams;
        if itiMean/res ~= floor(itiMean/res)
            error('The fixed ITI amount should be a multiple of the resolution.');
        end
    case 'uniform'
        itiMin = itiParams(1);
        itiMax = itiParams(2);
        if itiMin/res ~= floor(itiMin/res) || itiMax/res ~= floor(itiMax/res)
            error('The ITI maximum and minimums should be a multiple of the resolution.');
        end
        itiMean = mean(itiParams);
    case 'exponential'
        itiMin = itiParams(1);
        itiMax = itiParams(2);
        itiRange = itiMax - itiMin;
        
        itiMinResAdj = itiMin - res/2;
        itiRangeResAdj = itiRange + res;
        % When we round the ITIs into the temporal resolution (TR),
        % when we round values around X to X, there need to be values
        % down to X-TR/2 and up to X+TR/2, but this is false with
        % itiMin and itiMax, so we'll extend itiRange before rounding.
        
        if itiMin/res ~= floor(itiMin/res) || itiMax/res ~= floor(itiMax/res)
            error('The ITI maximum and minimums should be a multiple of the resolution.');
        end
            
        % True distribution of ITIs:
        itiMeanOfDist = lambda - itiRangeResAdj*(exp(itiRangeResAdj/lambda)-1)^(-1) + itiMinResAdj;
        itiMean = itiMeanOfDist;
end

% Determine # of trials
taskMean = condDurs * condProp'; % (roughly) average task duration
trialMean = itiMean + taskMean;

%% Compute number of trials per condition
% Assign particular integer values to frequencies of different trials
nRunsMax = max(nRuns);

% Trial and duration parameters
totalDur = runDur * nRunsMax; % totalDur can be a range [T1 T2] or just [T]
nTrials = floor(totalDur / trialMean);

% Find factor that makes condProp be a vector of whole numbers 
% (this gives the minimum number of trials per condition)
f = 1 / min(condProp);
F = 0;
while 1
    F = F + 1;
    condPropSet = (condProp*f) * F;
    cpsRec = round(condPropSet*10e10)/10e10; % round to 10th decimal point to eliminate tiny computer error
    if all(floor(cpsRec) == cpsRec)
        break;
    end
end
condPropSet = cpsRec;

% Find total number of trials per set
nTrialsPSet = sum(condPropSet);

%% Randomize trial sequence
% Abbreviations
sh.condID = 1; sh.condDur = 2; sh.cumDur = 3;

nSetsInitial = ceil(nTrials / nTrialsPSet) + 1;
trialSeqWDur = nan(nSetsInitial*nTrialsPSet, 3);

for s = 1:nSetsInitial
    % Shuffle sequence within the set
    trialSeq4Set = nan(1, nTrialsPSet);
    bIdx = 1; eIdx = condPropSet(1);
    for c = 1:nConds
        trialSeq4Set(bIdx:eIdx) = condIDs(c) .* ones(1, condPropSet(c));
        if c == nConds, break; end % so that c = nConds+1 isn't computed (outside array)
        bIdx = eIdx + 1;
        eIdx = eIdx + condPropSet(c+1);
    end
    trialSeq4Set = trialSeq4Set(randperm(nTrialsPSet));
    
    % Allocate to trialSeqWDur
    trialSeqWDur((s-1)*nTrialsPSet+1 : s*nTrialsPSet, sh.condID) = trialSeq4Set;
end

% Obtain duration of each trial type
trialSeqWDur(:,sh.condDur) = condDurs(trialSeqWDur(:,sh.condID)) + itiMean;
trialSeqWDur(:,sh.cumDur) = cumsum(trialSeqWDur(:,sh.condDur)); % cumulative durations

% Find indices in trialSeqWDur where to break the sequence into runs
lastIdxOfRun = nan(1, nRunsMax); % crucial that lastIdxOfRun's length does not exceed nRunsMax
for n = 1:nRunsMax
    if n == 1
        temp = trialSeqWDur(:,sh.cumDur) <= runDur;
    else
        temp = trialSeqWDur(:,sh.cumDur) <= runDur + trialSeqWDur(lastIdxOfRun(n-1), sh.cumDur);
    end
    lastIdxOfRun(n) = find(temp, 1, 'last');
end

% Recompute how many sets we can fit in based on our segmentation plan
% and how many runs that equals
nSetsAsFRun = floor(lastIdxOfRun ./ nTrialsPSet);
[nSets, nRuns4MaxSets] = max(nSetsAsFRun);
    % find the number of highest number of runs with whole sets
if length(nRuns) == 2
    if nRuns4MaxSets < nRuns(1)
        nRuns = nRuns(1);
    else
        nRuns = nRuns4MaxSets;
    end
end

% If addExtraTrials is OFF, eliminate extra trials
if ~addExtraTrials
    lastIdxOfRun(nRuns4MaxSets:end) = nSets * nTrialsPSet;
        % make what remains a whole set, to be deleted later anyways
end

% Produce errors or warnings if warranted
if nSets == 0
    error('Trials with requested probabilities of conditions cannot be fit even once. Try increasing number of runs, run duration, or changing condition ratios.')
end
fullSetsDur = (nSets * nTrialsPSet * trialMean);
fullSetPercent = fullSetsDur/(nRuns*runDur)*100;
if fullSetPercent < 85 && ~addExtraTrials
    warning('Only %2.0f%% of the available time contains desired condition proportions. Try adjusting condition probabilities, changing run duration, total duration, or trial timings.', fullSetPercent);
end
allDurUsed = trialSeqWDur(lastIdxOfRun(end), sh.cumDur);
allDurPercent = allDurUsed/(nRuns*runDur)*100;
if allDurPercent < 90 && addExtraTrials
    warning('Only %2.0f%% of the available total time is used. Try changing run duration, total duration or trial timings.', allDurPercent);
end

trialSeqAxRuns = cell(1, nRuns); % ax = across
% Segment trial sequence into runs
for n = 1:nRuns
    if n == 1
        idxRange = 1:lastIdxOfRun(n);
    else
        idxRange = lastIdxOfRun(n-1)+1:lastIdxOfRun(n);
    end
    trialSeq = trialSeqWDur(idxRange, sh.condID);
    trialSeqAxRuns{n} = trialSeq(randperm(length(trialSeq)));
        % Shuffle sequence within the run. This is useful to eliminate
        % predictability when more than one set can fit into a run
end

%% ITI order
itiSeqAxRuns = cell(1, nRuns);
for n = 1:nRuns
    nTrials = length(trialSeqAxRuns{n});
    
    if strcmpi(itiModel, 'fixed')
        itiSeq = itiMean * ones(1, nTrials);
    else
        % Make an itiSequence and then adjust it to match the nominal mean
        i = 0;
        while 1
            if strcmpi(itiModel, 'uniform')
                itiSeq = rand(1, nTrials) * (itiMax-itiMin+res) + (itiMin-res/2);
                itiSeq = round(itiSeq / res) * res;
            elseif strcmpi(itiModel, 'exponential')            
                % Generate distribution
                randvals = rand(1,nTrials);
                R = randvals*(1-exp(-itiRangeResAdj/lambda));
                itiSeqRes = -log(1-R)*lambda + itiMinResAdj;
                    % Lambda is not optimized
                itiSeq = round(itiSeqRes/res)*res;
            end  
    
            % Adjust some ITIs to make the real and requested means match
            % Difference between current and goal
            ch = 0;
            succ = 1; % success
            itiDiff = sum(itiSeq) - itiMean*nTrials;
            while (abs(itiDiff) > res) || (mean(itiSeq) > itiMean)
                changeableIdx = find(itiSeq ~= itiMin & itiSeq ~= itiMax); 
                    % the ITI we change shouldn't be min or max
                randIdx = changeableIdx(randi(length(changeableIdx)));
                itiSeq(randIdx) = itiSeq(randIdx) - sign(itiDiff)*res;
                ch = ch+1;
                    % This part of the code is based on Neurodesign's _fix_iti method in the
                    % 'generate' class ('generate.py').

                itiDiff = sum(itiSeq) - itiMean*nTrials; 
                if ch > nTrials
                    succ = 0; % if we have to adjust ITIs nTrial times, maybe we're going too far and we'll start over
                    break;
                end     
            end

            if succ
                break;
            end

            % Keep track of number of while loops
            i = i + 1;
            if i > 5, error('Going through while loop too many times. We can''t get nominal and actual mean to match.'); end
        end
    end
    
    itiSeqAxRuns{n} = itiSeq;
end

%% Sanity check for duration
for n = 1:nRuns
    usedDur = sum(itiSeqAxRuns{n}) + condDurs(trialSeqAxRuns{n});
    if usedDur > runDur
        error('The total run duration is exceeding requested amount. Time to troubleshoot.');
    end
end

%% Assemble runOutlines
runOutlinesAxRuns = cell(1, nRuns);

for n = 1:nRuns
    % Gather parameters of this run
    trialSeq = trialSeqAxRuns{n};
    itiSeq = itiSeqAxRuns{n};
    nEvents = sum(epochsPCond(trialSeqAxRuns{n}));
    nTrials = length(trialSeqAxRuns{n});

    runOutline = getRunOutlineOfHeight(nEvents);

    if isempty(trialSeq)
        % If addExtraTrials is off and preciseCondProbs is on, sometimes a
        % run will have no trials at all.
        runOutline = getRunOutlineOfHeight(0);
        
    else % if trialSeq actually contains events
        bIdx = 1; eIdx = epochsPCond(trialSeq(1));
        for t = 1:nTrials
            % Add trial
            slice = [epochIDs{trialSeq(t)}; ones(1, epochsPCond(trialSeq(t)));...
                taskDurs{trialSeq(t)}; zeros(1, epochsPCond(trialSeq(t)))];
            runOutline{bIdx:eIdx, :} = slice';

            % Add ITI
            runOutline{eIdx+1, :} = [0 0 itiSeq(t) 0];

            % Update indices
            if t == nTrials, break; end
            bIdx = eIdx+2; 
            eIdx = bIdx + epochsPCond(trialSeq(t+1)) - 1;
        end

        % Compute onset times
        runOutline.Onset = [0; cumsum(runOutline.Duration(1:end-1))];

        % Remove ITIs from runOutlines
        runOutline(runOutline.ID == 0, :) = [];
        
        % Add row with the ID of the last condition
        % runOutline(end+1, :) = [maxEpochID 0 0 0];
        % This causes issues later on with zero columns during trial-wise
        % regression.
    end
    
    % Assign to runOutlines (across runs)
    runOutlinesAxRuns{n} = runOutline;
end

function runOutline = getRunOutlineOfHeight(height)
    runOutline = table('Size', [height, 4],...
        'VariableNames', {'ID', 'Activity', 'Duration', 'Onset'},...
        'VariableTypes', {'double', 'double', 'double', 'double'});
end

end