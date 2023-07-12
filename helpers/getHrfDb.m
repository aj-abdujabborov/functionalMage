function [filtHrfDb, filtHrfCorr] = getHrfDb(R, nK, TR, hrfLibrary, dbLoc, hrfLen)
%GETHRFDB Find which HRF from a library best fits the BOLD data
% filtHrfDb = GETHRFDB(rng2Get, nGet, TR, hrfLibrary, dbLoc)
%
% <filtHrfDb> is a matrix of nGet HRFs with TR <TR>.
% <filtHrfCorr> is a matrix where each value is the correlation between
% the HRF in that index and the canonical.
%
% <R> is a range indicating the correlation with the canonical HRF
%   that you want. Numbers are Pearson's r and should be between 0 and 1.
% <nGet> is the number of HRFs to output.
% <hrfLibrary> is either 'simtb' or 'kk'.
% <dbLoc> is the location where the HRF database is / should be stored.
%
% Written by AJ.
% Sreenivasan Lab, 2022.

%% Checks
R = round(sort(R), 2);
assert(all(R >= 0), 'Minimum correlation value is 0.'); % even though the real minimum is -1
assert(all(R <= 1), 'Maximum correlation value is 1.'); % even though the real minimum is -1

%% Load existing database
pullHrfs = 1;

if exist(dbLoc, 'file')
    % Load file
    load(dbLoc, 'pars')
    parsAlgn = all([pars.TR == TR, strcmpi(pars.hrfLibrary, hrfLibrary), pars.hrfLen == hrfLen]);

    % See if it needs any modifications
    if parsAlgn
        load(dbLoc, 'baseKern', 'pCorrList');
        pullHrfs = 0;
    end
end

%% Pull HRFs
nBaseKern = 30000;

if pullHrfs
    % Get canonical HRF
    canHrf = getHRF('can', TR, 'hrfLen', hrfLen);
    
    % Pull
    baseKern = [];
    h = nBaseKern;
    while h > 0
        % Make HRF
        baseKern(:,h) = getHRF('var', TR, 'library', hrfLibrary, 'hrfLen', hrfLen); % matlab falsely showing inefficiency
        baseKern(:,h) = baseKern(:,h) / max(baseKern(:,h)); % normalize
        if ~any(isnan(baseKern(:,h)))
            h = h - 1;
        end
    end

    % Compute correlations
    pCorrList = nan(1, nBaseKern);
    for h = 1:nBaseKern
        % Get correlation value
        temp = corrcoef(baseKern(:,h), canHrf);
        pCorrList(h) = abs(temp(1,2));
    end

    % Save
    pars.TR = TR;
    pars.hrfLibrary = hrfLibrary;
    pars.hrfLen = hrfLen;
    save(dbLoc, 'baseKern', 'pCorrList', 'pars');
end

%% Divide by range
bin = 0.05;

% Determine some values
nLims = min(nK, ceil(range(R)/bin)); % number of ranges
lims = linspace(R(1), R(2), nLims+1); % the ranges
nKPerLimFloor = floor(nK/nLims); % least amount of kernels per range
nKPerLim = repmat(nKPerLimFloor, 1, nLims);

% Add extras
nKExtra = nK - nKPerLimFloor * nLims;
extraInds = randperm(nLims, nKExtra); % randomly add extra amount to some ranges
nKPerLim(extraInds) = nKPerLim(extraInds) + 1;

% Pick indices
inds2PickAll = cell(1, nLims);
for l = 1:nLims
    currRng = lims(l:l+1);
    indsInRng = (pCorrList > currRng(1)) & (pCorrList <= currRng(2));
    indsInRng = find(indsInRng);
    if nKPerLim(l) > length(indsInRng), error("There aren't enough HRFs in the database to produce this N HRfs in this range."); end
    inds2Pick = randperm(length(indsInRng), nKPerLim(l));
    inds2Pick = indsInRng(inds2Pick);
    inds2PickAll{l} = inds2Pick;
end

% Pick kernels
inds2PickAll = cat(2, inds2PickAll{:});
filtHrfDb = baseKern(:,inds2PickAll);
filtHrfCorr = pCorrList(:,inds2PickAll);