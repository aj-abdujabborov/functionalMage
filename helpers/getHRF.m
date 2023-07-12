function [hrf,hrfParams] = getHRF(type, TR, varargin)
%GETHRF Get a canonical or variable HRF from different libraries
%   [hrf, params] = GETHRF(type, TR, [parameters]). If type is 'can', get
%   the double gamma canonical HRF. If it is 'var', get a variable HRF.
%   Type can also be a vector of the 7 HRF parameters. TR is in seconds.
%   
%   Parameters:
%   Set 'Normalization' to 1 to normalize the HRFs with respect to the
%   canonical HRF (not by the max of the HRF you're getting!). Default is
%   1.
%   Set 'Library' to 'simtb' or 'kk' (Kendrick Kay) to select the library
%   from which you are picking the variable HRF. Default is 'simtb'.
%   Note: We make use of SimTB's HRF functions.

%   Written by AJ.

%% Parse inputs
p = inputParser;
validInput = @(x) any(strcmpi(x, {'var', 'can'})) || (isvector(x) && length(x) == 7);
addRequired(p, 'type', validInput);
addParameter(p, 'normalize', 1, @(x) length(x) == 1 & ismember(x, [0 1]));
addParameter(p, 'library', 'simtb', @(x) any(strcmpi(x, {'simtb', 'kk'})));
addParameter(p, 'hrfLen', [], @(x) numel(x) == 1);
parse(p, type, varargin{:});

% parameters
whichLibrary = p.Results.library;
toNormalize = p.Results.normalize;
hrfLen = p.Results.hrfLen;

%% Get HRF
% Canonical HRF parameters
canPars = [6 16 1 1 6 0 hrfLen];
canHrf = simtb_spm_hrf(TR, canPars);

for a = 1:10 % give 10 attempts to get a proper HRF, because sometimes we get all NaNs
    if strcmpi(type, 'can')
        hrf = canHrf;
        break; % HRF will be intact, so break for loop and move on
    end

    if length(type) == 7
        hrfParams = type;
    elseif strcmpi(whichLibrary, 'simtb') % then "type" MUST be 'var', so check for which library
        hrfParams = simtb_get_default_params(1);
    else
        kkparams = load('kkParams.mat'); kkparams = kkparams.params;
        hrfParams = kkparams(randi(size(kkparams, 1)), :);
    end

    if ~isempty(hrfLen)
        hrfParams(end) = hrfLen;
    end
    hrf = simtb_spm_hrf(TR, hrfParams);

    if all(~isnan(hrf)), break; end % if we have a non-NaN output, then we're done
    if a == 10, error('Somehow we''re repeatedly getting NaN filled HRFs.'); end
end

% Normalization
if toNormalize
    factor = 1/max(simtb_spm_hrf(TR, canPars));
    hrf = hrf .* factor;
end