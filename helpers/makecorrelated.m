function [good] = makecorrelated(mat, corrVal)
%MAKECORRELATED Make uncorrelated vectors within a matrix correlated
%   Use Cholensky decomposition to do it

%% Go
% Checks
if corrVal == 0, good = mat; return; end
if abs(corrVal) > 1, error('Absolute of correlation value exceeds 1.'); end
if abs(corrVal) > 0.99, corrVal = sign(corrVal) * 0.99; end

% Basics
nCol = size(mat, 2);

% Memorize
matMins = min(mat, [], 1);
matMaxs = max(mat, [], 1);
matRng = matMaxs - matMins;

% Do
corrMat = eye(nCol);
corrMat(corrMat == 0) = corrVal;
patProps = chol(corrMat); % how to sum the two uncorrBase to create vectors with corrMat
bad = mat * patProps;

% Re-range the output to the original
badMins = min(bad, [], 1);
badMaxs = max(bad, [], 1);
badRng = badMaxs - badMins;
good = ( (bad-badMins)./badRng ) .* matRng + matMins;