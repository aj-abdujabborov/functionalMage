classdef fm_hrf
%FM_HRF Get and cache HRFs 
%
% Static methods
%   > hrfs = fm_hrf.getAllNsdHrfs(TR) will return a matrix of HRFs
%     extracted from the Natural Scenes Dataset where each column is an HRF
%   > hrf  = fm_hrf.getNsdHrf(TR) will return a random NSD HRF
%   > hrf  = fm_hrf.getCanonicalHrf(TR) will return the double gamma HRF
%     from SPM
%   > hrf  = fm_hrf.getSimTbHrf(TR) will return a random SimTB HRF
%   > hrfs = fm_hrf.getDerivativeHrfs(TR, numDerivatives) will return the
%     double-gamma canonical HRF along with its first or first & second
%     derivatives
%
% Part of package funkyMage. November 2023.
% https://github.com/aj-abdujabborov/funkyMage

    properties (Constant = true)
        canonicalParams = [6 16 1 1 6 0 40];
    end

    methods (Static = true)
        function basisData = getAllNsdHrfs(TR)
            fm_hrf.checkTRValidity(TR);

            load('NaturalScenesDatasetHRFs.mat', 'params');
            for i = 20:-1:1
                basisData(:,i) = simtb_spm_hrf(TR, params(i,:));
            end
            basisData = basisData ./ max(basisData, [], 1);
        end

        function basisData = getNsdHrf(TR)
            fm_hrf.checkTRValidity(TR);

            load('NaturalScenesDatasetHRFs.mat', 'params');
            basisData = simtb_spm_hrf(TR, params(randi(20),:));
            basisData = basisData ./ max(basisData);
        end

        function basisData = getCanonicalHrf(TR)
            fm_hrf.checkTRValidity(TR);

            basisData = simtb_spm_hrf(TR, fm_hrf.canonicalParams);
            basisData = basisData ./ max(basisData);
        end

        function basisData = getSimTbHrf(TR)
            fm_hrf.checkTRValidity(TR);

            % Code from SimTB package for MATLAB
            % https://pubmed.ncbi.nlm.nih.gov/22178299/
            params = nan(1, 7);
            params(1) = 4+3*abs(randn(1));     % delay of response (relative to onset)
            params(2) = 12+3*abs(randn(1));    % delay of undershoot (relative to onset)
            params(3) = 1+0.2*randn(1);        % dispersion of response
            params(4) = 1+0.2*randn(1);        % dispersion of undershoot
            params(5) = 2+3*abs(randn(1));     % ratio of response to undershoot
            params(6) = 0;                     % onset (seconds)
            params(7) = 40;                    % length of kernel (seconds)

            basisData = simtb_spm_hrf(TR, params);
            basisData = basisData ./ max(basisData);
        end

        function basisData = getDerivativeHrfs(TR, numDerivatives)
            arguments
                TR (1,1) {mustBePositive};
                numDerivatives (1,1) {mustBeNonnegative, mustBeLessThan(numDerivatives, 3)};
            end
            
            [basisData, p] = simtb_spm_hrf(TR, fm_hrf.canonicalParams);
            
            % Based on SPM12 code
            if numDerivatives >= 1
                dp = 1;
                p(6) = p(6) + dp;
                D = (basisData(:,1) - simtb_spm_hrf(TR,p))/dp;
                basisData = [basisData D(:)];
                p(6) = p(6) - dp;
            end
            
            if numDerivatives == 2
                dp = 0.01;
                p(3) = p(3) + dp;
                D = (basisData(:,1) - simtb_spm_hrf(TR,p))/dp;
                basisData = [basisData D(:)];
            end

            basisData = fm_hrf.orthonormalize(basisData);
        end

        function basisData = getHrfsCache(TR, numHrfsToGet, cacheLocation, opts)
            arguments
                TR (1,1) {mustBePositive};
                numHrfsToGet (1,1) {mustBePositive};
                cacheLocation {mustBeText};
                opts {mustBeA(opts, 'struct')};
            end

            opts  = fillDefaultParameters(opts);
            cache = makeOrLoadCache(opts);

            numSubRngs            = length(cache.hrfsSubRngDivided);
            numHrfsToGetPerSubRng = computeNumHrfsToGetPerSubRng(numSubRngs);
            basisData             = getHrfsForEachSubRng(cache, numHrfsToGetPerSubRng);

            function opts = fillDefaultParameters(opts)
                if ~exist('opts', 'var')
                    opts = [];
                end
                if ~isfield(opts, 'numHrfsToStore')
                    opts.numHrfsToStore = 30000;
                end
                if ~isfield(opts, 'library')
                    opts.library = 'SimTB';
                end
                if ~isfield(opts, 'correlationWithCanonical')
                    opts.correlationWithCanonical = [];
                end
                if ~isfield(opts, 'binWidth')
                    opts.binWidth = 0.05;
                end
            end

            function cache = makeOrLoadCache(opts)
                makeNew = true;
                if exist(cacheLocation, 'file')
                    tmp = load(cacheLocation);
                    if ~isfield(tmp, 'cache')
                        error("Existing file was not created by fm_hrf. Delete or move file.");
                    end

                    cache = tmp.cache;
                    if isequal(cache.opts, opts)
                        makeNew = false;
                    end
                end
    
                if makeNew
                    cache = fm_hrf.getHrfsWithCorrelation(TR, opts);
                    save(cacheLocation, 'cache');
                end
            end

            function numHrfsToGetPerSubRng = computeNumHrfsToGetPerSubRng(numSubRngs)
                numHrfsToGetPerSubRng = floor(numHrfsToGet/numSubRngs);
                numHrfsToGetPerSubRng = repmat(numHrfsToGetPerSubRng, 1, numSubRngs);
                numRemaining = numHrfsToGet - sum(numHrfsToGetPerSubRng);
                subRngsIdxToIncrease = randperm(numSubRngs, numRemaining);
                numHrfsToGetPerSubRng(subRngsIdxToIncrease) = numHrfsToGetPerSubRng(subRngsIdxToIncrease) + 1;
            end

            function basisData = getHrfsForEachSubRng(cache, numHrfsToGetPerSubRng)
                for i = length(numHrfsToGetPerSubRng):-1:1
                    numHrfsStoredInCurrRng = width(cache.hrfsSubRngDivided{i});
                    if numHrfsStoredInCurrRng < numHrfsToGetPerSubRng(i)
                        error("Not enough HRFs in range of hrfsCorrelationWithCanonical.\n" + ...
                              "Try increasing the floor of the range or leave parameter empty.");
                    end
                    idx = randperm(numHrfsStoredInCurrRng, numHrfsToGetPerSubRng(i));
                    basisData{i} = cache.hrfsSubRngDivided{i}(:,idx);
                end
                basisData = cat(2, basisData{:}); 
            end
        end

        function cache = getHrfsWithCorrelation(TR, opts)
            arguments
                TR (1,1) {mustBeFinite};
                opts {mustBeA(opts, 'struct')};
            end

            cache.opts = opts;

            % Shortened variable
            corrRange = opts.correlationWithCanonical;

            hrfsPool = gatherHrfs();
            if isempty(corrRange)
                cache.hrfsSubRngDivided = {hrfsPool};
                cache.corrValuesSubRngDivided = [];
                return; 
            end

            corrValues = computeCorrelationWithCanonical(hrfsPool);
            cache.hrfsSubRngDivided = subdivideHrfsIntoBins(hrfsPool, corrValues);

            function hrfsPool = gatherHrfs()
                switch lower(opts.library)
                    case 'simtb'
                        hrfFunction = @fm_hrf.getSimTbHrf;
                    case 'nsd'
                        hrfFunction = @fm_hrf.getNsdHrf;
                end
            
                hrfsPool = [];
                h = opts.numHrfsToStore;
                while h > 0
                    hrfsPool(:,h) = hrfFunction(TR); %#ok<AGROW>
                    if ~any(isnan(hrfsPool(:,h)))
                        h = h - 1;
                    end
                end
            end

            function corrValues = computeCorrelationWithCanonical(hrfsPool)
                canonicalHrf = fm_hrf.getCanonicalHrf(TR);
                for h = opts.numHrfsToStore:-1:1
                    tmp = corrcoef(hrfsPool(:,h), canonicalHrf);
                    corrValues(h) = abs(tmp(1,2));
                end
            end

            function hrfsSubRngDivided = subdivideHrfsIntoBins(hrfsPool, corrValues)
                numSubRngs = min(opts.numHrfsToStore, ceil(range(corrRange)/opts.binWidth));
                subRngs = linspace(corrRange(1), corrRange(2), numSubRngs+1);
    
                hrfsSubRngDivided = {};
                for j = numSubRngs:-1:1
                    currRng = subRngs(j:j+1);
                    indsInRng = (corrValues > currRng(1)) & (corrValues <= currRng(2));
                    hrfsSubRngDivided{j} = hrfsPool(:,indsInRng);
                end
            end
        end
    end

    methods (Static = true, Access = protected)
        function [Q, R] = orthonormalize(X)
            % Modified Gram-Schmidt.
            % G. W. Stewart, "Matrix Algorithms, Volume 1", SIAM, 1998.
            % https://blogs.mathworks.com/cleve/2016/07/25/compare-gram-schmidt-and-householder-orthogonalization-algorithms/
            [n,p] = size(X);
            Q = zeros(n,p);
            R = zeros(p,p);
            for k = 1:p
                Q(:,k) = X(:,k);
                for i = 1:k-1
                    R(i,k) = Q(:,i)'*Q(:,k);
                    Q(:,k) = Q(:,k) - R(i,k)*Q(:,i);
                end
                R(k,k) = norm(Q(:,k))';
                Q(:,k) = Q(:,k)/R(k,k);
            end
        end

        function checkTRValidity(TR)
            assert(numel(TR) == 1 && TR > 0, "TR should be a single value larger than 0");
        end
    end
end