% We're not going to have a general fm_analysis function. The separate glm
% and mvpa steps of that will be done in functionalMage. 


classdef fm_analysisInfo
    properties
        taskTable;
        newDurationsAndOnsets;
        HRFs;
        glmAndLabelsStucture;
    end

    % SHOULD WE DO GLM ANALYSIS FOR MVPA SEPARATELY FOR EACH RUN? would that be
% conceptually and implenetationally eaiser / clearer? 
% fm_data could store BOLD data, the fitts, the betas, SSRegr, etc.

    methods
        function set.newDurationsAndOnsets()
            if newDurationsAndOnsets ~= taskTable
                redoTheEventList;
            end

            % this has to be able to be treated as a new analysis method.
            % in other words, you need to be able to hold multiple analysis
            % methods. actually, maybe you can have multiple fm_analysises
            % and iterate over them and do your analysis.

            % redoTheEventList with new durations and onsets

            % check that the size of redoTheEventList matches 

            % allow to specify + and - in the string to indicate relative
            %   to onset time
            % think about system where pre/appending + or - to the numbers
            % indicates how to subtract. but this might be quite a niche
            % thing that not many people need to use... ppl can definitely
            % figure out how to adjust the stuff the way they want. so.
            % could be fine.
        end
    end

    methods
        function prepareLSU()
            % prepare design matrix for GLM
            % prepare eventList for classification
            % these are the only two things we need to execute LSU
        end

        function prepareLSS1()


        end


        function prepareLSS2()


        end

        function prepareLSA()


        end
    end

    methods
        function prepareCanonical()

        end

        function prepareLibrary()

        end

        function prepareDerivative()

        end

        function prepareCustom()
 
        end

        function prepareNoConvolve()

        end
    end
end



classdef fm_glm
    properties
        data = fm_data;

        glm = 'OLS'; % or LSS
        hrf = 'Canonical'; % or Library or Custom or Derivative
    end

    properties (Access = protected)
        C_inv_XtX_C = nan(1, nHRFs);
        designMatrixRanks = nan(1, nHRFs);

        SSRegr = nan(1, nVoxels);
        betas = nan(1, nVoxels);

        designMatrices = nan(1, nHRFs);
    end

    properties (Dependent = true)
        fitTimeSeries;
            % check if saveFitTimeSeries is enabled. if so, then return.
            % check if saveDesignMatrices is enabled. if so, then return
            %   designMatrix * betas;
            % produce error.
    end

    methods
        function fm_glm()
            % convert the data types to (single).
        end
    end

    methods % GLM methods
        function computeLSU()
            % do not do any collapsing of the AnalysisIDs in eventList
            % 

        end
    end

    methods (Static = true)
        function [betas, rsqs, boldDesignMatrix, SSRegr] = computeOLS(data, eventList, HRF)
            % single HRF
            % convolution, make design matrix
            % OLS
            % get betas, rsqs, boldDesignMatrix, SSRegr
        end

        function computeOLSDerivative(data, ~)
            % the important thing is that the inputs and outputs to these
            % functions must be the same as inputs to computeLSS /
            % computeOLS etc.

            % and this function should be placedi n the same place in the
            % sequence of steps to get to the end.
        end

        function performTTestDerivative(data, ~)
            % the important thing is that the inputs and outputs are the
            % same as performTTest. it should be possible for the user to
            % have written this function and fed a handle to it.
        end

        function [betas, rsqs, boldDesignMatrix, SSRegr] = computeLSS(data, eventList, HRF)
            % uses computeOLS to do LSS
            % just make one event stand out, and do it

            % TODO: run a test to see if it's better to convolve all events
            % separately first and then add columns together rather than do
            % convolution multiple times (might not save a whole lot of
            % time, this, since main time consumer is matrix inversion).
        end

        function [SSTotal] = computeSSTotal(data)

        end

        function TStat = performTTest()

        end


    end
end


classdef fm_mvpa
    properties
        informationMetric = true;
        classificationAccuracyMetric = true;
    end

    methods

    end

end

classdef fm_data
    properties
        data;
        TR;
        numRuns;
        numVoxels;
        numTRs;
        IDs;
    end

    properties
        numRows;
    end

    % numTRs is a dependent property that just outputs numRows
    % numVoxels is an actual property because this data is always
    %   anticipated to have numVoxels on the column side

    methods
        function save2disk()

        end

    end
end


