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

classdef fm_mvpa
    properties
        informationMetric = true;
        classificationAccuracyMetric = true;
    end

    methods

    end

end



