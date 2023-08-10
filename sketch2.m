% We're not going to have a general fm_analysis function. The separate glm
% and mvpa steps of that will be done in functionalMage. 


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