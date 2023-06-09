
classdef fm_analysis < matlab.mixin.Copyable
    properties

    end

    methods
        % each analysis (e.g., LSA, Can) needs to have its own object to
        % keep things simple.
        
        function makeClassificationLabels()

        end
        function makeAwesomer()
            
        end

        function result = addDimension(obj, string)
            result = string;
            disp("Dope results");
        end
    end
end