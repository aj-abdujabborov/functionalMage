
classdef fm_trialsTable < table
    methods
        function r = vertcat(obj1, obj2)
            for i = 1:size(obj2, 2)
                obj2(1,i) = {obj2(1,i)};
            end
            r = vertcat@table(obj1, obj2);
        end
    end
end