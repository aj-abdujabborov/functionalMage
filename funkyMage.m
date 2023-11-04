classdef funkyMage < matlab.mixin.Copyable
    methods (Static = true)
        function output = getCacheDirectory()
            output = fullfile(fileparts(which('funkyMage.m')), 'cache/');
            if ~exist(output, 'dir')
                mkdir(output);
            end
        end
    end
end