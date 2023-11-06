classdef funkyMage < matlab.mixin.Copyable
    methods (Static = true)
        function output = getCacheDirectory()
            output = fullfile(fileparts(which('funkyMage.m')), 'cache/');
            if ~exist(output, 'dir')
                mkdir(output);
            end
        end
        
        function addToPath()
            location = fileparts(which('funkyMage.m'));
            if isempty(location)
                error('Cannot find funkyMage.m');
            end
            addpath(genpath(location));
        end
    end
end