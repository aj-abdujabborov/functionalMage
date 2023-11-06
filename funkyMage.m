classdef funkyMage < matlab.mixin.Copyable
%FUNKYMAGE Mix of basic methods for the toolbox
% 
% Static methods
%   > funkyMage.addToPath() will add all funkyMage files to path
%   > funkyMage.getCacheDirectory() will get the directory into which data
%     can be cached
%
% Part of package funkyMage. November 2023.
% https://github.com/aj-abdujabborov/funkyMage

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