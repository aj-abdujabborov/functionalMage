
classdef fm_simulation < matlab.mixin.Copyable
    properties
        trialsTable;

        TR
        nRuns

        conditions

        boldData

        noise
    end
 
    methods
        function obj = fm_simulation()
            obj.makeTrialsTable();
        end

        function generate()
            % produce ground-truth neural patterns
            % make runOutlines
            % createNTS
            % convolve with HRF
            % add noise
        end

        function saveToNifti(sampleNii)
            % take sample nifti and replace it with current data
            % or you could try to generate new data from scratch.
        end

        function generateNoise()
            fm_Noise(noiseParams, noiseData);
        end

        function saveToMat(sampleNii)

        end

    end

    methods (Access = protected)
        function makeTrialsTable(obj)
            columnNamesTypes = [["Probability", "double"]; ...
			            ["NeuralPatternIDs", "int16"]; ...
			            ["AnalysisIDs", "int16"]; ...
			            ["TaskDurations", "double"]; ...
			            ["TaskOnsets", "string"]];
            obj.trialsTable = table('Size', [0, size(columnNamesTypes,1)],... 
	            'VariableNames', columnNamesTypes(:,1),...
	            'VariableTypes', columnNamesTypes(:,2));
        end
    end
    
end