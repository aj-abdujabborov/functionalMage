function convolvedResult = convolveByColumn(timeCourses, kernels, varargin)
%CONVOLVEBYCOLUMN Convolve corresponding columns of two matrices
% Y = CONVOLVEBYCOLUMN(timeCourses, kernels) convolves the columns of
% timeCourse and kernels. Y is the same size as timeCourses.
%
% Y = CONVOLVEBYCOLUMN(timeCourses, kernels, kernelIdx) convolves the
% columns of timeCourses with the columns of kernels as indicated by
% kernelIdx, which is a row vector whose index indicates the column of
% timeCourses and value indicates the column of kernels to convolve
% together.
%
% Y = CONVOLVEBYCOLUMN(A, B, [], 'trimOutput', false) turns off trimming
% post-convolution.

%% Perform
p = inputParser;
addOptional(p, 'kernelIdx', []);
addParameter(p, 'trimOutput', true, @(x) numel(x) == 1);
parse(p, varargin{:});

kernelIdx = p.Results.kernelIdx;
trimOutput = p.Results.trimOutput;

nTimeCourseRows = height(timeCourses);
nTimeCourseCols = width(timeCourses);
nKernelsRows = height(kernels);
nKernelsCols = width(kernels);

convolveCorrespondingColumns = isempty(kernelIdx);
if convolveCorrespondingColumns
    assert(nTimeCourseCols == nKernelsCols,...
        'timeCourses and kernels should have the same number of columns'); 
    kernelIdx = 1:nTimeCourseRows;
else
    assert(isrow(kernelIdx) && width(kernelIdx) == nTimeCourseCols,...
        'The number of columns in kernelIdx and timeCourse should be the same.');
end

convolvedResult = zeros(nTimeCourseRows+nKernelsRows-1, nTimeCourseCols);
for i = 1:nKernelsCols
    convolvedResult(:, kernelIdx == i) = conv2(timeCourses(:, kernelIdx == i), kernels(:,i));
end

if trimOutput
    convolvedResult(nTimeCourseRows+1:end,:) = [];
end
