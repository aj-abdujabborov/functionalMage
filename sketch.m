% IDEAS
% --> Every fM object can have only one design indicated (if you want to
% compare fixed design with partial design, you need to make two objects).
% --> Inherit "MVPA" and "Univariate" from a functionalMage class.

% functionalMage needs to be split into modules, so that it's easier to
% maintain. If everything is thrown into one muddhole class, things could
% get very messy. Split into "simulate" and "analyze" sections.
% Put a lot of basic methods into the parent functionalMage class, so that
% the derived MVPA and Univarite options are simpler to program.

fm = functionalMage();
fm.taskTable = [fm.taskTable;...
    {3 "3 3 3" "0 6 9" "1 1 1", "A B B", "1 2 3" "1 2 3"};...
    {1 "3" "0" "1" "A", "1" "1"};...
    {1 "3 3" "0 3" "1 1", "A B", "1 2" "1 2"}
    ];
fm.simulate();

% 3	3 3 3	0 6 9	1 1 1	1 2 3	A B B	1 2 3
% 1	3	0	1 1 1	1	A	1
% 1	3 3	0 6	1 1 1	1 2	 A B	1 2
% 1	3 3	0 6	1 1 1	1 3	 1 B	1 3
% 3	4 3 3	0 6 9	1 1 1	1 4 5	A B C	1 2 3



























%% 
delete(fm); clear fm;
fm = functionalMage();

T = fm.simulation.trialsTable;
T = [T;...
    {1, [1 2 3], ["A", "B"], [3 3 3], nan};...
    {1, [4 5 6], ["C", "D"], [3 3 3], nan}];

%%
fm = functionalMage('MVPA');

disp(fm.task.getConditionExample());
disp(fm.task.getConditionFull());
cond1 = struct('epochIDs', [1 2 3], ...
    'analysisIDs', ["", "1A", ""], ...
    'taskDurs', [3 9 3], ...
    'taskOnset', [0 3 12], ...
    'condProp', 4);
cond2 = struct( ...
    'epochIDs', [4 5 6], ...
    'analysisIDs', ["", "1B", ""], ...
    'taskDurs', [3 9 3], ...
    'taskOnset', [0 3 12], ...
    'condProp', 4);
fm.task.addCondition(cond1, 'Visual');
fm.task.addCondition(cond2, 'Auditory');

analysis = struct('GLM', 'LSA',...
    'HRF', 'SPM12');
fm.addAnalysis('SPM12', analysis);
fm.addAnalysis('Derivative', struct(...
    'GLM', 'LSA', ...
    'HRF', 'Derivative'));

dataProperties = struct('TR', 1, ...
    'Duration', 400, ...
    'Voxels', 50, ...
    'Noise', 10, ...
    'NumRuns', 5);
fm.data.setProperties(dataProperties);
disp(fm.data.getProperties());

fm.addMetrics('classification accuracy', 'information');

fm.run('subjects', 10);
results = fm.getResults();

% results is a structure. It has 'Analysis1' and 'Analysis2'.
% 'Analysis1'. 
% Analysis1.ClassificationAccuracy = [0.9 ...];
% Analysis1.Information = 0.8;

%% Using functionalMage for closer to my study
fm = functionalMage('MVPA');

fm.simulation.addCondition("");

fm.analysis.addDimension("GLM", ["LSA", "LSS"]);
fm.analysis.addDimension("HRF", ["SPM12", "KK"]);

results = fm.getResults();
% results is a structure.
% results.dimensional has size 2 x 2 representing the dimensions
% results.Anlaysis1 contains another analysis
% in other words, you can submit dimensions AND individual anlayses, which
% are output differently.

fmCopy = copy(fm); % to copy; fmCopy = fm will just copy the handle
delete(fm); % to delete; clear fm will just delete the handle

%%
stuff();
function stuff()
    x = 5;

    function stuff2()
        x = x + 1;
    end

    stuff2();

    disp(x);
end