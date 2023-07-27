%% ideas
% --> Every fM object can have only one design indicated (if you want to
% compare fixed design with partial design, you need to make two objects).
% --> Inherit "MVPA" and "Univariate" from a functionalMage class.
% --> functionalMage() should be an optinal class. it should not be
% necessary to use it.

%% testing
parentPath = fileparts(which("functionalMage.m"));
addpath(genpath(parentPath));

fm = functionalMage();
fm.taskTable = [fm.taskTable;...
    {3 "3 3 3" "0 6 9" "1 0.3 0.15", "1 2 3", "1 2 3" "0 1 3"};...
    {1 "3" "0" "1" "1", "1" "0"};...
    {1 "3 3" "0 3" "0.11 1", "1 2", "1 2" "0 1"}
    ];
fm.simProperties.numRuns = 2;
fm.simulate();

aInfo = fm_analysisInfo(fm.taskTable, fm.simulation.eventList);
% aInfo.taskTable = fm.simulation.taskTable; % functionalMage should do this automatically
% aInfo.eventList = fm.simulation.eventList; % functionalMage should do this automatically
aInfo.newDurationsAndOnsets = []; % something here
aInfo.HRFs = 1; % convolution has no effect // also allow to set it [].
aInfo.Name = "Call this something you understand";
fm.analysisList.add(aInfo);

fm.analyze(); % does both GLM and the MVPA afterwards
fm.results; % results maybe a structure. or should it be a matrix? 
% or should it be a matrix with fm_data or fm_results type stuff?

%% testing without functionalMage class

taskTable = fm_taskTable();
taskTable.content = [taskTable.content;...
    {1, "3 1 3", "0 3 4", "1 1 1", "1 2 3", "1 1 2", "0 0 1"};...
    {1, "3 1 3", "0 3 4", "1 1 1", "1 2 4", "1 1 3", "0 0 1"};...
    ];

simProperties = fm_simulationProperties();
simProperties.TR = 0.5;

dm = fm_designMatrix(taskTable, simProperties);
dm.runwiseAnalysisIDs = [1 2];
glmLSA = dm.glmLSA;

mvpaLSA = dm.mvpaLSA;

simulation = fm_simulation(dm.neuralPatternIDEventList, simProperties);
simulation.generate();

analysisList = fm_analysisList;
glm.framework = "LSS";
glm.hrf = "Library";
glm.basis = randn(1,5);
analysisList(1) = glm;
analysisList.add('glm.framework', 'LSS',...
    'glm.hrf', 'Library',...
    'glm.basis', randn(10,5),...
    'mvpa.metric', 'information');
analysisList.display();

cellArrayOfAnalysisResults = analysisList.perform();

glm = fm_glm(simulation.boldTimeSeries, dm.glmLSA);
glm.hrf = "nsd+canonical";
glm.execute();
glm.framework = "LSS";
glm.basis = randn(1,5); % [some HRF matrix or HRF parameters]


betas = glm.betas;
results = fm_mvpa(dm.mvpaLSA, betas);

% NEXT STEP:
% write out the code to generate the GLM table and the MVPA table. (We
% don't want to throw in all analysis paramters into fm_analysisInfo. That
% will be muddy to the user.) If GLM and MVPA tables can be made
% semi-independently, we should go for a separated out approach. Maybe
% fm_glm will hold all the GLM parameters, fm_mvpa will hold all the MVPA
% parameters and we don't need a different class. Ideally, taskTable and
% the pure version of eventList will make it into neither of these classes.
% If there's some way to simplify all of what's going on with
% RegressionID/ClassificationGroups/NeuralPatterns into a simpler structure
% (or somehow make all this information *inherent* to the outputs that
% fm_glm / fm_mvpa produce) then this would reduce the complexity of this
% whole thing to the user. 

% sketch:
% aInfo = fm_analysisInfo(fm.taskTable, fm.simulation.eventList);
% aInfo.newDurationsAndOnsets = [];
% aInfo.Name = "Call this something you understand";
% glmDesignMatrix = aIinfo.glmDesignMatrix;
% 
% fm_glm(glmDesignMatrix, simulation.boldTimeSeries);

%% inputting multiple analysis methods possibility
glmMethods = ["LSA", "LSS"];
hrfMethods = ["SPM12", "NSD"];

for i = 1:length(glmMethods)
    for j = 1:length(hrfMethods)
        fm.analysisList.add(...
            "GLM", glmMethods(i),...
            "HRF", hrfMethods(j),...
            "newDurationsAndOnsets", fm.taskTable(:,["Onsets", "Durations"]));
    end
end

%% reminders
fmCopy = copy(fm); % to copy; fmCopy = fm will just copy the handle
delete(fm); % to delete; clear fm will just delete the handle


%%
TheNumbers = [-10 3 2];
old = [-10 2 2 3 2];
new = [-1 4 4 3 2];
[found, idx] = ismember(TheNumbers, old);
replaced_numbers = TheNumbers;
replaced_numbers(found) = new(idx(found));