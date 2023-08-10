%% ideas
% --> Every fM object can have only one design indicated (if you want to
% compare fixed design with partial design, you need to make two objects).
% --> Inherit "MVPA" and "Univariate" from a functionalMage class.
% --> functionalMage() should be an optinal class. it should not be
% necessary to use it.

%% testing
parentPath = fileparts(which("functionalMage.m"));
addpath(genpath(parentPath));

%% functionalMage class
fm = functionalMage();
fm.taskTable = [fm.taskTable;...
    {3 "3 3 3" "0 6 9" "1 0.3 0.15", "1 2 3", "1 2 3" "0 1 3"};...
    {1 "3" "0" "1" "1", "1" "0"};...
    {1 "3 3" "0 3" "0.11 1", "1 2", "1 2" "0 1"}
    ];
fm.simProperties.numRuns = 2;
fm.simulate();

fm.analyze(); % does both GLM and the MVPA afterwards
fm.results; % results maybe a structure. or should it be a matrix? 
% or should it be a matrix with fm_data or fm_results type stuff?

%% testing without functionalMage class
taskTable = fm_taskTable();
taskTable.content = [taskTable.content;...
    {1, "3 1 3", "0 3 4", "1 1 1", "1 5 3", "1  9 5", "2 0 1"};...
    {1, "3 1 3", "0 3 4", "1 1 1", "1 5 4", "2 10 6", "2 0 1"};...
    {1, "3 1 3", "0 3 4", "1 1 1", "2 5 3", "3 11 7", "2 0 1"};...
    {1, "3 1 3", "0 3 4", "1 1 1", "2 5 4", "4 12 8", "2 0 1"};...
    ];
% proportion, duration, onsets, activity, neuralPatternIDs, analysisIDs, classificaiton groups

simProperties = fm_simulationProperties();
simProperties.TR = 1;
simProperties.neuralFluctuationAmount = 0;

dm = fm_designMatrix(taskTable, simProperties);
dm.runwiseAnalysisIDs = [1];

simulation = fm_simulation(dm.neuralPatternIDEventList, simProperties);
simulation.generate();

glm = fm_glm(simulation.boldTimeSeries, dm.getGlmLSA());
glm.hrf = "nsd+canonical";
glm.framework = "findbesthrfs";
glm.execute();

glm.hrfsIdx = glm.results.hrfsIdx;
glm.framework = "ols";
glm.execute();

gtParams = glm.getGroundTruthInBetasForm(simulation.totalNeuralActivityPerEvent);

mvpa = fm_mvpa([glm.results.betas], dm.getMvpaLSA());
claAcc = mvpa.getClassificationAccuracy()
pattSim = mvpa.getPatternSimilarity(gtParams)

%%
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