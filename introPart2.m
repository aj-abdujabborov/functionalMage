%% 
% We'll repeat the same thing we did in introPart1.m but see a few more
% aspects of funkyMage.

%% Hypothetical study
% Compare lilac perception to hydrangea perception
% Task: 2s visual stimulus, 1s response
% Conditions: (1) see lilac, (2) see hydrangea
% Analysis:
% * Identify which method among the following methods is best suited for
%   the study:
%       Dimension 1: Least Squares All vs Least Squares Unitary
%       Dimension 2: Canonical HRF vs using library-fitting with the
%       Natural Scenes Dataset HRFs
% * Use classification accuracy of the stimulus period as metric for which
%   method is better. In other words, let's design the simulation such that
%   higher classification values will indicate better "method performance".

%% Task
funkyMage.addToPath();

% The obvious way to use classification accuracy as a metric for which
% method will provide the most accurate activity estimates would be to
% do what we did in introPart1.m:
task = fm_task();
task.content = [
    {1, "2 1", "", "", "1 2", "1 2", "1 0"};...
    {1, "2 1", "", "", "3 4", "3 4", "1 0"}
    ];

% .. and then compare the stimulus period of one condition against another.
% This is problematic, however. Because the fMRI signal is temporally
% protracted, in trying to estimate the stimulus activity, a method could
% actually estimate response activity, and classification accuracy would
% still be high since response NeuralPatternIDs covary with stimulus
% NeuralPatternIDs i.e., stimulus activity is different between conditions,
% but so is response activity. We cannot know if this is happening with the
% above task design.

% The solution is to inependently vary NeuralPatternIDs (the neural
% activity patterns) between the stimulus and response phases. Thus, the
% more the activity from the response phase leaks into the _estimates_ of
% activity from the stimulus phase, the more orthogonal the classifier
% labels of the stimulus phase become with respect to the information
% contained in the estimates of stimulus activity, pulling classifier
% accuracy toward chance.

% Independently varying conditions between two states across the two task
% phases means doing this:
% NeuralPatternIDs    AnalysisIDs
% ----------------    -----------
%      "1, 1"           "1, 2"       
%      "1, 2"           "3, 4"       
%      "2, 1"           "5, 6"       
%      "2, 2"           "7, 8"        

% Rather than write out 4 rows of the task table, we can simply write one
% row and use the expansion feature of fm_task. This feature will take
% whatever we write in the the fm_task object and vary the NeuralPatternIDs
% independently between two states while mapping out the values in
% AnalysisID to unique values for each new condition. To use this feature,
% we need to set 'expand' to true *before* we set our content:

task = fm_task();
task.expand = true;
task.content = {1, "2 1", "", "", "1 2", "1 2", "1 0"};

% Voila. Now see the result:
disp(task.content)

%% Sequence
seq = fm_sequence(task);
seq.itiModel = 'exponential';
seq.itiParams = [2 5]; % the range of our truncated exponential ITI distribution is 2s to 5s
seq.numRuns = 5;
seq.runDuration = 400;
seq.go();

if seq.occupiedPercentage < 95
    warning('We have a problem: less than 95% of the run contains trials');
end

%% Design matrix
dm = fm_designMatrix(task, seq.eventList);
dm.runwiseAnalysisIDs = [];
% If we wanted any of the AnalysisIDs to be treated as one regressor within
% each run, we could put that AnalysisID into the above property. Right
% now, it should be empty.

% From the previous intro, we know that dm.getNeuralPatternIDEventList()
% will give us what we need to simulate our data. Let's check what exactly
% it's giving us by visualizing it.

% Let's check what class seq.eventList is
class(dm.getNeuralPatternIDEventList()) 
    % this show the output is of class fm_eventList

% The fact that dm.getNeuralPatternIDEventList() returns a vector of
% fm_eventList objects means that we can use the methods that are on
% fm_eventList (which again, you can check with "help fm_eventList").

% Let's visualize
neurIDEventList = dm.getNeuralPatternIDEventList();
visual2D = neurIDEventList(1).computeDesignMatrix(seq.TR);
imagesc(visual2D)
xlabel('Neural Pattern IDs')
ylabel('Time (s)')

% We're picking the first run from the output and calling the
% 'computeDesignMatrix' method while supplying the TR into it. This gives
% us 'visual2D', which we're supplying to the MATLAB 'imagesc' plotting
% function.

% There are two columns in the figure, and these two columns simply
% represent the NeuralPatternID that's occuring at each time point. This
% neural pattern info is useful for the simulation. Hence we use it there.

%% Simulation
sim = fm_simulation(dm.getNeuralPatternIDEventList());
sim.TR = seq.TR;

% There is a MATLAB package SimTB that can generate variable HRFs. We're
% specifying that we want a random HRF from SimTB for each voxel.
sim.hrfLibrary = "SimTB";
sim.go();

% SimTB: https://github.com/calhounlab/simtb

%% Get best-fitting HRFs
% We can supply multiple HRFs to fm_glm and have it figure out which HRF
% best fits each voxel. In this case, we will use a predefined set of HRFs
% from the Natural Scenes Dataset (NSD) study, referenced here:
% + Allen, E. J., St-Yves, G., Wu, Y., Breedlove, J. L., Prince, J. S.,
%   Dowdle, L. T., Nau, M., Caron, B., Pestilli, F., Charest, I., Hutchinson,
%   J. B., Naselaris, T., & Kay, K. (2022). A massive 7T fMRI dataset to
%   bridge cognitive neuroscience and artificial intelligence. Nature
%   Neuroscience, 25(1), 116–126.

glm = fm_glm(sim.boldTimeSeries, dm.getGlmLSA());

% We can ask fm_glm to find the best-fitting HRF for each voxel using the
% 'findbesethrfs' parameter.
glm.framework = 'findbesthrfs';
glm.hrf = 'nsd'; % this says to use the NSD HRFs
glm.go();

% Now let's save glm.hrfsIdx, which contains the index of the best-fitting
% HRF for each voxel
bestHrfsIdx = glm.hrfsIdx;

% Make an fm_mvpa object for the next section
mvpa = fm_mvpa();

%% Analyze
% Let's add MVPA-light to our path
addpath('~/Downloads/MVPA-Light-master/startup')
startup_MVPA_Light

% We're about to get our results. Let's preallocate a structure to store
% our classification accuracy results for each method.
caResults = [];
caResults.LSA.canonical = nan;
caResults.LSA.nsd = nan;
caResults.LSS.canonical = nan;
caResults.LSS.nsd = nan;

% We'll get the results in the same order
% ------- LSA -------
% We already have an fm_glm object. Let's just change some properties.
glm.eventList = dm.getGlmLSA();
glm.framework = "ols";
mvpa.labelsAndGroups = dm.getMvpaLSA();

% - Canonical -
glm.hrf = "canonical";
glm.hrfsIdx = [];
glm.go();

mvpa.betas = [glm.results.betas];
caResults.LSA.canonical = mvpa.getClassificationAccuracy();

% - Library NSD -
glm.hrf = 'nsd';
glm.hrfsIdx = bestHrfsIdx;
glm.go();

mvpa.betas = [glm.results.betas];
caResults.LSA.nsd = mvpa.getClassificationAccuracy();

% We didn't need to get dm.getMvpaLSA() again because it's the same
% regardless of which HRF we choose

% ------- LSS -------
glm.eventList = dm.getGlmLSS2(); % update the design matrix
glm.framework = 'lss';
mvpa.labelsAndGroups = dm.getMvpaLSS(); % now we update the labels and groups

% - Canonical -
glm.hrf = 'canonical';
glm.hrfsIdx = [];
glm.go();

mvpa.betas = [glm.results.betas];
caResults.LSS.canonical = mvpa.getClassificationAccuracy();

% - Library NSD -
glm.hrf = 'nsd';
glm.hrfsIdx = bestHrfsIdx;
glm.go();

mvpa.betas = [glm.results.betas];
caResults.LSS.nsd = mvpa.getClassificationAccuracy();

% Now we can examine each of the 4 values, but one subject isn't enough for
% speculation.

%% Mutli-subject
% Let's do the entire process from above but for multiple subjects while
% keeping the comments to a minimum. This can take a while because of LSS.

% Basics
funkyMage.addToPath();
addpath('~/Downloads/MVPA-Light-master/startup'); startup_MVPA_Light;

nSubjects = 10;

methodDim1Names = ["LSA", "LSS"];
methodDim2Names = ["Canonical", "NSD"];

% Prep
task = fm_task();
task.expand = true;
task.content = {1, "2 1", "", "", "1 2", "1 2", "1 0"};

seq = fm_sequence(task);
seq.itiModel = 'exponential';
seq.itiParams = [2 5];
seq.numRuns = 5;
seq.runDuration = 400;

sim = fm_simulation();
sim.TR = seq.TR;
sim.hrfLibrary = "SimTB";

glm = fm_glm();
mvpa = fm_mvpa();

% Go
for s = nSubjects:-1:1
    seq.go();
    dm = fm_designMatrix(task, seq.eventList);
    sim.neurIDEventList = dm.getNeuralPatternIDEventList();
    sim.go();

    glm = fm_glm(sim.boldTimeSeries, dm.getGlmLSA());
    glm.framework = 'findbesthrfs';
    glm.hrf = 'nsd';
    glm.go();
    bestHrfsIdx = glm.hrfsIdx;

    for m1 = methodDim1Names
        if m1 == "LSA"
            glm.eventList = dm.getGlmLSA();
            glm.framework = "ols";
            mvpa.labelsAndGroups = dm.getMvpaLSA();
        else
            glm.eventList = dm.getGlmLSS2();
            glm.framework = 'lss';
            mvpa.labelsAndGroups = dm.getMvpaLSS();
        end

        for m2 = methodDim2Names
            if m2 == "Canonical"
                glm.hrf = "canonical";
                glm.hrfsIdx = [];
            else
                glm.hrf = 'nsd';
                glm.hrfsIdx = bestHrfsIdx;
            end

            glm.go();
            mvpa.betas = [glm.results.betas];
            caResults(s).(m1+m2) = mvpa.getClassificationAccuracy();
        end
    end

    fprintf('Subject %d out of %d done\n', nSubjects-s+1, nSubjects);
end

% Plot results
clf;
data = [caResults.LSACanonical; caResults.LSANSD;...
        caResults.LSSCanonical; caResults.LSSNSD];
bar(1:4,mean(data, 2));
hold on
er = errorbar(1:4, mean(data, 2), std(data, [], 2)./sqrt(nSubjects));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
set(gca, 'XTickLabels', {'LSA+Canonical', 'LSA+NSD', 'LSS+Canonical', 'LSS+NSD'})
xlabel('Method')
ylabel('Classification accuracy')
plot(gca().XLim, [0.5 0.5], 'k--')
ylim([-0.4 1])