%% Intro to MATLAB classes
% Classes are a way to bundle related *data* and computational *processes*
% into a single object by providing a template for how that object should
% behave. funkyMage is built on MATLAB classes. It provides you the
% templates, and you invoke them to create object instances.

% For example, there is an fm_simulation class for simulating fMRI data.
% You can get an object of that class simply by using the class name:
sim = fm_simulation();

% It might look like you're calling a function and receiving its output
% here, but you're not. You're creating an object based on a pre-built
% blueprint. Try:
disp(sim)

% You can see a variety of simulation-related "properties". These
% properties are "the data". Let's change one:
sim.numVoxels = 100; 
disp(sim) % you can see the change reflected here

% Let's try changing it to -1
sim.numVoxels = -1;

% This returns an error because the class has a predefined rule that
% numVoxels cannot be negative. If we try to add a new property, we will
% fail for the same reason:
sim.numSignificantResults = 99; % gives error

% Beyond these simple benefits, the major benefit of classes is that they
% also contain prebuilt functions called "methods". These are the
% "processes". These methods may use the properties you have supplied or
% the built-in defaults.

% So I can call the "go" method to simulate data after I supply all the
% necessary properties.
sim.TR = 0.5;
sim.noiseSD = 3;
sim.go(); % for now, gives error

% For now there's a missing property, but we will fix this issue in the
% following demo. Anyhow, once sim.go() is run, we can access the resulting
% data through the 'boldTimeSeries' property:
sim.boldTimeSeries()

% funkyMage documentation categorizes properties into "input properties"
% and "output properties", where "input properties" are used when any class
% method is called, and output properties are what the methods produce.
% (Note that this distinction is for convenience. There are no inherent
% input vs output properties in classes.)

% You can get more detail by using the help function, but that's for later.
help fm_simulation

%% Hypothetical study
% Compare lilac perception to hydrangea perception
% Task: 2s visual stimulus, 1s response
% Conditions: (1) see lilac, (2) see hydrangea
% Analysis: classify neural activity from the stimulus period

%% Task
% We will input the above information into a table. We can do this using
% the fm_task class. Let's get an fm_task object:
task = fm_task();

% Now we should fill in the task details and some core simulation and
% analysis parameters. We do this on task.content, which is a table
task.content % empty
summary(task.content) % get names of all the fields expected of us

% Let's fill the details. We will have a cell array for each condition,
% each each cell array will consist of 1 numeric value followed by 6
% strings (characters surrounded by double quotes).
task.content = [
    {1, "2 1", "", "", "1 2", "1 2", "1 0"};...
    {1, "2 1", "", "", "3 4", "3 4", "1 0"}
    ];

% Above we are saying: For the first condition:
% * Make its proportion of occurrences be 1, relative to other conditions
% * Make its durations 2s and 1 s
% * Calculate its onset times times automatically (since we left the third
%   string blank)
% * Set neural activation level to default of 100% (since we left fourth
%   string blank)
% * Give the two events different neural activation patterns (if we'd set
%   "1 1" they'd have the same pattern)
% * Treat the two events as analytically distinct (roughly meaning
%   different regressors)
% * Classify the first event (non-zero number), but do not classify the
%   second one

% For the second condition, things are the same except we're providing
% different NeuralPatternIDs, since the stimulus and (therefore) response
% are different from the first condition; and different AnalysisIDs, since
% the neural activities are different.

% Now let's see what we've made
disp(task.content)

% You can see the information we filled + the pre-built class functionality
% filled in the empty strings.

%% Sequence
% Now, since we have our design, we need to generate a random sequence of
% trials. For this we make an fm_sequence object
seq = fm_sequence(task);

% Note that the class needs a non-empty fm_task object if we want to use
% its methods, and it allows us to supply it the moment we create the
% fm_sequence object. You can see the different ways you can make a class
% object by doing "help fm_className" and looking at the "Constructors"
% section. Let's change some properties from their default.
seq.numRuns = 5;
seq.runDuration = 336;
seq.itiModel = 'fixed';
seq.itiParams = 5; % set a 5-second fixed ITI

% Let's call the "go" method, which generates the sequence
seq.go();

% Now let's see the first run
disp(seq.eventList(1).content)

% seq.eventList is a vector of fm_eventList objects. Just like numbers can
% be in an array, objects of the same class can be in an array

%% Design matrix
% fm_designMatrix is a one-stop shop for design-matrix related stuff. It
% will give us information to simulate data, perform a General Linear
% Model, and do an MVPA analysis.

% Let's see how we can call the fm_designMatrix object by looking at the
% "Constructors" section of:
help fm_designMatrix

% So there are several ways, but I'll use the following
dm = fm_designMatrix(task, seq.eventList);
% This object will become useful in all the following sections.

%% Simulation
% Now to simulate our data, we use the fm_simulation class, which requires
% that we supply it our sequence. We have to do this using the
% fm_designMarix object rather than the fm_sequence object, however. So
% let's get the sequence:
evList = dm.getNeuralPatternIDEventList();

% Now make the fm_simulation object and run it
sim = fm_simulation(evList);
sim.TR = seq.TR; % pass whatever the TR was in the fm_sequence object
sim.noiseSD = 3;
sim.go(); % first time might take a few seconds

% Plot the second voxel of the first run
clf;
plot(sim.boldTimeSeries(1).data(:,2))
xlabel('Time (s)')
ylabel('BOLD signal')
hold on; % do not close the plot just yet

%% General Linear Model
% For the GLM, we need our data and a model. Let's make an fm_glm object
% and first supply the data.
glm = fm_glm();
glm.fmriData = sim.boldTimeSeries;

% Now, as for the model we'll use, let's pick a Least Square All model,
% which will provide a neural activity estimate for every event. We can get
% this model using our fm_designMatrix object.
glm.eventList = dm.getGlmLSA();

% Specify a few other parameters
glm.hrf = "canonical";
glm.framework = "OLS"; % ordinary least squares
glm.saveFits = true;
glm.go();

% Let's plot the fit time series of the second voxel of the first run above
% the simulated data
plot(glm.results(1).fits.data(:,2))
legend('Data', 'Fit')

%% Multivoxel Pattern Analysis
% funkyMage uses the MVPA-Light toolbox to perform classification. So first
% download it to the 'Downloads' folder from here and extract it:
% 'https://github.com/treder/MVPA-Light'

% Assuming you're on Mac, add it to path:
addpath('~/Downloads/MVPA-Light-master/startup')
startup_MVPA_Light;

% Now, we make an fm_mvpa object and supply it the parameter estimates from
% all the runs. We can get that by:
betas = [glm.results.betas];

% Let's also get the labels for the betas. Once again, we use our
% fm_designMatrix object
labelsAndGroups = dm.getMvpaLSA();

% Let's make the object
mvpa = fm_mvpa(betas, labelsAndGroups);

% Now get classification accuracy
ca = mvpa.getClassificationAccuracy();
disp(ca)

msg = "";
if ca <= 0.5, msg = "not "; end
fprintf("Classification accuracy is %shigher than chance", msg);

% Because in the fm_task section, we only had one unique value under
% "ClassificationGroup", that value being 1, we only get 1 value
