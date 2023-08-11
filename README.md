# FunctionalMage
_FunctionalMage_ is a MATLAB toolbox for simulating and analyzing functional magnetic resonance imaging (fMRI) data. With it, you can easily specify complex task designs, simulate data with predefined defaults, and obtain parameter estimates with a General Linear Model. Afterwards, you can compare multiple designs and methods with each other to assess which will best fit your real-world experiments, whether they implement multivariate pattern analyses (MVPA) or univariate analyses.

> FunctionalMage is in development. Certain capabilities, guides and tutorials are not yet available.

# Guide
### How to find help
Help can be found using ```help classOrFunctionName ``` in the MATLAB Command Window. A general guide is contained here in the `README.md` file. See also `example1.m` and `example2.m`.

### Available designs
This toolbox was written to be capable of generating complex designs. This is done through the fm_taskTable class, through which you can specify multiple task conditions, and within each condition, the proportion of that condition (how often it occurs in the experiment sequence) and durations & onsets of the epochs / segments of that condition.

### Available analyses
You can implement several types of analyses with this package.

* With respect to general linear models, you can easily implement Ordinary Least Squares, [Least Squares Separate](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3408801/) and [Least Squares Unitary](https://pubmed.ncbi.nlm.nih.gov/30352691/) analyses.
* With respect to hemodynamic models (HRFs), data can be analyzed with the canonical (double gamma) model, the canonical + temporal/dispersion derivatives model, and the [Library](https://pubmed.ncbi.nlm.nih.gov/34916659/) approach, which, for each voxel, estimates betas with an HRF from a set of predefined HRFs that best fits the voxel's data. The Library method can be used with HRFs created from the [Natural Scenes Dataset](https://pubmed.ncbi.nlm.nih.gov/34916659/) (NSD) or with customly supplied HRFs.
* _FunctionalMage_ can also be used to analyze real (not simulated) data. However, all data pre-processing (removing drift, smoothing, masking, etc.) should be done beforehand.

### Classes and objects
*FunctionalMage* utilizes [MATLAB classes](https://www.mathworks.com/help/matlab/matlab_oop/create-a-simple-class.html) to conceptually and implementationally divide the different components of data simulation and analysis. Classes are a way to group related **attributes** and **methods** into the same object, which is accessed through a variable. We create objects by using their template defintion. For example, let's make an `fm_data` object.

```
fmriData = fm_data();
```
Now, through our variable `fmriData`, we can **set** and **get** the attributes of our object. For example, we can set its `data` and `TR` fields.

```
fmriData.data = randn(100, 2);
fmriData.TR   = 2;
```
These fields have been pre-defined to exist in the fm\_data.m file, and the object automatically checks that our input is valid. Now let's make and add another fm\_data object into `fmriData(2)`.

```
fmriData(2).data = randn(150, 2);
fmriData(2).TR = 2;
```

We can think of these two objects as containing two runs of data. Now, we can call a function that is *embedded* into the fm\_data object.

```
concatenatedFmriData = fmriData.cat();
```

The `cat` method takes the attributes of fmriData objects, runs the procedure of concatenating them together, and returns a new fm\_data object. This is a simple example, but it shows how powerful classes can be. They allow code to mirror conceptual distinctions. This can make code easier to use, understand and scale.

### Simple Example
Here,  we will design a simple, two-condition, two-epoch task, simulate data, and anlayze it with the HRFs generated from the NSD and the canonical approach.

```
fm = functionalMage(); % instantiate functionalMage object
fm.taskTable = [fm.taskTable; ... % define task
    {1, "3 3", "", "", "1 3", "1 3", "0 1"};...
    {2, "3 3", "", "", "2 4", "1 4", "0 2"}];
fm.designMatrix.runwiseAnalysisIDs = [1];
fm.analysisList.add('Library', glm.hrf', 'nsd+canonical', 'glm.framework', 'OLS'); % add analysis parameters
fm.analysisList.add('Canonical', 'glm.hrf', 'canonical', 'glm.framework', 'OLS');
classAcc = fm.getClassificationAccuracy();
doesLibraryHaveHigherClassifAcc = classAcc.Library > classAcc.Canonical
```

### Classes list
You can use the functionalMage.m class as the sole interface through which the entire simulation and analysis process is done. However, using the individual classes directly, such as in example2.m, might be more effective in certain cases.

This is a list of classes in _FunctionalMage_.

| Class | Purpose |
| ----- | ------- |
| fm_taskTable | Specify task timing, ground-truth and analysis information |
| fm_eventList | Store design matrix info, including timing of each event |
| fm_designMatrix | Take an fm_taskTable object and return various design matrices |
| fm_simulationProperties | Contains simulation properties with predefined defaults |
| fm_simulation | Simulate fMRI data |
| fm_glm | Analyze fMRI data with OLS/LSS |
| fm_mvpa | Classify estimates or correlate ground-truth with estimates|


### Dependencies
_FunctionalMage_ was tested on MATLAB R2020a and later. 

It may require the [mvpa-light](https://github.com/treder/MVPA-Light) toolbox for classification analyses and the [CoSMoMVPA](https://github.com/CoSMoMVPA/CoSMoMVPA) toolbox for loading and saving NIfTI files.

## Author
FunctionalMage was written by AJ Abdujabborov at the [Sreenivasan Lab](https://www.sreenivasanlab.org) in 2023. Inquiries can be sent to [aj [dot] abdujabborov [at] nyu [dot] edu]().