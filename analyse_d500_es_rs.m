%% analyse_d500_es_rs.m
% ES_RS DMS500 Experiment File (DM500 Experiment)
% Updated Dec 19, 2018 - updated to include both RSVP300 and DMS500 analysis
% Updated March 27, 2020 - updated to only look at DMS500
% Updated May 12, 2023 - updated to run on HOMER_B550 (PC)
%% TASK DESCRIPTION:
% This is a DMS task...

[d500specs.fList,d500specs.pList] = matlab.codetools.requiredFilesAndProducts('analyse_d500_es_rs.m');

clc; clear; close all;
dbstop if error;
global exptdata

if ispc
    rootdir='C:\Users\ahbel\OneDrive - King''s College London\MATLAB';
else % MAC
    rootdir='~/Documents/MATLAB/';
end
addpath(userpath);
addpath(genpath([rootdir,filesep,'ephys_projects']));
addpath(genpath([rootdir,filesep,'ephys_projects',filesep,'code_d500_es_rs']));
addpath(genpath([rootdir,filesep,'Common_Functions']));
addpath(genpath([rootdir,filesep,'MonkeyLogic'])); % may need to update

ephys_analysis_defaults;
exptdata.lastModified=datetime('today');
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

% Modifiable defaults and directories (Change These)
exptdata.analysisName='D500_ES-RS_Study'; % used for savenames, figures, etc. (pick whatever you want; will be used for filenames)
exptdata.projectdir=[exptdata.analysisdir,exptdata.analysisName,filesep]; mkdir(exptdata.projectdir); %
exptdata.figuredir500=[exptdata.projectdir,'figures',filesep,'d500_es_rs',filesep]; mkdir(exptdata.figuredir500);
mkdir([exptdata.figuredir500,'matlabFigFiles']); mkdir([exptdata.figuredir500,'neuronPrintouts']);
diary([exptdata.projectdir,lower(exptdata.analysisName),'_',exptdata.analysisName,'.txt']);

% Preprocessing Parameters
exptdata.behav_samplingFreq=1000;
exptdata.spike_samplingFreq=1000;
exptdata.LFP_samplingFreq=1000;
exptdata.gaussian_kernel=10; % gaussian sd = 10ms

% Analysis Parameters
exptdata.xrange_psths=-250:500; % window surrounding stimulus onset (in ms)
exptdata.reprocess=0; % recreate trialised files even if they already exist (preprocess) (includes prepping megaMatrices)
exptdata.reviewNeurons=0; % determines whether to generate figure for each neuron
exptdata.behav_reprocess=0; % set to one if you want to regenerate output files for existing files

% Save EXPTDATA structure
save([exptdata.analysisdir,filesep,exptdata.analysisName,'_exptdata.mat'],'exptdata');

%% STAGE 1: Compile and analyse behavioural data from both monkeys 
% (RUN THIS BEFORE NEURAL PREPROCESSING; does not require access to external drive)

exptdata.paradigm={'DMS500'};
do_ephys_compileBehavData; % should only have to run once
disp('Press a key to continue...'); pause;
do_d500_es_rs_behaviouralAnalysis1; % behavioural analysis that uses output of do_ephys_compileBehavData and includes regression analysis
disp('Press a key to continue...'); pause;

%% STAGE 2: (Pre)Process neural data for each monkey
do_ephys_processAllnevFiles; % run only once
disp('Press a key to continue...'); pause; 
 
exptdata.allMonkeyNames={'Vortex','Vulcan'};
for mm=1:length(exptdata.allMonkeyNames)
    exptdata.monkeyname=exptdata.allMonkeyNames{mm}; % monkey name

    % Preprocess! % requires access to raw data files and processed NEV files directory (external Drive)
    ephys_d500_es_rs_preprocessNeuralData;
    close all
    ephys_d500_es_rs_prepareMegaMatrices;
end

%% STAGE 3: Quality Assurance - doesn't require access to NEV files (external drive)
%do_d500_es_rs_export2Excel; % Up to date as of Feb 25, 2019 - only needs to be
% run if there are changes (Feb 25, 2019 - Once I mark TE neurons, I should run this again)
do_d500_es_rs_sessionQA; % Up to date as of Feb 18, 2019 - only needs to be run if there are changes

%% STAGE 4: Analyse neuronal data!
% May 30, 2023 - what is difference between these three analysis steps?)
do_d500_es_rs_models
do_d500_es_rs_regressionAnalysis; % incomplete and needs to be updated
do_d500_es_rs_univariateAnalysis1;

%% STAGE 5: Produce Figures
do_d500_es_rs_makePrettyFigures;



diary off

%% STAGE ??? LFP, CSD Data

