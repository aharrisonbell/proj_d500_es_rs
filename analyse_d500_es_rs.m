%% analyse_d500_es_rs.m
% ES_RS DMS500 Experiment File (DM500 Experiment)
% Updated Dec 19, 2018 - updated to include both RSVP300 and DMS500 analysis
% Updated March 27, 2020 - updated to only look at DMS500
% Updated May 12, 2023 - updated to run on HOMER_B550 (PC)
% Updated Jan 17, 2024 - updated for use by KCL students
% Updated Mar 14, 2024 - updated more
% Updated Apr 14, 2026 - updated even more (during Dadbatical)
clc; clear; close all;
program_version='1.1.2026-04-17';


%% Notes
% Uses: James Tursa (2024). MTIMESX - Fast Matrix Multiply with Multi-Dimensional Support (https://www.mathworks.com/matlabcentral/fileexchange/25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support), MATLAB Central File Exchange. Retrieved March 13, 2024.
% (13/03/2024 - no longer uses MTIMESX - replaced by MATLAB pagemtimes.m)

%% TASK DESCRIPTION:
% This is a DMS task...

%% SETUP DEFAULTS
dbstop if error;
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

if ispc
    %rootdir='C:\Users\ahbel\OneDrive - King''s College London\ephysProjects';
else % MAC
    rootdir='~/Library/CloudStorage/SynologyDrive-ephysProjects/';
    coderoot='~/MATLAB Drive/Projects';
end

addpath(genpath([coderoot,filesep,'commonProjectFunctions']));
addpath(genpath([coderoot,filesep,'proj_d500_es_rs_sua']));

ephys_analysis_defaults; % generates and loads config file
exptdata.lastModified=datetime('today');
exptdata.analysisName = 'd500_esrs'; % used for savenames, figures, etc. (pick whatever you want; will be used for filenames)

[d500specs.fList,d500specs.pList] = matlab.codetools.requiredFilesAndProducts('analyse_d500_es_rs.m');


% Preprocessing Parameters (most established by ephys_analysis_defaults.m)
exptdata.behav_samplingFreq=1000;
exptdata.spike_samplingFreq=1000;
exptdata.LFP_samplingFreq=1000;
exptdata.gaussian_kernel=10; % gaussian sd = 10ms

% Analysis Parameters 
exptdata.xrange_psths=-250:500; % window surrounding stimulus onset (in ms)
exptdata.reprocess=0; % recreate trialised files even if they already exist (preprocess) (includes prepping megaMatrices)
exptdata.reviewNeurons=0; % determines whether to generate figure for each neuron
exptdata.behav_reprocess=0; % set to one if you want to regenerate output files for existing files


qualityCutoff=2; % select only neurons that have at least this grade

% Analysis Parameters (many established by generate_td500_config.m)
exptdata.xrange_psths=-250:500; % window surrounding stimulus onset (in ms)
%exptdata.reprocess=1; % recreate trialised files even if they already exist (preprocess) (includes prepping megaMatrices)
%exptdata.reviewNeurons=1; % determines whether to generate figure for each neuron
%exptdata.behav_reprocess=1; % set to one if you want to regenerate output files for existing files
exptdata.figuredir=[exptdata.projectdir,'figures',filesep,'d500_es_rs',filesep]; mkdir(exptdata.figuredir);
mkdir([exptdata.figuredir,'matlabFigFiles']); mkdir([exptdata.figuredir,'neuronPrintouts']);


% Save EXPTDATA structure
save([exptdata.projectdir,filesep,exptdata.analysisName,'_exptdata.mat'],'exptdata');
diary([exptdata.projectdir,lower(exptdata.analysisName),'_',exptdata.analysisName,'.txt']);

%% --------------------------------------------------------------------------------------------------------------------
% STAGE 1: Compile data from both monkeys 
disp('+------------------------------------------------------------------------------------------------------+')
disp('| d500_esrs_paper.m - Main analysis program for EXPECTATION SUPPRESSION/REPETITION SUPPRESSION datasets|')
disp('+------------------------------------------------------------------------------------------------------+')
disp(['Program Version: ',program_version])
disp(['Project Directory: ', exptdata.projectdir])
disp(['Dataset Directory: ', exptdata.datalocation])
disp(['Figure Directory:  ', exptdata.figuredir])

%% STAGE 1: Compile and analyse behavioural data from both monkeys 
% (RUN THIS BEFORE NEURAL PREPROCESSING; does not require access to external drive)
exptdata.paradigm={'DMS500'};
do_ephys_compileBehavData; % should only have to run once
disp('Press a key to continue...'); pause;
do_d500_es_rs_behaviouralAnalysis1

%% STAGE 2: (Pre)Process neural data for each monkey
do_ephys_processAllnevFiles; % run only once
disp('Press a key to continue...'); pause; 
 
exptdata.allMonkeyNames={'Vulcan', 'Vortex'}; % 'Vortex',
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

