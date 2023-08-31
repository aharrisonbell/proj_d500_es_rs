% do_es_rs_regressionAnalysis1(vargin)
% by AHB, started June 11, 2020, based on do_d500_es_rs_univariateAnalysis1.m
clc
vers_ephys_es_rs='1.0; June 11, 2020';
% 1.0 - original version (June 11, 2020) - written to be self-starting
% 1.1 - not, this is not currently functional

%%
dbstop if error;
if ispc
    rootdir='C:\Users\Andrew H Bell\iCloudDrive\Documents\MATLAB';
else % MAC
    rootdir='~/Documents/MATLAB/';
end

global exptdata timeEpochs500 epochNum epochNames evdata sampleRate
ephys_analysis_defaults;
exptdata.analysisName='D500_ES-RS_Study';
addpath(userpath);
addpath(genpath([rootdir,filesep,'ephys_projects']));
addpath(genpath([rootdir,filesep,'ephys_projects',filesep,'code_d500_es_rs']));
addpath(genpath([rootdir,filesep,'Common_Functions']));
addpath(genpath([rootdir,filesep,'MonkeyLogic'])); % may need to update
addpath(genpath([rootdir,filesep,'spm12']));

exptdata.projectdir=[exptdata.analysisdir,exptdata.analysisName,filesep]; mkdir(exptdata.projectdir); %
exptdata.regenerateDists=0;
exptdata.lastModified=date;
warning('off','MATLAB:MKDIR:DirectoryExists');

% =====================================================================================================================
% megamatrices AND behaviouralAnalysis must be complete before running this program
% What this program will do:

% Data structure (<BRAIN-AREA>_megaMatrix<TASK-CODE>NoScrubs.mat)
%  1) monkey number
%  2) session number (unique within animal but not across animal)
%  3) neuron number
%  4) new neuron number (unique across ALL sessions)
%  5) trial number
%  6) trial outcome (0=correct)
%  7) block number
%  8) condition number
%  9) reaction time
% 10) total number of neurons per session

% PARADIGM SPECIFIC STUFF:
% 11) 1st stim presented (1=stimA; 2=stimB)
% 12) 2nd stim presented (1=stimA; 2=stimB)
% 13) expectation (1=expect repeat; 2=expect alternation; 3=no expectation)
% 14) actual (1=repeat, -1=alternation)
% 15) active (1) or passive (0)
% 16) nothing
% 17) nothing
% 18) task number (300,500,600)
% 19) new block number (600 only)
% 20) automatic neuron classification (1=excitatory visual, -1 suppressed visual, 0=nonresponsive)
% 21) manual Classification (1=excitatory, -1=suppressed, 0=non-responsive)
% 22) quality Classification (0=shit, problematic, or non-responsive, 1=meh, 2=ok, 3=awesome)
% 23) stimA (task alternates between stim A and B)
% 24) stimB
% 25) distractor1
% 26) distractor2
% ... lots of spike density stuff

% =====================================================================================================================
exptdata.projectdir=[exptdata.analysisdir,exptdata.analysisName,filesep];
clc
cprintf('yellow','*====================================*\n')
cprintf('yellow','| do_d500_es_rs_regressionAnalysis.m |\n')
cprintf('yellow','*====================================*\n')
cprintf(['Version: ',vers_ephys_es_rs,'\n'])
disp(['Data location:    ',exptdata.datalocation]);
disp(['Processed NEV dir: ',exptdata.processedDatadir]);
disp(['Study name:       ',exptdata.analysisName]);
disp(['Project Directory: ',exptdata.projectdir]);

% =====================================================================================================================
%% 1.0 Load Behavioural and Neurophysiological Data
fprintf('\nLoading megamatrices and behavioural data...')
load([exptdata.projectdir,'DMS500_behavAnalysis.mat'],'behavData','evdata'); % contains BOTH monkeys behavioural data (this includes markers)

temp1=load([exptdata.projectdir,'Vortex_',exptdata.analysisName,'_megaMatrix_synced.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
temp2=load([exptdata.projectdir,'Vulcan_',exptdata.analysisName,'_megaMatrix_synced.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
V4_megaMatrix500NoScrubs=[temp1.V4_megaMatrix500NoScrubs; temp2.V4_megaMatrix500NoScrubs];
TE_megaMatrix500NoScrubs=[temp1.TE_megaMatrix500NoScrubs; temp2.TE_megaMatrix500NoScrubs];
clear temp*
sessions = unique(evdata.session);
fprintf('done.\n')

%% 2.0 Set up variables
% Data structures (DMS500): (-250:500 ms) (exptdata.xrange_psths)
timeEpochs500 = [...
    21 771;... aligned on FP on (1)
    772 1522;... aligned on 1stStim on (2)
    1523 2273;... aligned on 2ndStim on (3)
    2274 3024;... aligned on choice array (4)
    3025 3775;... aligned on correct choice (5)
    3776 4526;... aligned on incorrect choice (6)
    4528 5277;... aligned on reward (7)
    772 3024]; % the BIG three (8)
epochNames={'FPonset','1stStim','2ndStim','ChoiceArrayOnset','CorrectChoices','IncorrectChoices','Reward','1st/2nd/Choice'};
sampleRate=10; % test every X timepoints
epochNum=8; % limiting it to the more interesting epochs for the moment

%% 3.0 Correlation Matrix EVs
% Check correlation between explanatory variables 
% specify which variables to correlate here:
tmp_ev=[...
    evdata.active, ...
    evdata.actual, ...
    evdata.choseRepeat, ...
    evdata.cor, ...
    evdata.delta_p_rep, ...
    evdata.expect, ...
    evdata.prep, ...
    evdata.stim1, ...
    evdata.stim2,...
    evdata.stim1_faceNonFace, ...
    evdata.stim2_faceNonFace, ...
    evdata.prevtrial_n1, ...
    evdata.prevtrial_n2, ...
    evdata.prevtrial_n3, ...
    ];
paperData.corrcoeff.evLabels={'active','actual','choseRepeat','cor','delta-p-rep','expect',...
    'prep','stim1','stim2','stim1-face?','stim2-face?','prevtrial-n1','prevtrial-n2','prevtrial-n3'};
[paperData.corrcoeff.coeffs, paperData.corrcoeff.pvalues]=corrcoef(tmp_ev);
figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial');
clabel = arrayfun(@(x){sprintf('%0.2f',x)}, paperData.corrcoeff.pvalues); 
heatmaps(paperData.corrcoeff.coeffs, paperData.corrcoeff.evLabels, paperData.corrcoeff.evLabels, clabel,'ShowAllTicks',1,'Colormap','money');
title('EV Correlation Matrix', 'FontSize', 14); % set title

% Ameya Deoras (2020). Customizable Heat Maps (https://www.mathworks.com/matlabcentral/fileexchange/24253-customizable-heat-maps), 
% MATLAB Central File Exchange. Retrieved April 16, 2020.

% April 16, 2020 plotmatrix(tmp_ev) - this reveals that there is a problem with this
% correlation approach.  most EVs are binary, which means Pearsons,
% Spearmans, etc. coefficients are not appropriate.  need something
% different.
% Based on Googling - possible solutions include:
% tetrachoric correlation (polychoric)
% Matthijs J. Warrens (2008) called "Similarity Coefficients for Binary Data", from Leiden University,
% SPhi is sometimes attributed to Yule (1912)

% Can try Hamann's Similarity Coefficient (and, I believe, significance is
% tested against a Chi-square distribution)

% Moving on for now...


%% 4.0 Model Parameters
% Here is where to specify which regressors will be included
brainArea=1;
if brainArea==1
    currentData=V4_megaMatrix500NoScrubs;
    currentLabel='V4';
else
    currentData=TE_megaMatrix500NoScrubs;
    currentLabel='TE';
end
disp(['Analysing ',currentLabel,'...'])



modelData.modelName=['model1-',currentLabel,'-active+actual'];
modelData.test_evs=[...
    zscore(evdata.active), ... % active or pasive
    zscore(evdata.actual) ...
    zscore(evdata.delta_p_rep), ...
    zscore(evdata.prep), ...
    zscore(evdata.trialNumber)
    ];
modelData.test_ev_labels={'active','actual', 'pe', 'p-rep', 'trialNum'};
modelData.regressType='linear';


%% 4.1 Prep Regressors (using SPM_ORTH)
modelData.test_evs_new=modelData.test_evs*0;
tmp_new_order=1:length(modelData.test_ev_labels);
for e=1:length(modelData.test_ev_labels) % scroll through each ev
    disp(['...processing regressor ', num2str(tmp_new_order(end))])
    tmp_orth = spm_orth(modelData.test_evs(:,tmp_new_order));
    modelData.test_evs_new(:,tmp_new_order(end)) = tmp_orth(:,end);
    tmp_new_order=circshift(tmp_new_order, 1);
end
figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.2 0.2 0.4 0.3]); set(gca,'FontName','Arial');
subplot(1, 2, 1)
h=heatmap(corrcoef(modelData.test_evs));
h.XDisplayLabels = modelData.test_ev_labels;
h.YDisplayLabels = modelData.test_ev_labels;
%h.title = 'Before Variance Partitioning';
clf
figure
h=heatmap(corrcoef(modelData.test_evs_new));
h.XDisplayLabels = modelData.test_ev_labels;
h.YDisplayLabels = modelData.test_ev_labels;
%h.title = 'After Variance Partitioning';






%% 5.0 General Regression Function

% For testing purposes, we will stick with one brain area
% for brainArea=1:2 % run on each brain area
clear currentData currentLabel




% run shuffle first
if ~exist([exptdata.projectdir,'d500_nullDist_',modelData(m).brainArea(brainArea).modelName,'.mat'],'file') || exptdata.regenerateDists==1
    cprintf('_magenta',['\n\nGenerating null distribution for model ',num2str(m),'...\n'])
    ndist_thresholds=ephys_multipleRegression_shuffle(...
        currentData, ...
        currentLabel, ...
        modelData(m).brainArea(brainArea).test_evs,...
        modelData(m).brainArea(brainArea).test_ev_labels,...
        modelData(m).brainArea(brainArea).modelName,...
        modelData(m).brainArea(brainArea).regressType);
    modelData(m).brainArea(brainArea).ndist_thresholds=ndist_thresholds;
    save([exptdata.projectdir,'d500_modelData.mat'],'modelData')
else
    cprintf('_magenta',['\n\nLoading null distribution for model ',num2str(m),'...\n'])
    load([exptdata.projectdir,'d500_nullDist_',modelData(m).brainArea(1).modelName,'.mat'],'ndist_thresholds')
    modelData(m).brainArea(brainArea).ndist_thresholds=ndist_thresholds;
end

% now run real regression
if ~isfield(modelData(m).brainArea(brainArea),'modelOutput') || isempty(modelData(m).brainArea(brainArea).modelOutput) || exptdata.regenerateDists==1
    cprintf('_Magenta',['\nRegressing actual data for model ',num2str(m),'...\n'])
    modelData(m).brainArea(brainArea).modelOutput=ephys_multipleRegression_spikes(...
        currentData, ...
        currentLabel, ...
        modelData(m).brainArea(brainArea).test_evs,...
        modelData(m).brainArea(brainArea).test_ev_labels,...
        modelData(m).brainArea(brainArea).modelName,...
        modelData(m).brainArea(brainArea).ndist_thresholds,...
        modelData(m).brainArea(brainArea).regressType);
    save([exptdata.projectdir,'d500_modelData.mat'],'modelData')
end

clear currentData currentLabel

save([exptdata.projectdir,'d500_modelData.mat'],'modelData')
clear modelData
fprintf('/n/n/n/n/n/n/n')


return

% HERE IS THE BIT OF CODE TO ADJUST THE MULTIREG ANALYSIS %%
% Possible regressors to include:
% evdata.trialNumber
% evdata.animal
% evdata.stim1 (stimA=1; stimB=-1)
% evdata.stim2 (stimA=1; stimB=-1)
% evdata.session: session number
% evdata.actual: (1=repeat, -1=alternation)
% evdata.cor: (correct=1, incorrect=0)
% evdata.active (1) or passive (0) block
% evdata.reactiontime: (in ms) ** not super reliable **
% evdata.expect: expectation (1=expect repeat; 2=expect alternation; % 3=no expectation) ('truth')
% evdata.prep: prob repetition (what the animal expects, modelled)
% evdata.smooth_p_rep: smoothed
% evdata.correctRight = whether correct answer was on right (1) or left (-1)
% evdata.congruent_n1 (n minus 1) = correct answer is same for trial n and n-1
% evdata.congruent_n1 (n minus 2) = correct answer is same for trial n and n-2
% evdata.congruent_n3 (n minus 3) = correct answer is same for trial n and n-3
% evdata.stimIDs
% evdata.stim1_faceNonFace %% SHOULDN'T USE THIS FOR ANYTHING BUT FACE-SELECTIVITY
% evdata.stim2_faceNonFace %% SHOULDN'T USE THIS FOR ANYTHING BUT FACE-SELECTIVITY
% evdata.choseRepeat





