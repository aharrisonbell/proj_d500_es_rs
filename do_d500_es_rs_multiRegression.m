% do_d500_es_rs_multiRegression(vargin)
% AHB, started May 20, 2020
clc; close all
vers_ephys_es_rs='1.0; May 20, 2020';
% 1.0 - original version (May 20, 2018)
%
tic
global exptdata timeEpochs500 epochNum epochNames evdata sampleRate
ephys_analysis_defaults;
exptdata.lastModified=date;
warning('off','MATLAB:MKDIR:DirectoryExists');

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
exptdata.xrange_psths=-250:500;     % window surrounding stimulus onset (in ms)
exptdata.reprocess=0;               % recreate trialised files even if they already exist (preprocess) (includes prepping megaMatrices)
exptdata.regenerateDists=1;
exptdata.regenerateNullDists=1;
exptdata.reviewNeurons=0;           % determines whether to generate figure for each neuron
exptdata.behav_reprocess=0;         % set to one if you want to regenerate output files for existing files
exptdata.paradigm={'DMS500'};
addpath(userpath);
exptdata.analysisName='D500_ES-RS_Study';
%addpath(genpath('~/Documents/MATLAB/spm12'));
%addpath(genpath('C:\Users\Andrew H Bell\iCloudDrive\Documents\MATLAB\spm12'));

% =====================================================================================================================
exptdata.projectdir=[exptdata.analysisdir,exptdata.analysisName,filesep];
clc
fprintf('<strong>*=================================*</strong>\n')
fprintf('<strong>| do_d500_es_rs_multiRegression.m |</strong>\n')
fprintf('<strong>*=================================*</strong>\n')
fprintf(['Version: ',vers_ephys_es_rs,'\n'])
disp(['Data location:    ',exptdata.datalocation]);
disp(['Processed NEV dir: ',exptdata.processedDatadir]);
disp(['Study name:       ',exptdata.analysisName]);
disp(['Project Directory: ',exptdata.projectdir]);

% =====================================================================================================================
% 1. Load Behavioural and Neurophysiological Data
fprintf('\n<strong>Loading megamatrices and behavioural data...</strong>')
load([exptdata.projectdir,'DMS500_behavAnalysis.mat'],'behavData','evdata'); % contains BOTH monkeys behavioural data (this includes markers)

temp1=load([exptdata.projectdir,'Vortex_',exptdata.analysisName,'_megaMatrix_synced.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
temp2=load([exptdata.projectdir,'Vulcan_',exptdata.analysisName,'_megaMatrix_synced.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
V4_megaMatrix500NoScrubs=[temp1.V4_megaMatrix500NoScrubs; temp2.V4_megaMatrix500NoScrubs];
TE_megaMatrix500NoScrubs=[temp1.TE_megaMatrix500NoScrubs; temp2.TE_megaMatrix500NoScrubs];
clear temp*
sessions = unique(evdata.session);
fprintf('done.\n')

% 2. Set up variables
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


%% 3. FEATURE NORMALISATION (Z-SCORING) - Normalise (most important) EVs (should be done BEFORE spm_orth)
evdata.norm_active        =zscore(evdata.active);
evdata.norm_actual        =zscore(evdata.actual);
evdata.norm_choseRepeat   =zscore(evdata.choseRepeat);
evdata.norm_cor           =zscore(evdata.cor);
evdata.norm_delta_p_rep   =zscore(evdata.delta_p_rep);
evdata.norm_expect        =zscore(evdata.expect);
evdata.norm_prep          =zscore(evdata.prep);
evdata.norm_stim1         =zscore(evdata.stim1);
evdata.norm_stim2         =zscore(evdata.stim2);
evdata.norm_stim1_faceNonFace =zscore(evdata.stim1_faceNonFace);
evdata.norm_stim2_faceNonFace =zscore(evdata.stim2_faceNonFace);

%% Part 1 - Correlation Matrix EVs
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
    ];

tmp_norm_ev=[...
    evdata.norm_active, ...
    evdata.norm_actual, ...
    evdata.norm_choseRepeat, ...
    evdata.norm_cor, ...
    evdata.norm_delta_p_rep, ...
    evdata.norm_expect, ...
    evdata.norm_prep, ...
    evdata.norm_stim1, ...
    evdata.norm_stim2,...
    evdata.norm_stim1_faceNonFace, ...
    evdata.norm_stim2_faceNonFace, ...
    ];
paperData.corrcoeff.evLabels={'active','actual','choseRepeat','cor','delta-p-rep','expect',...
    'prep','stim1','stim2','stim1-face?','stim2-face?'};
[paperData.corrcoeff.coeffs, paperData.corrcoeff.pvalues]=corrcoef(tmp_ev);
[paperData.corrcoeff.norm_coeffs, paperData.corrcoeff.norm_pvalues]=corrcoef(tmp_norm_ev);
figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial');
clf;
subplot(2, 2, 1)
heatmap(paperData.corrcoeff.evLabels, paperData.corrcoeff.evLabels, paperData.corrcoeff.coeffs, ...
    'Title','EV Correlation Matrix (BEFORE)',...
    'ColorLimits',[-0.15 1.00]);  colormap(jet); 

subplot(2, 2, 2)
heatmap(paperData.corrcoeff.evLabels, paperData.corrcoeff.evLabels, paperData.corrcoeff.norm_coeffs, ...
    'Title','EV Correlation Matrix (BEFORE) (norm''d EVs)',...
    'ColorLimits',[-0.15 1.00]);  colormap(jet); 

subplot(2, 2, 3)
heatmap(paperData.corrcoeff.evLabels, paperData.corrcoeff.evLabels, spm_orth(paperData.corrcoeff.coeffs), ...
    'Title','EV Correlation Matrix (AFTER)',...
    'ColorLimits',[-0.15 1.00]);  colormap(jet); 

subplot(2, 2, 4)
heatmap(paperData.corrcoeff.evLabels, paperData.corrcoeff.evLabels, spm_orth(paperData.corrcoeff.norm_coeffs), ...
    'Title','EV Correlation Matrix (AFTER) (norm''d EVs)',...
    'ColorLimits',[-0.15 1.00]);  colormap(jet); 

savefig([exptdata.figuredir500,'d500_EVcorrelation_heatmap.fig'])
jpgfigname=[exptdata.figuredir500,'d500_EVcorrelation_heatmap.jpg'];
print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure

figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial');
subplot(1,2,1)
imagesc(paperData.corrcoeff.norm_coeffs); axis square
set(gca,'XTick',1:14,'XTickLabels',paperData.corrcoeff.evLabels); xtickangle(45);
set(gca,'YTick',1:14,'YTickLabels',paperData.corrcoeff.evLabels)
title('EVs before Orthogonalisation (normed)')
clim([-0.15 1]); colorbar; colormap(jet)
subplot(1,2,2)
imagesc(spm_orth(paperData.corrcoeff.norm_coeffs)); axis square
set(gca,'XTick',1:14,'XTickLabels',paperData.corrcoeff.evLabels); xtickangle(45);
set(gca,'YTick',1:14,'YTickLabels',paperData.corrcoeff.evLabels)
clim([-0.15 1]); colorbar; colormap(jet)
title('After Orthogonalisation')
savefig([exptdata.figuredir500,'d500_EVcorrelationBOTH.fig'])

% Ameya Deoras (2020). Customizable Heat Maps (https://www.mathworks.com/matlabcentral/fileexchange/24253-customizable-heat-maps), 
% MATLAB Central File Exchange. Retrieved April 16, 2020. [this appears to be a hack of an existing matlab function!]

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



%% FEATURE ENGINEERING

evdata.norm_unsigned_pred_err = abs(evdata.norm_delta_p_rep);


%% REGRESSION
for brainArea=1:1 % run on each brain area
    clear currentData currentLabel
    if brainArea==1
        currentData=V4_megaMatrix500NoScrubs;
        currentLabel='V4';
    else
        currentData=TE_megaMatrix500NoScrubs;
        currentLabel='TE';
    end
    disp(['Analysing ',currentLabel,'...'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Model:
    % July 4, 2020 - 
    % Here are the possible regressors:
    %      trialNumber: [6139×1 double]
    %      animal: [6139×1 double]
    %      actual: [6139×1 double]
    %      session: [6139×1 double]
    %      cor: [6139×1 double]
    %      active: [6139×1 double]
    %      stim1: [6139×1 double]
    %      stim2: [6139×1 double]
    %      reactiontime: [6139×1 double]
    %      expect: [6139×1 double]
    %      correctRight: [6139×1 double]
    %      stimIDs: [6139×4 double]
    %      stim1_faceNonFace: [6139×1 double]
    %      stim2_faceNonFace: [6139×1 double]
    %      choseRepeat: [6139×1 double]
    %      congruent_n1: [6139×1 double]
    %      congruent_n2: [6139×1 double]
    %      congruent_n3: [6139×1 double]
    %      prevtrial_n1: [6139×1 double]
    %      prevtrial_n2: [6139×1 double]
    %      prevtrial_n3: [6139×1 double]
    %      prep: [6139×1 double]
    %      smooth_p_rep: [1×6139 double]
    %      delta_p_rep: [6139×1 double]
    %
    % The "important" ones are:
    % - active (1, -1) (Active vs. Passive)
    % - actual (1, -1) (Active vs. Passive)
    % - cor (1, -1) (Active vs. Passive)
    % - expect (1, -1) (Active vs. Passive)
    % - prep (1, -1) (Active vs. Passive)
    % - delta_p_rep (1, -1) (Active vs. Passive)
    % ...????
    
    modelDate = '19-Jul-2020';
    
    % BASIC MODEL LOOKING AT EXPECTATION
    %=====================================================================
    modelData(1).brainArea(brainArea).modelName=['model1-',currentLabel,'-',modelDate,'-expectationModel'];
    modelData(1).brainArea(brainArea).test_evs=[...
        zscore(evdata.norm_active .* evdata.norm_actual), ...
        evdata.norm_prep, ...
        evdata.norm_delta_p_rep, ...
        evdata.norm_unsigned_pred_err ...
        ];
    modelData(1).brainArea(brainArea).test_ev_labels={'actual*active','p-rep','deltaP','unsigned-PE'};
    modelData(1).brainArea(brainArea).regressType='linear_resid';
    %=====================================================================
    % EXPANDED MODEL
    % This model has a rank deficiency issue with "ACTUAL"
    % - model5 will try same things in different order
    %=====================================================================
    modelData(2).brainArea(brainArea).modelName=['model2-',currentLabel,'-',modelDate,'-expandedModel'];
    modelData(2).brainArea(brainArea).test_evs=[...
        zscore(evdata.norm_active .* evdata.norm_actual), ...
        evdata.norm_actual, ...
        evdata.norm_prep, ...
        evdata.norm_delta_p_rep, ...
        evdata.norm_stim1_faceNonFace, ...
        evdata.norm_stim2_faceNonFace ...
        ];
    modelData(2).brainArea(brainArea).test_ev_labels={'actual*active','actual','p-rep','deltaP','stim1face','stim2face'};
    modelData(2).brainArea(brainArea).regressType='linear_resid';
    %=====================================================================
    % REPLICATE REP SUPPRESSION      
    %=====================================================================
    modelData(3).brainArea(brainArea).modelName=['model3-',currentLabel,'-',modelDate,'-repetitionSuppression'];
    modelData(3).brainArea(brainArea).test_evs=[...
        zscore(evdata.norm_active .* evdata.norm_actual), ...
        zscore(evdata.norm_active), ...
        zscore(evdata.norm_actual), ...
        zscore(evdata.norm_choseRepeat), ...
        ];
    modelData(3).brainArea(brainArea).test_ev_labels={'actual*active','active','actual','choseRepeat'};
    modelData(3).brainArea(brainArea).regressType='linear_resid';
    %=====================================================================
    % STIM FEATURE MODEL
    %=====================================================================
    modelData(4).brainArea(brainArea).modelName=['model4-',currentLabel,'-',modelDate,'-stimFeatures'];
    modelData(4).brainArea(brainArea).test_evs=[...
        evdata.trialNumber, ...
        evdata.norm_stim1, ...
        evdata.norm_stim2, ...
        evdata.norm_stim1_faceNonFace, ...
        evdata.norm_stim2_faceNonFace ...
        ];
    modelData(4).brainArea(brainArea).test_ev_labels={'trialNumber','stim1', 'stim2', 'stim1face','stim2face'};
    modelData(4).brainArea(brainArea).regressType='linear_resid';
    %=====================================================================

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    for m=1:numel(modelData) % scroll through each model
        % run shuffle first
        if ~exist([exptdata.projectdir,'d500_nullDist_',modelData(m).brainArea(brainArea).modelName,'.mat'],'file') || exptdata.regenerateNullDists==1
            cprintf('_magenta',['\n\nGenerating null distribution for model ',num2str(m),'...\n'])
            ndist_thresholds=ephys_multipleRegression_shuffle(...
                currentData, ...
                currentLabel, ...        
                modelData(m).brainArea(brainArea).test_evs,...
                modelData(m).brainArea(brainArea).test_ev_labels,...
                modelData(m).brainArea(brainArea).modelName,...
                modelData(m).brainArea(brainArea).regressType);
            modelData(m).brainArea(brainArea).ndist_thresholds=ndist_thresholds; %#ok<*SAGROW>
        else
            cprintf('_magenta',['\n\nLoading null distribution for model ',num2str(m),'...\n'])
            load([exptdata.projectdir,'d500_nullDist_',modelData(m).brainArea(brainArea).modelName,'.mat'],'ndist_thresholds');
            modelData(m).brainArea(brainArea).ndist_thresholds=ndist_thresholds;
        end
        
        % now run real regression
        if ~exist([exptdata.projectdir,'d500_',modelData(m).brainArea(brainArea).modelName,'.mat'],'file') || exptdata.regenerateDists==1
            cprintf('_Magenta',['\nRegressing actual data for model ',num2str(m),'...\n'])
            modelData(m).brainArea(brainArea).modelOutput=ephys_multipleRegression_spikes(...
                currentData, ...
                currentLabel, ...
                modelData(m).brainArea(brainArea).test_evs,...
                modelData(m).brainArea(brainArea).test_ev_labels,...
                modelData(m).brainArea(brainArea).modelName,...
                modelData(m).brainArea(brainArea).ndist_thresholds,...
                modelData(m).brainArea(brainArea).regressType);
            model=  modelData(m).brainArea(brainArea).modelOutput;
            save([exptdata.projectdir,'d500_',modelData(m).brainArea(brainArea).modelName,'.mat'],'model','ndist_thresholds')
            clear model
        else
            cprintf('_magenta',['\n\nLoading printouts for model ',num2str(m),'...\n'])
            openfig([exptdata.figuredir500,'d500_modelPrintout_',modelData(m).brainArea(brainArea).modelName,'.fig']);
        end
        
        %% Multivariate Approach (Figure 4B in Bell, Summerfield, et al., 2016)
        % The goal here is to:
        % 1) randomly split the data into two halves
        % 2) repeat the regression analysis
        % 3) correlate the Beta coefficients across the two halves for all
        % time points:
%          
%         if ~exist([exptdata.projectdir,'d500_',modelData(m).brainArea(brainArea).modelName,'_randomSplit.mat'],'file') || exptdata.regenerateDists==1
%             cprintf('_Magenta',['\nRegressing actual data for model ',num2str(m),'...\n'])
%             modelData(m).brainArea(brainArea).modelOutput=ephys_multipleRegression_split(...
%                 currentData, ...
%                 currentLabel, ...
%                 modelData(m).brainArea(brainArea).test_evs,...
%                 modelData(m).brainArea(brainArea).test_ev_labels,...
%                 modelData(m).brainArea(brainArea).modelName,...
%                 modelData(m).brainArea(brainArea).ndist_thresholds,...
%                 modelData(m).brainArea(brainArea).regressType, ...
%                 2);
%             model=  modelData(m).brainArea(brainArea).modelOutput;
%             save([exptdata.projectdir,'d500_',modelData(m).brainArea(brainArea).modelName,'_randomSplit.mat'],'model','ndist_thresholds');
%             clear model
%         else
%             cprintf('_magenta',['\n\nLoading printouts for model ',num2str(m),'...\n'])
%             openfig([exptdata.figuredir500,'d500_modelPrintoutC_',modelData(m).brainArea(brainArea).modelName,'.fig']);
%         end
%         
        
        
        
        
    end
    clear currentData currentLabel
end
save([exptdata.projectdir,'d500_modelData.mat'],'modelData')
clear modelData
% modelData is very very large (10x largest spikedata matrix). Should fix
% this up so that it isn't so large. 
toc
fprintf('\n\n\n\n\n\n\n')
