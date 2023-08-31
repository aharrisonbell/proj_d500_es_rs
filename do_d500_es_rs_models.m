% do_d500_es_rs_models(vargin)
% AHB, started July 24, 2020, based on do_d500_es_rs_multipleRegression
% This version is essentially a cleaned-up version of the above, with an
% effort to:
% a) prepare all figures for the paper
% b) follow all of Yinan's advice, including vectorisation
% BEHAV DATA STRUCTURE:
%  1) monkey number
%  2) session number
%  5) trial number
%  6) trial outcome (0=correct)
%  7) block number
%  8) condition number
%  9) reaction time
% 10) total number of neurons per session


%% TO DO (Sept 10, 2020)
% 1. create plots that compare the values from spm_orth 'normal', 'reverse', 'wiki'
% 2. create plots that compares the regression results from the three methods above
% 3. maybe consider annotating the figures (though not as important)
% 4. create spike density functions


clc; close all; clearvars
vers_ephys_es_rs='1.1; Sept 12, 2020';
% 1.0 - original version (July 20, 2018), added vectorisation
% 1.1 - updated figures, added functionality to compare residualisation methods
% 1.2 - updated to run on modern hardware and PC (HOMER_B550)

tic
ephys_analysis_defaults;
exptdata.lastModified=date;
warning('off','MATLAB:MKDIR:DirectoryExists');

% Modifiable defaults and directories (Change These)
exptdata.analysisName='D500_ES-RS_Study'; % used for savenames, figures, etc. (pick whatever you want; will be used for filenames)
projectDir=[exptdata.analysisdir,exptdata.analysisName,filesep]; mkdir(projectDir); %
figureDir=[projectDir,'figures',filesep]; mkdir(figureDir);
mkdir([figureDir,'matlabFigFiles']); mkdir([figureDir,'neuronPrintouts']);
diary([projectDir,lower(exptdata.analysisName),'_',exptdata.analysisName,'.txt']);

% Analysis Parameters
exptdata.xrange_psths=-250:500;     % window surrounding stimulus onset (in ms)
exptdata.reprocess=0;               % recreate trialised files even if they already exist (preprocess) (includes prepping megaMatrices)
regenerateDists=1;         % generate model distributions (0 if already exist and just want to produce figures)
regenerateNullDists=1;     % generate null distributions (0 if already exist and just want to produce figures) (shouldn't have to be done again unless model changes)
addpath(userpath);

%addpath(genpath('~/Documents/MATLAB/spm12'));
%addpath(genpath('C:\Users\Andrew H Bell\iCloudDrive\Documents\MATLAB\spm12'));

% =====================================================================================================================
clc
fprintf('<strong>*========================*</strong>\n')
fprintf('<strong>| do_d500_es_rs_models.m |</strong>\n')
fprintf('<strong>*========================*</strong>\n')
fprintf(['Version:              ',vers_ephys_es_rs,'\n'])
disp(['Study name:           ',exptdata.analysisName]);
disp(['Data location:        ',exptdata.datalocation]);
disp(['Processed NEV dir:    ',exptdata.processedDatadir]);
disp(['Project Directory:    ',projectDir]);
disp(['Figure Directory:     ',figureDir]);

% =====================================================================================================================
%% 1. Load Behavioural Data (Neuro Data loaded at a later stage)
fprintf('\n<strong>Loading behavioural data...</strong>')
load([projectDir,'DMS500_behavAnalysis.mat'],'evdata','behavData'); % contains BOTH monkeys behavioural data (this includes markers)
fprintf('done.\n')

%% 2. Set up variables
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

%% 3. Feature Engineering and Normalisation (Z-SCORING) - Normalise (most important) EVs (should be done BEFORE spm_orth)
% evdata.active            => don't normalize (1=active, -1=passive)
% evdata.actual            => don't normalize (1=repeat, -1=alternate)
% evdata.choseRepeat       => don't normalize (1=chose repeat, 0=did not
% evdata.cor               => don't normalize (1=correct, 0=not correct)
% evdata.expect            => don't normalize (1=expect repeat, 0=none, -1=expect alternate)
% evdata.stim1             => don't normalize (1=stim A, -1=stim B)
% evdata.stim2             => don't normalize (1=stim A, -1=stim B)
% evdata.stim1_faceNonFace => don't normalize (1=face, -1=nonface)
% evdata.stim2_faceNonFace => don't normalize (1=face, -1=nonface)

evdata.unsigned_pe      = abs(evdata.delta_p_rep);
evdata.norm_prep        = normalize(evdata.prep,'range',[-1 1]); % looks good
evdata.norm_unsigned_pe = normalize(evdata.delta_p_rep, 'range', [-1 1]);
evdata.norm_signed_pe   = zscore(evdata.unsigned_pe);
evdata.norm_trialNumber = normalize(evdata.trialNumber,'range',[-1 1]); % looks good

%%% 3.1 - Feature Engineering
% 3.1.1 - Trial number within block
tmp_tr_within = 1; tmp_sess=1; tmp_block=behavData(1,7); % initialise starting values
for tr = 1:size(behavData, 1)
    if behavData(tr, 2) ~= tmp_sess % new session?
        tmp_sess=behavData(tr, 2);
        tmp_block = behavData(tr, 7);
        tmp_tr_within = 1;
    elseif behavData(tr, 7) ~= tmp_block % new block?
        tmp_block = behavData(tr, 7);
        tmp_tr_within = 1;
    elseif behavData(tr, 7) == tmp_block % trial within block
        tmp_tr_within = tmp_tr_within+1;
    end
    evdata.trial_num_within_block(tr) = tmp_tr_within;
end
evdata.trial_num_within_block = evdata.trial_num_within_block';
evdata.norm_trial_num_within_block = normalize(evdata.trial_num_within_block, 'range', [-1 1] );

% 3.1.2 - Interaction Features
%%% Both of these features are set up such that a POSITIVE value means an
%%% expectation was MET, and a NEGATIVE VALUE 
evdata.norm_prep_x_expect = evdata.norm_prep .* evdata.expect;
evdata.actual_x_expect = evdata.actual.* evdata.expect;         

%%% 3.2 - Plot Regressors 
h=figure; set(gcf, 'Units', 'Normalized', 'OuterPosition',[0.1 0.1 0.9 0.8],'NumberTitle','Off','Name','Regressors after normalization');
subplot(2, 2, [1 2])
boxplot([evdata.prep evdata.norm_prep evdata.delta_p_rep evdata.unsigned_pe evdata.norm_signed_pe ...
    evdata.norm_unsigned_pe evdata.norm_prep_x_expect evdata.actual_x_expect, evdata.norm_trialNumber],...
    'notch','on','labels', {'p(rep)', 'z(p(rep))', 'predErr','unsignedPredErr','normPredErr',...
    'normUnPredErr','normPrepxExpect','ActualXExpect','TrialNum'})
title('regressors')
subplot(2, 2, 3)
histogram(evdata.actual_x_expect)
title('actual*expect')
subplot(2, 2, 4)
histogram(evdata.norm_prep_x_expect)
title('norm-prep*expect')
savefig(h,[figureDir,'d500_generalEDA_newFeatures.fig']);
jpgfigname=[figureDir,'d500_generalEDA_newFeatures.jpg'];
print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure

%%% July 24, 2020 - What this figure shows is that there is little point in
%%% normalizing prediction error because it is essentially ALREADY normalized.
%%% So, we can scrap norm_unsigned_pe and norm_signed_pe and just keep the
%%% delta_p and unsigned_pe. 
%%% July 24, 2020 - These look ok. It actually makes sense that there would
%%% be fewer values below ZERO (i.e., expectations that were NOT met)

plot_range = [1500 2000];
h=figure; set(gcf,'Units','Normalized','NumberTitle','Off','Name','Snapshot of New Regressors');  set(gcf,'Position',[0.1 0.1 0.9 0.8]);
subplot(4, 4, [1 2 3]); hold on
plot(zscore(evdata.trialNumber), 'k-', 'LineWidth',1)
plot(zscore(evdata.trial_num_within_block), 'b-', 'LineWidth',1.5)
plot(evdata.expect,'r-', 'LineWidth',2)
xlim(plot_range); ylabel('Normalised Value')
legend('trialNum','withinBlock','expect')
title('SAMPLE - Comparison of Trial Number within block and within session')

subplot(4, 4, 4); hold on
boxplot([zscore(evdata.trialNumber) zscore(evdata.trial_num_within_block) evdata.expect],...
    'notch','on','labels', {'trialNum', 'withinBlock','expect'})

subplot(4, 4, [5 6 7]); hold on
p1=plot(evdata.actual, 'k-', 'LineWidth', 2);
p2=plot(evdata.actual_x_expect*0.5, 'r-', 'LineWidth', 0.75);
p1.Color(4)=0.50; % p2.Color(4)=0.50;
xlim(plot_range); ylabel('Normalised Value')
legend('Actual','Actual*Expect (scaled)')

subplot(4, 4, 8); hold on
histogram(evdata.actual,'DisplayStyle','Stairs')
histogram(evdata.actual_x_expect)
%boxplot([evdata.actual evdata.actual_x_expect],...
% 'notch','on','labels', {'actual', 'actual*expect'})

subplot(4, 4, [9 10 11 13 14 15]); hold on
p1=plot(evdata.norm_prep_x_expect, 'k-', 'LineWidth', 2);
p2=plot(evdata.prep, 'r-', 'LineWidth', 2);
p3=plot(evdata.norm_prep, 'b-', 'LineWidth', 2);
p4=plot(evdata.actual_x_expect, 'm-', 'LineWidth', 1);
xlim(plot_range); ylabel('Normalised Value')
legend('Expect * z(p-rep)','prep','norm-prep','actual*expect')

subplot(4, 4, [12 16]); hold on
boxplot([evdata.norm_prep_x_expect evdata.prep evdata.norm_prep evdata.actual_x_expect], 'notch', 'on', 'labels', ...
    {'n(prep*expect)', 'p(rep)', 'n(p(rep))', 'actual*expect'})

savefig(h,[figureDir,'d500_generalEDA_newFeatures_snapshot.fig']);
jpgfigname=[figureDir,'d500_generalEDA_newFeatures_snapshot.jpg'];
print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure

%% 4. Model Definition
%%% July 24, 2020 - What EV's are we left with currently?
% Check correlation between explanatory variables 
% specify which variables to correlate here:

%%% Possible Regressors
%    evdata.active, ...                      % (1=active, -1=passive)
%    evdata.actual, ...                      % (1=repeat, -1=alternate)
%    evdata.choseRepeat, ...                 % (1=chose repeat, 0=did not
%    evdata.cor, ...                         % (1=correct, 0=not correct)
%    evdata.delta_p_rep, ...                 % continuous between -1 and 1
%    evdata.expect, ...                      % (1=expect repeat, 0=none, -1=expect alternate)
%    evdata.norm_prep, ...                   % mean-normalized to -1 and 1     
%    evdata.stim1, ...                       % (1=stim A, -1=stim B)
%    evdata.stim2,...                        % (1=stim A, -1=stim B)
%    evdata.stim1_faceNonFace, ...           % (1=face, -1=nonface)
%    evdata.stim2_faceNonFace, ...           % (1=face, -1=nonface)
%    evdata.unsigned_pe, ...                 % continuous between 0-1
%    evdata.norm_trial_num_within_block, ... % mean-normalized to -1 and 1
%    evdata.norm_prep_x_expect, ...          % continuous (-1 to 1) based on outcome POSITIVE=expected, NEGATIVE=unexpected
%    evdata.actual_x_expect, ...             % categorical based on outcome -1=unexpected, 0=no expect, 1=expected
%    evdata.norm_trialNumber, ...

%%% July 25, 2020 - Construct models:

modelsuffix='_resid_reverse';

%%%------------------------------------------------------------------
% 1) Examine Prediction Signals
modelData(1).modelName='Model1-ES-aXe+prXe+pr+pe+upe';
modelData(1).ev_suffix=modelsuffix;
modelData(1).evs=[...
        evdata.actual_x_expect, ...
        evdata.norm_prep_x_expect, ...
        evdata.norm_prep, ...
        evdata.delta_p_rep, ...
        evdata.unsigned_pe ...
        ];
modelData(1).ev_labels={...
    'actual*expect', 'p(rep)*expect', 'p(rep)','PE','uPE' ...
    };

%%%------------------------------------------------------------------
% 2) Replicate Repetition Suppression
modelData(2).modelName='Model2-RS-atv+atu+cr';
modelData(2).ev_suffix=modelsuffix;
modelData(2).evs=[...
        evdata.active, ...
        evdata.actual, ...
        evdata.choseRepeat, ...
        ];
modelData(2).ev_labels={...
    'active', 'actual', 'choseRepeat'...
    };

%%%------------------------------------------------------------------
% 3) Stimulus Feature Model
modelData(3).modelName='Model3-Stim-s1f+s2f';
modelData(3).ev_suffix=modelsuffix;
modelData(3).evs=[...
    evdata.stim1_faceNonFace, ...
    evdata.stim2_faceNonFace, ...
    ];
modelData(3).ev_labels={...
    'stim1face?', 'stim2face?'...
    };

%%%------------------------------------------------------------------
% 4) Stimulus Feature Model
modelData(4).modelName='Model4-Act+actXexp+PE+prep';
modelData(4).ev_suffix=modelsuffix;
modelData(4).evs=[...
    evdata.actual, ...
    evdata.actual_x_expect, ...
    evdata.delta_p_rep, ...                       % (1=stim A, -1=stim B)
    evdata.norm_prep,...                        % (1=stim A, -1=stim B)
    ];
modelData(4).ev_labels={...
    'actual', 'act*exp', 'PE', 'nPrep'...
    };

%% 5. Examine EVs, Correlation and Residualisation
for m = 1:numel(modelData)
    %%% Initial CorrPlot
    fprintf(['\n<strong>Examining EVs in </strong>',modelData(m).modelName,' ...']);
    figure; set(gcf,'Units','Normalized', 'NumberTitle', 'off','Name', ['Model ',num2str(m),' - Correlation between regressors before/after Residualisation']);  set(gcf,'Position',[0.1 0.1 0.8 0.8]);
    subplot(2,2,1)
    corrplot(modelData(m).evs, 'varnames', modelData(m).ev_labels)
    title({modelData(m).modelName, 'Before Residualisation'}, 'FontSize', 14)
    axis square
    
    %%% Residualisation
    subplot(2,2,2) % forward
    modelData(m).evs_resid = ephys_my_spm_orth(modelData(m).evs);
    corrplot(modelData(m).evs_resid, 'varnames', modelData(m).ev_labels)
    title({modelData(m).modelName, 'Residualisation (forward)'}, 'FontSize', 14)
    
    subplot(2,2,3) % reverse
    modelData(m).evs_resid_reverse = ephys_my_spm_orth_reverse(modelData(m).evs);
    corrplot(modelData(m).evs_resid_reverse, 'varnames', modelData(m).ev_labels)
    title({modelData(m).modelName, 'Residualisation (reverse)'}, 'FontSize', 14)
    
    subplot(2,2,4) % wiki
    modelData(m).evs_resid_wiki = ephys_my_gram_schmidt(modelData(m).evs);
    corrplot(modelData(m).evs_resid_wiki, 'varnames', modelData(m).ev_labels)
    title({modelData(m).modelName, 'Residualisation (wiki algorithm)'}, 'FontSize', 14)
    
    savefig([figureDir,'d500_',modelData(m).modelName,'EDA_corrplot_before+after.fig'])
    jpgfigname=[figureDir,'d500_',modelData(m).modelName,'EDA_corrplot_before+after.jpg'];
    print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure

    
%     
%     %%% Examine effect of Residualisation on individual regressors
%     tmp_numEV = size(modelData(m).evs, 2);
%     numRows = ceil(tmp_numEV/3);
%     figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]);
%     for e = 1:size(modelData(m).evs, 2)
%         ax1 = subplot(numRows, 3, e); hold on
%         scatter(ax1,modelData(m).evs(:,e), modelData(m).evs_resid(:,e))
%         scatter(ax1,modelData(m).evs(:,e), modelData(m).evs_resid_wiki(:,e),'ro')
%         lsline(ax1)
%         axis square
%         title([modelData(m).modelName, modelData(m).ev_labels(e)], 'FontSize', 14, 'FontWeight', 'Bold')
%         xlabel('Original', 'FontSize', 12); ylabel('Residualised', 'FontSize', 12);
%     end
%     savefig([figureDir,'d500_',modelData(m).modelName,'_EDA_ev_corrBeforeAfter.fig'])
%     jpgfigname=[figureDir,'d500_',modelData(m).modelName,'_EDA_ev_corrBeforeAfter.jpg'];
%     print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
    
    %%% Within Model EV Correlation
    figure; set(gcf,'Units','Normalized','NumberTitle', 'Off', 'Name', ['Model #',num2str(m),' - Correlations between EVs']);  set(gcf,'Position',[0.1 0.1 0.8 0.8]);
    [modelData(m).corr_coeffs, modelData(m).corr_pvalues] = corrcoef(modelData(m).evs);
    [modelData(m).corr_coeffs_resid, modelData(m).corr_pvalues_resid] = corrcoef(modelData(m).evs_resid);
    [modelData(m).corr_coeffs_resid_reverse, modelData(m).corr_pvalues_resid_reverse] = corrcoef(modelData(m).evs_resid_reverse);
    [modelData(m).corr_coeffs_resid_wiki, modelData(m).corr_pvalues_resid_wiki] = corrcoef(modelData(m).evs_resid_wiki);

    subplot(2, 2, 1)
    heatmap(modelData(m).ev_labels, modelData(m).ev_labels, modelData(m).corr_coeffs, ...
        'Title',[modelData(m).modelName, ' - EV Correlation Matrix (BEFORE)'],...
        'ColorLimits',[-0.15 1.00]);  colormap(jet);
    subplot(2, 2, 2)
    heatmap(modelData(m).ev_labels, modelData(m).ev_labels, modelData(m).corr_coeffs_resid, ...
        'Title','EV Correlation Matrix (AFTER RESIDUALISATION)',...
        'ColorLimits',[-0.15 1.00]);  colormap(jet);
    subplot(2, 2, 3)
    heatmap(modelData(m).ev_labels, modelData(m).ev_labels, modelData(m).corr_coeffs_resid_reverse, ...
        'Title','EV Correlation Matrix (AFTER WIKI RESIDUALISATION)',...
        'ColorLimits',[-0.15 1.00]);  colormap(jet);
    subplot(2, 2, 4)
    heatmap(modelData(m).ev_labels, modelData(m).ev_labels, modelData(m).corr_coeffs_resid_wiki, ...
        'Title','EV Correlation Matrix (AFTER WIKI RESIDUALISATION)',...
        'ColorLimits',[-0.15 1.00]);  colormap(jet);
    
    savefig([figureDir,'d500_',modelData(m).modelName,'_EDA_ev_heatmap.fig'])
    jpgfigname=[figureDir,'d500_',modelData(m).modelName,'_EDA_ev_heatmap.jpg'];
    print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
    fprintf('done.');
end
disp(' ')
%%% Notes about how to evaluate the correlation strength (statistical significance)
%%% between EVs:
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

%% 6. Perform Regression
% Hardcode some variables to allow parallel processing
regressType = 'linear_vect';
all_trialNumbers = evdata.trialNumber;
all_sessions = evdata.session;
timeWindow = [timeEpochs500(epochNum, 1) timeEpochs500(epochNum, 2)];

for ba=1:2 % run on each brain area
    
    clear currentData currentLabel
    if ba==1
        temp1=load([projectDir,'Vortex_',exptdata.analysisName,'_megaMatrix_synced.mat'],'V4_megaMatrix500NoScrubs');
        temp2=load([projectDir,'Vulcan_',exptdata.analysisName,'_megaMatrix_synced.mat'],'V4_megaMatrix500NoScrubs');
        currentData=[temp1.V4_megaMatrix500NoScrubs; temp2.V4_megaMatrix500NoScrubs];
        currentLabel='V4';
        clear temp*
    else
        temp1=load([projectDir,'Vortex_',exptdata.analysisName,'_megaMatrix_synced.mat'],'TE_megaMatrix500NoScrubs');
        temp2=load([projectDir,'Vulcan_',exptdata.analysisName,'_megaMatrix_synced.mat'],'TE_megaMatrix500NoScrubs');
        currentData=[temp1.TE_megaMatrix500NoScrubs; temp2.TE_megaMatrix500NoScrubs];
        currentLabel='TE';
        clear temp*
    end
    disp(['<strong>Analysing ',currentLabel,'... </strong>'])
   
    for m=1:numel(modelData) % scroll through each model
        modelData(m).regressType = regressType;
        
        %% 6.1 Generate Null Distribution First
        if ~exist([projectDir,'d500_nullDist_',modelData(m).modelName,modelData(m).ev_suffix,'_',currentLabel,'.mat'],'file') || regenerateNullDists==1
            fprintf('<strong>..Generating null distribution for: </strong>\n')
            ndist_thresholds=ephys_multipleRegression_shuffle(...
                projectDir, ...
                figureDir, ...
                timeWindow, ...
                sampleRate, ...
                currentData, ...
                currentLabel, ...
                normalize(modelData(m).(['evs', modelData(m).ev_suffix])), ...
                modelData(m).ev_labels,...
                [modelData(m).modelName,modelData(m).ev_suffix],...
                modelData(m).regressType, ...
                all_trialNumbers, ...
                all_sessions);
            modelData(m).brainArea(ba).ndist_thresholds=ndist_thresholds; %#ok<*SAGROW>
        else
            fprintf('\n<strong>..Loading null distribution: </strong>\n')
            load([projectDir,'d500_nullDist_',modelData(m).modelName,modelData(m).ev_suffix,'_',currentLabel,'.mat'],'ndist_thresholds');
            modelData(m).brainArea(ba).ndist_thresholds=ndist_thresholds;
        end
        
        %% 6.2 Regression
        if ~exist([projectDir,'d500_',modelData(m).modelName,modelData(m).ev_suffix,'.mat'],'file') || regenerateDists==1
            fprintf('\n<strong>..Regressing actual data for: </strong>\n')
            modelData(m).brainArea(ba).modelOutput=ephys_multipleRegression_spikes(...
                projectDir, ...
                figureDir, ...
                timeWindow, ...
                sampleRate, ...
                currentData, ...
                currentLabel, ...
                normalize(modelData(m).(['evs', modelData(m).ev_suffix])), ...
                modelData(m).ev_labels,...
                [modelData(m).modelName,modelData(m).ev_suffix],...
                ndist_thresholds,...
                modelData(m).regressType, ...
                all_trialNumbers, ...
                all_sessions);
            model=  modelData(m).brainArea(ba).modelOutput;
            save([projectDir,'d500_',modelData(m).modelName,modelData(m).ev_suffix,'_',currentLabel,'.mat'],'model','ndist_thresholds')
            clear model
        else
            fprintf(['\n<strong>..Loading printouts for </strong>',num2str(m),'...\n'])
            openfig([figureDir,'d500_modelPrintout_wiki',modelData(m).modelName,modelData(m).ev_suffix,'_', currentLabel,'.fig']);
        end
        
%         %% 6.3 Multivariate Approach (Figure 4B in Bell, Summerfield, et al., 2016)
%         % The goal here is to:
%         % 1) randomly split the data into two halves
%         % 2) repeat the regression analysis
%         % 3) correlate the Beta coefficients across the two halves for all
%         % time points:
%         
%         if ~exist([projectDir,'d500_',modelData(m).modelName,'_randomSplit.mat'],'file') || exptdata.regenerateDists==1
%             fprintf('\n..<strong>Conducting multivariate analysis on: </strong>\n')
%             modelData(m).brainArea(ba).modelOutput=ephys_multipleRegression_split(...
%                 projectDir, ...
%                 figureDir, ...
%                 timeWindow, ...
%                 sampleRate, ...
%                 currentData, ...
%                 currentLabel, ...
%                 modelData(m).evs_resid,... % modelData(m).evs_resid
%                 modelData(m).ev_labels,...
%                 modelData(m).modelName,...
%                 ndist_thresholds,...
%                 modelData(m).regressType, ...
%                 all_trialNumbers, ...
%                 all_sessions, ...
%                 2);
%             model=  modelData(m).brainArea(ba).modelOutput;
%             save([projectDir,'d500_',modelData(m).modelName,'_',currentLabel,'_randomSplit.mat'],'model','ndist_thresholds');
%             clear model
%         else
%             cprintf('_magenta',['\n\nLoading printouts for model ',num2str(m),'...\n'])
%             openfig([figureDir,'d500_',test_modelName,'_',brainLabel,'_splitEV.jpg']);
%             openfig([figureDir,'d500_',test_modelName,'_',brainLabel,'_corrMatrix.fig']);
%         end
    end
    clear currentData currentLabel
end
save([projectDir,'d500_modelData.mat'],'modelData')
clear modelData

% modelData is very very large (10x largest spikedata matrix). Should fix
% this up so that it isn't so large. 
toc
fprintf('\n\n')


%% 7. Plot Spike Density Functions
do_d500_es_rs_makePrettyFigures;
% Want the following figures:
% Expected vs. Unexpected vs. Neutral Repeat vs. Alternation
% for individual example neuron + population?


