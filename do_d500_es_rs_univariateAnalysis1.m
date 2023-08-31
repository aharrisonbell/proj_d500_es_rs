% do_es_rs_univariateAnalysis1(vargin)
% by CS, AHB, started Nov 6, 2018
clc
vers_ephys_es_rs='1.6; May 19, 2020';
% 1.0 - original version (May 3, 2018)
% 1.1 - added model variables (Nov 6, 2018)
% 1.2 - added model variables (Dec 17, 2018)
% 1.3 - quality of life additions (Jan 15, 2019)
% 1.4 - new analysis, including stim ID, defaults to load both monkeys (Feb 25, 2019)
% 1.5 - updating (April 7, 2020)
% 1.6 - added possibilty of spm_orthogonalisation (May 19, 2020)
% 1.7 - ????? split model output to separate mat files.

%%
 
global exptdata timeEpochs500 epochNum epochNames evdata currentSpikeData currentBrainLabel sampleRate
load([exptdata.analysisdir,filesep,exptdata.analysisName,'_exptdata.mat'],'exptdata');
addpath(userpath);
exptdata.analysisName='D500_ES-RS_Study';
exptdata.regenerateDists=0;

addpath(genpath('C:\Users\Andrew H Bell\iCloudDrive\Documents\MATLAB\spm12'));


% =====================================================================================================================
% megamatrices AND behaviouralAnalysis must be complete before running this program
% What this program will do:

% Note - from Feb 23, 2018
% 1. need to select visual neurons only (?).
% 2. need to add average SPDEN functions of visual neurons, possibly
%    subselecting for stimulus selectivity
%    Check Vogels' work to see if they subselect
% 3. still need to debug multiple regression
% 4. when creating figures, normalise 2ndPres to 1stPres (of same stimuli)


% Need to add:
% regressor for active/passive
% maybe filter out just visually responsive neurons
% decoding analysis

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
%
%
%
%
%
... lots of spike density stuff

% =====================================================================================================================
exptdata.projectdir=[exptdata.analysisdir,exptdata.analysisName,filesep];
clc
cprintf('yellow','*=====================================*\n')
cprintf('yellow','| do_d500_es_rs_univariateAnalysis1.m |\n')
cprintf('yellow','*=====================================*\n')
cprintf(['Version: ',vers_ephys_es_rs,'\n'])
disp(['Data location:    ',exptdata.datalocation]);
disp(['Processed NEV dir: ',exptdata.processedDatadir]);
disp(['Study name:       ',exptdata.analysisName]);
disp(['Project Directory: ',exptdata.projectdir]);

% =====================================================================================================================
%% 1. Load Behavioural and Neurophysiological Data
fprintf('\nLoading megamatrices and behavioural data...')
load([exptdata.projectdir,'DMS500_behavAnalysis.mat'],'behavData','evdata'); % contains BOTH monkeys behavioural data (this includes markers)

temp1=load([exptdata.projectdir,'Vortex_',exptdata.analysisName,'_megaMatrix_synced.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
temp2=load([exptdata.projectdir,'Vulcan_',exptdata.analysisName,'_megaMatrix_synced.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
V4_megaMatrix500NoScrubs=[temp1.V4_megaMatrix500NoScrubs; temp2.V4_megaMatrix500NoScrubs];
TE_megaMatrix500NoScrubs=[temp1.TE_megaMatrix500NoScrubs; temp2.TE_megaMatrix500NoScrubs];
clear temp*
sessions = unique(evdata.session);
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

paperData=struct('rq1_demographics',[],'rq2_BetaStuff',struct('session',[]));

%% 3. Neuron "demographics"
% Identify visually responsive neurons as per Bell et al., 2016
% "Overall, 174 neurons (69%) responded to visual stimuli (comparing baseline versus post-stimulus activity in low-noise trials; p < 0.05)"
% 1) how many visual neurons (baseline vs. post-stim activity in neutral trials?)
% 2) face-responsive?
% 3) excitatory, suppressed, both?
% 4) neuron type?
% *FIGURE: show examples

cprintf('Red*','* Q1: Demographics\n')
for brainArea=1:2 % run on each brain area
    if brainArea==1
        currentSpikeData=V4_megaMatrix500NoScrubs;
        currentBrainLabel='V4';
    else
        currentSpikeData=TE_megaMatrix500NoScrubs;
        currentBrainLabel='TE';
    end

    [~,tmp_uniqueRows,~]=unique([currentSpikeData(:,2) currentSpikeData(:,3)],'rows'); % identify unique row with neuron classification info
    tmp_neuronDemos=currentSpikeData(tmp_uniqueRows,[1:4 10 20:22]);
   
    % APRIL 15, 2020 - CURRENTLY USING AUTOMATIC CLASSIFICATIONS
    paperData.rq1_demographics(brainArea).neurontype(1)=size(tmp_neuronDemos,1); % total number of neurons
    paperData.rq1_demographics(brainArea).neurontype(2)=length(find(tmp_neuronDemos(:,6)==1)); % visually "enhanced"
    paperData.rq1_demographics(brainArea).neurontype(3)=length(find(tmp_neuronDemos(:,6)==-1)); % visually suppressed
    paperData.rq1_demographics(brainArea).neurontype(4)=length(find(tmp_neuronDemos(:,6)==0)); % non-responsive
    cprintf('Green',[' Total neurons (',currentBrainLabel,'): ',num2str(paperData.rq1_demographics(brainArea).neurontype(1)),'\n'])
    cprintf('Green',['    Excitatory (',currentBrainLabel,'): ',num2str(paperData.rq1_demographics(brainArea).neurontype(2)),'\n'])
    cprintf('Green',['    Suppressed (',currentBrainLabel,'): ',num2str(paperData.rq1_demographics(brainArea).neurontype(3)),'\n'])
    cprintf('Green',['Non-Responsive (',currentBrainLabel,'): ',num2str(paperData.rq1_demographics(brainArea).neurontype(4)),'\n\n'])
end
clear tmp*

%% 4. Correlation Matrix EVs
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

figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial');
[paperData.corrcoeff.coeffsOrth, paperData.corrcoeff.pvaluesOrth]=corrcoef(spm_orth(tmp_ev));
clabel = arrayfun(@(x){sprintf('%0.2f',x)}, paperData.corrcoeff.pvaluesOrth); 
heatmaps(paperData.corrcoeff.coeffsOrth, paperData.corrcoeff.evLabels, paperData.corrcoeff.evLabels, clabel,'ShowAllTicks',1,'Colormap','money');
title('EV Correlation Matrix (Orth)', 'FontSize', 14); % set title


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

%% 5. General Regression Function
% April 16, 2020 - I would like to set something up here that refers to a
% separate function.m file with the following abilities:
% - variable number of EVs
% - significance testing (as per Bell et al., 2016)
% - plot data
% - specify type of regression (stepwise, linear, logistic, etc.)
% - output pvalues, deviance, and residuals along with betas
% 
% if ~exist([exptdata.projectdir,'d500_modelData.mat'],'file')
%     disp('modelData.mat doesn''t exist. Creating new structure.')
%     modelData=struct('brainArea',struct('modelName',[],...
%         'test_evs',[],...
%         'test_ev_labels',[],...
%         'ndist_thresholds',[],...
%         'modelOutput',[]...
%         ));
% else
%     fprintf('Loading preexisting modelData.mat...')
%     load([exptdata.projectdir,'d500_modelData.mat'])
%     fprintf('done.\n')
% end

for brainArea=1:2 % run on each brain area
    clear currentData currentLabel
    if brainArea==1
        currentData=V4_megaMatrix500NoScrubs;
        currentLabel='V4';
    else
        currentData=TE_megaMatrix500NoScrubs;
        currentLabel='TE';
    end
    disp(['Analysing ',currentLabel,'...'])
    
    %evdata.active, ...
    %evdata.actual, ...
    %evdata.choseRepeat, ...
    %evdata.cor, ...
    %evdata.delta_p_rep, ...
    %evdata.expect, ...
    %evdata.prep, ...
    %evdata.stim1, ...
    %evdata.stim2,...
    %evdata.stim1_faceNonFace, ...
    %evdata.stim2_faceNonFace, ...
    %evdata.prevtrial_n1, ...
    %evdata.prevtrial_n2, ...
    %evdata.prevtrial_n3, ...
    
    %%%%%% models
    modelData(1).brainArea(brainArea).modelName=['model1-',currentBrainLabel,'-active+actual'];
    modelData(1).brainArea(brainArea).test_evs=[...
        evdata.active, ...
        evdata.actual ...
        ];
    modelData(1).brainArea(brainArea).test_ev_labels={'active','actual'};
    modelData(1).brainArea(brainArea).regressType='linear';
    
    %%%%%%%
    modelData(2).brainArea(brainArea).modelName=['model2-',currentBrainLabel,'-act+actual+Prep+deltaPrep'];
    modelData(2).brainArea(brainArea).test_evs=[...
        evdata.active, ...
        evdata.actual, ...
        evdata.prep, ...
        evdata.delta_p_rep ...
        ];
    modelData(2).brainArea(brainArea).test_ev_labels={'active','actual','prep','delta_p_rep'};
    modelData(2).brainArea(brainArea).regressType='linear';
    
    %%%%%%%
    modelData(3).brainArea(brainArea).modelName=['model2-',currentBrainLabel,'-act+actual+Prep+deltaPrep'];
    modelData(3).brainArea(brainArea).test_evs=[...
        evdata.active, ...
        evdata.actual, ...
        evdata.prep, ...
        evdata.delta_p_rep ...
        ];
    modelData(3).brainArea(brainArea).test_ev_labels={'active','actual','prep','delta_p_rep'};
    modelData(3).brainArea(brainArea).regressType='linear_orth';
    
%     %%%%%%
%     modelData(3).brainArea(brainArea).modelName=['model3-',currentBrainLabel,'-Prep+deltaPrep'];
%     modelData(3).brainArea(brainArea).test_evs=[...
%         evdata.prep, ...
%         evdata.delta_p_rep ...
%         ];
%     modelData(3).brainArea(brainArea).test_ev_labels={'Prep','deltaPrep'};
%     modelData(3).brainArea(brainArea).regressType='linear';
%     
%     %%%%%%%
%     modelData(4).brainArea(brainArea).modelName=['model4-',currentBrainLabel,'-stim1+face+stim2+face'];
%     modelData(4).brainArea(brainArea).test_evs=[...
%         evdata.stim1_faceNonFace, ...
%         evdata.stim2_faceNonFace, ...
%         evdata.stim1, ...
%         evdata.stim2...
%         ];
%     modelData(4).brainArea(brainArea).test_ev_labels={'stim1_face','stim2_face','stim1','stim2'};
%     modelData(4).brainArea(brainArea).regressType='linear';
%     
%     %%%%%%% Stepwise regression
%     modelData(5).brainArea(brainArea).modelName=['model5-',currentBrainLabel,'-stepwise-huge'];
%     modelData(5).brainArea(brainArea).test_evs=[...
%         evdata.active, ...
%         evdata.actual, ...
%         evdata.expect, ...
%         evdata.prep, ...
%         evdata.delta_p_rep, ...
%         evdata.trialNumber, ...
%         evdata.cor, ...
%         evdata.stim1_faceNonFace, ...
%         evdata.stim2_faceNonFace ,...
%         evdata.choseRepeat...
%         ];
%     modelData(5).brainArea(brainArea).test_ev_labels={'active','actual','expect','prep','delta-p-rep','trialNum','correct?','stim1face?','stim2face?','choserepeat?'};
%     modelData(5).brainArea(brainArea).regressType='stepwise';
%     
%     %%%%%% Stepwise regression
%     modelData(6).brainArea(brainArea).modelName=['model6-',currentBrainLabel,'-stepwise-notSohuge'];
%     modelData(6).brainArea(brainArea).test_evs=[...
%         evdata.actual, ...
%         evdata.expect, ...
%         evdata.prep, ...
%         evdata.delta_p_rep, ...
%         evdata.stim1_faceNonFace, ...
%         evdata.stim2_faceNonFace ,...
%         evdata.choseRepeat...
%         ];
%     modelData(6).brainArea(brainArea).test_ev_labels={'actual','expect','prep','delta-p-rep','stim1face?','stim2face?','choserepeat?'};
%     modelData(6).brainArea(brainArea).regressType='stepwise';
%     
%     
%     %%%%%% Stepwise regression
%     modelData(7).brainArea(brainArea).modelName=['model7-',currentBrainLabel,'-stepwise-reallyNotSoHuge'];
%     modelData(7).brainArea(brainArea).test_evs=[...
%         evdata.actual, ...
%         evdata.expect, ...
%         evdata.prep, ...
%         evdata.delta_p_rep ...
%         ];
%     modelData(7).brainArea(brainArea).test_ev_labels={'actual','expect','prep','delta-p-rep'};
%     modelData(7).brainArea(brainArea).regressType='stepwise';
    
    %for m=1:numel(modelData) % scroll through each model
    for m=1:numel(modelData) % scroll through each model
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
    end
    clear currentData currentLabel
end
save([exptdata.projectdir,'d500_modelData.mat'],'modelData')
clear modelData
% modelData is very very large (10x largest spikedata matrix). Should fix
% this up so that it isn't so large. 

fprintf('/n/n/n/n/n/n/n')


return






%% 6. Research #2 - REPETITION SUPPRESSION (passive)
% Model the effect of repeating a face stimulus y=b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9
% b0=intercept
% b1=repeat or alternate
% b4=trial number
% b5=prevtrial-1
% b6=prevtrial-2
% b7=prevtrial-3

paperData.rq2_BetaStuff=[];

for brainArea=1:2 % run on each brain area
    if brainArea==1
        currentSpikeData=V4_megaMatrix500NoScrubs;
        currentBrainLabel='V4';
    else
        currentSpikeData=TE_megaMatrix500NoScrubs;
        currentBrainLabel='TE';
    end
    cprintf('Red*',['* Q2: RepSuppression (',currentBrainLabel,') around ',epochNames{epochNum},' (session: '])
   
    % scroll through each session
    for ss=1:length(sessions)
        clear tmp_indxB tmp_indxN tmp_predictrz tmp_numNeuronsSession tmp_tt tmp_newTimeRange
        fprintf([num2str(sessions(ss)),'.'])
        tmp_indxB=find(evdata.session == sessions(ss) & evdata.active==-1);
        
        % Create set of EVs for each session
        tmp_predictrz = [ ...
            evdata.actual(tmp_indxB) ...
            evdata.trialNumber(tmp_indxB) ...
            evdata.stim1(tmp_indxB) ...
            evdata.stim2(tmp_indxB) ...
            %evdata.prevtrial_n1(tmp_indxB) ...
            %evdata.prevtrial_n2(tmp_indxB) ...
            %evdata.prevtrial_n3(tmp_indxB) ...
            ];
        ev_labels={'rep_v_alt','trialNum','stim1','stim2'};
        
        tmp_numNeuronsSession=max(currentSpikeData(currentSpikeData(:,2)==sessions(ss),10)); % how many neurons for this session
        for nn=1:tmp_numNeuronsSession
            tmp_indxN=find(currentSpikeData(:,2) == sessions(ss) & currentSpikeData(:,3) == nn & currentSpikeData(:,15)==0 & ismember(currentSpikeData(:,5),evdata.trialNumber(tmp_indxB)));
            %indxN=find(currentSpikeData(:,2) == sessions(ss) & currentSpikeData(:,3) == nn);
            
            if ~isempty(tmp_indxN)
                tmp_newTimeRange=timeEpochs500(epochNum,1):sampleRate:timeEpochs500(epochNum,2);
                if size(tmp_predictrz,1)~=length(tmp_indxN)
                    tmp_predictrz=tmp_predictrz(1:length(tmp_indxN),:); % fix needed in case behav and neuro trials don't match
                    fprintf('*PROB*')
                end
                for tt=1:length(tmp_newTimeRange)
                    paperData.rq2_BetaStuff(brainArea).session(ss).epoch(epochNum).betaWeights(nn,tt,:) = glmfit(tmp_predictrz, currentSpikeData(tmp_indxN,tmp_newTimeRange(tt)))';
                    %BetaStuff(brainArea).session(ss).epoch(epochNum).betaWeights(nn,tt,:)= glmfit(predictrz, ztransfnan_ahb(AllspikeData(indx,tt)))';
                end
            end
            clear tmp_indxN tt
            
        end
    end
    fprintf(') => done\n')
    clear tmp_* tt ss nn
end % brainArea

clear V4data TEdata
tmp_V4data=paperData.rq2_BetaStuff(1);
tmp_TEdata=paperData.rq2_BetaStuff(2);
epoch=8; tempData=[]; tempStuff=[]; %#ok<NASGU>
for sess=1:length(sessions)
    tempStuff=tmp_V4data.session(sess).epoch(epoch).betaWeights; % select V4
    for regressVar=1:length(ev_labels)+1
        tempData(sess,regressVar,:)=mean(squeeze(tempStuff(:,:,regressVar))); %#ok<SAGROW>
    end
end

% Figure 1
figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial')
subplot(2,3,1);  hold on 
plot(mean(squeeze(tempData(:,1,:))),'k-','LineWidth',1); % intercept
plot([750/sampleRate 750/sampleRate],[0 25],'k-','LineWidth',1.5)
plot([1500/sampleRate 1500/sampleRate],[0 25],'k-','LineWidth',1.5)
plot([250/sampleRate 250/sampleRate],[0 25],'k:','LineWidth',1.5)
plot([1000/sampleRate 1000/sampleRate],[0 25],'k:','LineWidth',1.5)
plot([1750/sampleRate 1750/sampleRate],[0 25],'k:','LineWidth',1.5)
fill([25 25 50 50],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
fill([100 100 125 125],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
fill([175 175 200 200],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
title({'Intercept',['Area V4, numSessions: ',num2str(length(sessions))]},'FontSize',16); ylim([0 20])
xlabel(['Time from Event (downsampled to x/',num2str(sampleRate),' ms, 1st Stim | 2nd Stim | Choice Array )']);
ylabel('Parameter Estimate')

subplot(2,3,[2 3]);  hold on 
h2=fillsteplot(squeeze(tempData(:,2,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'r'});
h3=fillsteplot(squeeze(tempData(:,3,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'b'});
h4=fillsteplot(squeeze(tempData(:,4,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'g'});
h5=fillsteplot(squeeze(tempData(:,5,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'m'});
%h6=fillsteplot(squeeze(tempData(:,6,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'c'});
plot([750/sampleRate 750/sampleRate],[-25 25],'k-','LineWidth',1.5)
plot([1500/sampleRate 1500/sampleRate],[-25 25],'k-','LineWidth',1.5)
plot([250/sampleRate 250/sampleRate],[-25 25],'k:','LineWidth',1.5)
plot([1000/sampleRate 1000/sampleRate],[-25 25],'k:','LineWidth',1.5)
plot([1750/sampleRate 1750/sampleRate],[-25 25],'k:','LineWidth',1.5)
fill([25 25 50 50],[-5 5 5 -5],'k','Edgecolor','none','FaceAlpha',0.1)
fill([100 100 125 125],[-5 5 5 -5],'k','Edgecolor','none','FaceAlpha',0.1)
fill([175 175 200 200],[-5 5 5 -5],'k','Edgecolor','none','FaceAlpha',0.1)
ylim([ -1.5 1.5]);
legend([h2 h3 h4 h5],ev_labels,'Location','south','NumColumns',numel(ev_labels))
xlabel(['Time from Event (downsampled to x/',num2str(sampleRate),' ms, 1st Stim | 2nd Stim | Choice Array )']);
ylabel('Parameter Estimate')
title({'Regressors of Interest',string(datetime)},'FontSize',16)

% TE
epoch=8; tempData=[]; tempStuff=[]; 
for sess=1:length(sessions)
    tempStuff=tmp_TEdata.session(sess).epoch(epoch).betaWeights; % select TE
    for regressVar=1:size(ev_labels,2)+1
        tempData(sess,regressVar,:)=mean(squeeze(tempStuff(:,:,regressVar))); %#ok<SAGROW>
    end
end

subplot(2,3,4);  hold on % Onset of 2nd Stim
plot(mean(squeeze(tempData(:,1,:))),'k-','LineWidth',1); % intercept
plot([750/sampleRate 750/sampleRate],[0 25],'k-','LineWidth',1.5)
plot([1500/sampleRate 1500/sampleRate],[0 25],'k-','LineWidth',1.5)
plot([250/sampleRate 250/sampleRate],[0 25],'k:','LineWidth',1.5)
plot([1000/sampleRate 1000/sampleRate],[0 25],'k:','LineWidth',1.5)
plot([1750/sampleRate 1750/sampleRate],[0 25],'k:','LineWidth',1.5)
fill([25 25 50 50],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
fill([100 100 125 125],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
fill([175 175 200 200],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
title({'Intercept',['Area TE, numSessions: ',num2str(length(sessions))]},'FontSize',16); ylim([0 20])
xlabel(['Time from Event (downsampled to x/',num2str(sampleRate),' ms, 1st Stim | 2nd Stim | Choice Array )']);
ylabel('Parameter Estimate')

subplot(2,3,[5 6]);  hold on % Onset of 1st Stim
h2=fillsteplot(squeeze(tempData(:,2,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'r'});
h3=fillsteplot(squeeze(tempData(:,3,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'b'});
h4=fillsteplot(squeeze(tempData(:,4,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'g'});
h5=fillsteplot(squeeze(tempData(:,5,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'m'});
%h6=fillsteplot(squeeze(tempData(:,6,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'c'});
plot([750/sampleRate 750/sampleRate],[-25 25],'k-','LineWidth',1.5)
plot([1500/sampleRate 1500/sampleRate],[-25 25],'k-','LineWidth',1.5)
plot([250/sampleRate 250/sampleRate],[-25 25],'k:','LineWidth',1.5)
plot([1000/sampleRate 1000/sampleRate],[-25 25],'k:','LineWidth',1.5)
plot([1750/sampleRate 1750/sampleRate],[-25 25],'k:','LineWidth',1.5)
fill([25 25 50 50],[-1.5 1.5 1.5 -1.5],'k','Edgecolor','none','FaceAlpha',0.1)
fill([100 100 125 125],[-1.5 1.5 1.5 -1.5],'k','Edgecolor','none','FaceAlpha',0.1)
fill([175 175 200 200],[-1.5 1.5 1.5 -1.5],'k','Edgecolor','none','FaceAlpha',0.1)
ylim([ -1.5 1.5]);
legend([h2 h3 h4 h5 ],ev_labels,'Location','south','NumColumns',numel(ev_labels))
xlabel(['Time from Event (downsampled to x/',num2str(sampleRate),' ms, 1st Stim | 2nd Stim | Choice Array )']);
ylabel('Parameter Estimate')
title('Regressors of Interest','FontSize',16)

clear tmp*

return

%%%--------------------------------------------------




%% RQ2.5 - Repetition Suppression in ACTIVE case



%% RESEARCH QUESTION #3 - EXPECTATION SUPPRESSION (both)
% Model the effect of repeating a face stimulus y=b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9
% b0=intercept
% b1=active or passive
% b2=repeat or alternate
% b3=p(rep) * active/passive
% b4=p(rep) * repeat/alternate
% b5=repeat/alternate * active/passive

% Bell et al: y = b0 + b1stimulus + b2p(face)RL + b3deltap(face)RL + b4prevchoice1 * stimulus + b5trial

cprintf('Red*','* Q3: Expectation Suppression\n')
paperData.rq3_BetaStuff=[];

for brainArea=1:2 % run on each brain area
    if brainArea==1
        currentSpikeData=V4_megaMatrix500NoScrubs;
        currentBrainLabel='V4';
    else
        currentSpikeData=TE_megaMatrix500NoScrubs;
        currentBrainLabel='TE';
    end
    fprintf('Analysing Session (session: ')
   
    % scroll through each session
    for ss=1:length(sessions)
        clear tmp_indxB tmp_indxN tmp_predictrz tmp_numNeuronsSession tmp_tt tmp_newTimeRange
        fprintf([num2str(sessions(ss)),'.'])
        tmp_indxB=find(evdata.session == sessions(ss));
        
        % Create set of EVs for each session
        % b1=active or passive
        % b2=repeat or alternate
        % b3=p(rep) * active/passive
        % b4=p(rep) * repeat/alternate
        % b5=repeat/alternate * active/passive
        
        tmp_predictrz = [ ...
            %evdata.prep(tmp_indxB) evdata.actual(tmp_indxB) .* evdata.prep(tmp_indxB) evdata.delta_p_rep(tmp_indxB)]; % model 1
            evdata.active(tmp_indxB) evdata.actual(tmp_indxB).* evdata.prep(tmp_indxB) evdata.actual(tmp_indxB)]; % model 2
            %evdata.actual(tmp_indxB) ...
            %evdata.active(tmp_indxB) .* evdata.actual(tmp_indxB) ...
            %evdata.prep(tmp_indxB).* evdata.active(tmp_indxB) ...
            %evdata.prep(tmp_indxB).* evdata.actual(tmp_indxB) ...
            %evdata.stim1(tmp_indxB) ...
            %evdata.stim2(tmp_indxB) ...
            %evdata.prevtrial_n1(tmp_indxB) ...
            %evdata.prevtrial_n2(tmp_indxB) ...
            %evdata.prevtrial_n3(tmp_indxB) ...
            
        ev_labels={'active','actual*prep','actual'};
        tmp_numNeuronsSession=max(currentSpikeData(currentSpikeData(:,2)==sessions(ss),10)); % how many neurons for this session
        for nn=1:tmp_numNeuronsSession
            tmp_indxN=find(currentSpikeData(:,2) == sessions(ss) & currentSpikeData(:,3) == nn ...
                & ismember(currentSpikeData(:,5),evdata.trialNumber(tmp_indxB)));
            if ~isempty(tmp_indxN)
                tmp_newTimeRange=timeEpochs500(epochNum,1):sampleRate:timeEpochs500(epochNum,2);
                if size(tmp_predictrz,1)~=length(tmp_indxN)
                    tmp_predictrz=tmp_predictrz(1:length(tmp_indxN),:); % fix needed in case behav and neuro trials don't match
                    fprintf('*PROB*')
                end
                for tt=1:length(tmp_newTimeRange)
                    [paperData.rq3_BetaStuff(brainArea).session(ss).epoch(epochNum).betaWeights(nn,tt,:), ...
                        paperData.rq3_BetaStuff(brainArea).session(ss).epoch(epochNum).modelDev(nn,tt)] = ...
                        glmfit(tmp_predictrz, currentSpikeData(tmp_indxN,tmp_newTimeRange(tt)));
                end
            end
            clear tmp_indxN tt 
        end
    end
    fprintf(') => done\n')
    clear tmp_* tt ss nn
end % brainArea


% [b,dev] = glmfit(...)returns dev, the deviance of the fit at the solution vector. The deviance is a generalization of the residual sum of squares. It is possible to perform an analysis of deviance to compare several models, each a subset of the other, and to test whether the model with more terms is significantly better than the model with fewer terms.


%stepwisefit

clear V4data TEdata
tmp_V4data=paperData.rq3_BetaStuff(1);
tmp_TEdata=paperData.rq3_BetaStuff(2);
epoch=8; tempData=[]; tempStuff=[]; %#ok<NASGU>
for sess=1:length(sessions)
    tempStuff=tmp_V4data.session(sess).epoch(epoch).betaWeights; % select V4
    for regressVar=1:length(ev_labels)+1
        tempData(sess,regressVar,:)=mean(squeeze(tempStuff(:,:,regressVar))); %#ok<SAGROW>
        
    end
end

% Figure 1
figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial')
subplot(2,3,1);  hold on 
plot(mean(squeeze(tempData(:,1,:))),'k-','LineWidth',1); % intercept
plot([750/sampleRate 750/sampleRate],[0 25],'k-','LineWidth',1.5)
plot([1500/sampleRate 1500/sampleRate],[0 25],'k-','LineWidth',1.5)
plot([250/sampleRate 250/sampleRate],[0 25],'k:','LineWidth',1.5)
plot([1000/sampleRate 1000/sampleRate],[0 25],'k:','LineWidth',1.5)
plot([1750/sampleRate 1750/sampleRate],[0 25],'k:','LineWidth',1.5)
fill([25 25 50 50],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
fill([100 100 125 125],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
fill([175 175 200 200],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
title({'Intercept',['Area V4, numSessions: ',num2str(length(sessions))]},'FontSize',16); ylim([0 20])
xlabel(['Time from Event (downsampled to x/',num2str(sampleRate),' ms, 1st Stim | 2nd Stim | Choice Array )']);
ylabel('Parameter Estimate')

subplot(2,3,[2 3]);  hold on 
h2=fillsteplot(squeeze(tempData(:,2,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'r'});
h3=fillsteplot(squeeze(tempData(:,3,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'b'});
h4=fillsteplot(squeeze(tempData(:,4,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'g'});
%h5=fillsteplot(squeeze(tempData(:,5,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'m'});
%h6=fillsteplot(squeeze(tempData(:,6,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'c'});
%h7=fillsteplot(squeeze(tempData(:,5,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'y'});
%h8=fillsteplot(squeeze(tempData(:,6,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'k'});
plot([750/sampleRate 750/sampleRate],[-25 25],'k-','LineWidth',1.5)
plot([1500/sampleRate 1500/sampleRate],[-25 25],'k-','LineWidth',1.5)
plot([250/sampleRate 250/sampleRate],[-25 25],'k:','LineWidth',1.5)
plot([1000/sampleRate 1000/sampleRate],[-25 25],'k:','LineWidth',1.5)
plot([1750/sampleRate 1750/sampleRate],[-25 25],'k:','LineWidth',1.5)
fill([25 25 50 50],[-5 5 5 -5],'k','Edgecolor','none','FaceAlpha',0.1)
fill([100 100 125 125],[-5 5 5 -5],'k','Edgecolor','none','FaceAlpha',0.1)
fill([175 175 200 200],[-5 5 5 -5],'k','Edgecolor','none','FaceAlpha',0.1)
ylim([ -2 2]);
legend([h2 h3 h4],ev_labels,'Location','south','NumColumns',numel(ev_labels))
xlabel(['Time from Event (downsampled to x/',num2str(sampleRate),' ms, 1st Stim | 2nd Stim | Choice Array )']);
ylabel('Parameter Estimate')
title({'Regressors of Interest',string(datetime)},'FontSize',16)

% TE
epoch=8; tempData=[]; tempStuff=[]; 
for sess=1:length(sessions)
    tempStuff=tmp_TEdata.session(sess).epoch(epoch).betaWeights; % select TE
    for regressVar=1:size(ev_labels,2)+1
        tempData(sess,regressVar,:)=mean(squeeze(tempStuff(:,:,regressVar))); %#ok<SAGROW>
    end
end

subplot(2,3,4);  hold on % Onset of 2nd Stim
plot(mean(squeeze(tempData(:,1,:))),'k-','LineWidth',1); % intercept
plot([750/sampleRate 750/sampleRate],[0 25],'k-','LineWidth',1.5)
plot([1500/sampleRate 1500/sampleRate],[0 25],'k-','LineWidth',1.5)
plot([250/sampleRate 250/sampleRate],[0 25],'k:','LineWidth',1.5)
plot([1000/sampleRate 1000/sampleRate],[0 25],'k:','LineWidth',1.5)
plot([1750/sampleRate 1750/sampleRate],[0 25],'k:','LineWidth',1.5)
fill([25 25 50 50],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
fill([100 100 125 125],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
fill([175 175 200 200],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
title({'Intercept',['Area TE, numSessions: ',num2str(length(sessions))]},'FontSize',16); ylim([0 20])
xlabel(['Time from Event (downsampled to x/',num2str(sampleRate),' ms, 1st Stim | 2nd Stim | Choice Array )']);
ylabel('Parameter Estimate')

subplot(2,3,[5 6]);  hold on % 
h2=fillsteplot(squeeze(tempData(:,2,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'r'});
h3=fillsteplot(squeeze(tempData(:,3,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'b'});
h4=fillsteplot(squeeze(tempData(:,4,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'g'});
%h5=fillsteplot(squeeze(tempData(:,5,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'m'});
%h6=fillsteplot(squeeze(tempData(:,6,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'c'});
%h7=fillsteplot(squeeze(tempData(:,5,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'y'});
%h8=fillsteplot(squeeze(tempData(:,6,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'k'});
plot([750/sampleRate 750/sampleRate],[-25 25],'k-','LineWidth',1.5)
plot([1500/sampleRate 1500/sampleRate],[-25 25],'k-','LineWidth',1.5)
plot([250/sampleRate 250/sampleRate],[-25 25],'k:','LineWidth',1.5)
plot([1000/sampleRate 1000/sampleRate],[-25 25],'k:','LineWidth',1.5)
plot([1750/sampleRate 1750/sampleRate],[-25 25],'k:','LineWidth',1.5)
fill([25 25 50 50],[-5 5 5 -5],'k','Edgecolor','none','FaceAlpha',0.1)
fill([100 100 125 125],[-5 5 5 -5],'k','Edgecolor','none','FaceAlpha',0.1)
fill([175 175 200 200],[-5 5 5 -5],'k','Edgecolor','none','FaceAlpha',0.1)
ylim([ -1.5 1.5]);
legend([h2 h3 h4],ev_labels,'Location','south','NumColumns',numel(ev_labels))
xlabel(['Time from Event (downsampled to x/',num2str(sampleRate),' ms, 1st Stim | 2nd Stim | Choice Array )']);
ylabel('Parameter Estimate')
title({'Regressors of Interest',string(datetime)},'FontSize',16)
















    
    for sess=1:length(sessions)
        clear tmp_indxB tmp_indxN tmp_numNeuronsSession tmp_sessRegressors
    
    fprintf([num2str(sess),'.'])
    
    tmp_indxB=find(regressors.session == sessions(sess));
    tmp_numNeuronsSession=max(data(data(:,2)==sessions(sess),10)); % how many neurons in this unit (this code is a bit of a hack but does the trick)
    
    % Create session-specific regressors
    tmp_sessRegressors=[];
    for nr=1:length(predictor_fields)
        tmp_sessRegressors=[tmp_sessRegressors, ...
            getfield(regressors,predictor_fields{1},{tmp_indxB})]; %#ok<AGROW>
    end
    
    % Run GLM on each neuron within session
    for nn=1:tmp_numNeuronsSession
        tmp_indxN=find(data(:,2) == sessions(sess) & data(:,3) == nn & ismember(data(:,5),regressors.trialNumber(tmp_indxB)), 1);
        if ~isempty(tmp_indxN)
            tmp_newTimeRange=timeEpochs500(epochNum,1):sampleRate:timeEpochs500(epochNum,2);
            if size(tmp_sessRegressors,1)~=length(tmp_indxN)
                tmp_sessRegressors=tmp_sessRegressors(1:length(tmp_indxN),:); % fix needed in case behav and neuro trials don't match
                fprintf('**numTr~=numEVs**')
            end
            for tt=1:length(tmp_newTimeRange)
                betas.session(sess).epoch(epochNum).betaWeights(nn,tt,:) = glmfit(tmp_sessRegressors, data(tmp_indxN,tmp_newTimeRange(tt)))';
                %betas.session(sess).epoch(epochNum).betaWeights(nn,tt,:)= glmfit(tmp_sessRegressors, ztransfnan_ahb(AllspikeData(indx,tt)))';
            end
        end
    end
    fprintf(') => done\n')
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
epochNum=8; tempData=[]; tempStuff=[]; %#ok<NASGU>
for sess=1:length(sessions)
    tempStuff=tmp_V4data.session(sess).epoch(epoch).betaWeights; % select V4
    for regressVar=1:size(tmp_predictrz,2)+1
        tempData(sess,regressVar,:)=mean(squeeze(tempStuff(:,:,regressVar))); %#ok<SAGROW>
    end
end


















%% 3. Perform Multiple Regression
clear BetaStuff
BetaStuff=struct('session',[]); warning off; close all
sampleRate=10; % test every X timepoints
epochNum=8; % limiting it to the more interesting epochs for the moment

for brainArea=1:2
    if brainArea==1
        AllspikeData=V4_megaMatrix500NoScrubs;
        currentBrainLabel='V4';
    else
        AllspikeData=TE_megaMatrix500NoScrubs;
        currentBrainLabel='TE';
    end
    fprintf(['* Analysing neural activity (',currentBrainLabel,') around ',epochNames{epochNum},' (session: '])
    
    %% THIS IS WHERE TO CHANGE EVs
    tmp_betas=ephys_multipleRegression_spikes(AllspikeData,data,{...
        'stim1_faceNonFace',...
        'stim2_faceNonFace',...
        'actual',...
        'smooth_p_rep',...
        'active'});
    ev_labels={'stim1','stim2','rep*act','smooth_p_rep','Active?'};
    %%
    
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
    BetaStuff(brainArea)=tmp_betas;

end % brainArea
clear tmp_*


tmp_V4data=BetaStuff(1);
tmp_TEdata=BetaStuff(2);

% END OF MULTIPLE REGRESSION

%% 4. Plot Data
epoch=8; tempData=[]; tempStuff=[]; %#ok<NASGU>
for sess=1:length(sessions)
    tempStuff=tmp_V4data.session(sess).epoch(epoch).betaWeights; % select V4
    for regressVar=1:size(tmp_predictrz,2)+1
        tempData(sess,regressVar,:)=mean(squeeze(tempStuff(:,:,regressVar))); %#ok<SAGROW>
    end
end

% Figure 1
figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial')
subplot(2,3,1);  hold on % Onset of 2nd Stim
plot(mean(squeeze(tempData(:,1,:))),'k-','LineWidth',1); % intercept
plot([750/sampleRate 750/sampleRate],[0 25],'k-','LineWidth',1.5)
plot([1500/sampleRate 1500/sampleRate],[0 25],'k-','LineWidth',1.5)
plot([250/sampleRate 250/sampleRate],[0 25],'k:','LineWidth',1.5)
plot([1000/sampleRate 1000/sampleRate],[0 25],'k:','LineWidth',1.5)
plot([1750/sampleRate 1750/sampleRate],[0 25],'k:','LineWidth',1.5)
fill([25 25 50 50],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
fill([100 100 125 125],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
fill([175 175 200 200],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
title({'Intercept',['Area V4, numSessions: ',num2str(length(sessions))]},'FontSize',16); ylim([0 20])
xlabel(['Time from Event (downsampled to x/',num2str(sampleRate),' ms, 1st Stim | 2nd Stim | Choice Array )']);
ylabel('Parameter Estimate')

subplot(2,3,[2 3]);  hold on % Onset of 1st Stim
h2=fillsteplot(squeeze(tempData(:,2,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'r'});
h3=fillsteplot(squeeze(tempData(:,3,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'b'});
h4=fillsteplot(squeeze(tempData(:,4,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'g'});
h5=fillsteplot(squeeze(tempData(:,5,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'m'});
h6=fillsteplot(squeeze(tempData(:,6,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'c'});
plot([750/sampleRate 750/sampleRate],[-25 25],'k-','LineWidth',1.5)
plot([1500/sampleRate 1500/sampleRate],[-25 25],'k-','LineWidth',1.5)
plot([250/sampleRate 250/sampleRate],[-25 25],'k:','LineWidth',1.5)
plot([1000/sampleRate 1000/sampleRate],[-25 25],'k:','LineWidth',1.5)
plot([1750/sampleRate 1750/sampleRate],[-25 25],'k:','LineWidth',1.5)
fill([25 25 50 50],[-1.5 1.5 1.5 -1.5],'k','Edgecolor','none','FaceAlpha',0.1)
fill([100 100 125 125],[-1.5 1.5 1.5 -1.5],'k','Edgecolor','none','FaceAlpha',0.1)
fill([175 175 200 200],[-1.5 1.5 1.5 -1.5],'k','Edgecolor','none','FaceAlpha',0.1)
ylim([ -1.5 1.5]);
legend([h2 h3 h4 h5 h6],ev_labels,'Location','south','NumColumns',numel(ev_labels))
xlabel(['Time from Event (downsampled to x/',num2str(sampleRate),' ms, 1st Stim | 2nd Stim | Choice Array )']);
ylabel('Parameter Estimate')
title({'Regressors of Interest',string(datetime)},'FontSize',16)

% TE
epoch=8; tempData=[]; tempStuff=[]; 
for sess=1:length(sessions)
    tempStuff=tmp_TEdata.session(sess).epoch(epoch).betaWeights; % select TE
    for regressVar=1:size(tmp_predictrz,2)+1
        tempData(sess,regressVar,:)=mean(squeeze(tempStuff(:,:,regressVar))); %#ok<SAGROW>
    end
end

subplot(2,3,4);  hold on % Onset of 2nd Stim
plot(mean(squeeze(tempData(:,1,:))),'k-','LineWidth',1); % intercept
plot([750/sampleRate 750/sampleRate],[0 25],'k-','LineWidth',1.5)
plot([1500/sampleRate 1500/sampleRate],[0 25],'k-','LineWidth',1.5)
plot([250/sampleRate 250/sampleRate],[0 25],'k:','LineWidth',1.5)
plot([1000/sampleRate 1000/sampleRate],[0 25],'k:','LineWidth',1.5)
plot([1750/sampleRate 1750/sampleRate],[0 25],'k:','LineWidth',1.5)
fill([25 25 50 50],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
fill([100 100 125 125],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
fill([175 175 200 200],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
title({'Intercept',['Area TE, numSessions: ',num2str(length(sessions))]},'FontSize',16); ylim([0 20])
xlabel(['Time from Event (downsampled to x/',num2str(sampleRate),' ms, 1st Stim | 2nd Stim | Choice Array )']);
ylabel('Parameter Estimate')

subplot(2,3,[5 6]);  hold on % Onset of 1st Stim
h2=fillsteplot(squeeze(tempData(:,2,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'r'});
h3=fillsteplot(squeeze(tempData(:,3,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'b'});
h4=fillsteplot(squeeze(tempData(:,4,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'g'});
h5=fillsteplot(squeeze(tempData(:,5,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'m'});
h6=fillsteplot(squeeze(tempData(:,6,:)),2,'-','v1',0:size(squeeze(tempData(:,2,:)),2),{'c'});
plot([750/sampleRate 750/sampleRate],[-25 25],'k-','LineWidth',1.5)
plot([1500/sampleRate 1500/sampleRate],[-25 25],'k-','LineWidth',1.5)
plot([250/sampleRate 250/sampleRate],[-25 25],'k:','LineWidth',1.5)
plot([1000/sampleRate 1000/sampleRate],[-25 25],'k:','LineWidth',1.5)
plot([1750/sampleRate 1750/sampleRate],[-25 25],'k:','LineWidth',1.5)
fill([25 25 50 50],[-1.5 1.5 1.5 -1.5],'k','Edgecolor','none','FaceAlpha',0.1)
fill([100 100 125 125],[-1.5 1.5 1.5 -1.5],'k','Edgecolor','none','FaceAlpha',0.1)
fill([175 175 200 200],[-1.5 1.5 1.5 -1.5],'k','Edgecolor','none','FaceAlpha',0.1)
ylim([ -1.5 1.5]);
legend([h2 h3 h4 h5 h6],ev_labels,'Location','south','NumColumns',numel(ev_labels))
xlabel(['Time from Event (downsampled to x/',num2str(sampleRate),' ms, 1st Stim | 2nd Stim | Choice Array )']);
ylabel('Parameter Estimate')
title('Regressors of Interest','FontSize',16)


% =================================================================================================================================================
% =================================================================================================================================================
% =================================================================================================================================================
% =================================================================================================================================================
% =================================================================================================================================================
% =================================================================================================================================================
% =================================================================================================================================================
% =================================================================================================================================================
% =================================================================================================================================================


%!feb 26, 2019
%!panel 1,2 - plot intercept and then regressors (simple ones, like "active"
%!panel 3,4 - plot intercept and regressors (stimulus ones - figure out the face/non-face thing)
%!panel 5,6 - prep/smooth-prep
%!next look at prep*expect






