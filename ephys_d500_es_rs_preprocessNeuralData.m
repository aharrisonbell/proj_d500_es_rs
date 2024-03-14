% ephys_d500_es_rs_preprocessNeuralData
% AHB, Mar 4, 2018, based on ephys_preprocessRSVP300, part of ES_RS study
vers_ephys_dms='1.1; Dec 20, 2018';
% 1.0 - original version Mar 04, 2018
% 1.1 - added QA figure + fixed spike density function scaling
% This program preprocesses data obtained using standard ephys tasks


[d500specs.fList,d500specs.pList] = matlab.codetools.requiredFilesAndProducts('ephys_d500_es_rs_preprocessNeuralData.m');

clear datafiles ddir df dd numtrials 
clc
global exptdata
fprintf('*=========================================*\n')
fprintf('| ephys_d500_es_rs_preprocessNeuralData.m |\n')
fprintf('*=========================================*\n')
fprintf(['Version: ',vers_ephys_dms,'\n'])
disp(['Data location:     ',exptdata.datalocation]);
disp(['Processed NEV dir: ',exptdata.processedDatadir]);
disp(['Project Directory: ',exptdata.projectdir]);
disp(['Monkey name:       ',exptdata.monkeyname]);
disp(['Study name:        ',exptdata.analysisName]);

disp(' ')
%% 1. Search directory for datafiles
ddir=dir([exptdata.processedDatadir,filesep,exptdata.monkeyname,'*DMS500*','mat']); % search for any files starting with MONKEYNAME
if isempty(ddir)
    disp('No DMS500 datafiles found')
    return
end
for dd=1:numel(ddir), datafiles(dd)={ddir(dd).name}; end %#ok<*SAGROW>

%% 2. Set task-specific variables
%exptdata.datafiles_behavSamplingRate=exptdata.dms500.datafiles_behavSamplingRate;
%exptdata.datafiles_spikeSamplingRate=exptdata.dms500.datafiles_spikeSamplingRate;

% Scroll through each datafile
for df=1:numel(datafiles)
    clear numtrials behavMarkers trialData behavData diffConditions numV4neurons numTEpneurons tempname
    clear V4_trial* V4_cond* TEp_trial* TEp_cond*
    tempname=char(datafiles(df));
    % Check if file has already been preprocessed
    if exist([exptdata.analysisdir,exptdata.analysisName,filesep,tempname(1:end-13),'sessionData.mat'],'file') && exptdata.reprocess~=1
        fprintf(['*** Already processed ',tempname,'. Skipping...\n'])
        continue
    end
    
    %% 3. Load datafile (must already have been run through EPHYS_PROCESSNEVFILE)
    cprintf('*blue',['\nAnalysing ',tempname,' (',num2str(df),'/',num2str(numel(datafiles)),')...\n'])
    load([exptdata.processedDatadir,tempname]);
    if isempty(cell2mat(fileInfo(5))) %exptdata.datafiles_behavSamplingRate(df)=cell2mat(fileInfo(5));
        exptdata.datafiles_behavSamplingRate(df)=1000;
    else
        exptdata.datafiles_behavSamplingRate(df)=cell2mat(fileInfo(5));
    end
    numtrials=min([length(trialData) length(behavData.TrialNumber)]); % TOTAL NUMBER OF TRIALS (shitty compromise because sometimes shit not good)
    numV4neurons=trialData(1).numV4neurons;
    try %% catch code for old method of naming TEp vs. TeO
        numTEpneurons=trialData(1).numTeOneurons;
    catch
        numTEpneurons=trialData(1).numTEpneurons;
    end
    
    %%******************************************************************************************************************************
    %%******************************************************************************************************************************
    %% 3. Analyse Behaviour
    % Reaction Time
    behavData.newReactionTime=nan(numtrials,1);
    for tt=1:numtrials
        clear tempCodes tempTime tempChoice*
        tempCodes(:,1)=cell2mat(behavData.CodeNumbers(tt));
        tempCodes(:,2)=cell2mat(behavData.CodeTimes(tt));
        if behavData.TrialError(tt)==0
            tempChoiceOnIndex=find(tempCodes==exptdata.lsn_ML.goSignal,1,'first');
            tempChoiceMadeIndex=find(tempCodes==exptdata.lsn_ML.eyeOnTarget,1,'first');
            behavData.newReactionTime(tt,1)=tempCodes(tempChoiceMadeIndex,2)-tempCodes(tempChoiceOnIndex,2);
        elseif behavData.TrialError(tt)==6
            tempChoiceOnIndex=find(tempCodes(:,1)==exptdata.lsn_ML.goSignal,1,'first');
            tempChoiceMadeIndex=find(ismember(tempCodes(:,1),[exptdata.lsn_ML.eyeIncorrectTarget exptdata.lsn_ML.eyeIncorrectTarget2]),1,'first');
            behavData.newReactionTime(tt,1)=tempCodes(tempChoiceMadeIndex,2)-tempCodes(tempChoiceOnIndex,2);
        end
        behavData.newReactionTime(tt,2)=ztransfnan_ahb(behavData.newReactionTime(tt,1));
    end
    behavData.newReactionTime(:,2)=ztransfnan_ahb(behavData.newReactionTime(:,1));
    
    % if number of trials in ML doesn't match BlackRock...
    behavData.BlockNumber=behavData.BlockNumber(1:length(behavData.newReactionTime));
    
    % Generate quick figure
    figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial')
    subplot(1,2,1) % reaction time
    bardata=nan(6,2);
    bardata(1,1)=nanmean(behavData.newReactionTime(behavData.TrialError==0 & ismember(behavData.BlockNumber,1) & ismember(behavData.ConditionNumber,exptdata.dms500.expectRepeat_actualRepeat),2)); %#ok<*NANMEAN>
    bardata(1,2)=nanmean(behavData.newReactionTime(behavData.TrialError==0 & ismember(behavData.BlockNumber,1) & ismember(behavData.ConditionNumber,exptdata.dms500.expectRepeat_actualAltern),2));
    bardata(2,1)=nanmean(behavData.newReactionTime(behavData.TrialError==0 & ismember(behavData.BlockNumber,2) & ismember(behavData.ConditionNumber,exptdata.dms500.expectAltern_actualRepeat),2));
    bardata(2,2)=nanmean(behavData.newReactionTime(behavData.TrialError==0 & ismember(behavData.BlockNumber,2) & ismember(behavData.ConditionNumber,exptdata.dms500.expectAltern_actualAltern),2));
    bardata(3,1)=nanmean(behavData.newReactionTime(behavData.TrialError==0 & ismember(behavData.BlockNumber,3) & ismember(behavData.ConditionNumber,exptdata.dms500.expectEither_actualRepeat),2));
    bardata(3,2)=nanmean(behavData.newReactionTime(behavData.TrialError==0 & ismember(behavData.BlockNumber,3) & ismember(behavData.ConditionNumber,exptdata.dms500.expectEither_actualAltern),2));
    bardata(4,1)=nanmean(behavData.newReactionTime(behavData.TrialError==0 & ismember(behavData.BlockNumber,4) & ismember(behavData.ConditionNumber,exptdata.dms500.expectRepeat_actualRepeat),2));
    bardata(4,2)=nanmean(behavData.newReactionTime(behavData.TrialError==0 & ismember(behavData.BlockNumber,4) & ismember(behavData.ConditionNumber,exptdata.dms500.expectRepeat_actualAltern),2));
    bardata(5,1)=nanmean(behavData.newReactionTime(behavData.TrialError==0 & ismember(behavData.BlockNumber,5) & ismember(behavData.ConditionNumber,exptdata.dms500.expectAltern_actualRepeat),2));
    bardata(5,2)=nanmean(behavData.newReactionTime(behavData.TrialError==0 & ismember(behavData.BlockNumber,5) & ismember(behavData.ConditionNumber,exptdata.dms500.expectAltern_actualAltern),2));
    bardata(6,1)=nanmean(behavData.newReactionTime(behavData.TrialError==0 & ismember(behavData.BlockNumber,6) & ismember(behavData.ConditionNumber,exptdata.dms500.expectEither_actualRepeat),2));
    bardata(6,2)=nanmean(behavData.newReactionTime(behavData.TrialError==0 & ismember(behavData.BlockNumber,6) & ismember(behavData.ConditionNumber,exptdata.dms500.expectEither_actualAltern),2));
    bar(bardata,'group')
    text(1,.9,'ActiveTask','FontSize',12,'FontWeight','Bold'); text(4,.9,'PassiveTask','FontSize',12,'FontWeight','Bold')
    xticklabels({'ExpRep','ExpAlt','NoExp','ExpRep','ExpAlt','NoExp'}); legend('Actual Repeat','Actual Altern')
    ylabel('Z-score (ReactionTime)'); ylim([0 1])
    title(['Behavioural Results - ',tempname],'Interpreter','None','FontWeight','Bold','FontSize',14);
    subplot(1,2,2) % choice probability
    bardata=nan(6,2);
    bardata(1,1)=length(find(behavData.TrialError==0 & ismember(behavData.BlockNumber,1) & ismember(behavData.ConditionNumber,exptdata.dms500.expectRepeat_actualRepeat))) / length(find(behavData.TrialError~=4 & ismember(behavData.BlockNumber,1) & ismember(behavData.ConditionNumber,exptdata.dms500.expectRepeat_actualRepeat))); % choose repeat
    bardata(1,2)=length(find(behavData.TrialError==0 & ismember(behavData.BlockNumber,1) & ismember(behavData.ConditionNumber,exptdata.dms500.expectRepeat_actualAltern))) / length(find(behavData.TrialError~=4 & ismember(behavData.BlockNumber,1) & ismember(behavData.ConditionNumber,exptdata.dms500.expectRepeat_actualAltern))); % choose repeat
    bardata(2,1)=length(find(behavData.TrialError==0 & ismember(behavData.BlockNumber,2) & ismember(behavData.ConditionNumber,exptdata.dms500.expectAltern_actualRepeat))) / length(find(behavData.TrialError~=4 & ismember(behavData.BlockNumber,2) & ismember(behavData.ConditionNumber,exptdata.dms500.expectAltern_actualRepeat))); % choose repeat
    bardata(2,2)=length(find(behavData.TrialError==0 & ismember(behavData.BlockNumber,2) & ismember(behavData.ConditionNumber,exptdata.dms500.expectAltern_actualAltern))) / length(find(behavData.TrialError~=4 & ismember(behavData.BlockNumber,2) & ismember(behavData.ConditionNumber,exptdata.dms500.expectAltern_actualAltern))); % choose repeat
    bardata(3,1)=length(find(behavData.TrialError==0 & ismember(behavData.BlockNumber,3) & ismember(behavData.ConditionNumber,exptdata.dms500.expectEither_actualRepeat))) / length(find(behavData.TrialError~=4 & ismember(behavData.BlockNumber,3) & ismember(behavData.ConditionNumber,exptdata.dms500.expectEither_actualRepeat))); % choose repeat
    bardata(3,2)=length(find(behavData.TrialError==0 & ismember(behavData.BlockNumber,3) & ismember(behavData.ConditionNumber,exptdata.dms500.expectEither_actualAltern))) / length(find(behavData.TrialError~=4 & ismember(behavData.BlockNumber,3) & ismember(behavData.ConditionNumber,exptdata.dms500.expectEither_actualAltern))); % choose repeat
    bardata(4,1)=length(find(behavData.TrialError==0 & ismember(behavData.BlockNumber,4) & ismember(behavData.ConditionNumber,exptdata.dms500.expectRepeat_actualRepeat))) / length(find(behavData.TrialError~=4 & ismember(behavData.BlockNumber,4) & ismember(behavData.ConditionNumber,exptdata.dms500.expectRepeat_actualRepeat))); % choose repeat
    bardata(4,2)=length(find(behavData.TrialError==0 & ismember(behavData.BlockNumber,4) & ismember(behavData.ConditionNumber,exptdata.dms500.expectRepeat_actualAltern))) / length(find(behavData.TrialError~=4 & ismember(behavData.BlockNumber,4) & ismember(behavData.ConditionNumber,exptdata.dms500.expectRepeat_actualAltern))); % choose repeat
    bardata(5,1)=length(find(behavData.TrialError==0 & ismember(behavData.BlockNumber,5) & ismember(behavData.ConditionNumber,exptdata.dms500.expectAltern_actualRepeat))) / length(find(behavData.TrialError~=4 & ismember(behavData.BlockNumber,5) & ismember(behavData.ConditionNumber,exptdata.dms500.expectAltern_actualRepeat))); % choose repeat
    bardata(5,2)=length(find(behavData.TrialError==0 & ismember(behavData.BlockNumber,5) & ismember(behavData.ConditionNumber,exptdata.dms500.expectAltern_actualAltern))) / length(find(behavData.TrialError~=4 & ismember(behavData.BlockNumber,5) & ismember(behavData.ConditionNumber,exptdata.dms500.expectAltern_actualAltern))); % choose repeat
    bardata(6,1)=length(find(behavData.TrialError==0 & ismember(behavData.BlockNumber,6) & ismember(behavData.ConditionNumber,exptdata.dms500.expectEither_actualRepeat))) / length(find(behavData.TrialError~=4 & ismember(behavData.BlockNumber,6) & ismember(behavData.ConditionNumber,exptdata.dms500.expectEither_actualRepeat))); % choose repeat
    bardata(6,2)=length(find(behavData.TrialError==0 & ismember(behavData.BlockNumber,6) & ismember(behavData.ConditionNumber,exptdata.dms500.expectEither_actualAltern))) / length(find(behavData.TrialError~=4 & ismember(behavData.BlockNumber,6) & ismember(behavData.ConditionNumber,exptdata.dms500.expectEither_actualAltern))); % choose repeat
    bar(bardata,'group')
    text(1,1.1,'ActiveTask','FontSize',12,'FontWeight','Bold'); text(4,1.1,'PassiveTask','FontSize',12,'FontWeight','Bold')
    xticklabels({'ExpRep','ExpAlt','NoExp','ExpRep','ExpAlt','NoExp'}); legend('Actual Repeat','Actual Altern')
    ylabel('Percent Correct (%)'); ylim([0.5 1.2])
        
    %%******************************************************************************************************************************
    %%******************************************************************************************************************************
    %% 4. Neurophysiology Data Analysis
    diffConditions=unique(behavData.ConditionNumber);
    totalConditions=1:192;
    % Preallocate matrices
    V4_trialSpden_fpON=nan(numtrials,numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_trialSpden_1stStimON=nan(numtrials,numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_trialSpden_2ndStimON=nan(numtrials,numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_trialSpden_ChoiceArrayON=nan(numtrials,numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_trialSpden_ResponseGOOD=nan(numtrials,numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_trialSpden_ResponseBAD=nan(numtrials,numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_trialSpden_reward=nan(numtrials,numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_trialSpike_fpON=nan(numtrials,numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_trialSpike_1stStimON=nan(numtrials,numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_trialSpike_2ndStimON=nan(numtrials,numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_trialSpike_ChoiceArrayON=nan(numtrials,numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_trialSpike_ResponseGOOD=nan(numtrials,numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_trialSpike_ResponseBAD=nan(numtrials,numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_trialSpike_reward=nan(numtrials,numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    
    V4_condSpden_fpON=nan(length(totalConditions),numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_condSpike_fpON=nan(length(totalConditions),numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_condSpden_1stStimON=nan(length(totalConditions),numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_condSpike_1stStimON=nan(length(totalConditions),numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_condSpden_2ndStimON=nan(length(totalConditions),numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_condSpike_2ndStimON=nan(length(totalConditions),numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_condSpden_ChoiceArrayON=nan(length(totalConditions),numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_condSpike_ChoiceArrayON=nan(length(totalConditions),numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_condSpden_ResponseGOOD=nan(length(totalConditions),numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_condSpike_ResponseGOOD=nan(length(totalConditions),numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_condSpden_ResponseBAD=nan(length(totalConditions),numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_condSpike_ResponseBAD=nan(length(totalConditions),numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_condSpden_reward=nan(length(totalConditions),numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    V4_condSpike_reward=nan(length(totalConditions),numV4neurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    
    TEp_trialSpden_fpON=nan(numtrials,numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_trialSpden_1stStimON=nan(numtrials,numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_trialSpden_2ndStimON=nan(numtrials,numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_trialSpden_ChoiceArrayON=nan(numtrials,numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_trialSpden_ResponseGOOD=nan(numtrials,numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_trialSpden_ResponseBAD=nan(numtrials,numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_trialSpden_reward=nan(numtrials,numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_trialSpike_fpON=nan(numtrials,numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_trialSpike_1stStimON=nan(numtrials,numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_trialSpike_2ndStimON=nan(numtrials,numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_trialSpike_ChoiceArrayON=nan(numtrials,numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_trialSpike_ResponseGOOD=nan(numtrials,numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_trialSpike_ResponseBAD=nan(numtrials,numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_trialSpike_reward=nan(numtrials,numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    
    TEp_condSpden_fpON=nan(length(totalConditions),numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_condSpike_fpON=nan(length(totalConditions),numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_condSpden_1stStimON=nan(length(totalConditions),numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_condSpike_1stStimON=nan(length(totalConditions),numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_condSpden_2ndStimON=nan(length(totalConditions),numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_condSpike_2ndStimON=nan(length(totalConditions),numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_condSpden_ChoiceArrayON=nan(length(totalConditions),numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_condSpike_ChoiceArrayON=nan(length(totalConditions),numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_condSpden_ResponseGOOD=nan(length(totalConditions),numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_condSpike_ResponseGOOD=nan(length(totalConditions),numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_condSpden_ResponseBAD=nan(length(totalConditions),numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_condSpike_ResponseBAD=nan(length(totalConditions),numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_condSpden_reward=nan(length(totalConditions),numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    TEp_condSpike_reward=nan(length(totalConditions),numTEpneurons,length(exptdata.xrange_psths)); % going to pack everything into one big matrix
    
    tic
   
    %% 4.1 Create spike density functions for each trial
    fprintf(' * Creating spike density functions for each trial... (be patient) . ')
    h = waitbar(0,'Preparing','Name','Creating spike density functions (trials)...');
    for tt=1:numtrials
        waitbar(tt/numtrials,h,['Trial No. ',num2str(tt)])
        clear trialIndex tempSpikeData
        % identify CORRECT trials with that condition number
        tempSpikeData=trialData(tt);
        [~,V4_trialSpike_fpON(tt,:,:),~,V4_trialSpden_fpON(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'V4spikeData',exptdata.lsn_ML.eyeOnFP);
        [~,V4_trialSpike_1stStimON(tt,:,:),~,V4_trialSpden_1stStimON(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'V4spikeData',exptdata.lsn_ML.stim1stPresON);
        [~,V4_trialSpike_2ndStimON(tt,:,:),~,V4_trialSpden_2ndStimON(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'V4spikeData',exptdata.lsn_ML.stim2ndPresON);
        [~,V4_trialSpike_ChoiceArrayON(tt,:,:),~,V4_trialSpden_ChoiceArrayON(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'V4spikeData',exptdata.lsn_ML.goSignal);
        [~,V4_trialSpike_ResponseGOOD(tt,:,:),~,V4_trialSpden_ResponseGOOD(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'V4spikeData',exptdata.lsn_ML.eyeOnTarget);
        [~,V4_trialSpike_ResponseBAD(tt,:,:),~,V4_trialSpden_ResponseBAD(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'V4spikeData',exptdata.lsn_ML.eyeIncorrectTarget);
        [~,V4_trialSpike_reward(tt,:,:),~,V4_trialSpden_reward(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'V4spikeData',exptdata.lsn_ML.Reward);
        try % this is necessary to control for differences in naming scheme in UTAH data
            [~,TEp_trialSpike_fpON(tt,:,:),~,TEp_trialSpden_fpON(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TEpspikeData',exptdata.lsn_ML.eyeOnFP);
            [~,TEp_trialSpike_1stStimON(tt,:,:),~,TEp_trialSpden_1stStimON(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TEpspikeData',exptdata.lsn_ML.stim1stPresON);
            [~,TEp_trialSpike_2ndStimON(tt,:,:),~,TEp_trialSpden_2ndStimON(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TEpspikeData',exptdata.lsn_ML.stim2ndPresON);
            [~,TEp_trialSpike_ChoiceArrayON(tt,:,:),~,TEp_trialSpden_ChoiceArrayON(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TEpspikeData',exptdata.lsn_ML.goSignal);
            [~,TEp_trialSpike_ResponseGOOD(tt,:,:),~,TEp_trialSpden_ResponseGOOD(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TEpspikeData',exptdata.lsn_ML.eyeOnTarget);
            [~,TEp_trialSpike_ResponseBAD(tt,:,:),~,TEp_trialSpden_ResponseBAD(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TEpspikeData',exptdata.lsn_ML.eyeIncorrectTarget);
            [~,TEp_trialSpike_reward(tt,:,:),~,TEp_trialSpden_reward(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TEpspikeData',exptdata.lsn_ML.Reward);
        catch
            [~,TEp_trialSpike_fpON(tt,:,:),~,TEp_trialSpden_fpON(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TeOspikeData',exptdata.lsn_ML.eyeOnFP);
            [~,TEp_trialSpike_1stStimON(tt,:,:),~,TEp_trialSpden_1stStimON(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TeOspikeData',exptdata.lsn_ML.stim1stPresON);
            [~,TEp_trialSpike_2ndStimON(tt,:,:),~,TEp_trialSpden_2ndStimON(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TeOspikeData',exptdata.lsn_ML.stim2ndPresON);
            [~,TEp_trialSpike_ChoiceArrayON(tt,:,:),~,TEp_trialSpden_ChoiceArrayON(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TeOspikeData',exptdata.lsn_ML.goSignal);
            [~,TEp_trialSpike_ResponseGOOD(tt,:,:),~,TEp_trialSpden_ResponseGOOD(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TeOspikeData',exptdata.lsn_ML.eyeOnTarget);
            [~,TEp_trialSpike_ResponseBAD(tt,:,:),~,TEp_trialSpden_ResponseBAD(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TeOspikeData',exptdata.lsn_ML.eyeIncorrectTarget);
            [~,TEp_trialSpike_reward(tt,:,:),~,TEp_trialSpden_reward(tt,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TeOspikeData',exptdata.lsn_ML.Reward);
        end
    end
    fprintf('Done.\n'); close(h);
    
    %% 4.2 Create average spike density functions for each condition (CORRECT TRIALS ONLY)
    fprintf(' * Creating spike density functions for each condition... (be patient) . ')
    h = waitbar(0,'Preparing','Name','Creating spike density functions (conditions)...');
    for cond=1:length(totalConditions)
        waitbar(cond/length(totalConditions),h,['Condition No. ',num2str(cond)])
        clear trialIndex tempSpikeData
        % identify CORRECT trials with that condition number
        trialIndex=find(behavData.TrialError==0 & behavData.ConditionNumber==cond);
        if ~isempty(trialIndex)
            tempSpikeData=trialData(trialIndex);
            [~,V4_condSpike_fpON(cond,:,:),~,V4_condSpden_fpON(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'V4spikeData',exptdata.lsn_ML.eyeOnFP);
            [~,V4_condSpike_1stStimON(cond,:,:),~,V4_condSpden_1stStimON(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'V4spikeData',exptdata.lsn_ML.stim1stPresON);
            [~,V4_condSpike_2ndStimON(cond,:,:),~,V4_condSpden_2ndStimON(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'V4spikeData',exptdata.lsn_ML.stim2ndPresON);
            [~,V4_condSpike_ChoiceArrayON(cond,:,:),~,V4_condSpden_ChoiceArrayON(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'V4spikeData',exptdata.lsn_ML.goSignal);
            [~,V4_condSpike_ResponseGOOD(cond,:,:),~,V4_condSpden_ResponseGOOD(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'V4spikeData',exptdata.lsn_ML.eyeOnTarget);
            [~,V4_condSpike_ResponseBAD(cond,:,:),~,V4_condSpden_ResponseBAD(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'V4spikeData',exptdata.lsn_ML.eyeIncorrectTarget);
            [~,V4_condSpike_reward(cond,:,:),~,V4_condSpden_reward(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'V4spikeData',exptdata.lsn_ML.Reward);
            try
                [~,TEp_condSpike_fpON(cond,:,:),~,TEp_condSpden_fpON(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TEpspikeData',exptdata.lsn_ML.eyeOnFP);
                [~,TEp_condSpike_1stStimON(cond,:,:),~,TEp_condSpden_1stStimON(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TEpspikeData',exptdata.lsn_ML.stim1stPresON);
                [~,TEp_condSpike_2ndStimON(cond,:,:),~,TEp_condSpden_2ndStimON(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TEpspikeData',exptdata.lsn_ML.stim2ndPresON);
                [~,TEp_condSpike_ChoiceArrayON(cond,:,:),~,TEp_condSpden_ChoiceArrayON(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TEpspikeData',exptdata.lsn_ML.goSignal);
                [~,TEp_condSpike_ResponseGOOD(cond,:,:),~,TEp_condSpden_ResponseGOOD(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TEpspikeData',exptdata.lsn_ML.eyeOnTarget);
                [~,TEp_condSpike_ResponseBAD(cond,:,:),~,TEp_condSpden_ResponseBAD(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TEpspikeData',exptdata.lsn_ML.eyeIncorrectTarget);
                [~,TEp_condSpike_reward(cond,:,:),~,TEp_condSpden_reward(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TEpspikeData',exptdata.lsn_ML.Reward);
            catch
                [~,TEp_condSpike_fpON(cond,:,:),~,TEp_condSpden_fpON(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TeOspikeData',exptdata.lsn_ML.eyeOnFP);
                [~,TEp_condSpike_1stStimON(cond,:,:),~,TEp_condSpden_1stStimON(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TeOspikeData',exptdata.lsn_ML.stim1stPresON);
                [~,TEp_condSpike_2ndStimON(cond,:,:),~,TEp_condSpden_2ndStimON(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TeOspikeData',exptdata.lsn_ML.stim2ndPresON);
                [~,TEp_condSpike_ChoiceArrayON(cond,:,:),~,TEp_condSpden_ChoiceArrayON(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TeOspikeData',exptdata.lsn_ML.goSignal);
                [~,TEp_condSpike_ResponseGOOD(cond,:,:),~,TEp_condSpden_ResponseGOOD(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TeOspikeData',exptdata.lsn_ML.eyeOnTarget);
                [~,TEp_condSpike_ResponseBAD(cond,:,:),~,TEp_condSpden_ResponseBAD(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TeOspikeData',exptdata.lsn_ML.eyeIncorrectTarget);
                [~,TEp_condSpike_reward(cond,:,:),~,TEp_condSpden_reward(cond,:,:)]=ephys_sortSpikeTrains(tempSpikeData,'TeOspikeData',exptdata.lsn_ML.Reward);
            end
        end
    end
    fprintf('Done.\n'); close(h);
    toc
    
    %% 4.3 Export Data for population-level analysis
    % V4 data
    % Create empty matrix for this session
    V4_bigMatrix=nan(numV4neurons*length(trialData),(length(exptdata.xrange_psths)*7)+20);
    % Paste in trial values
    if strcmp(exptdata.monkeyname,'Vortex')
        V4_bigMatrix(:,1)=1; % monkey number
    else
        V4_bigMatrix(:,1)=2; % monkey number
    end
    V4_bigMatrix(:,2)=df; % session number
    V4_bigMatrix(:,3:9)=createTrial_UnitVectors(trialData,numV4neurons,behavData);
    V4_bigMatrix(:,10)=numV4neurons;
    trialGap=1:numtrials:numV4neurons*numtrials;
    windowGap=21:length(exptdata.xrange_psths):length(exptdata.xrange_psths)*7;
    for un=1:numV4neurons
        V4_bigMatrix(trialGap(un):trialGap(un)+numtrials-1,windowGap(1):windowGap(1)+length(exptdata.xrange_psths)-1)=squeeze(V4_trialSpden_fpON(:,un,:));
        V4_bigMatrix(trialGap(un):trialGap(un)+numtrials-1,windowGap(2):windowGap(2)+length(exptdata.xrange_psths)-1)=squeeze(V4_trialSpden_1stStimON(:,un,:));
        V4_bigMatrix(trialGap(un):trialGap(un)+numtrials-1,windowGap(3):windowGap(3)+length(exptdata.xrange_psths)-1)=squeeze(V4_trialSpden_2ndStimON(:,un,:));
        V4_bigMatrix(trialGap(un):trialGap(un)+numtrials-1,windowGap(4):windowGap(4)+length(exptdata.xrange_psths)-1)=squeeze(V4_trialSpden_ChoiceArrayON(:,un,:));
        V4_bigMatrix(trialGap(un):trialGap(un)+numtrials-1,windowGap(5):windowGap(5)+length(exptdata.xrange_psths)-1)=squeeze(V4_trialSpden_ResponseGOOD(:,un,:));
        V4_bigMatrix(trialGap(un):trialGap(un)+numtrials-1,windowGap(6):windowGap(6)+length(exptdata.xrange_psths)-1)=squeeze(V4_trialSpden_ResponseBAD(:,un,:));
        V4_bigMatrix(trialGap(un):trialGap(un)+numtrials-1,windowGap(7):windowGap(7)+length(exptdata.xrange_psths)-1)=squeeze(V4_trialSpden_reward(:,un,:));
    end
    
    % TEp data
    % Create empty matrix for this session
    TEp_bigMatrix=nan(numTEpneurons*length(trialData),(length(exptdata.xrange_psths)*7)+20);
    % Paste in trial values
    if strcmp(exptdata.monkeyname,'Vortex')
        TEp_bigMatrix(:,1)=1; % monkey number
    else
        TEp_bigMatrix(:,1)=2; % monkey number
    end
    TEp_bigMatrix(:,2)=df; % session number
    TEp_bigMatrix(:,3:9)=createTrial_UnitVectors(trialData,numTEpneurons,behavData);
    TEp_bigMatrix(:,10)=numTEpneurons;
    trialGap=1:numtrials:numTEpneurons*numtrials;
    windowGap=21:length(exptdata.xrange_psths):length(exptdata.xrange_psths)*7;
    for un=1:numTEpneurons
        TEp_bigMatrix(trialGap(un):trialGap(un)+numtrials-1,windowGap(1):windowGap(1)+length(exptdata.xrange_psths)-1)=squeeze(TEp_trialSpden_fpON(:,un,:));
        TEp_bigMatrix(trialGap(un):trialGap(un)+numtrials-1,windowGap(2):windowGap(2)+length(exptdata.xrange_psths)-1)=squeeze(TEp_trialSpden_1stStimON(:,un,:));
        TEp_bigMatrix(trialGap(un):trialGap(un)+numtrials-1,windowGap(3):windowGap(3)+length(exptdata.xrange_psths)-1)=squeeze(TEp_trialSpden_2ndStimON(:,un,:));
        TEp_bigMatrix(trialGap(un):trialGap(un)+numtrials-1,windowGap(4):windowGap(4)+length(exptdata.xrange_psths)-1)=squeeze(TEp_trialSpden_ChoiceArrayON(:,un,:));
        TEp_bigMatrix(trialGap(un):trialGap(un)+numtrials-1,windowGap(5):windowGap(5)+length(exptdata.xrange_psths)-1)=squeeze(TEp_trialSpden_ResponseGOOD(:,un,:));
        TEp_bigMatrix(trialGap(un):trialGap(un)+numtrials-1,windowGap(6):windowGap(6)+length(exptdata.xrange_psths)-1)=squeeze(TEp_trialSpden_ResponseBAD(:,un,:));
        TEp_bigMatrix(trialGap(un):trialGap(un)+numtrials-1,windowGap(7):windowGap(7)+length(exptdata.xrange_psths)-1)=squeeze(TEp_trialSpden_reward(:,un,:));
    end
    save([exptdata.projectdir,tempname(1:end-13),'sessionData.mat'],'V4_bigMatrix','TEp_bigMatrix','-v7.3');
    
    
    
    %% 4.4 Export Data for population-level analysis (spikes)
    % V4 data
    % Create empty matrix for this session
    V4_bigMatrix_spikes=nan(numV4neurons*length(trialData),(length(exptdata.xrange_psths)*7)+20);
    % Paste in trial values
    V4_bigMatrix_spikes(:,1:20)=V4_bigMatrix(:,1:20);
    trialGap=1:numtrials:numV4neurons*numtrials;
    windowGap=21:length(exptdata.xrange_psths):length(exptdata.xrange_psths)*7;
    for un=1:numV4neurons
        V4_bigMatrix_spikes(trialGap(un):trialGap(un)+numtrials-1,windowGap(1):windowGap(1)+length(exptdata.xrange_psths)-1)=squeeze(V4_trialSpike_fpON(:,un,:));
        V4_bigMatrix_spikes(trialGap(un):trialGap(un)+numtrials-1,windowGap(2):windowGap(2)+length(exptdata.xrange_psths)-1)=squeeze(V4_trialSpike_1stStimON(:,un,:));
        V4_bigMatrix_spikes(trialGap(un):trialGap(un)+numtrials-1,windowGap(3):windowGap(3)+length(exptdata.xrange_psths)-1)=squeeze(V4_trialSpike_2ndStimON(:,un,:));
        V4_bigMatrix_spikes(trialGap(un):trialGap(un)+numtrials-1,windowGap(4):windowGap(4)+length(exptdata.xrange_psths)-1)=squeeze(V4_trialSpike_ChoiceArrayON(:,un,:));
        V4_bigMatrix_spikes(trialGap(un):trialGap(un)+numtrials-1,windowGap(5):windowGap(5)+length(exptdata.xrange_psths)-1)=squeeze(V4_trialSpike_ResponseGOOD(:,un,:));
        V4_bigMatrix_spikes(trialGap(un):trialGap(un)+numtrials-1,windowGap(6):windowGap(6)+length(exptdata.xrange_psths)-1)=squeeze(V4_trialSpike_ResponseBAD(:,un,:));
        V4_bigMatrix_spikes(trialGap(un):trialGap(un)+numtrials-1,windowGap(7):windowGap(7)+length(exptdata.xrange_psths)-1)=squeeze(V4_trialSpike_reward(:,un,:));
    end
    
    % TEp data
    % Create empty matrix for this session
    TEp_bigMatrix_spikes=nan(numTEpneurons*length(trialData),(length(exptdata.xrange_psths)*7)+20);
    % Paste in trial values
    TEp_bigMatrix_spikes(:,1:20)=TEp_bigMatrix(:,1:20);
    trialGap=1:numtrials:numTEpneurons*numtrials;
    windowGap=21:length(exptdata.xrange_psths):length(exptdata.xrange_psths)*7;
    for un=1:numTEpneurons
        TEp_bigMatrix_spikes(trialGap(un):trialGap(un)+numtrials-1,windowGap(1):windowGap(1)+length(exptdata.xrange_psths)-1)=squeeze(TEp_trialSpike_fpON(:,un,:));
        TEp_bigMatrix_spikes(trialGap(un):trialGap(un)+numtrials-1,windowGap(2):windowGap(2)+length(exptdata.xrange_psths)-1)=squeeze(TEp_trialSpike_1stStimON(:,un,:));
        TEp_bigMatrix_spikes(trialGap(un):trialGap(un)+numtrials-1,windowGap(3):windowGap(3)+length(exptdata.xrange_psths)-1)=squeeze(TEp_trialSpike_2ndStimON(:,un,:));
        TEp_bigMatrix_spikes(trialGap(un):trialGap(un)+numtrials-1,windowGap(4):windowGap(4)+length(exptdata.xrange_psths)-1)=squeeze(TEp_trialSpike_ChoiceArrayON(:,un,:));
        TEp_bigMatrix_spikes(trialGap(un):trialGap(un)+numtrials-1,windowGap(5):windowGap(5)+length(exptdata.xrange_psths)-1)=squeeze(TEp_trialSpike_ResponseGOOD(:,un,:));
        TEp_bigMatrix_spikes(trialGap(un):trialGap(un)+numtrials-1,windowGap(6):windowGap(6)+length(exptdata.xrange_psths)-1)=squeeze(TEp_trialSpike_ResponseBAD(:,un,:));
        TEp_bigMatrix_spikes(trialGap(un):trialGap(un)+numtrials-1,windowGap(7):windowGap(7)+length(exptdata.xrange_psths)-1)=squeeze(TEp_trialSpike_reward(:,un,:));
    end
    save([exptdata.projectdir,tempname(1:end-13),'sessionData_spikes.mat'],'V4_bigMatrix_spikes','TEp_bigMatrix_spikes','-v7.3');
    
    
    
    
    
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %% 5. Plot Average Spikes
    % V4 Neurons
    figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.05 0.1 0.9 0.8]); set(gca,'FontName','Arial')
    % StimA
    subplot(4,5,[1 2]) % V4 EXEMPLAR A
    hold on % stimulus onset
    a=nanmean(V4_condSpden_1stStimON(intersect(exptdata.dms500.stimA_1st,exptdata.dms500.activeCond),:,:),1); b=sum(a,1); c=sum(squeeze(b),1); % active stim1
    plot(exptdata.xrange_psths,c,'k-','LineWidth',0.5)
    a=nanmean(V4_condSpden_1stStimON(intersect(exptdata.dms500.stimA_1st,exptdata.dms500.passiveCond),:,:),1); b=sum(a,1); c=sum(squeeze(b),1); % passive stim1
    plot(exptdata.xrange_psths,c,'k-.','LineWidth',0.5)
    a=nanmean(V4_condSpden_2ndStimON(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimA_expectRepeat_actualRepeat),:,:),1); b=sum(a,1); c1=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c1,'r-','LineWidth',2)
    a=nanmean(V4_condSpden_2ndStimON(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimA_expectRepeat_actualRepeat),:,:),1); b=sum(a,1); c2=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c2,'r-.','LineWidth',2)
    a=nanmean(V4_condSpden_2ndStimON(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimA_expectAltern_actualRepeat),:,:),1); b=sum(a,1); c3=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c3,'b-','LineWidth',2)
    a=nanmean(V4_condSpden_2ndStimON(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimA_expectAltern_actualRepeat),:,:),1); b=sum(a,1); c4=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c4,'b-.','LineWidth',2)
    plot([0 0],[0 1],'k:'); grid on;
    legend('1stPres-Active','1stPres-Passive','2ndPresExpRep,ActRep-Active','2ndPresExpRep,ActRep-Passive','2ndPresExpAlt,ActRep-Active','2ndPresExpAlt,ActRep-Passive')
    title({tempname,'Responses to Repeated Exemplar A - V4'},'FontSize',14,'FontWeight','Bold','Interpreter','None');
    xlabel('Time From Stimulus Onset (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
    
    % StimB
    subplot(4,5,[6 7]) % V4 EXEMPLAR B
    hold on % stimulus onset
    a=nanmean(V4_condSpden_1stStimON(intersect(exptdata.dms500.stimB_1st,exptdata.dms500.activeCond),:,:),1); b=sum(a,1); c=sum(squeeze(b),1); % active stim1
    plot(exptdata.xrange_psths,c,'k-','LineWidth',0.5)
    a=nanmean(V4_condSpden_1stStimON(intersect(exptdata.dms500.stimB_1st,exptdata.dms500.passiveCond),:,:),1); b=sum(a,1); c=sum(squeeze(b),1); % passive stim1
    plot(exptdata.xrange_psths,c,'k-.','LineWidth',0.5)
    a=nanmean(V4_condSpden_2ndStimON(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimBstimB_expectRepeat_actualRepeat),:,:),1); b=sum(a,1); c1=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c1,'r-','LineWidth',2)
    a=nanmean(V4_condSpden_2ndStimON(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimBstimB_expectRepeat_actualRepeat),:,:),1); b=sum(a,1); c2=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c2,'r-.','LineWidth',2)
    a=nanmean(V4_condSpden_2ndStimON(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimBstimB_expectAltern_actualRepeat),:,:),1); b=sum(a,1); c3=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c3,'b-','LineWidth',2)
    a=nanmean(V4_condSpden_2ndStimON(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimBstimB_expectAltern_actualRepeat),:,:),1); b=sum(a,1); c4=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c4,'b-.','LineWidth',2)
    plot([0 0],[0 1],'k:'); grid on;
    legend('1stPres-Active','1stPres-Passive','2ndPresExpRep,ActRep-Active','2ndPresExpRep,ActRep-Passive','2ndPresExpAlt,ActRep-Active','2ndPresExpAlt,ActRep-Passive')
    title('Responses to Repeated Exemplar B - V4'); xlabel('Time From Stimulus Onset (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
    
    subplot(4,5,[11 12]) % TE EXEMPLAR A
    hold on % stimulus onset
    a=nanmean(TEp_condSpden_1stStimON(intersect(exptdata.dms500.stimA_1st,exptdata.dms500.activeCond),:,:),1); b=sum(a,1); c=sum(squeeze(b),1); % active stim1
    plot(exptdata.xrange_psths,c,'k-','LineWidth',0.5)
    a=nanmean(TEp_condSpden_1stStimON(intersect(exptdata.dms500.stimA_1st,exptdata.dms500.passiveCond),:,:),1); b=sum(a,1); c=sum(squeeze(b),1); % passive stim1
    plot(exptdata.xrange_psths,c,'k-.','LineWidth',0.5)
    a=nanmean(TEp_condSpden_2ndStimON(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimA_expectRepeat_actualRepeat),:,:),1); b=sum(a,1); c1=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c1,'r-','LineWidth',2)
    a=nanmean(TEp_condSpden_2ndStimON(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimA_expectRepeat_actualRepeat),:,:),1); b=sum(a,1); c2=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c2,'r-.','LineWidth',2)
    a=nanmean(TEp_condSpden_2ndStimON(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimA_expectAltern_actualRepeat),:,:),1); b=sum(a,1); c3=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c3,'b-','LineWidth',2)
    a=nanmean(TEp_condSpden_2ndStimON(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimA_expectAltern_actualRepeat),:,:),1); b=sum(a,1); c4=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c4,'b-.','LineWidth',2)
    plot([0 0],[0 1],'k:'); grid on; % ylim([0 0.5])
    legend('1stPres-Active','1stPres-Passive','2ndPresExpRep,ActRep-Active','2ndPresExpRep,ActRep-Passive','2ndPresExpAlt,ActRep-Active','2ndPresExpAlt,ActRep-Passive')
    title('Responses to Repeated Exemplar A - TE'); xlabel('Time From Stimulus Onset (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
    
    % StimB
    subplot(4,5,[16 17]) % TE EXEMPLAR B
    hold on % stimulus onset
    a=nanmean(TEp_condSpden_1stStimON(intersect(exptdata.dms500.stimB_1st,exptdata.dms500.activeCond),:,:),1); b=sum(a,1); c=sum(squeeze(b),1); % active stim1
    plot(exptdata.xrange_psths,c,'k-','LineWidth',0.5)
    a=nanmean(TEp_condSpden_1stStimON(intersect(exptdata.dms500.stimB_1st,exptdata.dms500.passiveCond),:,:),1); b=sum(a,1); c=sum(squeeze(b),1); % passive stim1
    plot(exptdata.xrange_psths,c,'k--','LineWidth',0.5)
    a=nanmean(TEp_condSpden_2ndStimON(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimBstimB_expectRepeat_actualRepeat),:,:),1); b=sum(a,1); c1=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c1,'r-','LineWidth',2)
    a=nanmean(TEp_condSpden_2ndStimON(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimBstimB_expectRepeat_actualRepeat),:,:),1); b=sum(a,1); c2=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c2,'r--','LineWidth',2)
    a=nanmean(TEp_condSpden_2ndStimON(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimBstimB_expectAltern_actualRepeat),:,:),1); b=sum(a,1); c3=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c3,'b-','LineWidth',2)
    a=nanmean(TEp_condSpden_2ndStimON(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimBstimB_expectAltern_actualRepeat),:,:),1); b=sum(a,1); c4=sum(squeeze(b),1);
    plot(exptdata.xrange_psths,c4,'b--','LineWidth',2)
    plot([0 0],[0 1],'k:'); grid on;
    legend('1stPres-Active','1stPres-Passive','2ndPresExpRep,ActRep-Active','2ndPresExpRep,ActRep-Passive','2ndPresExpAlt,ActRep-Active','2ndPresExpAlt,ActRep-Passive')
    title('Responses to Repeated Exemplar B - TE'); xlabel('Time From Stimulus Onset (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
    
    % All Neuron Panels
    subplot(4,5,[3 8]) % All neurons (V4)(aligned on 1stStimulusPresentation - allStimA)
    a=squeeze(nanmean(V4_condSpden_1stStimON(intersect(exptdata.dms500.stimA_1st,exptdata.dms500.activeCond),:,:),1));
    imagesc(a);
    xlabel('Time From Stimulus Onset (ms)'); ylabel('Neuron Number');
    title({'Responses of each V4 neuron in session','aligned on 1st Stim presentation'},'FontSize',8)
    cb = colorbar;
    cb.Label.String = 'Firing Rate (sp/s)';
    
    subplot(4,5,[4 9]) % All neurons (V4)(aligned on 1stStimulusPresentation - allStimA)
    a=squeeze(nanmean(V4_condSpden_2ndStimON(intersect(exptdata.dms500.stimA_1st,exptdata.dms500.activeCond),:,:),1));
    imagesc(a);
    xlabel('Time From Stimulus Onset (ms)'); ylabel('Neuron Number');
    title({'Responses of each V4 neuron in session','aligned on 2nd Stim presentation'},'FontSize',8)
    cb = colorbar;
    cb.Label.String = 'Firing Rate (sp/s)';
    
    subplot(4,5,[5 10]) % All neurons (V4)(aligned on 1stStimulusPresentation - allStimA)
    a=squeeze(nanmean(V4_trialSpden_ChoiceArrayON(intersect(exptdata.dms500.stimA_1st,exptdata.dms500.activeCond),:,:),1));
    imagesc(a);
    xlabel('Time From Stimulus Onset (ms)'); ylabel('Neuron Number');
    title({'Responses of each V4 neuron in session','aligned on Choice Array presentation'},'FontSize',8)
    cb = colorbar;
    cb.Label.String = 'Firing Rate (sp/s)';
    
    subplot(4,5,[13 18]) % All neurons (TEp)(aligned on 1stStimulusPresentation - allStimA)
    a=squeeze(nanmean(TEp_condSpden_1stStimON(intersect(exptdata.dms500.stimA_1st,exptdata.dms500.activeCond),:,:),1));
    imagesc(a);
    xlabel('Time From Stimulus Onset (ms)'); ylabel('Neuron Number');
    title({'Responses of each TEp neuron in session','aligned on 1st Stim presentation'},'FontSize',8)
    cb = colorbar;
    cb.Label.String = 'Firing Rate (sp/s)';
    
    subplot(4,5,[14 19]) % All neurons (TEp)(aligned on 1stStimulusPresentation - allStimA)
    a=squeeze(nanmean(TEp_condSpden_2ndStimON(intersect(exptdata.dms500.stimA_1st,exptdata.dms500.activeCond),:,:),1));
    imagesc(a);
    xlabel('Time From Stimulus Onset (ms)'); ylabel('Neuron Number');
    title({'Responses of each TEp neuron in session','aligned on 2nd Stim presentation'},'FontSize',8)
    cb = colorbar;
    cb.Label.String = 'Firing Rate (sp/s)';
    
    subplot(4,5,[15 20]) % All neurons (TEp)(aligned on choice - allStimA)
    a=squeeze(nanmean(TEp_trialSpden_ChoiceArrayON(intersect(exptdata.dms500.stimA_1st,exptdata.dms500.activeCond),:,:),1));
    imagesc(a);
    xlabel('Time From Stimulus Onset (ms)'); ylabel('Neuron Number');
    title({'Responses of each TEp neuron in session','aligned on Choice Array presentation'},'FontSize',8)
    cb = colorbar;
    cb.Label.String = 'Firing Rate (sp/s)';
    
    jpgfigname=[exptdata.figuredir500,tempname(1:end-13),'_summaryFigure1.jpg'];
    print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
    savefig([exptdata.figuredir500,'matlabFigFiles',filesep,tempname(1:end-13),'_summaryFigure1.fig'])
end