%% do_d500_es_rs_sessionQA;
% by AHB, started Jan 30, 2019
vers_ephys_es_rs='1.1; May 24, 2023';
% 1.0 - original version (Feb 4, 2019)
% 1.1 - modified to run on modern hardware/MATLAB version (HOMER_B550) (May 24, 2023)

%% SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clc; close all;
dbstop if error;
global exptdata
if ispc
    rootdir='C:\Users\ahbel\OneDrive - King''s College London\MATLAB';
else % MAC
    rootdir='~/Documents/MATLAB/';
end
addpath(userpath);
addpath(genpath([rootdir,filesep,'currentProjects',filesep,'proj_d500_es_rs']));
addpath(genpath([rootdir,filesep,'currentProjects',filesep,'commonProjectFunctions']));
addpath(genpath([rootdir,filesep,'Common_Functions']));
ephys_analysis_defaults;
exptdata.lastModified=datetime('today');
warning('off','MATLAB:MKDIR:DirectoryExists');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% COMMENT OUT THE ONES YOU DON'T WANT
%exptdata.monkeyname='Vortex';
%exptdata.monkeyname='Vulcan';
exptdata.monkeyname='both';

exptdata.reprocessQA=1;

%% Figure Selector
% Panel 1 - imagesc figures
% Panel 2 - average spikedensity functions
% Panel 3 - histogram of spike counts


%% =====================================================================================================================
% megamatrices AND behaviouralAnalysis must be complete before running this program
%% What this program will do:
% 1) Scroll through each session and generate imagesc
% 2)
% 3)


%% Data structure (<BRAIN-AREA>_megaMatrix<TASK-CODE>NoScrubs.mat)
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

% PARADIGM SPECIFIC STUFF (these are set by lsn_ephys_addTrialClassifers):
% 11) 1st stim presented (1=stimA; 2=stimB)
% 12) 2nd stim presented (1=stimA; 2=stimB)
% 13) expectation (1=expect repeat; 2=expect alternation; 3=no expectation)
% 14) actual (1=repeat, -1=alternation)
% 15) active (1) or passive (0)
% 16) right (1) or left (-1) for current trial
% 17) nothing
% 18) task number (300,500,600)
% 19) new block number (600 only)
% 20) nothing (WILL BECOME VISUAL RESPONSIVENESS WITH THIS FILE)
% 21...) lots of spike density stuff

%% =====================================================================================================================
exptdata.analysisName='D500_ES-RS_Study'; % used for savenames, figures, etc. (pick whatever you want; will be used for filenames)
exptdata.projectdir=[exptdata.analysisdir,exptdata.analysisName,filesep];

cprintf('*blue','*---------------------------*\n')
cprintf('*blue','| do_d500_es_rs_sessionQA.m |\n')
cprintf('*blue','*---------------------------*\n')
cprintf('*red',['Version: ',vers_ephys_es_rs,'\n'])
disp(['Data location:    ',exptdata.datalocation]);
disp(['Processed NEV dir: ',exptdata.processedDatadir]);
disp(['Project Directory: ',exptdata.projectdir]);
disp(['Monkey name:      ',exptdata.monkeyname]);
disp(['Study name:       ',exptdata.analysisName]);

%% =====================================================================================================================
% load spike density data
fprintf('\nLoading megamatrices and behavioural data...')
load([exptdata.projectdir,'DMS500_behavAnalysis.mat']); % contains BOTH monkeys behavioural data

%sessions = unique(V4_megaMatrix500NoScrubs(:,2));
if strcmp(exptdata.monkeyname,'Vortex')
    sessions = unique(evdata.session(evdata.animal==1)); % if Vortex
    load([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_megaMatrix.mat']);
elseif strcmp(exptdata.monkeyname,'Vulcan')
    sessions = unique(evdata.session(evdata.animal==2)); % if Vulcan
    load([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_megaMatrix.mat']);
else
    disp('combining monkey data')
    temp1=load([exptdata.projectdir,'Vortex_',exptdata.analysisName,'_megaMatrix.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
    temp2=load([exptdata.projectdir,'Vulcan_',exptdata.analysisName,'_megaMatrix.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
    V4_megaMatrix500NoScrubs=[temp1.V4_megaMatrix500NoScrubs; temp2.V4_megaMatrix500NoScrubs];
    TE_megaMatrix500NoScrubs=[temp1.TE_megaMatrix500NoScrubs; temp2.TE_megaMatrix500NoScrubs];
    clear temp*
    sessions = unique(evdata.session); %#ok<*NASGU>
end
fprintf('done.\n')


%% =====================================================================================================================
% load spiking data
fprintf('\nLoading spiketrain data...')
if strcmp(exptdata.monkeyname,'Vortex')
    sessions = unique(evdata.session(evdata.animal==1)); % if Vortex
    load([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_megaMatrix_spikes.mat']);
elseif strcmp(exptdata.monkeyname,'Vulcan')
    sessions = unique(evdata.session(evdata.animal==2)); % if Vulcan
    load([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_megaMatrix_spikes.mat']);
else
    disp('combining monkey data')
    temp1=load([exptdata.projectdir,'Vortex_',exptdata.analysisName,'_megaMatrix_spikes.mat'],'V4_megaMatrix500NoScrubs_spikes','TE_megaMatrix500NoScrubs_spikes');
    temp2=load([exptdata.projectdir,'Vulcan_',exptdata.analysisName,'_megaMatrix_spikes.mat'],'V4_megaMatrix500NoScrubs_spikes','TE_megaMatrix500NoScrubs_spikes');
    V4_megaMatrix500NoScrubs_spikes=[temp1.V4_megaMatrix500NoScrubs_spikes; temp2.V4_megaMatrix500NoScrubs_spikes];
    TE_megaMatrix500NoScrubs_spikes=[temp1.TE_megaMatrix500NoScrubs_spikes; temp2.TE_megaMatrix500NoScrubs_spikes];
    clear temp*
    sessions = unique(evdata.session);
end
fprintf('done.\n')


%% Begin Brain Loops
for brainArea=1:2
    if brainArea==1
        AllspikeData=V4_megaMatrix500NoScrubs;
        AllspikeData_spikes=V4_megaMatrix500NoScrubs_spikes;
        labelBrain='V4';       
    else
        AllspikeData=TE_megaMatrix500NoScrubs;
        AllspikeData_spikes=TE_megaMatrix500NoScrubs_spikes;
        labelBrain='TE';
    end
    
    %% Classify Neurons
    % This will perform a paired two-tailed t-test on mean(baseline) (-200-0 ms prior to
    % STIM 1 ONSET) compared to mean 50-250 ms after onset across all conditions
    % If baseline < sensory (p<0.05) => (x,20)=1 (visual neuron, excitatory)
    % If baseline > sensory (p<0.05) => (x,20)=-1 (visual neuron, suppressed)
    % If baseline = sensory (n.s.) => (x,20)=0 (non-responsive neuron)
    
    disp('Classifying neurons according to visual responsiveness...')
    for sess=1:length(sessions)
        clear indxN predictrz numNeuronsSession tt
        numNeuronsSession=max(AllspikeData(AllspikeData(:,2)==sessions(sess),10)); % how many neurons in this unit (this code is a bit of a hack but does the trick)
        tempVisualBad=0; tempVisualGood=0; tempVisualBoring=0;
        if numNeuronsSession > 0
            for nn=1:numNeuronsSession
                indxN=find(AllspikeData(:,2) == sessions(sess) & AllspikeData(:,3) == nn);  % row pointer for all trials for this neuron (nn) in this session (sess)
                baselineWindow=772+50:772+250;
                sensoryWindow =772+250+50:772+250+250; % (-250ms is start of window, so -250+300:-250+500 = +50 to +250ms)
                tempH=ttest(mean(AllspikeData(indxN,baselineWindow),2) , mean(AllspikeData(indxN,sensoryWindow),2));
                if tempH==1
                    if mean(mean(AllspikeData(indxN,baselineWindow),2)) > mean(mean(AllspikeData(indxN,sensoryWindow),2))
                        AllspikeData(indxN,20)=-1; % visual suppressed
                        tempVisualBad=tempVisualBad+1;
                    else
                        AllspikeData(indxN,20)=1; % visual excitatory
                        tempVisualGood=tempVisualGood+1;
                    end
                else
                    AllspikeData(indxN,20)=0;
                    tempVisualBoring=tempVisualBoring+1;
                end
            end
        end
        disp(['...Session: ',num2str(sess),' -- (E:',num2str(tempVisualGood),' / S:',num2str(tempVisualBad),' / nR:',num2str(tempVisualBoring),' / Total:',num2str(numNeuronsSession),')'])
    end
    clear temp*
   
    %FILE,ARRAY,SHEET,RANGE
    
    %% New version - One figure, One neuron
    numMonkeys=length(unique(AllspikeData(:,1)));
    for monk=1:numMonkeys
        sessionsMonk=unique(AllspikeData(AllspikeData(:,1)==monk,2));
        for sess=1:length(sessionsMonk)
            neuronsSession=unique(AllspikeData(AllspikeData(:,1)==monk & AllspikeData(:,2)==sessionsMonk(sess),3));
            for nn=1:length(neuronsSession)
                if exist([exptdata.QA_figureDir,'matlabFigFiles',filesep,exptdata.allMonkeyNames{monk},'_',labelBrain,'_QA_S',num2str(sessionsMonk(sess)),'_Neuron',num2str(nn),'_FullSummary.fig'],'file') && exptdata.reprocessQA~=1
                    fprintf(['*** Already processed ',exptdata.allMonkeyNames{monk},'_',labelBrain,'_QA_S',num2str(sessionsMonk(sess)),'_Neuron',num2str(nn),'. Skipping...\n'])
                    continue
                end
                
                h=figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial') % new figure
                
                %% Panel 1 - Average Spike Density Function
                subplot(5,1,1); hold on
                tempPointer=AllspikeData(:,1)==monk & AllspikeData(:,2)==sessionsMonk(sess) & AllspikeData(:,3)==neuronsSession((nn));
                imagesc(AllspikeData(tempPointer,[772:1522 1523:2273 2274:3024]),[0 ceil(max(max(AllspikeData(tempPointer,[772:1522 1523:2273 2274:3024]))))])
                plot([750 750],[0 sum(tempPointer)],'w-','LineWidth',3)
                plot([1500 1500],[0 sum(tempPointer)],'w-','LineWidth',3)
                plot([250 250],[0 sum(tempPointer)],'w:','LineWidth',2)
                plot([1000 1000],[0 sum(tempPointer)],'w:','LineWidth',2)
                plot([1750 1750],[0 sum(tempPointer)],'w:','LineWidth',2)
                ylabel('Trial Number'); xlabel('Time from Event ( ms, 1st Stim | 2nd Stim | Choice Array )');
                if mean(AllspikeData(tempPointer,20))>0
                    tempNeuronLabel='(Visual, Excitatory)';
                elseif mean(AllspikeData(tempPointer,20))<0
                    tempNeuronLabel='(Visual, Suppressed)';
                else
                    tempNeuronLabel='(Non-Responsive)';
                end
                title([exptdata.allMonkeyNames{monk},' (',labelBrain,'): Session No.',num2str(sessionsMonk(sess)),' / Neuron No.',num2str(neuronsSession(nn)),' ',tempNeuronLabel],'FontSize',14)
                c = colorbar; colormap jet
                c.Label.String = 'Firing Rate (sp/s)';
                c.Limits=[0 ceil(max(max(AllspikeData(tempPointer,[772:1522 1523:2273 2274:3024]))))];
                xlim([0 2250])
                                
                %% Panel 2
                subplot(5,1,2); hold on
                clear tempNeuronData tempDataPassive tempDataActive tempScale avgSpden h1 h2
                tempNeuronData=AllspikeData(tempPointer,:);
                [avgSpden, NavgSpden]=ephys_prepAverageSpden(tempNeuronData);
                tempScale=10;
                try
                    tempDataPassive=[avgSpden.passive.RR(:,772:3024);avgSpden.passive.RA(:,772:3024);avgSpden.passive.AR(:,772:3024);avgSpden.passive.AA(:,772:3024)];
                    h1=fillsteplot_ahb(tempDataPassive,1,'-',{'b'});
                    clear tempDataPassive
                    tempScale=max(max(nanmean([avgSpden.passive.RR(:,772:3024);avgSpden.passive.RA(:,772:3024);avgSpden.passive.AR(:,772:3024);avgSpden.passive.AA(:,772:3024)])));
                catch
                    h1=plot([0 2250],[10 10],'b--');
                    clear tempDataPassive
                end
                try
                    tempDataActive=[avgSpden.active.RR(:,772:3024);avgSpden.active.RA(:,772:3024);avgSpden.active.AR(:,772:3024);avgSpden.active.AA(:,772:3024)];
                    h2=fillsteplot_ahb(tempDataActive,1,'-',{'r'});
                    tempScale=max(max(nanmean(tempDataActive)));
                    clear tempDataActive
                catch
                    h2=plot([0 2250],[10 10],'r--');
                    clear tempDataActive
                end
                if tempScale<1
                    tempScale=10;
                end
                plot([750 750],[0 ceil(tempScale*1.5)],'k-','LineWidth',3)
                plot([1500 1500],[0 ceil(tempScale*1.5)],'k-','LineWidth',3)
                plot([250 250],[0 ceil(tempScale*1.5)],'k:','LineWidth',2)
                plot([1000 1000],[0 ceil(tempScale*1.5)],'k:','LineWidth',2)
                plot([1750 1750],[0 ceil(tempScale*1.5)],'k:','LineWidth',2)
                legend([h1 h2],'Passive','Active'); xlim([0 2250]); ylim([0 ceil(tempScale*1.5)]);
                ylabel('Average Firing Rate (sp/s)'); xlabel('Time from Event ( ms, 1st Stim | 2nd Stim | Choice Array )');
                title('Average Spike Density Functions','FontSize',14)
                
                %% Panel 3
                subplot(5,1,3); hold on
                tempNeuronData=AllspikeData(tempPointer,[773:1522 1524:2273 2274:3024]);
                bar(1:size(tempNeuronData,1),sum(tempNeuronData,2))
                ylabel({'Sum of Spike Density','(sp/s/trial)'}); xlabel('Trial Number');
                if mean(AllspikeData(tempPointer,20))>0
                    tempNeuronLabel='(Visual, Excitatory)';
                elseif mean(AllspikeData(tempPointer,20))<0
                    tempNeuronLabel='(Visual, Suppressed)';
                else
                    tempNeuronLabel='(Non-Responsive)';
                end
                title('Sum of Spike Density (ordered according to Trial Number','FontSize',14)
                                
                %% Panel 4
                subplot(5,1,4); hold on
                hist(sum(tempNeuronData,2),1:100:50000)
                ylabel('Trial Count'); xlabel('Sum of Spike Density/Trial');
                title('Histogram of Sum of Spike Density','FontSize',14)
                
                %% Panel 5
                subplot(5,1,5); hold on
                clear tempPointer
                tempPointer=AllspikeData_spikes(:,1)==monk & AllspikeData_spikes(:,2)==sessionsMonk(sess) & AllspikeData_spikes(:,3)==neuronsSession((nn));
                tempNeuronDataSpikes=AllspikeData_spikes(tempPointer,[773:1522 1524:2273 2274:3024]);
                hist(sum(tempNeuronDataSpikes,2),100)
                ylabel('Trial Count'); xlabel('Sum of Spikes/Trial');
                
                title('Histogram of Sum of Spikes','FontSize',14)
                clear tempNeuronLabel tempScale
                
                %% Save and Print
                jpgfigname=[exptdata.QA_figureDir,exptdata.allMonkeyNames{monk},'_',labelBrain,'_QA_S',num2str(sessionsMonk(sess)),'_Neuron',num2str(nn),'_FullSummary.jpg'];
                print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
                savefig([exptdata.QA_figureDir,'matlabFigFiles',filesep,exptdata.allMonkeyNames{monk},'_',labelBrain,'_QA_S',num2str(sessionsMonk(sess)),'_Neuron',num2str(nn),'_FullSummary.fig'])
                
                clear temp*
            end % end neuron
            close all
        end % end session
        close all % close all figures once per session
    end % end monkey
    toc
end % end brain area
return


























%% Make Pretty Pictures
% This section will be ad hoc (with the assumption that every panel will be repeated for each brain area)
% Data structure:
% 21 771;... aligned on FP on (1)
% 772 1522;... aligned on 1stStim on (2)
% 1523 2273;... aligned on 2ndStim on (3)
% 2274 3024;... aligned on choice array (4)
% 3025 3775;... aligned on correct choice (5)
% 3776 4526;... aligned on incorrect choice (6)
% 4528 5277]; % aligned on reward (7)

%% QA Figure 1 - imagesc, rough plots of all neural activity
if ~isempty(find(ismember(includeFigs,1)==1, 1))
    tic
    numMonkeys=length(unique(AllspikeData(:,1)));
    for monk=1:numMonkeys
        sessionsMonk=unique(AllspikeData(AllspikeData(:,1)==monk,2));
        for sess=1:length(sessionsMonk)
            neuronsSession=unique(AllspikeData(AllspikeData(:,1)==monk & AllspikeData(:,2)==sessionsMonk(sess),3));
            tempPanelCount=0; tempFigCount=1;
            h=figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial') % new figure
            for nn=1:length(neuronsSession)
                tempPanelCount=tempPanelCount+1;
                if tempPanelCount==6 % allow 10 plots per figure
                    jpgfigname=[exptdata.QA_figureDir,exptdata.allMonkeyNames{monk},'_',labelBrain,'_QA_S',num2str(sessionsMonk(sess)),'_F',num2str(tempFigCount),'_imagesc.jpg'];
                    print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
                    savefig([exptdata.QA_figureDir,'matlabFigFiles',filesep,exptdata.allMonkeyNames{monk},'_',labelBrain,'_QA_S',num2str(sessionsMonk(sess)),'_F',num2str(tempFigCount),'_imagesc.fig'])
                    h=figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial') % new figure
                    tempPanelCount=1; tempFigCount=tempFigCount+1;
                end
                subplot(5,1,tempPanelCount); hold on
                tempPointer=AllspikeData(:,1)==monk & AllspikeData(:,2)==sessionsMonk(sess) & AllspikeData(:,3)==neuronsSession((nn));
                imagesc(AllspikeData(tempPointer,[772:1522 1523:2273 2274:3024]),[0 ceil(max(max(AllspikeData(tempPointer,[772:1522 1523:2273 2274:3024]))))])
                plot([750 750],[0 sum(tempPointer)],'w-','LineWidth',3)
                plot([1500 1500],[0 sum(tempPointer)],'w-','LineWidth',3)
                plot([250 250],[0 sum(tempPointer)],'w:','LineWidth',2)
                plot([1000 1000],[0 sum(tempPointer)],'w:','LineWidth',2)
                plot([1750 1750],[0 sum(tempPointer)],'w:','LineWidth',2)
                ylabel('Trial Number'); xlabel('Time from Event ( ms, 1st Stim | 2nd Stim | Choice Array )');
                if mean(AllspikeData(tempPointer,20))>0
                    tempNeuronLabel='(Visual, Excitatory)';
                elseif mean(AllspikeData(tempPointer,20))<0
                    tempNeuronLabel='(Visual, Suppressed)';
                else
                    tempNeuronLabel='(Non-Responsive)';
                end
                title([exptdata.allMonkeyNames{monk},' (',labelBrain,'): Session No.',num2str(sessionsMonk(sess)),' / Neuron No.',num2str(neuronsSession(nn)),' ',tempNeuronLabel],'FontSize',14)
                c = colorbar; colormap jet
                c.Label.String = 'Firing Rate (sp/s)';
                c.Limits=[0 ceil(max(max(AllspikeData(tempPointer,[772:1522 1523:2273 2274:3024]))))];
                xlim([0 2250])
                clear tempNeuronLabel
            end
            close all
        end
        clear temp*
    end
    toc
end

%% QA Figure 2 - spike density functions
if ~isempty(find(ismember(includeFigs,2)==1, 1))
    tic
    numMonkeys=length(unique(AllspikeData(:,1)));
    for monk=1:numMonkeys
        sessionsMonk=unique(AllspikeData(AllspikeData(:,1)==monk,2));
        for sess=1:length(sessionsMonk)
            neuronsSession=unique(AllspikeData(AllspikeData(:,1)==monk & AllspikeData(:,2)==sessionsMonk(sess),3));
            tempPanelCount=0; tempFigCount=1;
            h=figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial') % new figure
            for nn=1:length(neuronsSession)
                tempPanelCount=tempPanelCount+1;
                if tempPanelCount==6 % allow 10 plots per figure
                    jpgfigname=[exptdata.QA_figureDir,exptdata.allMonkeyNames{monk},'_',labelBrain,'_QA_S',num2str(sessionsMonk(sess)),'_F',num2str(tempFigCount),'_spdens.jpg'];
                    print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
                    savefig([exptdata.QA_figureDir,'matlabFigFiles',filesep,exptdata.allMonkeyNames{monk},'_',labelBrain,'_QA_S',num2str(sessionsMonk(sess)),'_F',num2str(tempFigCount),'_spdens.fig'])
                    h=figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial') % new figure
                    tempPanelCount=1; tempFigCount=tempFigCount+1;
                end
                subplot(5,1,tempPanelCount); hold on
                tempPointer=AllspikeData(:,1)==monk & AllspikeData(:,2)==sessionsMonk(sess) & AllspikeData(:,3)==neuronsSession((nn));
                clear tempNeuronData tempDataPassive tempDataActive tempScale tempNeuronLabel avgSpden h1 h2
                tempNeuronData=AllspikeData(tempPointer,:);
                [avgSpden, NavgSpden]=lsn_ephys_prepAverageSpden(tempNeuronData);
                tempScale=10;
                try
                    tempDataPassive=[avgSpden.passive.RR(:,772:3024);avgSpden.passive.RA(:,772:3024);avgSpden.passive.AR(:,772:3024);avgSpden.passive.AA(:,772:3024)];
                    h1=fillsteplot_ahb(tempDataPassive,1,'-',{'b'});
                    clear tempDataPassive
                    tempScale=max(max(nanmean([avgSpden.passive.RR(:,772:3024);avgSpden.passive.RA(:,772:3024);avgSpden.passive.AR(:,772:3024);avgSpden.passive.AA(:,772:3024)])));
                catch
                    h1=plot([0 2250],[10 10],'b--');
                    clear tempDataPassive
                end
                try
                    tempDataActive=[avgSpden.active.RR(:,772:3024);avgSpden.active.RA(:,772:3024);avgSpden.active.AR(:,772:3024);avgSpden.active.AA(:,772:3024)];
                    h2=fillsteplot_ahb(tempDataActive,1,'-',{'r'});
                    tempScale=max(max(nanmean(tempDataActive)));
                    clear tempDataActive
                catch
                    h2=plot([0 2250],[10 10],'r--');
                    clear tempDataActive
                end
                if tempScale<1
                    tempScale=10;
                end
                plot([750 750],[0 ceil(tempScale*1.5)],'k-','LineWidth',3)
                plot([1500 1500],[0 ceil(tempScale*1.5)],'k-','LineWidth',3)
                plot([250 250],[0 ceil(tempScale*1.5)],'k:','LineWidth',2)
                plot([1000 1000],[0 ceil(tempScale*1.5)],'k:','LineWidth',2)
                plot([1750 1750],[0 ceil(tempScale*1.5)],'k:','LineWidth',2)
                legend([h1 h2],'Passive','Active'); xlim([0 2250]); ylim([0 ceil(tempScale*1.5)]);
                ylabel('Average Firing Rate (sp/s)'); xlabel('Time from Event ( ms, 1st Stim | 2nd Stim | Choice Array )');
                if mean(AllspikeData(tempPointer,20))>0
                    tempNeuronLabel='(Visual, Excitatory)';
                elseif mean(AllspikeData(tempPointer,20))<0
                    tempNeuronLabel='(Visual, Suppressed)';
                else
                    tempNeuronLabel='(Non-Responsive)';
                end
                title([exptdata.allMonkeyNames{monk},' (',labelBrain,'): Session No.',num2str(sessionsMonk(sess)),' / Neuron No.',num2str(neuronsSession(nn)),' ',tempNeuronLabel],'FontSize',14)
                clear tempNeuronLabel tempScale
            end % neuron
            close all
        end % session
        clear temp*
    end % monkey
    toc
end

%% QA Figure 3 - bar graphs
if ~isempty(find(ismember(includeFigs,3)==1, 1))
    tic
    numMonkeys=length(unique(AllspikeData(:,1)));
    for monk=1:numMonkeys
        sessionsMonk=unique(AllspikeData(AllspikeData(:,1)==monk,2));
        for sess=1:length(sessionsMonk)
            neuronsSession=unique(AllspikeData(AllspikeData(:,1)==monk & AllspikeData(:,2)==sessionsMonk(sess),3));
            tempPanelCount=0; tempFigCount=1;
            h=figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial') % new figure
            for nn=1:length(neuronsSession)
                tempPanelCount=tempPanelCount+1;
                if tempPanelCount==6 % allow 10 plots per figure
                    jpgfigname=[exptdata.QA_figureDir,exptdata.allMonkeyNames{monk},'_',labelBrain,'_QA_S',num2str(sessionsMonk(sess)),'_F',num2str(tempFigCount),'_barSumSpikes.jpg'];
                    print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
                    savefig([exptdata.QA_figureDir,'matlabFigFiles',filesep,exptdata.allMonkeyNames{monk},'_',labelBrain,'_QA_S',num2str(sessionsMonk(sess)),'_F',num2str(tempFigCount),'_barSumSpikes.fig'])
                    h=figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial') % new figure
                    tempPanelCount=1; tempFigCount=tempFigCount+1;
                end
                subplot(5,1,tempPanelCount); hold on
                clear tempNeuronData tempDataPassive tempDataActive tempScale tempNeuronLabel tempPointer
                tempPointer=AllspikeData(:,1)==monk & AllspikeData(:,2)==sessionsMonk(sess) & AllspikeData(:,3)==neuronsSession((nn));
                tempNeuronData=AllspikeData(tempPointer,[773:1522 1524:2273 2274:3024]);
                bar(1:size(tempNeuronData,1),sum(tempNeuronData,2))
                ylabel({'Sum of Spike Density','(sp/s/trial)'}); xlabel('Trial Number');
                if mean(AllspikeData(tempPointer,20))>0
                    tempNeuronLabel='(Visual, Excitatory)';
                elseif mean(AllspikeData(tempPointer,20))<0
                    tempNeuronLabel='(Visual, Suppressed)';
                else
                    tempNeuronLabel='(Non-Responsive)';
                end
                title([exptdata.allMonkeyNames{monk},' (',labelBrain,'): Session No.',num2str(sessionsMonk(sess)),' / Neuron No.',num2str(neuronsSession(nn)),' ',tempNeuronLabel],'FontSize',14)
                clear tempNeuronLabel tempScale
            end % neuron
            close all
        end % session
        clear temp*
    end % monkey
    toc
end

%% QA Figure 4 - bar graphs (per session)
if ~isempty(find(ismember(includeFigs,4)==1, 1))
    tic
    numMonkeys=length(unique(AllspikeData(:,1)));
    for monk=1:numMonkeys
        sessionsMonk=unique(AllspikeData(AllspikeData(:,1)==monk,2));
        for sess=1:length(sessionsMonk)
            neuronsSession=unique(AllspikeData(AllspikeData(:,1)==monk & AllspikeData(:,2)==sessionsMonk(sess),3));
            tempPanelCount=0; tempFigCount=1;
            h=figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial') % new figure
            for nn=1:length(neuronsSession)
                tempPanelCount=tempPanelCount+1;
                if tempPanelCount==6 % allow 10 plots per figure
                    jpgfigname=[exptdata.QA_figureDir,exptdata.allMonkeyNames{monk},'_',labelBrain,'_QA_S',num2str(sessionsMonk(sess)),'_F',num2str(tempFigCount),'_barSumSpikes.jpg'];
                    print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
                    savefig([exptdata.QA_figureDir,'matlabFigFiles',filesep,exptdata.allMonkeyNames{monk},'_',labelBrain,'_QA_S',num2str(sessionsMonk(sess)),'_F',num2str(tempFigCount),'_barSumSpikes.fig'])
                    h=figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial') % new figure
                    tempPanelCount=1; tempFigCount=tempFigCount+1;
                end
                subplot(5,1,tempPanelCount); hold on
                clear tempNeuronData tempDataPassive tempDataActive tempScale tempNeuronLabel tempPointer
                tempPointer=AllspikeData(:,1)==monk & AllspikeData(:,2)==sessionsMonk(sess) & AllspikeData(:,3)==neuronsSession((nn));
                tempNeuronData=AllspikeData(tempPointer,[773:1522 1524:2273 2274:3024]);
                hist(sum(tempNeuronData,2),1:100:50000)
                ylabel('Trial Count'); xlabel('Sum of Spike Density/Trial');
                if mean(AllspikeData(tempPointer,20))>0
                    tempNeuronLabel='(Visual, Excitatory)';
                elseif mean(AllspikeData(tempPointer,20))<0
                    tempNeuronLabel='(Visual, Suppressed)';
                else
                    tempNeuronLabel='(Non-Responsive)';
                end
                title([exptdata.allMonkeyNames{monk},' (',labelBrain,'): Session No.',num2str(sessionsMonk(sess)),' / Neuron No.',num2str(neuronsSession(nn)),' ',tempNeuronLabel],'FontSize',14)
                clear tempNeuronLabel tempScale
            end % neuron
            close all
        end % session
        clear temp*
    end % monkey
    toc
end

%% QA Figure 5 - bar graphs SPIKES (per session)
if ~isempty(find(ismember(includeFigs,5)==1, 1))
    tic
    numMonkeys=length(unique(AllspikeData_spikes(:,1)));
    for monk=1:numMonkeys
        sessionsMonk=unique(AllspikeData_spikes(AllspikeData_spikes(:,1)==monk,2));
        for sess=1:length(sessionsMonk)
            neuronsSession=unique(AllspikeData_spikes(AllspikeData_spikes(:,1)==monk & AllspikeData_spikes(:,2)==sessionsMonk(sess),3));
            tempPanelCount=0; tempFigCount=1;
            h=figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial') % new figure
            for nn=1:length(neuronsSession)
                tempPanelCount=tempPanelCount+1;
                if tempPanelCount==6 % allow 5 plots per figure
                    jpgfigname=[exptdata.QA_figureDir,exptdata.allMonkeyNames{monk},'_',labelBrain,'_QA_S',num2str(sessionsMonk(sess)),'_F',num2str(tempFigCount),'_barSumSpikeTrain.jpg'];
                    print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
                    savefig([exptdata.QA_figureDir,'matlabFigFiles',filesep,exptdata.allMonkeyNames{monk},'_',labelBrain,'_QA_S',num2str(sessionsMonk(sess)),'_F',num2str(tempFigCount),'_barSumSpikeTrain.fig'])
                    h=figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial') % new figure
                    tempPanelCount=1; tempFigCount=tempFigCount+1;
                end
                subplot(5,1,tempPanelCount); hold on
                clear tempNeuronData tempDataPassive tempDataActive tempScale tempNeuronLabel tempPointer
                tempPointer=AllspikeData_spikes(:,1)==monk & AllspikeData_spikes(:,2)==sessionsMonk(sess) & AllspikeData_spikes(:,3)==neuronsSession((nn));
                tempNeuronData=AllspikeData_spikes(tempPointer,[773:1522 1524:2273 2274:3024]);
                hist(sum(tempNeuronData,2),100)
                ylabel('Trial Count'); xlabel('Sum of Spikes/Trial');
                if mean(AllspikeData_spikes(tempPointer,20))>0
                    tempNeuronLabel='(Visual, Excitatory)';
                elseif mean(AllspikeData_spikes(tempPointer,20))<0
                    tempNeuronLabel='(Visual, Suppressed)';
                else
                    tempNeuronLabel='(Non-Responsive)';
                end
                title([exptdata.allMonkeyNames{monk},' (',labelBrain,'): Session No.',num2str(sessionsMonk(sess)),' / Neuron No.',num2str(neuronsSession(nn)),' ',tempNeuronLabel],'FontSize',14)
                clear tempNeuronLabel tempScale
            end % neuron
            close all
        end % session
        clear temp*
    end % monkey
    toc
end