% do_d500_es_rs_makePrettyFigures_500;
% by AHB, started Jan 30, 2019
vers_ephys_es_rs='1.0; Jan 30, 2019';
% 1.0 - original version (Jan 30, 2019)
% 1.1 - added ability to screen neurons according to quality and manual classifications (Feb 25, 2019)

global exptdata
addpath(userpath);
ephys_analysis_defaults;
load([exptdata.analysisdir,filesep,exptdata.analysisName,'_exptdata.mat'],'exptdata');
filterNeurons=0;

% COMMENT OUT THE ONES YOU DON'T WANT 
%exptdata.monkeyname='Vortex';
%exptdata.monkeyname='Vulcan';
exptdata.monkeyname='both';

% =====================================================================================================================
% megamatrices AND behaviouralAnalysis must be complete before running this program
% What this program will do:
% 1) Collect SpikeDensity Functions for all sessions, all neurons, all monkeys
% 2) Do quick identification of visual vs. non-visual neurons (two-tailed, to capture neurons that are suppressed)
% 3) generate some user-friendly pictures


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

% PARADIGM SPECIFIC STUFF (these are set by ephys_addTrialClassifers):
% 11) 1st stim presented (1=stimA; 2=stimB)
% 12) 2nd stim presented (1=stimA; 2=stimB)
% 13) expectation (1=expect repeat; 2=expect alternation; 3=no expectation)
% 14) actual (1=repeat, -1=alternation)
% 15) active (1) or passive (0)
% 16) right (1) or left (-1) for current trial
% 17) nothing
% 18) task number (300,500,600)
% 19) new block number (600 only)
% 20) Auto Classification
% 21) Manual Classification (1=excitatory, -1=suppressed, 0=non-responsive)
% 22) Quality Classification (0=shit, problematic, or non-responsive, 1=meh, 2=ok, 3=awesome)
% 23) stim1
% 24) stim2
% 25) stim3
% 26) stim4
% 27...) lots of spike density stuff

% =====================================================================================================================
exptdata.projectdir=[exptdata.analysisdir,filesep,exptdata.analysisName,filesep];
exptdata.figuredir=[exptdata.projectdir,'figures',filesep,'d500_es_rs',filesep]; mkdir(exptdata.figuredir);

clc
fprintf('<strong>*===================================*</strong>\n')
fprintf('<strong>| do_d500_es_rs_makePrettyFigures.m |</strong>\n')
fprintf('<strong>*===================================*</strong>\n')
fprintf(['Version:              ',vers_ephys_es_rs,'\n'])
disp(['Study name:           ',exptdata.analysisName]);
disp(['Data location:        ',exptdata.datalocation]);
disp(['Processed NEV dir:    ',exptdata.processedDatadir]);
disp(['Project Directory:    ',exptdata.projectdir]);
disp(['Figure Directory:     ',exptdata.figuredir]);

% =====================================================================================================================
%% 1. Load behavioural and neurophysiological data
fprintf('\nLoading megamatrices and behavioural data...')
load([exptdata.projectdir,'DMS500_behavAnalysis.mat']); % contains BOTH monkeys behavioural data

%sessions = unique(V4_megaMatrix500NoScrubs(:,2));
if strcmp(exptdata.monkeyname,'Vortex')
    sessions = unique(evdata.session(evdata.animal==1)); % if Vortex
    load([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_megaMatrix_synced.mat']);
elseif strcmp(exptdata.monkeyname,'Vulcan')
    sessions = unique(evdata.session(evdata.animal==2)); % if Vulcan
    load([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_megaMatrix_synced.mat']);
else
    disp('combining monkey data')
    temp1=load([exptdata.projectdir,'Vortex_',exptdata.analysisName,'_megaMatrix.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
    temp2=load([exptdata.projectdir,'Vulcan_',exptdata.analysisName,'_megaMatrix.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
    V4_megaMatrix500NoScrubs=[temp1.V4_megaMatrix500NoScrubs; temp2.V4_megaMatrix500NoScrubs];
    TE_megaMatrix500NoScrubs=[temp1.TE_megaMatrix500NoScrubs; temp2.TE_megaMatrix500NoScrubs];
    clear temp*
    sessions = unique(evdata.session);
end
fprintf('done.\n')

%% 2. Begin Brain Loops
for brainArea=1:2
    if brainArea==1
        AllspikeData=V4_megaMatrix500NoScrubs;
        AllspikeData_unfiltered=AllspikeData;
        labelBrain='V4';
    else
        AllspikeData=TE_megaMatrix500NoScrubs;
        labelBrain='TE';
        AllspikeData_unfiltered=AllspikeData;
    end
        
    %% 2.1 Classify Neurons
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
        disp(['...Session: ',num2str(sess),' -- (E:',num2str(tempVisualGood),' / S:',num2str(tempVisualBad),' / nR:',num2str(tempVisualBoring),' / Total:',num2str(numNeuronsSession),')'])
    end
    clear temp*
          
    %% 2.2 Make Pretty Pictures
    % This section will be ad hoc (with the assumption that every panel will be repeated for each brain area)
    % Data structure:
    % 21 771;... aligned on FP on (1)
    % 772 1522;... aligned on 1stStim on (2)
    % 1523 2273;... aligned on 2ndStim on (3)
    % 2274 3024;... aligned on choice array (4)
    % 3025 3775;... aligned on correct choice (5)
    % 3776 4526;... aligned on incorrect choice (6)
    % 4528 5277]; % aligned on reward (7)
     
    % ========================================================================================================
    % Select visual data only...
    VISUALspikeData=AllspikeData(AllspikeData(:,20)==1,:);

    % Add filtering options:
    if filterNeurons==1
        cprintf('*red','FILTERING NEURONS!!!\n')
        % Filter Visual Only (1=excitatory, -1=suppressed, 0=non-responsive)
        VISUALspikeData=VISUALspikeData(VISUALspikeData(:,21)==1,:); % visual excitatory
        %AllspikeData=AllspikeData(AllspikeData(:,21)~=0,:); % all visual
        
        % Quality (0=shit, problematic, or non-responsive, 1=meh, 2=ok, 3=awesome)
        VISUALspikeData=VISUALspikeData(VISUALspikeData(:,22)>1,:); % visual excitatory
    end
 
    [avgSpden, NavgSpden]=ephys_prepAverageSpden(VISUALspikeData);
    
   
    %% 2.2.1 Figure 1 - Raw Firing Rates
    figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial')
    subplot(2,2,1); hold on % passive
    h1=plot(-250:500,nanmean([avgSpden.passive.RR(:,772:1522);avgSpden.passive.RA(:,772:1522);avgSpden.passive.AR(:,772:1522);avgSpden.passive.AA(:,772:1522)],1),'k-');
    h2=plot(-250:500,nanmean(avgSpden.passive.RR(:,1523:2273),1),'r-');
    h3=plot(-250:500,nanmean(avgSpden.passive.RA(:,1523:2273),1),'b-');
    h4=plot(-250:500,nanmean(avgSpden.passive.AR(:,1523:2273),1),'g-');
    h5=plot(-250:500,nanmean(avgSpden.passive.AA(:,1523:2273),1),'m-');
    h6=plot(-250:500,nanmean(avgSpden.passive.ER(:,1523:2273),1),'c-');
    h7=plot(-250:500,nanmean(avgSpden.passive.EA(:,1523:2273),1),'y-');
    legend([h1 h2 h3 h4 h5 h6 h7],{'1stPres','ExpRep-Rep','ExpRep-Alt','ExpAlt-Rep','ExpAlt-Alt','ExpEit-Rep','ExpEit-Alt'},'Location','south','NumColumns',5)
    xlabel('Time from Stimulus Onset (ms)'); ylabel('Firing Rate (sp/s)')
    title({[labelBrain,' Visual Neurons Only (Passive)'],'Raw Firing Rates without SEM'})
    
    subplot(2,2,2); hold on % passive
    h1=fillsteplot_ahb([avgSpden.passive.RR(:,772:1522);avgSpden.passive.RA(:,772:1522);avgSpden.passive.AR(:,772:1522);avgSpden.passive.AA(:,772:1522)],1,'-',{'k'},-250:500);
    h2=fillsteplot_ahb(avgSpden.passive.RR(:,1523:2273),1,'-',{'r'},-250:500);
    h3=fillsteplot_ahb(avgSpden.passive.RA(:,1523:2273),1,'-',{'b'},-250:500);
    h4=fillsteplot_ahb(avgSpden.passive.AR(:,1523:2273),1,'-',{'g'},-250:500);
    h5=fillsteplot_ahb(avgSpden.passive.AA(:,1523:2273),1,'-',{'m'},-250:500);
    legend([h1 h2 h3 h4 h5],{'1stPres','ExpRep-Rep','ExpRep-Alt','ExpAlt-Rep','ExpAlt-Alt'},'Location','south','NumColumns',5)
    xlabel('Time from Stimulus Onset (ms)'); ylabel('Firing Rate (sp/s)')
    title({[labelBrain,' Visual Neurons Only (Passive)'],'Raw Firing Rates with SEM'})
    
    subplot(2,2,3); hold on % passive
    h1=plot(-250:500,nanmean([avgSpden.active.RR(:,772:1522);avgSpden.active.RA(:,772:1522);avgSpden.active.AR(:,772:1522);avgSpden.active.AA(:,772:1522)],1),'k-');
    h2=plot(-250:500,nanmean(avgSpden.active.RR(:,1523:2273),1),'r-');
    h3=plot(-250:500,nanmean(avgSpden.active.RA(:,1523:2273),1),'b-');
    h4=plot(-250:500,nanmean(avgSpden.active.AR(:,1523:2273),1),'g-');
    h5=plot(-250:500,nanmean(avgSpden.active.AA(:,1523:2273),1),'m-');
    legend([h1 h2 h3 h4 h5],{'1stPres','ExpRep-Rep','ExpRep-Alt','ExpAlt-Rep','ExpAlt-Alt'},'Location','south','NumColumns',5)
    xlabel('Time from Stimulus Onset (ms)'); ylabel('Firing Rate (sp/s)')
    title({[labelBrain,' Visual Neurons Only (Active)'],'Raw Firing Rates without SEM'})
    
    subplot(2,2,4); hold on % passive
    h1=fillsteplot_ahb([avgSpden.active.RR(:,772:1522);avgSpden.active.RA(:,772:1522);avgSpden.active.AR(:,772:1522);avgSpden.active.AA(:,772:1522)],1,'-',{'k'},-250:500);
    h2=fillsteplot_ahb(avgSpden.active.RR(:,1523:2273),1,'-',{'r'},-250:500);
    h3=fillsteplot_ahb(avgSpden.active.RA(:,1523:2273),1,'-',{'b'},-250:500);
    h4=fillsteplot_ahb(avgSpden.active.AR(:,1523:2273),1,'-',{'g'},-250:500);
    h5=fillsteplot_ahb(avgSpden.active.AA(:,1523:2273),1,'-',{'m'},-250:500);
    legend([h1 h2 h3 h4 h5],{'1stPres','ExpRep-Rep','ExpRep-Alt','ExpAlt-Rep','ExpAlt-Alt'},'Location','south','NumColumns',5)
    xlabel('Time from Stimulus Onset (ms)'); ylabel('Firing Rate (sp/s)')
    title({[labelBrain,' Visual Neurons Only (Active)'],'Raw Firing Rates with SEM'})
    savefig([exptdata.figuredir500,'Raw_',labelBrain,'.fig'])
        
    %% 2.2.2 Figure 2 - Raw Firing Rates
    % SOMEWHERE THERE IS A BAD NEURON. For the time being, trimming it out using a cludgy fix
    cprintf('_green','\n*** THERE IS A PROBLEM WITH ONE NEURON SOMEWHERE IN THE MIDDLE ***\n')
    NavgSpden.passive.RR(218:222,:)=[];
    NavgSpden.passive.RA(218:222,:)=[];
    NavgSpden.passive.AR(218:222,:)=[];
    NavgSpden.passive.AA(218:222,:)=[];
    NavgSpden.passive.ER([203 218:222],:)=[];
    NavgSpden.passive.EA([203 218:222],:)=[];
    NavgSpden.active.RR([203 218:222],:)=[];
    NavgSpden.active.RA([203 218:222],:)=[];
    NavgSpden.active.AR([203 218:222],:)=[];
    NavgSpden.active.AA([203 218:222],:)=[];
    NavgSpden.active.ER([203 218:222],:)=[];
    NavgSpden.active.EA([203 218:222],:)=[];
    
    figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial')
    subplot(2,2,1); hold on % passive
    h1=plot(-250:500,nanmean([NavgSpden.passive.RR(:,772:1522);NavgSpden.passive.RA(:,772:1522);NavgSpden.passive.AR(:,772:1522);NavgSpden.passive.AA(:,772:1522)],1),'k-');
    h2=plot(-250:500,nanmean(NavgSpden.passive.RR(:,1523:2273),1),'r-');
    h3=plot(-250:500,nanmean(NavgSpden.passive.RA(:,1523:2273),1),'b-');
    h4=plot(-250:500,nanmean(NavgSpden.passive.AR(:,1523:2273),1),'g-');
    h5=plot(-250:500,nanmean(NavgSpden.passive.AA(:,1523:2273),1),'m-');
    h6=plot(-250:500,nanmean(NavgSpden.passive.ER(:,1523:2273),1),'c-');
    h7=plot(-250:500,nanmean(NavgSpden.passive.EA(:,1523:2273),1),'y-');
    legend([h1 h2 h3 h4 h5 h6 h7],{'1stPres','ExpRep-Rep','ExpRep-Alt','ExpAlt-Rep','ExpAlt-Alt','ExpEit-Rep','ExpEit-Alt'},'Location','south','NumColumns',5)
    xlabel('Time from Stimulus Onset (ms)'); ylabel('Normalised Firing Rate (Normed to 1st Pres RR)')
    title({[labelBrain,' Visual Neurons Only (Passive)'],'Normed Firing Rates without SEM'})
    
    subplot(2,2,2); hold on % passive
    h1=fillsteplot_ahb([NavgSpden.passive.RR(:,772:1522);NavgSpden.passive.RA(:,772:1522);NavgSpden.passive.AR(:,772:1522);NavgSpden.passive.AA(:,772:1522)],1,'-',{'k'},-250:500);
    h2=fillsteplot_ahb(NavgSpden.passive.RR(:,1523:2273),1,'-',{'r'},-250:500);
    h3=fillsteplot_ahb(NavgSpden.passive.RA(:,1523:2273),1,'-',{'b'},-250:500);
    h4=fillsteplot_ahb(NavgSpden.passive.AR(:,1523:2273),1,'-',{'g'},-250:500);
    h5=fillsteplot_ahb(NavgSpden.passive.AA(:,1523:2273),1,'-',{'m'},-250:500);
    legend([h1 h2 h3 h4 h5],{'1stPres','ExpRep-Rep','ExpRep-Alt','ExpAlt-Rep','ExpAlt-Alt'},'Location','south','NumColumns',5)
    xlabel('Time from Stimulus Onset (ms)'); ylabel('Normalised Firing Rate (Normed to 1st Pres RR)')
    title({[labelBrain,' Visual Neurons Only (Passive)'],'Normed Firing Rates with SEM'})
    
    subplot(2,2,3); hold on % passive
    h1=plot(-250:500,nanmean([NavgSpden.active.RR(:,772:1522);NavgSpden.active.RA(:,772:1522);NavgSpden.active.AR(:,772:1522);NavgSpden.active.AA(:,772:1522)],1),'k-');
    h2=plot(-250:500,nanmean(NavgSpden.active.RR(:,1523:2273),1),'r-');
    h3=plot(-250:500,nanmean(NavgSpden.active.RA(:,1523:2273),1),'b-');
    h4=plot(-250:500,nanmean(NavgSpden.active.AR(:,1523:2273),1),'g-');
    h5=plot(-250:500,nanmean(NavgSpden.active.AA(:,1523:2273),1),'m-');
    legend([h1 h2 h3 h4 h5],{'1stPres','ExpRep-Rep','ExpRep-Alt','ExpAlt-Rep','ExpAlt-Alt'},'Location','south','NumColumns',5)
    xlabel('Time from Stimulus Onset (ms)'); ylabel('Normalised Firing Rate (Normed to 1st Pres RR)')
    title({[labelBrain,' Visual Neurons Only (Active)'],'Normed Firing Rates without SEM'})
    
    subplot(2,2,4); hold on % passive
    h1=fillsteplot_ahb([NavgSpden.active.RR(:,772:1522);NavgSpden.active.RA(:,772:1522);NavgSpden.active.AR(:,772:1522);NavgSpden.active.AA(:,772:1522)],1,'-',{'k'},-250:500);
    h2=fillsteplot_ahb(NavgSpden.active.RR(:,1523:2273),1,'-',{'r'},-250:500);
    h3=fillsteplot_ahb(NavgSpden.active.RA(:,1523:2273),1,'-',{'b'},-250:500);
    h4=fillsteplot_ahb(NavgSpden.active.AR(:,1523:2273),1,'-',{'g'},-250:500);
    h5=fillsteplot_ahb(NavgSpden.active.AA(:,1523:2273),1,'-',{'m'},-250:500);
    legend([h1 h2 h3 h4 h5],{'1stPres','ExpRep-Rep','ExpRep-Alt','ExpAlt-Rep','ExpAlt-Alt'},'Location','south','NumColumns',5)
    xlabel('Time from Stimulus Onset (ms)'); ylabel('Normalised Firing Rate (Normed to 1st Pres RR)')
    title({[labelBrain,' Visual Neurons Only (Active)'],'Normed Firing Rates with SEM'})
    savefig([exptdata.figuredir500,'Normalised_',labelBrain,'.fig'])
    
    %% 2.2.3 Figure 3 - Raw Firing Rates (Active vs. Passive)
    figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.5 0.8]); set(gca,'FontName','Arial')
    subplot(3,1,1); hold on % passive vs. active expect repeat, actual repeat
    h1=plot(nanmean(avgSpden.passive.RR(:,772:3024),1),'r:','LineWidth',2);
    h2=plot(nanmean(avgSpden.active.RR(:,772:3024),1),'r-','LineWidth',2);
    h3=plot(nanmean(avgSpden.passive.RA(:,772:3024),1),'b:','LineWidth',2);
    h4=plot(nanmean(avgSpden.active.RA(:,772:3024),1),'b-','LineWidth',2);
    tempScale=ceil(max(nanmean(avgSpden.active.RR(:,772:3024),1)));
    plot([750 750],[0 ceil(tempScale*1.5)],'k-','LineWidth',1.5)
    plot([1500 1500],[0 ceil(tempScale*1.5)],'k-','LineWidth',1.5)
    plot([250 250],[0 ceil(tempScale*1.5)],'k:','LineWidth',1.5)
    plot([1000 1000],[0 ceil(tempScale*1.5)],'k:','LineWidth',1.5)
    plot([1750 1750],[0 ceil(tempScale*1.5)],'k:','LineWidth',1.5)
    legend([h1 h2 h3 h4],'expR-actR(P)','expR-actR(A)','expR-actA(P)','expR-actA(A)'); xlim([0 2250]); ylim([0 ceil(tempScale*1.5)]);
    ylabel('Average Firing Rate (sp/s)'); xlabel('Time from Event ( ms, 1st Stim | 2nd Stim | Choice Array )');
    title({[labelBrain,' Visual Neurons Only (Active vs. Passive)'],'Expect Repeat/Actual Repeat|Expect Repeat/Actual Altern'})
    
    subplot(3,1,2); hold on % passive vs. active expect altern, actual repeat
    h1=plot(nanmean(avgSpden.passive.AA(:,772:3024),1),'r:','LineWidth',2);
    h2=plot(nanmean(avgSpden.active.AA(:,772:3024),1),'r-','LineWidth',2);
    h3=plot(nanmean(avgSpden.passive.AR(:,772:3024),1),'b:','LineWidth',2);
    h4=plot(nanmean(avgSpden.active.AR(:,772:3024),1),'b-','LineWidth',2);
    tempScale=ceil(max(nanmean(avgSpden.active.RR(:,772:3024),1)));
    plot([750 750],[0 ceil(tempScale*1.5)],'k-','LineWidth',1.5)
    plot([1500 1500],[0 ceil(tempScale*1.5)],'k-','LineWidth',1.5)
    plot([250 250],[0 ceil(tempScale*1.5)],'k:','LineWidth',1.5)
    plot([1000 1000],[0 ceil(tempScale*1.5)],'k:','LineWidth',1.5)
    plot([1750 1750],[0 ceil(tempScale*1.5)],'k:','LineWidth',1.5)
    legend([h1 h2 h3 h4],'expA-actA(P)','expA-actA(A)','expA-actR(P)','expA-actR(A)'); xlim([0 2250]); ylim([0 ceil(tempScale*1.5)]);
    ylabel('Average Firing Rate (sp/s)'); xlabel('Time from Event ( ms, 1st Stim | 2nd Stim | Choice Array )');
    title({[labelBrain,' Visual Neurons Only (Active vs. Passive)'],'Expect Altern/Actual Altern|Expect Altern/Actual Repeat'})
    
    subplot(3,1,3); hold on % passive vs. active expect altern, actual repeat
    h1=plot(nanmean(avgSpden.passive.ER(:,772:3024),1),'k:','LineWidth',2);
    h2=plot(nanmean(avgSpden.active.ER(:,772:3024),1),'k-','LineWidth',2);
    h3=plot(nanmean(avgSpden.passive.EA(:,772:3024),1),'g:','LineWidth',2);
    h4=plot(nanmean(avgSpden.active.EA(:,772:3024),1),'g-','LineWidth',2);
    tempScale=ceil(max(nanmean(avgSpden.active.RR(:,772:3024),1)));
    plot([750 750],[0 ceil(tempScale*1.5)],'k-','LineWidth',1.5)
    plot([1500 1500],[0 ceil(tempScale*1.5)],'k-','LineWidth',1.5)
    plot([250 250],[0 ceil(tempScale*1.5)],'k:','LineWidth',1.5)
    plot([1000 1000],[0 ceil(tempScale*1.5)],'k:','LineWidth',1.5)
    plot([1750 1750],[0 ceil(tempScale*1.5)],'k:','LineWidth',1.5)
    legend([h1 h2 h3 h4],'expE-actR(P)','expE-actR(A)','expE-actA(P)','expE-actA(A)'); xlim([0 2250]); ylim([0 ceil(tempScale*1.1)]);
    ylabel('Average Firing Rate (sp/s)'); xlabel('Time from Event ( ms, 1st Stim | 2nd Stim | Choice Array )');
    title({[labelBrain,' Visual Neurons Only (Active vs. Passive)'],'No Expect/Actual Repeat|No Expect/Actual Altern'})
    savefig([exptdata.figuredir500,'RawActiveVPassive_',labelBrain,'.fig'])
    
    %% 2.2.4 Figure 4 - Norm Firing Rates (Active vs. Passive)
    figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.5 0.8]); set(gca,'FontName','Arial')
    subplot(3,1,1); hold on % passive vs. active expect repeat, actual repeat
    h1=plot(nanmean(NavgSpden.passive.RR(:,772:3024),1),'r:','LineWidth',2);
    h2=plot(nanmean(NavgSpden.active.RR(:,772:3024),1),'r-','LineWidth',2);
    h3=plot(nanmean(NavgSpden.passive.RA(:,772:3024),1),'b:','LineWidth',2);
    h4=plot(nanmean(NavgSpden.active.RA(:,772:3024),1),'b-','LineWidth',2);
    tempScale=ceil(max(nanmean(NavgSpden.active.RR(:,772:3024),1)));
    plot([750 750],[0 round(tempScale*1.1,1)],'k-','LineWidth',1.5)
    plot([1500 1500],[0 round(tempScale*1.1,1)],'k-','LineWidth',1.5)
    plot([250 250],[0 round(tempScale*1.1,1)],'k:','LineWidth',1.5)
    plot([1000 1000],[0 round(tempScale*1.1,1)],'k:','LineWidth',1.5)
    plot([1750 1750],[0 round(tempScale*1.1,1)],'k:','LineWidth',1.5)
    legend([h1 h2 h3 h4],'expR-actR(P)','expR-actR(A)','expR-actA(P)','expR-actA(A)'); xlim([0 2250]); ylim([0 round(tempScale*1.1,1)]);
    ylabel('Average Firing Rate (sp/s)'); xlabel('Time from Event ( ms, 1st Stim | 2nd Stim | Choice Array )');
    title({[labelBrain,' Visual Neurons Only (Active vs. Passive)'],'Expect Repeat/Actual Repeat|Expect Repeat/Actual Altern'})
    
    subplot(3,1,2); hold on % passive vs. active expect altern, actual repeat
    h1=plot(nanmean(NavgSpden.passive.AA(:,772:3024),1),'r:','LineWidth',2);
    h2=plot(nanmean(NavgSpden.active.AA(:,772:3024),1),'r-','LineWidth',2);
    h3=plot(nanmean(NavgSpden.passive.AR(:,772:3024),1),'b:','LineWidth',2);
    h4=plot(nanmean(NavgSpden.active.AR(:,772:3024),1),'b-','LineWidth',2);
    tempScale=ceil(max(nanmean(NavgSpden.active.RR(:,772:3024),1)));
    plot([750 750],[0 round(tempScale*1.1,1)],'k-','LineWidth',1.5)
    plot([1500 1500],[0 round(tempScale*1.1,1)],'k-','LineWidth',1.5)
    plot([250 250],[0 round(tempScale*1.1,1)],'k:','LineWidth',1.5)
    plot([1000 1000],[0 round(tempScale*1.1,1)],'k:','LineWidth',1.5)
    plot([1750 1750],[0 round(tempScale*1.1,1)],'k:','LineWidth',1.5)
    legend([h1 h2 h3 h4],'expA-actA(P)','expA-actA(A)','expA-actR(P)','expA-actR(A)'); xlim([0 2250]); ylim([0 round(tempScale*1.1,1)]);
    ylabel('Average Firing Rate (sp/s)'); xlabel('Time from Event ( ms, 1st Stim | 2nd Stim | Choice Array )');
    title({[labelBrain,' Visual Neurons Only (Active vs. Passive)'],'Expect Altern/Actual Altern|Expect Altern/Actual Repeat'})
    
    subplot(3,1,3); hold on % passive vs. active expect altern, actual repeat
    h1=plot(nanmean(NavgSpden.passive.ER(:,772:3024),1),'k:','LineWidth',2);
    h2=plot(nanmean(NavgSpden.active.ER(:,772:3024),1),'k-','LineWidth',2);
    h3=plot(nanmean(NavgSpden.passive.EA(:,772:3024),1),'g:','LineWidth',2);
    h4=plot(nanmean(NavgSpden.active.EA(:,772:3024),1),'g-','LineWidth',2);
    tempScale=ceil(max(nanmean(NavgSpden.active.RR(:,772:3024),1)));
    plot([750 750],[0 round(tempScale*1.1,1)],'k-','LineWidth',1.5)
    plot([1500 1500],[0 round(tempScale*1.1,1)],'k-','LineWidth',1.5)
    plot([250 250],[0 round(tempScale*1.1,1)],'k:','LineWidth',1.5)
    plot([1000 1000],[0 round(tempScale*1.1,1)],'k:','LineWidth',1.5)
    plot([1750 1750],[0 round(tempScale*1.1,1)],'k:','LineWidth',1.5)
    legend([h1 h2 h3 h4],'expA-actA(P)','expA-actA(A)','expA-actR(P)','expA-actR(A)'); xlim([0 2250]); ylim([0 round(tempScale*1.1,1)]);
    ylabel('Average Firing Rate (sp/s)'); xlabel('Time from Event ( ms, 1st Stim | 2nd Stim | Choice Array )');
    title({[labelBrain,' Visual Neurons Only (Active vs. Passive)'],'Expect Altern/Actual Altern|Expect Altern/Actual Repeat'})
    savefig([exptdata.figuredir500,'NormActiveVPassive_',labelBrain,'.fig'])    
end 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

