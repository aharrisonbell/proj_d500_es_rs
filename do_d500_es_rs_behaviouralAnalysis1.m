%% do_es_rs_behaviouralAnalysis1
% by CS, AHB, started Nov 6, 2018
vers_ephys_es_rs='1.3; May 15, 2023';
% 1.0 - original version (Nov 6, 2018)
% 1.1 - added stimulus ID features (face ID, etc.) (Feb 27, 2019)
% 1.2 - updated with new analyses, etc. (March 29, 2020)
% 1.3 - updated to run on current version and on PC (HOMER_B550) (May 15, 2023)

[d500specs.fList,d500specs.pList] = matlab.codetools.requiredFilesAndProducts('do_d500_es_rs_behaviouralAnalysis1.m');

global exptdata
close all;
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
exptdata.analysisName='D500_ES-RS_Study'; % used for savenames, figures, etc. (pick whatever you want; will be used for filenames)
exptdata.projectdir=[exptdata.analysisdir,exptdata.analysisName,filesep]; 
exptdata.figuredir500=[exptdata.projectdir,'figures',filesep,'d500_es_rs',filesep]; mkdir(exptdata.figuredir500);

%% Set Analysis Parameters
d500specs.rt_cutoff=500; % eliminate RTs longer than 500 ms
d500specs.upperLim=.66;
d500specs.lowerLim=.33;

%% SCRATCH PAD (March 30, 2020)
% - 

%% x. INTRODUCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
clc;
fprintf('<strong>*======================================*</strong>\n')
fprintf('<strong>| do_d500_es_rs_behaviouralAnalysis1.m |</strong>\n')
fprintf('<strong>*======================================*</strong>\n')
fprintf(['Version: ',vers_ephys_es_rs,'\n'])
disp(['Study name:             ',exptdata.analysisName]);
disp(['ePhys Data location:    ',exptdata.datalocation]);
disp(['Processed NEV dir:      ',exptdata.processedDatadir]);
disp(['Project Directory:      ',exptdata.projectdir]);
disp(' ');
disp(' ');

%% 0. Load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fprintf('* Loading behavioural data (from all behavioural sessions)...')
load([exptdata.projectdir,'DMS500_behavData.mat'],'MegaBehav','SessionNames','trialCounts'); 
behavData=MegaBehav; dms500sessionNames=SessionNames;
clear MegaBehav


% Data structure:
%  1) monkey number
%  2) session number
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
% 16) right/left
% 17) nothing
% 18) task number (300,500,600)
% 19) new block number (600 only)
% 20,21,22) save these for neurons
% 23) stimA (task alternates between stim A and B)
% 24) stimB
% 25) distractor1
% 26) distractor2

%% 1. Plot Trial Counts / QA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.2 0.1 0.4 0.8]); set(gca,'FontName','Arial')
subplot(1,2,1)
plot(behavData(:,5),behavData(:,2),'k.')
xlabel('Number of Trials'); ylabel('Session Number'); title('Behavioural Sessions: BEFORE FIX')
fprintf('done.\n')

% Trim the fat (get rid of shit trials: those that are not either correct or incorrect choice)
behavData=behavData(ismember(behavData(:,6),[0 6]),:);
% behavData=behavData(behavData(:,9)>d500specs.rt_cutoff,:);

subplot(1,2,2)
plot(behavData(:,5),behavData(:,2),'r.')
xlabel('Number of Trials'); ylabel('Session Number'); title('Behavioural Sessions: AFTER FIX')
mtit('Filtering out Bad Trials','FontSize',16,'color',[0 0 1],'xoff',-.1,'yoff',.025);

for ss=1:length(unique(behavData(:,2))) % trialCounts(:,3) = GOOD trials
    trialCounts(ss,3)=length(find(behavData(:,2)==ss)); %#ok<*SAGROW>
end

% Simple "Demographics"
behavSummary(1,1)=length(unique(behavData(behavData(:,1)==1,2))); % number of sessions
behavSummary(1,2)=length(unique(behavData(behavData(:,1)==2,2)));
behavSummary(2,1)=length(behavData(behavData(:,1)==1,1)); % number of trials
behavSummary(2,2)=length(behavData(behavData(:,1)==2,1));
behavSummary(3,1)=length(behavData(behavData(:,1)==1 & behavData(:,7)>3,1));
behavSummary(3,2)=length(behavData(behavData(:,1)==2 & behavData(:,7)>3,1));
behavSummary(4,1)=length(behavData(behavData(:,1)==1 & behavData(:,7)>3 & behavData(:,6)==0,1));
behavSummary(4,2)=length(behavData(behavData(:,1)==2 & behavData(:,7)>3 & behavData(:,6)==0,1));
behavSummary(5,1)=length(behavData(behavData(:,1)==1 & behavData(:,7)<4,1));
behavSummary(5,2)=length(behavData(behavData(:,1)==2 & behavData(:,7)<4,1));
behavSummary(6,1)=length(behavData(behavData(:,1)==1 & behavData(:,7)<4 & behavData(:,6)==0,1));
behavSummary(6,2)=length(behavData(behavData(:,1)==2 & behavData(:,7)<4 & behavData(:,6)==0,1));

disp('-------------------------------------------------------------------------------------------')
fprintf(['Total Number of Sessions (M1): ',num2str(behavSummary(1,1)),'\n'])
fprintf(['Total Number of Sessions (M2): ',num2str(behavSummary(1,2)),'\n'])
fprintf(['Total Number of Trials (M1): ',num2str(behavSummary(2,1)),'\n'])
fprintf(['Total Number of Trials (M2): ',num2str(behavSummary(2,2)),'\n'])

disp('-------------------------------------------------------------------------------------------')
fprintf(['Total Number of Passive Trials (M1): ',...
    num2str(behavSummary(3,1)),'\n'])
fprintf(['Total Number of Passive Trials (M2): ',...
    num2str(behavSummary(3,2)),'\n'])

disp(['Total Number of Passive Trials (Correct) (M1): ',...
    num2str(behavSummary(4,1)),' (',num2str(behavSummary(4,1)/behavSummary(3,1)*100),'%)'])
disp(['Total Number of Passive Trials (Correct) (M2): ',...
    num2str(behavSummary(4,2)),' (',num2str(behavSummary(4,2)/behavSummary(3,2)*100),'%)'])
disp('-------------------------------------------------------------------------------------------')
fprintf(['Total Number of Active Trials (M1): ',...
    num2str(behavSummary(5,1)),'\n'])
fprintf(['Total Number of Active Trials (M2): ',...
    num2str(behavSummary(5,2)),'\n'])

disp(['Total Number of Active Trials (Correct) (M1): ',...
    num2str(behavSummary(6,1)),' (',num2str(behavSummary(6,1)/behavSummary(5,1)*100),'%)'])
disp(['Total Number of Active Trials (Correct) (M2): ',...
    num2str(behavSummary(6,2)),' (',num2str(behavSummary(6,2)/behavSummary(5,2)*100),'%)'])




%% 2. Create Regressor Data Structure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clear evdata
evdata=[];
evdata.trialNumber = behavData(:,5);              % paste trial number (note, will not be unique for each session)
evdata.animal = behavData(:,1);                   % paste animal number
evdata.actual = logic2sign(behavData(:,14)==1);   % actual repeat or alternation
evdata.session = behavData(:,2);                  % session
evdata.cor = double(behavData(:,6)==0);           % correct or incorrect choice
evdata.active = logic2sign(behavData(:,15)==1);   % active or passive blocks
evdata.stim1 = logic2sign(behavData(:,11) ==1);   % stim on first presentation (1 if A, 2 if B)
evdata.stim2 = logic2sign(behavData(:,12) ==1);   % stim on second presentation (1 if A, 2 if B)
evdata.reactiontime = behavData(:,9);             % reaction time
evdata.expect = behavData(:,13);                  % expectation (1,2,3)
evdata.correctRight = behavData(:,16);            % whether correct answer was on right (1) or left (-1)
evdata.date_lastModified=datetime('today');

% relabel evdata.expect (1=1, 2=-1, 3=0)
evdata.expect(evdata.expect==1)=1;
evdata.expect(evdata.expect==2)=-1;                  % expectation (1,2,3)
evdata.expect(evdata.expect==3)=0;                  % expectation (1,2,3)

% Paste in Stim IDs (hard-coded in exptdata)
evdata.stimIDs = [evdata.animal evdata.animal evdata.animal evdata.animal]; % create blanks
for ss=1:length(unique(evdata.session))
    tempPointer=find(evdata.session==ss);
    evdata.stimIDs(tempPointer,:)=repmat(exptdata.dms500.stimIDs(ss,2:5),length(tempPointer),1);
end

% Solve for face/non-face
evdata.stim1_faceNonFace=evdata.trialNumber*0; % create empty vectors
evdata.stim2_faceNonFace=evdata.trialNumber*0; % create empty vectors
evdata.choseRepeat=evdata.trialNumber*0; % create empty vectors
for tr=1:length(evdata.trialNumber)
    % is first stim presentation a face? (stim1)?
    if evdata.stim1(tr)==1 && ismember(evdata.stimIDs(tr,1),exptdata.dms500.stimFaces) % if stim1 is "A" and "A" is a face?
        evdata.stim1_faceNonFace(tr)=1;
    elseif evdata.stim1(tr)==-1 && ismember(evdata.stimIDs(tr,2),exptdata.dms500.stimFaces) % if stim is "B" and "B" is a face?
        evdata.stim1_faceNonFace(tr)=1;
    else
        evdata.stim1_faceNonFace(tr)=-1;
    end

    % is second stim presentation a face? (stim2)?
    if evdata.stim2(tr)==1 && ismember(evdata.stimIDs(tr,1),exptdata.dms500.stimFaces) % if stim1 is "A" and "A" is a face?
        evdata.stim2_faceNonFace(tr)=1;
    elseif evdata.stim2(tr)==-1 && ismember(evdata.stimIDs(tr,2),exptdata.dms500.stimFaces) % if sti1m is "B" and "B" is a face?
        evdata.stim2_faceNonFace(tr)=1;
    else
        evdata.stim2_faceNonFace(tr)=-1;
    end
    
    % chose repeated stimulus(1=yes, 0=no)
    evdata.choseRepeat(tr)=0; % default to no
    if evdata.cor(tr)==1 && evdata.actual(tr)==1 % if actual = repeat and correct
        evdata.choseRepeat(tr)=1;
    end
    if evdata.cor(tr)==0 && evdata.actual(tr)==-1 % if actual = alternate and incorrect
        evdata.choseRepeat(tr)=1;
    end
end
clear tr

% this variable determines whether the correct answer on the current trial
% was on the same side as the previous trial
evdata.congruent_n1=zeros(size(behavData,1),1);    % initialise variable
for tt=2:size(behavData,1)
    if evdata.correctRight(tt)==evdata.correctRight(tt-1)
        evdata.congruent_n1(tt)=1;
    end
end
evdata.congruent_n2=zeros(size(behavData,1),1);    % initialise variable
for tt=3:size(behavData,1)
    if evdata.correctRight(tt)==evdata.correctRight(tt-2)
        evdata.congruent_n2(tt)=1;
    end
end
evdata.congruent_n3=zeros(size(behavData,1),1);    % initialise variable
for tt=4:size(behavData,1)
    if evdata.correctRight(tt)==evdata.correctRight(tt-3)
        evdata.congruent_n3(tt)=1;
    end
end
clear tt

% this variable determines whether previous STIM2 is the same as current STIM1
evdata.prevtrial_n1=zeros(size(behavData,1),1);    % initialise variable
evdata.prevtrial_n2=zeros(size(behavData,1),1);    % initialise variable
evdata.prevtrial_n3=zeros(size(behavData,1),1);    % initialise variable


for tt=4:size(behavData,1)
    if evdata.stim2(tt-1)==evdata.stim1(tt)
        evdata.prevtrial_n1(tt)=1;
    end
    if evdata.stim2(tt-2)==evdata.stim1(tt)
        evdata.prevtrial_n2(tt)=1;
    end
    if evdata.stim2(tt-3)==evdata.stim1(tt)
        evdata.prevtrial_n3(tt)=1;
    end
end
clear tt


%% 3. Calculate latent expectation of repeat ("p_rep", using delta rule) and plot each session ~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

evdata.prep=[]; evdata.smooth_p_rep=[];
clear Betas*;
d500specs.alpha = 0.1; % learning rate

%
%Use MLE to claculate top d500specs.alpha
%

h=figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial')
sessmat = unique(evdata.session);
panelCount=0; indexCounter=0; figCount=1;
for s = 1:length(sessmat)
    temp_indx = find(evdata.session == sessmat(s));
    panelCount=panelCount+1;
    if panelCount<4
        subplot(3,1,panelCount); hold on
    else
        savefig(h,[exptdata.figuredir500,'d500_es_rs_sessionBehaviour500_',num2str(figCount),'.fig']);
        jpgfigname=[exptdata.figuredir500,'d500_es_rs_sessionBehaviour500_',num2str(figCount),'.jpg'];
        print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
        
        h=figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial')
        panelCount=1; figCount=figCount+1;
        subplot(3,1,panelCount); hold on
    end
        
    %plot(data.actual(temp_indx));
    clear log_p_rep delta_p_rep
    p_rep = 0.5;
    for t = 1:length(temp_indx)
        delta_p_rep = (evdata.actual(temp_indx(t))==1)-p_rep;
        p_rep = p_rep + d500specs.alpha*delta_p_rep;
        log_p_rep(t) = p_rep; 
        evdata.prep(temp_indx(t),1) = p_rep;
        evdata.delta_p_rep(temp_indx(t),1) = delta_p_rep;
    end
    h1=plot(log_p_rep,'k-'); % plot p_rep
    ylim([0 1]);
    
    % calculate session choice probability (modelled)
    tempChoice(s,1)=length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)>d500specs.upperLim & evdata.choseRepeat(temp_indx)==1)) / length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)>d500specs.upperLim));
    tempChoice(s,2)=length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)<d500specs.lowerLim & evdata.choseRepeat(temp_indx)==1)) / length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)<d500specs.lowerLim));
    tempChoice(s,3)=length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)<d500specs.upperLim & evdata.prep(temp_indx)>d500specs.lowerLim & evdata.choseRepeat(temp_indx)==1)) / length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)<d500specs.upperLim & evdata.prep(temp_indx)>d500specs.lowerLim));
    tempChoice=tempChoice*100;
       
    % by generative prob 
    tempChoice(s,1)=length(find(evdata.active(temp_indx)==1 & evdata.expect(temp_indx)==1 & evdata.choseRepeat(temp_indx)==1)) / length(find(evdata.active(temp_indx)==1 & evdata.expect(temp_indx)==1));
    tempChoice(s,2)=length(find(evdata.active(temp_indx)==1 & evdata.expect(temp_indx)==-1 & evdata.choseRepeat(temp_indx)==1)) / length(find(evdata.active(temp_indx)==1 & evdata.expect(temp_indx)==-1));
    tempChoice(s,3)=length(find(evdata.active(temp_indx)==1 & evdata.expect(temp_indx)==0 & evdata.choseRepeat(temp_indx)==1)) / length(find(evdata.active(temp_indx)==1 & evdata.expect(temp_indx)==0));
    tempChoice=tempChoice*100;

    % confirm generative probability by block
    tempG(s,1)=length(find(evdata.active(temp_indx)==1 & evdata.expect(temp_indx)==1 & evdata.actual(temp_indx)==1)) / length(find(evdata.active(temp_indx)==1 & evdata.expect(temp_indx)==1));
    tempG(s,2)=length(find(evdata.active(temp_indx)==1 & evdata.expect(temp_indx)==-1 & evdata.actual(temp_indx)==1)) / length(find(evdata.active(temp_indx)==1 & evdata.expect(temp_indx)==-1));
    tempG(s,3)=length(find(evdata.active(temp_indx)==1 & evdata.expect(temp_indx)==0 & evdata.actual(temp_indx)==1)) / length(find(evdata.active(temp_indx)==1 & evdata.expect(temp_indx)==0));
    tempG(s,4)=length(find(evdata.active(temp_indx)~=1 & evdata.expect(temp_indx)==1 & evdata.actual(temp_indx)==1)) / length(find(evdata.active(temp_indx)~=1 & evdata.expect(temp_indx)==1));
    tempG(s,5)=length(find(evdata.active(temp_indx)~=1 & evdata.expect(temp_indx)==-1 & evdata.actual(temp_indx)==1)) / length(find(evdata.active(temp_indx)~=1 & evdata.expect(temp_indx)==-1));
    tempG(s,6)=length(find(evdata.active(temp_indx)~=1 & evdata.expect(temp_indx)==0 & evdata.actual(temp_indx)==1)) / length(find(evdata.active(temp_indx)~=1 & evdata.expect(temp_indx)==0));
    tempG=tempG*100;

    % Detailed Graph
    xlabel('Trial')
    ylabel({'Expectation','(1=Repeat/0=Alternation'})
    title({[dms500sessionNames{s},' / Session #',num2str(s),' Total Trials: ',num2str(trialCounts(s,2)),'/Included Trials: ',num2str(trialCounts(s,3))],...
        ['GP (Act): ',num2str(tempG(s,1),'%2.1f'),'/',num2str(tempG(s,2),'%2.1f'),'/',num2str(tempG(s,3),'%2.1f'),' // ',...
        'GP (Pas): ',num2str(tempG(s,4),'%2.1f'),'/',num2str(tempG(s,5),'%2.1f'),'/',num2str(tempG(s,6),'%2.1f')],...
        ['CP: ',num2str(tempChoice(s,1),'%2.1f'),'/',num2str(tempChoice(s,2),'%2.1f'),'/',num2str(tempChoice(s,3),'%2.1f')]},...
        'FontSize',14,'Interpreter','None');    
    h2=plot(evdata.active(temp_indx),'ms'); % active vs. passive
    
    
    temp_indx1 = find(evdata.session == sessmat(s) & evdata.expect == 1 );
    h3=plot(temp_indx1-indexCounter,ones(length(temp_indx1),1)*0.75,'c.','LineWidth',2);
    temp_indx1 = find(evdata.session == sessmat(s) & evdata.expect == -1 );
    plot(temp_indx1-indexCounter,ones(length(temp_indx1),1)*0.25,'c.','LineWidth',2);
    temp_indx1 = find(evdata.session == sessmat(s) & evdata.expect == 0 );
    plot(temp_indx1-indexCounter,ones(length(temp_indx1),1)*0.50,'c.','LineWidth',2);
    
    % trial outcome markers
    temp_indx1 = find(evdata.session == sessmat(s) & evdata.cor == 1 & evdata.actual == 1);
    h4=plot(temp_indx1-indexCounter,ones(length(temp_indx1),1)*0.85,'gd','MarkerFaceColor',[0 .7 .25]);
    
    temp_indx1 = find(evdata.session == sessmat(s) & evdata.cor == 0 & evdata.actual == -1);
    plot(temp_indx1-indexCounter,ones(length(temp_indx1),1)*0.15,'gd','MarkerFaceColor',[1 0 .25]);
    
    temp_indx1 = find(evdata.session == sessmat(s) & evdata.cor == 0 & evdata.actual == 1);
    h6=plot(temp_indx1-indexCounter,ones(length(temp_indx1),1)*0.85,'gd','MarkerFaceColor',[1 0 .25]);
    
    temp_indx1 = find(evdata.session == sessmat(s) & evdata.cor == 1 & evdata.actual == -1);
    h7=plot(temp_indx1-indexCounter,ones(length(temp_indx1),1)*0.15,'gd','MarkerFaceColor',[0 .7 .25]);
    
    clear temp_indx1 
    % Plot floating average
    floatAverageKernel=5; % average across 5 trials
    tmp_smoothData=mySmooth(log_p_rep',floatAverageKernel,1,'backward');
    h5=plot(tmp_smoothData,'-','LineWidth',3,'Color',[.5 .7 .9]) ;
    evdata.smooth_p_rep(temp_indx)=tmp_smoothData';
    indexCounter=indexCounter+length(temp_indx);
    if panelCount==1
        legend([h1 h5 h3 h2 h4],{'p-rep(est)','p-rep(estSmooth)','p-rep(actual)','active','chooseR'},'Location','south','NumColumns',5)
    end
end
clear h* tmp* temp* panelCount


%% 4. Multiple Regression and Plot Figures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear Betas_*
sessmat = unique(evdata.session);
for s = 1:length(sessmat) % run GLM on each session
    temp_indx = find(evdata.session == sessmat(s));
    
    %**** Possible regressors to include: ****
    % btw: glmfit and fitglm are essentially the same thing
    
    % data.stim1 (stimA=1; stimB=-1)
    % data.stim2 (stimA=1; stimB=-1)
    % data.session: session number
    % data.actual: (ACTUAL trial type: 1=repeat, -1=alternation) 
    % data.cor: (correct=1, incorrect=0)
    % data.active (1) or passive (-1) block
    % data.reactiontime: (in ms)
    % data.expect: expectation (1=expect repeat; -1=expect alternation; 0=no expectation)
    % data.choseRepeat: whether they chose repeated stimulus
    % data.correctRight = whether correct answer was on right (1) or left (-1)
    % data.congruent_n1 (n minus 1) = correct side is same for trial n and n-1
    % data.congruent_n1 (n minus 2) = correct side is same for trial n and n-2
    % data.congruent_n3 (n minus 3) = correct side is same for trial n and n-3
    
    %** glm on CORRECT OR INCORRECT
    tempregressors = [evdata.actual(temp_indx) evdata.expect(temp_indx).*evdata.active(temp_indx) evdata.prep(temp_indx) evdata.prep(temp_indx).*evdata.active(temp_indx) evdata.congruent_n1(temp_indx) evdata.congruent_n2(temp_indx) evdata.congruent_n3(temp_indx)];
    xtick1labels={'intercept','actual','expect*active','p_rep','p_rep*active','congruent-n1','congruent-n2','congruent-n3'};
    Betas_cor(s,:) = glmfit(tempregressors,evdata.cor(temp_indx));  
    clear tempregressors
    
    %** glm on CHOSE REPEAT
    tempregressors = [evdata.actual(temp_indx) evdata.active(temp_indx) evdata.prep(temp_indx) evdata.active(temp_indx).*evdata.prep(temp_indx)];
    xtick2labels={'intercept','actual','active','p_rep','active*prep'};
    Betas_choice(s,:) = glmfit(tempregressors,evdata.choseRepeat(temp_indx),'binomial','probit');  
    clear tempregressors
    
    %** stepwise glm on CHOSE REPEAT
    tempregressors = [evdata.expect(temp_indx) evdata.active(temp_indx) evdata.prep(temp_indx) evdata.prep(temp_indx).*evdata.active(temp_indx) evdata.congruent_n1(temp_indx) evdata.congruent_n2(temp_indx) evdata.congruent_n3(temp_indx)];
    xtick4labels={'intercept','expect','active','prep','active*prep','congruent-n1','congruent-n2','congruent-n3'};
    [Betas_choiceStep(s,:),~,~,Betas_choiceStepInModel(s,:),~,~,~] = stepwisefit(tempregressors,evdata.choseRepeat(temp_indx));  
    clear tempregressors
    
    %** glm on REACTION TIME
    tempregressors = [evdata.actual(temp_indx) evdata.expect(temp_indx) evdata.active(temp_indx) evdata.prep(temp_indx) evdata.active(temp_indx).*evdata.prep(temp_indx)];
    xtick3labels={'intercept','actual','expect','active','prep','active*prep'};
    Betas_rt(s,:) = glmfit(tempregressors,evdata.reactiontime(temp_indx));  
    clear tempregressors
    
    %----------------
    % Calculate SDT measures for each session
  
    temp_hit(1)=length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)>d500specs.upperLim & evdata.actual(temp_indx)==1 & evdata.choseRepeat(temp_indx)==1))/...
        length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)>d500specs.upperLim & evdata.actual(temp_indx)==1));
    temp_hit(2)=length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)<d500specs.upperLim & evdata.prep(temp_indx)>d500specs.lowerLim & evdata.actual(temp_indx)==1 & evdata.choseRepeat(temp_indx)==1))/...
        length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)<d500specs.upperLim & evdata.prep(temp_indx)>d500specs.lowerLim & evdata.actual(temp_indx)==1));
    temp_hit(3)=length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)<d500specs.lowerLim & evdata.actual(temp_indx)==1 & evdata.choseRepeat(temp_indx)==1))/...
        length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)<d500specs.lowerLim & evdata.actual(temp_indx)==1));
    temp_hit(4)=length(find(evdata.active(temp_indx)~=1 & evdata.prep(temp_indx)>d500specs.upperLim & evdata.actual(temp_indx)==1 & evdata.choseRepeat(temp_indx)==1))/...
        length(find(evdata.active(temp_indx)~=1 & evdata.prep(temp_indx)>d500specs.upperLim & evdata.actual(temp_indx)==1));
    temp_hit(5)=length(find(evdata.active(temp_indx)~=1 & evdata.prep(temp_indx)<d500specs.upperLim & evdata.prep(temp_indx)>d500specs.lowerLim & evdata.actual(temp_indx)==1 & evdata.choseRepeat(temp_indx)==1))/...
        length(find(evdata.active(temp_indx)~=1 & evdata.prep(temp_indx)<d500specs.upperLim & evdata.prep(temp_indx)>d500specs.lowerLim & evdata.actual(temp_indx)==1));
    temp_hit(6)=length(find(evdata.active(temp_indx)~=1 & evdata.prep(temp_indx)<d500specs.lowerLim & evdata.actual(temp_indx)==1 & evdata.choseRepeat(temp_indx)==1))/...
        length(find(evdata.active(temp_indx)~=1 & evdata.prep(temp_indx)<d500specs.lowerLim & evdata.actual(temp_indx)==1));
    
    
    temp_falseAlarm(1)=length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)>d500specs.upperLim & evdata.actual(temp_indx)~=1 & evdata.choseRepeat(temp_indx)==1))/...
        length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)>d500specs.upperLim & evdata.actual(temp_indx)~=1));    
    temp_falseAlarm(2)=length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)<d500specs.upperLim & evdata.prep(temp_indx)>d500specs.lowerLim & evdata.actual(temp_indx)~=1 & evdata.choseRepeat(temp_indx)==1))/...
        length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)<d500specs.upperLim & evdata.prep(temp_indx)>d500specs.lowerLim & evdata.actual(temp_indx)~=1));    
    temp_falseAlarm(3)=length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)<d500specs.lowerLim & evdata.actual(temp_indx)~=1 & evdata.choseRepeat(temp_indx)==1))/...
        length(find(evdata.active(temp_indx)==1 & evdata.prep(temp_indx)<d500specs.lowerLim & evdata.actual(temp_indx)~=1));
    temp_falseAlarm(4)=length(find(evdata.active(temp_indx)~=1 & evdata.prep(temp_indx)>d500specs.upperLim & evdata.actual(temp_indx)~=1 & evdata.choseRepeat(temp_indx)==1))/...
        length(find(evdata.active(temp_indx)~=1 & evdata.prep(temp_indx)>d500specs.upperLim & evdata.actual(temp_indx)~=1));    
    temp_falseAlarm(5)=length(find(evdata.active(temp_indx)~=1 & evdata.prep(temp_indx)<d500specs.upperLim & evdata.prep(temp_indx)>d500specs.lowerLim & evdata.actual(temp_indx)~=1 & evdata.choseRepeat(temp_indx)==1))/...
        length(find(evdata.active(temp_indx)~=1 & evdata.prep(temp_indx)<d500specs.upperLim & evdata.prep(temp_indx)>d500specs.lowerLim & evdata.actual(temp_indx)~=1));    
    temp_falseAlarm(6)=length(find(evdata.active(temp_indx)~=1 & evdata.prep(temp_indx)<d500specs.lowerLim & evdata.actual(temp_indx)~=1 & evdata.choseRepeat(temp_indx)==1))/...
        length(find(evdata.active(temp_indx)~=1 & evdata.prep(temp_indx)<d500specs.lowerLim & evdata.actual(temp_indx)~=1));
    
    [sdt(s,1,1), sdt(s,1,2)]=dprime_simple(temp_hit(1),temp_falseAlarm(1)); % dprime & crit
    [sdt(s,2,1), sdt(s,2,2)]=dprime_simple(temp_hit(2),temp_falseAlarm(2)); % dprime & crit
    [sdt(s,3,1), sdt(s,3,2)]=dprime_simple(temp_hit(3),temp_falseAlarm(3)); % dprime & crit
    [sdt(s,4,1), sdt(s,4,2)]=dprime_simple(temp_hit(4),temp_falseAlarm(4)); % dprime & crit
    [sdt(s,5,1), sdt(s,5,2)]=dprime_simple(temp_hit(5),temp_falseAlarm(5)); % dprime & crit
    [sdt(s,6,1), sdt(s,6,2)]=dprime_simple(temp_hit(6),temp_falseAlarm(6)); % dprime & crit
  
    clear temp*
end
h=figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial');
subplot(5,1,1); hold on
errorbar(mean(Betas_cor),sem(Betas_cor)); xlim([0 length(mean(Betas_cor))+1]);
plot([0 length(mean(Betas_cor))],[0 0],'k:'); ylabel('Parameter Estimate')
for x=1:length(mean(Betas_cor))
   [P,~]=signrank(Betas_cor(:,x),Betas_cor(:,x)*0);
   text(x,0.35,['P=',num2str(P)],'FontSize',14)
end
set(gca,'XTick',1:length(mean(Betas_cor)),'XTickLabels',xtick1labels)
title('DMS500 - behaviour (correct trials)','FontSize',16)

subplot(5,1,2); hold on
errorbar(mean(Betas_choice),sem(Betas_choice)); xlim([0 length(mean(Betas_choice))+1]);
plot([0 length(mean(Betas_choice))],[0 0],'k:'); ylabel('Parameter Estimate')
for x=1:length(mean(Betas_choice))
   [P,~]=signrank(Betas_choice(:,x),Betas_choice(:,x)*0);
   text(x,0.35,['P=',num2str(P)],'FontSize',14)
end
set(gca,'XTick',1:length(mean(Betas_choice)),'XTickLabels',xtick2labels)
title('DMS500 - behaviour (chose repeat)','FontSize',16)

subplot(5,1,3); hold on
bar(sum(Betas_choiceStepInModel)/size(Betas_choiceStepInModel,1))
errorbar(mean(Betas_choiceStep),sem(Betas_choiceStep)); xlim([0 length(mean(Betas_choiceStep))+1]);
plot([0 length(mean(Betas_choiceStep))],[0 0],'k:'); ylabel('Parameter Estimate')
for x=1:length(mean(Betas_choiceStep))
   [P,~]=signrank(Betas_choiceStep(:,x),Betas_choiceStep(:,x)*0);
   text(x,0.15,['P=',num2str(P)],'FontSize',14)
end
set(gca,'XTick',1:length(mean(Betas_choiceStep)),'XTickLabels',xtick4labels(2:end))
title('DMS500 - behaviour (chose repeat STEP)','FontSize',16)





subplot(5,1,4); hold on
errorbar(mean(Betas_rt),sem(Betas_rt)); xlim([0 length(mean(Betas_rt))+1]);
plot([0 length(mean(Betas_rt))],[0 0],'k:'); ylabel('Parameter Estimate')
for x=1:length(mean(Betas_rt))
    [P,~]=signrank(Betas_rt(:,x),Betas_rt(:,x)*0);
   text(x,15,['P=',num2str(P)],'FontSize',14)
end
set(gca,'XTick',1:length(mean(Betas_rt)),'XTickLabels',xtick3labels)
title('DMS500 - behaviour (reactiontime)','FontSize',16)

subplot(5,2,9); hold on % sdt
temp_dprime=squeeze(sdt(:,:,1));
temp_dprime(isinf(temp_dprime))=nan;
temp_crit=squeeze(sdt(:,:,2));
temp_crit(isinf(temp_crit))=nan;
errorbar(nanmean(temp_dprime),nanstd(temp_dprime)/sqrt(length(sessmat))); xlim([0 length(nanmean(temp_dprime))+1]); 
ylabel('std units')
for x=1:length(mean(temp_dprime))
    [P,~]=signrank(temp_dprime(:,x),temp_dprime(:,x)*0);
   text(x,.15,['P=',num2str(P)],'FontSize',8)
end
set(gca,'XTick',1:length(mean(temp_dprime)),'XTickLabels',{'AexpR' 'AnoExp' 'AexpA' 'PexpR' 'PnoExp' 'PexpA'})
title('DMS500 - dprime','FontSize',16)

subplot(5,2,10); hold on % sdt
errorbar(nanmean(temp_crit),nanstd(temp_crit)/sqrt(length(sessmat))); xlim([0 length(nanmean(temp_crit))+1]); %#ok<*NANSTD>
ylabel('std units')
for x=1:length(mean(temp_crit))
    [P,~]=signrank(temp_crit(:,x),temp_crit(:,x)*0);
   text(x,.15,['P=',num2str(P)],'FontSize',8)
end
set(gca,'XTick',1:length(mean(Betas_rt)),'XTickLabels',{'AexpR' 'AnoExp' 'AexpA' 'PexpR' 'PnoExp' 'PexpA'})
title('DMS500 - criterion','FontSize',16)

savefig(h,[exptdata.figuredir500,'d500_behaviour_regressionAnalysis.fig']);
jpgfigname=[exptdata.figuredir500,'d500_behaviour_regressionAnalysis.jpg'];
print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure

clear xtick* x tt t tr s P p_rep temp*

%% 5. Histogram to identify anticipatory movements/saccades? ~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bins=0:10:580;
h=figure; hold on
subplot(2,1,1) % active
temp_indx=find(evdata.active==1);
histogram(evdata.reactiontime(temp_indx),bins,'FaceColor','r')
xlabel('Saccadic Reaction Time (ms)'); ylabel('Trial Counts')
title('Evaluating Anticipatory Behaviour','FontSize',16)
subplot(2,1,2); hold on
corProb=zeros(length(bins),1);
for tempB=1:length(bins)-1
    temp_trialCount(tempB,1)=length(find(evdata.active==1 & evdata.reactiontime>bins(tempB) & evdata.reactiontime<=bins(tempB+1) & evdata.cor==1));
    temp_trialCount(tempB,2)=length(find(evdata.active==1 & evdata.reactiontime>bins(tempB) & evdata.reactiontime<=bins(tempB+1)));
    corProb(tempB)=(temp_trialCount(tempB,1)/temp_trialCount(tempB,2))*100;
end
plot(bins,corProb,'k-','LineWidth',2)
xlabel('Saccadic Reaction Time (ms)'); ylabel('% Correct')

savefig(h,[exptdata.figuredir500,'d500_behaviour_histogram_rts.fig']);
jpgfigname=[exptdata.figuredir500,'d500_behaviour_histogram_rts.jpg'];
print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure

%% 6. Bar Graph of Choice Probability  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% behavProbs      |     Active         |   Passive              
%                 | Repeat | Alternate | Repeat | Alternate |
%   expect repeat |        |           |
%   expect altern |        |           |
%       no expect |        |           |
%

% behavProbsM          |     Active         |   Passive 
%                      | Repeat | Alternate | Repeat | Alternate |
%   expect repeat >.65 |        |           |
%   expect altern <.35 |        |           |
%       no expect      |        |           |
%

behavProbs=[];
behavProbs(1,1)=length(find(evdata.active==1 & evdata.expect==1 & evdata.choseRepeat==1));
behavProbs(1,2)=length(find(evdata.active==1 & evdata.expect==1 & evdata.choseRepeat==0));
behavProbs(2,1)=length(find(evdata.active==1 & evdata.expect==0 & evdata.choseRepeat==1));
behavProbs(2,2)=length(find(evdata.active==1 & evdata.expect==0 & evdata.choseRepeat==0));
behavProbs(3,1)=length(find(evdata.active==1 & evdata.expect==-1 & evdata.choseRepeat==1));
behavProbs(3,2)=length(find(evdata.active==1 & evdata.expect==-1 & evdata.choseRepeat==0));
behavProbs(4,1)=length(find(evdata.active~=1 & evdata.expect==1 & evdata.choseRepeat==1));
behavProbs(4,2)=length(find(evdata.active~=1 & evdata.expect==1 & evdata.choseRepeat==0));
behavProbs(5,1)=length(find(evdata.active~=1 & evdata.expect==0 & evdata.choseRepeat==1));
behavProbs(5,2)=length(find(evdata.active~=1 & evdata.expect==0 & evdata.choseRepeat==0));
behavProbs(6,1)=length(find(evdata.active~=1 & evdata.expect==-1 & evdata.choseRepeat==1));
behavProbs(6,2)=length(find(evdata.active~=1 & evdata.expect==-1 & evdata.choseRepeat==0));
behavProbs(:,3)=behavProbs(:,1)./sum(behavProbs(:,1:2),2)*100;
% sdt
behavProbs(:,3)=behavProbs(:,1)./sum(behavProbs(:,1:2),2)*100;
behavProbs(:,3)=behavProbs(:,1)./sum(behavProbs(:,1:2),2)*100;
behavProbs(:,3)=behavProbs(:,1)./sum(behavProbs(:,1:2),2)*100;

%------------

d500specs.upperLim=.66;
d500specs.lowerLim=.33;
behavProbsM=[];
behavProbsM(1,1)=length(find(evdata.active==1 & evdata.prep>d500specs.upperLim & evdata.choseRepeat==1));
behavProbsM(1,2)=length(find(evdata.active==1 & evdata.prep>d500specs.upperLim & evdata.choseRepeat==0));
behavProbsM(2,1)=length(find(evdata.active==1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.choseRepeat==1));
behavProbsM(2,2)=length(find(evdata.active==1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.choseRepeat==0));
behavProbsM(3,1)=length(find(evdata.active==1 & evdata.prep<d500specs.lowerLim & evdata.choseRepeat==1));
behavProbsM(3,2)=length(find(evdata.active==1 & evdata.prep<d500specs.lowerLim & evdata.choseRepeat==0));
behavProbsM(4,1)=length(find(evdata.active~=1 & evdata.prep>d500specs.upperLim & evdata.choseRepeat==1));
behavProbsM(4,2)=length(find(evdata.active~=1 & evdata.prep>d500specs.upperLim & evdata.choseRepeat==0));
behavProbsM(5,1)=length(find(evdata.active~=1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.choseRepeat==1));
behavProbsM(5,2)=length(find(evdata.active~=1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.choseRepeat==0));
behavProbsM(6,1)=length(find(evdata.active~=1 & evdata.prep<d500specs.lowerLim & evdata.choseRepeat==1));
behavProbsM(6,2)=length(find(evdata.active~=1 & evdata.prep<d500specs.lowerLim & evdata.choseRepeat==0));
behavProbsM(:,3)=behavProbsM(:,1)./sum(behavProbsM(:,1:2),2)*100;

%------------


h=figure; hold on
subplot(2,2,1) % 
bar(behavProbs(:,1:2),'group')
legend('chooseR','chooseA')
set(gca,'XTick',1:6,'XTickLabels',{'AexpR' 'AnoExp' 'AexpA' 'PexpR' 'PnoExp' 'PexpA'})
ylabel('Trial Counts','FontSize',12)
title('Generative Probability')

subplot(2,2,2) % 
bar(behavProbsM(:,3))
ylabel('% Chose Repeat','FontSize',12); ylim([0 100])
set(gca,'XTick',1:6,'XTickLabels',{'AexpR' 'AnoExp' 'AexpA' 'PexpR' 'PnoExp' 'PexpA'})

subplot(2,2,3) % 
bar(behavProbsM(:,1:2),'group')
legend('chooseR','chooseA')
set(gca,'XTick',1:6,'XTickLabels',{'AexpR' 'AnoExp' 'AexpA' 'PexpR' 'PnoExp' 'PexpA'})
ylabel('Trial Counts','FontSize',12)
title('Modelled Associative Strength Cutoffs')

subplot(2,2,4) % 
bar(behavProbsM(:,3))
ylabel('% Chose Repeat','FontSize',12); ylim([0 100])
set(gca,'XTick',1:6,'XTickLabels',{'AexpR' 'AnoExp' 'AexpA' 'PexpR' 'PnoExp' 'PexpA'})

savefig(h,[exptdata.figuredir500,'d500_behaviour_choiceProbability.fig']);
jpgfigname=[exptdata.figuredir500,'d500_behaviour_choiceProbability.jpg'];
print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure


%% 7. Signal Detection Theory Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sdtData_modeled=[];
% 1 hits                (actual rep - chose rep)
% 2 false alarms        (actual alt - chose rep)
% 3 misses 1-H          (actual rep - chose alt)
% 4 correct reject 1-F  (actual alt - chose alt)
% d-prime = [z(Hits) - z(False Alarms)]
% criterion

sdtData_modeled(1,1)=length(find(evdata.active==1 & evdata.prep>d500specs.upperLim & evdata.actual==1 & evdata.choseRepeat==1))  / length(find(evdata.active==1 & evdata.prep>d500specs.upperLim & evdata.actual==1))  ;
sdtData_modeled(1,2)=length(find(evdata.active==1 & evdata.prep>d500specs.upperLim & evdata.actual==-1 & evdata.choseRepeat==1)) / length(find(evdata.active==1 & evdata.prep>d500specs.upperLim & evdata.actual==-1))  ;
sdtData_modeled(1,3)=length(find(evdata.active==1 & evdata.prep>d500specs.upperLim & evdata.actual==1 & evdata.choseRepeat~=1))  / length(find(evdata.active==1 & evdata.prep>d500specs.upperLim & evdata.actual==1))  ;
sdtData_modeled(1,4)=length(find(evdata.active==1 & evdata.prep>d500specs.upperLim & evdata.actual==-1 & evdata.choseRepeat~=1)) / length(find(evdata.active==1 & evdata.prep>d500specs.upperLim & evdata.actual==-1))  ;
[sdtData_modeled(1,5),sdtData_modeled(1,6)]=dprime_simple(sdtData_modeled(1,1),sdtData_modeled(1,2)); % 

sdtData_modeled(2,1)=length(find(evdata.active==1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==1 & evdata.choseRepeat==1))  / length(find(evdata.active==1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==1))  ;
sdtData_modeled(2,2)=length(find(evdata.active==1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==-1 & evdata.choseRepeat==1)) / length(find(evdata.active==1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==-1))  ;
sdtData_modeled(2,3)=length(find(evdata.active==1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==1 & evdata.choseRepeat~=1))  / length(find(evdata.active==1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==1))  ;
sdtData_modeled(2,4)=length(find(evdata.active==1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==-1 & evdata.choseRepeat~=1)) / length(find(evdata.active==1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==-1))  ;
[sdtData_modeled(2,5),sdtData_modeled(2,6)]=dprime_simple(sdtData_modeled(2,1),sdtData_modeled(2,2)); % 

sdtData_modeled(3,1)=length(find(evdata.active==1 & evdata.prep<d500specs.lowerLim & evdata.actual==1 & evdata.choseRepeat==1))  / length(find(evdata.active==1 & evdata.prep<d500specs.lowerLim & evdata.actual==1))  ;
sdtData_modeled(3,2)=length(find(evdata.active==1 & evdata.prep<d500specs.lowerLim & evdata.actual==-1 & evdata.choseRepeat==1)) / length(find(evdata.active==1 & evdata.prep<d500specs.lowerLim & evdata.actual==-1))  ;
sdtData_modeled(3,3)=length(find(evdata.active==1 & evdata.prep<d500specs.lowerLim & evdata.actual==1 & evdata.choseRepeat~=1))  / length(find(evdata.active==1 & evdata.prep<d500specs.lowerLim & evdata.actual==1))  ;
sdtData_modeled(3,4)=length(find(evdata.active==1 & evdata.prep<d500specs.lowerLim & evdata.actual==-1 & evdata.choseRepeat~=1)) / length(find(evdata.active==1 & evdata.prep<d500specs.lowerLim & evdata.actual==-1))  ;
[sdtData_modeled(3,5),sdtData_modeled(3,6)]=dprime_simple(sdtData_modeled(3,1),sdtData_modeled(3,2)); % 


sdtData_modeled(4,1)=length(find(evdata.active~=1 & evdata.prep>d500specs.upperLim & evdata.actual==1 & evdata.choseRepeat==1))  / length(find(evdata.active~=1 & evdata.prep>d500specs.upperLim & evdata.actual==1))  ;
sdtData_modeled(4,2)=length(find(evdata.active~=1 & evdata.prep>d500specs.upperLim & evdata.actual==-1 & evdata.choseRepeat==1)) / length(find(evdata.active~=1 & evdata.prep>d500specs.upperLim & evdata.actual==-1))  ;
sdtData_modeled(4,3)=length(find(evdata.active~=1 & evdata.prep>d500specs.upperLim & evdata.actual==1 & evdata.choseRepeat~=1))  / length(find(evdata.active~=1 & evdata.prep>d500specs.upperLim & evdata.actual==1))  ;
sdtData_modeled(4,4)=length(find(evdata.active~=1 & evdata.prep>d500specs.upperLim & evdata.actual==-1 & evdata.choseRepeat~=1)) / length(find(evdata.active~=1 & evdata.prep>d500specs.upperLim & evdata.actual==-1))  ;
[sdtData_modeled(4,5),sdtData_modeled(4,6)]=dprime(sdtData_modeled(4,1),sdtData_modeled(4,2)); % 

sdtData_modeled(5,1)=length(find(evdata.active~=1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==1 & evdata.choseRepeat==1))  / length(find(evdata.active~=1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==1))  ;
sdtData_modeled(5,2)=length(find(evdata.active~=1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==-1 & evdata.choseRepeat==1)) / length(find(evdata.active~=1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==-1))  ;
sdtData_modeled(5,3)=length(find(evdata.active~=1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==1 & evdata.choseRepeat~=1))  / length(find(evdata.active~=1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==1))  ;
sdtData_modeled(5,4)=length(find(evdata.active~=1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==-1 & evdata.choseRepeat~=1)) / length(find(evdata.active~=1 & evdata.prep<d500specs.upperLim & evdata.prep>d500specs.lowerLim & evdata.actual==-1))  ;
[sdtData_modeled(5,5),sdtData_modeled(5,6)]=dprime(sdtData_modeled(5,1),sdtData_modeled(5,2)); % 

sdtData_modeled(6,1)=length(find(evdata.active~=1 & evdata.prep<d500specs.lowerLim & evdata.actual==1 & evdata.choseRepeat==1))  / length(find(evdata.active~=1 & evdata.prep<d500specs.lowerLim & evdata.actual==1))  ;
sdtData_modeled(6,2)=length(find(evdata.active~=1 & evdata.prep<d500specs.lowerLim & evdata.actual==-1 & evdata.choseRepeat==1)) / length(find(evdata.active~=1 & evdata.prep<d500specs.lowerLim & evdata.actual==-1))  ;
sdtData_modeled(6,3)=length(find(evdata.active~=1 & evdata.prep<d500specs.lowerLim & evdata.actual==1 & evdata.choseRepeat~=1))  / length(find(evdata.active~=1 & evdata.prep<d500specs.lowerLim & evdata.actual==1))  ;
sdtData_modeled(6,4)=length(find(evdata.active~=1 & evdata.prep<d500specs.lowerLim & evdata.actual==-1 & evdata.choseRepeat~=1)) / length(find(evdata.active~=1 & evdata.prep<d500specs.lowerLim & evdata.actual==-1))  ;
[sdtData_modeled(6,5),sdtData_modeled(6,6)]=dprime(sdtData_modeled(6,1),sdtData_modeled(6,2)); % 

h=figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Courier');
subplot(2,2,1) % hits / false alarms
bar(sdtData_modeled(1:3,1:2),'group')
set(gca,'XTick',1:3,'XTickLabels',{'expR','noExp','expA'},'FontSize',10)
ylabel('p(chooseR)','FontSize',12); legend('Hits','FAs')
title('Active: SDT Hits and False Alarms','FontSize',14,'FontWeight','Bold')

subplot(2,2,2) % hits / false alarms
bar(sdtData_modeled(1:3,5:6),'group')
set(gca,'XTick',1:3,'XTickLabels',{'expR','noExp','expA'},'FontSize',10)
ylabel('std units','FontSize',12); legend('Hits','FAs')
title('Active: d''prime & crit','FontSize',14,'FontWeight','Bold')

subplot(2,2,3) % hits / false alarms
bar(sdtData_modeled(4:6,1:2),'group')
set(gca,'XTick',1:3,'XTickLabels',{'expR','noExp','expA'},'FontSize',10)
ylabel('p(chooseR)','FontSize',12); legend('Hits','FAs')
title('Passive: SDT Hits and False Alarms','FontSize',14,'FontWeight','Bold')

subplot(2,2,4) % hits / false alarms
bar(sdtData_modeled(4:6,5:6),'group')
set(gca,'XTick',1:3,'XTickLabels',{'expR','noExp','expA'},'FontSize',10)
ylabel('std units','FontSize',12); legend('Hits','FAs')
title('Passive: d''prime & crit','FontSize',14,'FontWeight','Bold')

savefig(h,[exptdata.figuredir500,'d500_behaviour_sdt.fig']);
jpgfigname=[exptdata.figuredir500,'d500_behaviour_sdt.jpg'];
print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure

clear h* temp* tmp* ss t tr trialCount

%% SUMMARY
% Useful analyses include the SDT and regression analyses. They (mostly)
% show that there is a behavioural difference 
% They mostly show that there is indeed a behavioural effect:
% 1) p(repeat) influences choice probability
% 2) bias shifts with p(repeat)
% 3) active vs. passive shows, not surprisingly, different behaviours.
% All in all - this mostly confirms the monkeys were influenced by the task
% parameters.

d500specs.dateModified=datetime;
save([exptdata.projectdir,'DMS500_behavAnalysis.mat'],'behav*','Beta*','SessionNames','trialCounts','d500specs','evdata','sdt*','vers_ephys_es_rs'); 

disp('Finished!')
return


