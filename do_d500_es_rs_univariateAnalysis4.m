%% do_es_rs_univariateAnalysis1(vargin)
% by CS, AHB, started Nov 6, 2018
clear all
cd /Users/aharrisonbell/ePhys_Studies/ES-RS_Study/

vers_ephys_es_rs='1.4; Jan 30, 2019';
% 1.0 - original version (May 3, 2018)
% 1.1 - added model variables (Nov 6, 2018)
% 1.2 - added model variables (Dec 17, 2018)
% 1.3 - quality of life additions (Jan 15, 2019)

global exptdata
addpath(userpath);
lsn_ephys_analysis_defaults;

exptdata.analysisdir='/Users/aharrisonbell/ePhys_Studies/ES-RS_Study/'; % specifies overall directory for study (e.g., paper)


%% COMMENT OUT THE ONES YOU DON'T WANT
exptdata.monkeyname='both';
%exptdata.monkeyname='Vulcan';
%exptdata.monkeyname='BOTH MONKEYS';
%%
%% =====================================================================================================================
% megamatrices AND behaviouralAnalysis must be complete before running this program
%% What this program will do:

%% Note - from Feb 23, 2018
% 1. need to select visual neurons only (?).
% 2. need to add average SPDEN functions of visual neurons, possibly
%    subselecting for stimulus selectivity
%    Check Vogels' work to see if they subselect
% 3. still need to debug multiple regression
% 4. when creating figures, normalise 2ndPres to 1stPres (of same stimuli)


%% Need to add:
% regressor for active/passive
% maybe filter out just visually responsive neurons
% decoding analysis

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
% 20) ... lots of spike density stuff


%% HERE IS THE BIT OF CODE TO ADJUST THE MULTIREG ANALYSIS %%
%%% Possible regressors to include:
% data.stim1 (stimA=1; stimB=-1)
% data.stim2 (stimA=1; stimB=-1)
% data.session: session number
% data.repeat: (1=repeat, -1=alternation)
% data.cor: (correct=1, incorrect=0)
% data.active (1) or passive (0) block
% data.reactiontime: (in ms) ** not super reliable **
% data.expect: expectation (1=expect repeat; 2=expect alternation; 3=no expectation)
% data.correctRight = whether correct answer was on right (1) or left (-1)
% data.congruent_n1 (n minus 1) = correct answer is same for trial n and n-1
% data.congruent_n1 (n minus 2) = correct answer is same for trial n and n-2
% data.congruent_n3 (n minus 3) = correct answer is same for trial n and n-3


%% =====================================================================================================================
exptdata.analysisName='ES-RS_Study'; % used for savenames, figures, etc. (pick whaTEver you want; will be used for filenames)
exptdata.projectdir=[exptdata.analysisdir];

cprintf('*blue','*--------------------------------*\n')
cprintf('*blue','| do_es_rs_univariateAnalysis1.m |\n')
cprintf('*blue','*--------------------------------*\n')
cprintf('*red',['Version: ',vers_ephys_es_rs,'\n'])
disp(['Data location:    ',exptdata.datalocation]);
disp(['Processed NEV dir: ',exptdata.processedDatadir]);
disp(['Analysis Dir:      ',exptdata.analysisdir]);
disp(['Monkey name:      ',exptdata.monkeyname]);
disp(['Study name:       ',exptdata.analysisName]);
disp(['Project Directory: ',exptdata.projectdir]);

%% =====================================================================================================================
% load data
fprintf('\nLoading megamatrices and behavioural data...')
load([exptdata.analysisdir,'DMS500_processedbehavData.mat']); % contains BOTH monkeys behavioural data

%sessions = unique(V4_megaMatrix500NoScrubs(:,2));
if strcmp(exptdata.monkeyname,'Vortex')
    sessions = unique(data.session(data.animal==1)); % if Vortex
    load([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_megaMatrix_DMS500.mat']);
elseif strcmp(exptdata.monkeyname,'Vulcan')
    sessions = unique(data.session(data.animal==2)); % if Vulcan
    load([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_megaMatrix_DMS500.mat']);
else
    disp('combining monkey data')
    temp1=load([exptdata.projectdir,'Vortex_',exptdata.analysisName,'_megaMatrix_DMS500.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
    temp2=load([exptdata.projectdir,'Vulcan_',exptdata.analysisName,'_megaMatrix_DMS500.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
    V4_megaMatrix500NoScrubs=[temp1.V4_megaMatrix500NoScrubs; temp2.V4_megaMatrix500NoScrubs];
    TE_megaMatrix500NoScrubs=[temp1.TE_megaMatrix500NoScrubs; temp2.TE_megaMatrix500NoScrubs];
    clear temp*
    sessions = unique(data.session);
end
fprintf('done.\n')




%% Set up variables
% Data structures (DMS500): (-250:500 ms) (exptdata.xrange_psths)
timeEpochs500 = [...
    21 771;... aligned on FP on (1)
    772 1522;... aligned on 1stStim on (2)
    1523 2273;... aligned on 2ndStim on (3)
    2274 3024;... aligned on choice array (4)
    3025 3775;... aligned on correct choice (5)
    3776 4526;... aligned on incorrect choice (6)
    4528 5277]; % aligned on reward (7)
epochNames={'FPonset','1stStim','2ndStim','ChoiceArrayOnset','CorrectChoices','IncorrectChoices','Reward'};


%% Prep Spiking Data
BetaStuff=struct('epoch',[]);
brainArea = 2; % 1 = V4; 2 = TE

if brainArea==1
    BehavData=V4_megaMatrix500NoScrubs(:,1:20);
    NeuralData=V4_megaMatrix500NoScrubs(:,[773:1522 1524:2273 2275:3024]);
    NeuralData = reshape(NeuralData,[size(NeuralData,1),5,450]);
    NeuralData = double(squeeze(mean(NeuralData,2)));
    labelBrain='V4 (Normalised Spike Density)';
else
    BehavData=TE_megaMatrix500NoScrubs(:,1:20);
    NeuralData=TE_megaMatrix500NoScrubs(:,[773:1522 1524:2273 2275:3024]);
    NeuralData = reshape(NeuralData,[size(NeuralData,1),5,450]);
    NeuralData = double(squeeze(mean(NeuralData,2)));
    labelBrain='TE (Normalised Spike Density)';
end

%return


%% main analysis loop

close all;
clear bigdata glmdata;
ncount = 0;
nbins = 11;
clear indx*

options = statset('glmfit');
options.MaxIter = 100;

BehavData(:,4) = BehavData(:,3) + (BehavData(:,2)*100);

neuronmat = unique(BehavData(:,4));
sessmat = 1:26;

for n = 1:length(sessmat)
    
    if mod(n,100)==0
    fprintf([num2str(n),'.'])
    end
    
    indx_1 =find(BehavData(:,2) == sessmat(n) & BehavData(:,3)==1);
    indx = find(BehavData(:,2) == sessmat(n));
    
    cor = (BehavData(indx_1,6)==0);
    active = (BehavData(indx_1,15)==1);
    if sum(isnan(active))>0
        keyboard
    end
    
    rep = (BehavData(indx_1,14)==1);
    exp = rep.*0;
    exp(BehavData(indx_1,13)==1) = 0.75;  % expect rep
    exp(BehavData(indx_1,13)==2) = 0.25;  % expect alternate
    exp(BehavData(indx_1,13)==3) = 0.5;  % neutral

    
    [Eq, EH, KL]=BayesVolatilityChris(rep,exp,1/40,0);
    
    n_neurons = length(indx)./length(indx_1);
    
    rep = repmat(rep',[1 n_neurons])';
    Eq = repmat(Eq,[1 n_neurons]);
    active = repmat(active',[1 n_neurons])';
    cor = repmat(cor',[1 n_neurons])';

    % ERPs
    bigdata(n,1,:,1,1) = mean(NeuralData(indx(rep==1 & Eq' > 0.5 & active == 1 & cor == 1),:)',2);
    bigdata(n,2,:,1,1) = mean(NeuralData(indx(rep==0 & Eq' > 0.5 & active == 1 & cor == 1),:)',2);
    bigdata(n,1,:,2,1) = mean(NeuralData(indx(rep==1 & Eq' < 0.5 & active == 1 & cor == 1),:)',2);
    bigdata(n,2,:,2,1) = mean(NeuralData(indx(rep==0 & Eq' < 0.5 & active == 1 & cor == 1),:)',2);
    
    bigdata(n,1,:,1,2) = mean(NeuralData(indx(rep==1 & Eq' > 0.5 & active == 0 & cor == 1),:)',2);
    bigdata(n,2,:,1,2) = mean(NeuralData(indx(rep==0 & Eq' > 0.5 & active == 0 & cor == 1),:)',2);
    bigdata(n,1,:,2,2) = mean(NeuralData(indx(rep==1 & Eq' < 0.5 & active == 0 & cor == 1),:)',2);
    bigdata(n,2,:,2,2) = mean(NeuralData(indx(rep==0 & Eq' < 0.5 & active == 0 & cor == 1),:)',2);
    
%     for t = 1:300;
%     b = robustfit(logodds(Eq(rep==1 & active == 1))',NeuralData(indx(rep==1 & active == 1),t)');
%     regdata(n,1,t,1) = b(2);
%     b = robustfit(logodds(Eq(rep==0 & active == 1))',NeuralData(indx(rep==0 & active == 1),t)');
%     regdata(n,2,t,1) = b(2);
%     b = robustfit(logodds(Eq(rep==1 & active == 0))',NeuralData(indx(rep==1 & active == 0),t)');
%     regdata(n,1,t,2) = b(2);
%     b = robustfit(logodds(Eq(rep==0 & active == 0))',NeuralData(indx(rep==0 & active == 0),t)');
%     regdata(n,2,t,2) = b(2);
%     end
%     
end

disp('done')

return


%%

bigdata(:,:,:,:,1) = nanny(bigdata(:,:,:,:,1));
bigdata(:,:,:,:,2) = nanny(bigdata(:,:,:,:,2));

mresp = squeeze(mean(mean(mean(bigdata(:,:,50:100,:),2),3),4));
[i sorty] = sort(mresp,'descend');
sort_bigdata = bigdata(sorty,:,:,:,:);

figure;
for i = 1:2
    for j = 1:2;subplot(2,2,j+((i-1)*2));
        imagesc(squeeze(sort_bigdata(:,i,:,j,1)));
        set(gca,'clim',[0 20]);
    end
end

%%
%
sortindx = 1:26;

tp = [225 250];
tp = [400 450];

yy = 20;

close all
figure('color',[1 1 1]);

subplot(2,2,1)
hold on;
patch([tp(2) tp(1) tp(1) tp(2)],[0 0 yy yy],[0.9 0.9 0.9],'edgecolor',[0.9 0.9 0.9]);
fillsteplot(sort_bigdata(sortindx,:,:,1,1));
set(gca,'FontSize',16);
xlabel('time');
ylabel('spike rate');
title('active,expect rep');

subplot(2,2,2);
hold on
patch([tp(2) tp(1) tp(1) tp(2)],[0 0 yy yy],[0.9 0.9 0.9],'edgecolor',[0.9 0.9 0.9]);
fillsteplot(sort_bigdata(sortindx,:,:,2,1));
set(gca,'FontSize',16);
xlabel('time');
ylabel('spike rate');
title('active,expect alt');

subplot(2,2,3);
hold on;
patch([tp(2) tp(1) tp(1) tp(2)],[0 0 yy yy],[0.9 0.9 0.9],'edgecolor',[0.9 0.9 0.9]);
fillsteplot(sort_bigdata(sortindx,:,:,1,2));
set(gca,'FontSize',16);
xlabel('time');
ylabel('spike rate');
title('passive,expect rep');

subplot(2,2,4);
hold on;
patch([tp(2) tp(1) tp(1) tp(2)],[0 0 yy yy],[0.9 0.9 0.9],'edgecolor',[0.9 0.9 0.9]);
fillsteplot(sort_bigdata(sortindx,:,:,2,2));
set(gca,'FontSize',16);
xlabel('time');
ylabel('spike rate');
title('passive, exppect alt');


tmp = meeze(sort_bigdata(:,:,tp(1):tp(2),:,:),3);
disp('X1 = rep, X2 = Eq, X3 = active');
myanova(tmp);

return

% 
% [tt(:,1) p] = masst(r_bigdata(:,1,:,1)-r_bigdata(:,2,:,1));
% [tt(:,2) p] = masst(r_bigdata(:,1,:,2)-r_bigdata(:,2,:,2));
% 
% RSindex(:,1) = mean(bigdata(:,1,1001:1400,1)-bigdata(:,2,1001:1400,1),3);
% RSindex(:,2) = mean(bigdata(:,1,1001:1400,2)-bigdata(:,2,1001:1400,2),3);
% 
% figure('color',[1 1 1]);
% subplot(1,2,1);
% plot(log(mresp),RSindex(:,1),'ro');
% 
% subplot(1,2,2);
% plot(log(mresp),RSindex(:,2),'ro');
% ylim
% %
% figure('color',[1 1 1]);
% hold on;
% plot(tt(:,1),'color','r','linewidth',3);
% plot(tt(:,2),'color',[0 0.8 0],'linewidth',3);
% t_thresh = tinv(0.001,size(bigdata,1));
% line([0 size(tt,1)],[t_thresh t_thresh],'color','k','linestyle','--');
% xlabel('time');
% ylabel('t-value');
% set(gca,'FontSize',16);