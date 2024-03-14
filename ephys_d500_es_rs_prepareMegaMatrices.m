% ephys_d500_es_rs_prepareMegaMatrices
% AHB, Feb 22, 2018
clear
vers_ephys_es_rs='2.2; May 15, 2023';
% 1.0 - original version (Feb 22, 2018)
% 1.1 - session numbers are reset in this program (June 5, 2018)
% 1.2 - split program so it deals with one task at a time (for memory management reasons) (June 5, 2018)
% 1.3 - fixed up some bugs re: behavioural data (Nov 6, 2018)
% 2.0 - DMS500 - now loads files according to order and list specified by output of compileBehavData (so that session numbers - and hopefully trial 
% numbers match across data modalities)
% 2.1 - added code to compile spike train data
% This program collects session data
% 2.2 - updated to run using modern functions and on PC (HOMER_B550) May 15, 2023
% exptdata.reprocess=1; % set to one if you want to regenerate megaMatrices

% Notes - Feb 23, 2018
% 1. need to select visual neurons only (?).
% 2. need to add average SPDEN functions of visual neurons, possibly
% subselecting for stimulus selectivity
% Check Vogels' work to see if they subselect
% 3. still need to debug multiple regression
% 4. when creating figures, normalise 2ndPres to 1stPres (of same stimuli)

% Analyse Behaviour
% Blocks:
% 1) Expect Repeat
% 2) Expect Alternation
% 3) No Repeat
% 4) Expect Repeat
% 5) Expect Alternation
% 6) No Repeat

% Trial Outcome Markers
% 0 Correct
% 1 No Response
% 2 Late Response
% 3 Break Fixation
% 4 No Fixation
% 5 Early Response
% 6 Incorrect Response
% 7 Lever Break
% 8 Ignored
% 9 Aborted

% choice probability
%behavData.TrialError 0=correct, 4=no fixation, 6=incorrect response
%behav
clc
global exptdata
fprintf('*-------------------------------------------*\n')
fprintf('| ephys_d500_es_rs_prepareMegaMatrices.m    |\n')
fprintf('*-------------------------------------------*\n')
fprintf(['Version: ',vers_ephys_es_rs,'\n'])
disp(['Data location:    ',exptdata.datalocation]);
disp(['Processed NEV dir: ',exptdata.processedDatadir]);
disp(['Project Directory: ',exptdata.projectdir]);
disp(' ')
disp(['Monkey name:      ',exptdata.monkeyname]);
disp(['Study name:       ',exptdata.analysisName]);


clear ddir dd df V4* TE* ll ss *megaMatrix* *header* *300 *500
if exist([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_megaMatrix.mat'],'file') && ...
        exist([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_megaMatrix_spikes.mat'],'file') && exptdata.reprocess~=1
    fprintf('\n*** 500-Series MegaMatrix already exist. Skipping...\n')
else
    %% 1. Load data and prepare matrices
    load([exptdata.projectdir,'DMS500_behavAnalysis.mat']);
    clear temp*
    fprintf('* Loading DMS500 data from individual sessions...\n')
    V4_megaMatrix500=[]; TE_megaMatrix500=[];
    V4_megaMatrix500_spikes=[]; TE_megaMatrix500_spikes=[];

    for df=1:length(SessionNames) % scroll through each file
        % Here is where the code could be switched to include both animals in one GIGANTIC matrix
        if contains(SessionNames(df),exptdata.monkeyname)
            clear temp*
            fprintf(['  ...Session Number ',num2str(df),': B:',SessionNames{df},'\t N: '])
            tempPointer=cellfun(@(s) ~isempty(strfind(SessionNames{df}, s)), exptdata.datafiles(:,2));
            tempFilenameN=exptdata.datafiles{tempPointer,1};
            tempFilenameP=[exptdata.monkeyname,'_',char(datetime(tempFilenameN(1:8),'InputFormat','ddMMyyyy','Format','yyyy-MM-dd')),'_DMS500_sessionData.mat'];
            fprintf([tempFilenameN,'\t P: ',tempFilenameP,'\n']);
            load([exptdata.projectdir,filesep,'sessionData',filesep,tempFilenameP]);
            tempFilenameP=[exptdata.monkeyname,'_',char(datetime(tempFilenameN(1:8),'InputFormat','ddMMyyyy','Format','yyyy-MM-dd')),'_DMS500_sessionData_spikes.mat'];
            load([exptdata.projectdir,filesep,'sessionData',filesep,tempFilenameP]);
            V4_bigMatrix(:,2)=df;
            TEp_bigMatrix(:,2)=df;
            V4_bigMatrix_spikes(:,2)=df;
            TEp_bigMatrix_spikes(:,2)=df;
            V4_megaMatrix500=[V4_megaMatrix500; single(V4_bigMatrix)]; %#ok<AGROW>
            TE_megaMatrix500=[TE_megaMatrix500; single(TEp_bigMatrix)]; %#ok<AGROW>
            V4_megaMatrix500_spikes=[V4_megaMatrix500_spikes; single(V4_bigMatrix_spikes)]; %#ok<AGROW>
            TE_megaMatrix500_spikes=[TE_megaMatrix500_spikes; single(TEp_bigMatrix_spikes)]; %#ok<AGROW>
            clear V4_bigMatrix* TEp_bigMatrix* temp*
        end
    end
    fprintf('* Done.\n')

    %% 2. Need to clear out NANs
    for ll=1:size(V4_megaMatrix500,1)
        if ~isempty(find(isnan(V4_megaMatrix500(ll,21:end)==1), 1)) % if row contains NaNs kill it (set condition number to 1000)
            V4_megaMatrix500(ll,8)=1000; disp(['* (V4-dms500) killing nan data row: ',num2str(ll)]) %#ok<*SAGROW>
        end
    end
    for ll=1:size(TE_megaMatrix500,1)
        if ~isempty(find(isnan(TE_megaMatrix500(ll,21:end)==1), 1)) % if row contains NaNs kill it (set condition number to 1000)
            TE_megaMatrix500(ll,8)=1000; disp(['* (TE-dms500) killing nan data row: ',num2str(ll)])
        end
    end
    for ll=1:size(V4_megaMatrix500_spikes,1)
        if ~isempty(find(isnan(V4_megaMatrix500_spikes(ll,21:end)==1), 1)) % if row contains NaNs kill it (set condition number to 1000)
            V4_megaMatrix500_spikes(ll,8)=1000; disp(['* (V4-dms500_spikes) killing nan data row: ',num2str(ll)])
        end
    end
    for ll=1:size(TE_megaMatrix500_spikes,1)
        if ~isempty(find(isnan(TE_megaMatrix500_spikes(ll,21:end)==1), 1)) % if row contains NaNs kill it (set condition number to 1000)
            TE_megaMatrix500_spikes(ll,8)=1000; disp(['* (TE-dms500_spikes) killing nan data row: ',num2str(ll)])
        end
    end

    %% 3. Trim the fat (select only correct and incorrect trials)
    fprintf('* Clearing out bad trials\n')
    V4_megaMatrix500NoScrubs=V4_megaMatrix500;
    V4_megaMatrix500NoScrubs=V4_megaMatrix500NoScrubs(ismember(V4_megaMatrix500NoScrubs(:,6),[0 6]) & V4_megaMatrix500NoScrubs(:,8)<193,:);
    TE_megaMatrix500NoScrubs=TE_megaMatrix500;
    TE_megaMatrix500NoScrubs=TE_megaMatrix500NoScrubs(ismember(TE_megaMatrix500NoScrubs(:,6),[0 6]) & TE_megaMatrix500NoScrubs(:,8)<193,:);

    V4_megaMatrix500NoScrubs_spikes=V4_megaMatrix500_spikes;
    V4_megaMatrix500NoScrubs_spikes=V4_megaMatrix500NoScrubs_spikes(ismember(V4_megaMatrix500NoScrubs_spikes(:,6),[0 6]) & V4_megaMatrix500NoScrubs_spikes(:,8)<193,:);
    TE_megaMatrix500NoScrubs_spikes=TE_megaMatrix500_spikes;
    TE_megaMatrix500NoScrubs_spikes=TE_megaMatrix500NoScrubs_spikes(ismember(TE_megaMatrix500NoScrubs_spikes(:,6),[0 6]) & TE_megaMatrix500NoScrubs_spikes(:,8)<193,:);

    clear V4_megaMatrix500 TE_megaMatrix500 V4_megaMatrix500_spikes TE_megaMatrix500_spikes

    %% 4. Add trial classifers
    fprintf('* Add trial classifers\n')
    V4header500=[]; TEheader500=[]; %#ok<NASGU>
    V4header500=ephys_addTrialClassifers(V4_megaMatrix500NoScrubs(:,1:20),'DMS500');
    TEheader500=ephys_addTrialClassifers(TE_megaMatrix500NoScrubs(:,1:20),'DMS500');
    numSessions500=unique(V4header500(:,2));

    % DMS500 - New Unit Numbers
    tic
    for ss=1:length(numSessions500)
        tempUnitPerSessionV4(ss)=length(unique(V4header500(V4header500(:,2)==ss,3)));
        tempUnitPerSessionTE(ss)=length(unique(TEheader500(TEheader500(:,2)==ss,3)));
    end
    toc
    for ss=1:length(numSessions500)
        if ss==1
            V4header500(V4header500(:,2)==ss,4)=V4header500(V4header500(:,2)==ss,3);
            TEheader500(TEheader500(:,2)==ss,4)=TEheader500(TEheader500(:,2)==ss,3);
        else
            V4header500(V4header500(:,2)==ss,4)=V4header500(V4header500(:,2)==ss,3)+sum(tempUnitPerSessionV4(1:ss-1));
            TEheader500(TEheader500(:,2)==ss,4)=TEheader500(TEheader500(:,2)==ss,3)+sum(tempUnitPerSessionTE(1:ss-1));
        end
    end
    V4_megaMatrix500NoScrubs(:,1:20)=V4header500;
    TE_megaMatrix500NoScrubs(:,1:20)=TEheader500;
    V4_megaMatrix500NoScrubs_spikes(:,1:20)=V4header500;
    TE_megaMatrix500NoScrubs_spikes(:,1:20)=TEheader500;

    clear temp*

    % Review Neurons (scroll through individual neuron outputs)
    % DMS500
    temp_windowGap500=21:length(exptdata.xrange_psths):length(exptdata.xrange_psths)*8; % we want the third
    if exptdata.reviewNeurons==1
        ephys_d500_es_rs_plotDMS500_neurons(V4_megaMatrix500NoScrubs,V4header500,temp_windowGap500,'dms500V4')
        ephys_d500_es_rs_plotDMS500_neurons(TE_megaMatrix500NoScrubs,TEheader500,temp_windowGap500,'dms500TE')
    end

    % Classify Neurons and generate univariate descriptive statistics (e.g., peak firing rates, etc.)
    fprintf('* Classify neurons\n')
    [V4_neuronData500,V4header_new500]=ephys_sortDMS500_neurons(V4_megaMatrix500NoScrubs,V4header500,temp_windowGap500,'V4',0);
    [TE_neuronData500,TEheader_new500]=ephys_sortDMS500_neurons(TE_megaMatrix500NoScrubs,TEheader500,temp_windowGap500,'TE',0);

    %



    % SHRINK MEGAMATRICES
    V4_megaMatrix500NoScrubs=V4_megaMatrix500NoScrubs(:,1:3024); % eliminate spden after choice array
    V4_megaMatrix500NoScrubs_spikes=V4_megaMatrix500NoScrubs_spikes(:,1:3024); % eliminate spden after choice array

    TE_megaMatrix500NoScrubs=TE_megaMatrix500NoScrubs(:,1:3024); % eliminate spden after choice array
    TE_megaMatrix500NoScrubs_spikes=TE_megaMatrix500NoScrubs_spikes(:,1:3024); % eliminate spden after choice array

    disp('* Saving MEGAMATRICES...')
    save([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_megaMatrix.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs',...
        'V4_neuronData500','V4header_new500','TE_neuronData500','TEheader_new500','-v7.3');
    save([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_megaMatrix_spikes.mat'],'V4_megaMatrix500NoScrubs_spikes','TE_megaMatrix500NoScrubs_spikes','-v7.3');
    fprintf('Done.\n')
    clear temp*
end


return

















