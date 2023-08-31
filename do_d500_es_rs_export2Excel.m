% do_es_rs_export2Excel;
% by AHB, started Feb 13, 2019
vers_ephys_es_rs='1.0; Feb 13, 2019';
% 1.0 - original version (Feb 13, 2019)
% 1.1 - (Feb 25, 2019) - added SYNCing FROM Excel
global exptdata
addpath(userpath);
ephys_analysis_defaults;
load([exptdata.analysisdir,filesep,exptdata.analysisName,'_exptdata.mat'],'exptdata');

exptdata.monkeyname='both';
syncDir='READ'; % 'WRITE','BOTH','READ'

%% =====================================================================================================================
% megamatrices AND behaviouralAnalysis must be complete before running this program
%% What this program will do:
% 1) Exports data to Excel Spreadsheet (ES_RS), including Session Info, Neuron Info
% 2) Retrieves data from Excel Spreadsheet - specifically Quality and Manual Classifications


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
% 20) Auto Classification
% 21) Manual Classification (1=excitatory, -1=suppressed, 0=non-responsive)
% 22) Quality Classification (0=shit, problematic, or non-responsive, 1=meh, 2=ok, 3=awesome)


% =====================================================================================================================
exptdata.projectdir=[exptdata.analysisdir,exptdata.analysisName,filesep];
clc
cprintf('*yellow','*-------------------------*\n')
cprintf('*yellow','| do_es_rs_export2excel.m |\n')
cprintf('*yellow','*-------------------------*\n')
cprintf('*red',['Version: ',vers_ephys_es_rs,'\n'])
disp(['Data location:    ',exptdata.datalocation]);
disp(['Processed NEV dir: ',exptdata.processedDatadir]);
disp(['Analysis Dir:      ',exptdata.analysisdir]);
disp(['Monkey name:      ',exptdata.monkeyname]);
disp(['Study name:       ',exptdata.analysisName]);
disp(['Project Directory: ',exptdata.projectdir]);


if strcmp(syncDir,'WRITE') || strcmp (syncDir,'BOTH')
    
    %% =====================================================================================================================
    % load spike density data
    fprintf('\nLoading megamatrices and behavioural data...\n')
    load([exptdata.analysisdir,'DMS500_behavAnalysis.mat']); % contains BOTH monkeys behavioural data
    
    if strcmp(exptdata.monkeyname,'Vortex')
        sessions = unique(data.session(data.animal==1)); % if Vortex
        load([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_megaMatrix.mat']);
    elseif strcmp(exptdata.monkeyname,'Vulcan')
        sessions = unique(data.session(data.animal==2)); % if Vulcan
        load([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_megaMatrix.mat']);
    else
        disp('...Combining monkey data...')
        temp1=load([exptdata.projectdir,'Vortex_',exptdata.analysisName,'_megaMatrix.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
        temp2=load([exptdata.projectdir,'Vulcan_',exptdata.analysisName,'_megaMatrix.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
        V4_megaMatrix500NoScrubs=[temp1.V4_megaMatrix500NoScrubs; temp2.V4_megaMatrix500NoScrubs];
        TE_megaMatrix500NoScrubs=[temp1.TE_megaMatrix500NoScrubs; temp2.TE_megaMatrix500NoScrubs];
        clear temp*
        sessions = unique(data.session);
    end
    fprintf('done.\n')
    
    %% Begin Brain Loops
    for brainArea=1:2
        if brainArea==1
            AllspikeData=V4_megaMatrix500NoScrubs;
            labelBrain='V4';
        else
            AllspikeData=TE_megaMatrix500NoScrubs;
            labelBrain='TE';
        end
        xlMatrix=[];
        disp(['Processing data from ',labelBrain])
        for sess=1:length(sessions)
            clear numNeuronsSession
            fprintf(['...Session ',num2str(sess),'/',num2str(length(sessions)),': '])
            numNeuronsSession=max(AllspikeData(AllspikeData(:,2)==sessions(sess),10)); % how many neurons in this unit (this code is a bit of a hack but does the trick)
            tempVisualBad=0; tempVisualGood=0; tempVisualBoring=0;
            
            for nn=1:numNeuronsSession
                fprintf([num2str(nn),'.'])
                % Determine Neuron Type
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
                tempBehav=AllspikeData(indxN,6);
                
                temp_xlMatrix(1) = mean(AllspikeData(indxN,1)); % animal
                temp_xlMatrix(2) = mean(AllspikeData(indxN,2)); % session number
                temp_xlMatrix(3) = nan;                         % add later (session name)
                temp_xlMatrix(4) = mean(AllspikeData(indxN,3)); % neuron session number
                temp_xlMatrix(5) = mean(AllspikeData(indxN,4)); % neuron cohort number
                temp_xlMatrix(6) = length(find(tempBehav==0))/length(tempBehav); % behavioural performance
                temp_xlMatrix(7) = mean(AllspikeData(indxN,20)); % neuron type (automatic)
                temp_xlMatrix(8) = nan; % neuron type (confirmed)
                temp_xlMatrix(9) = nan; % quality
                
                xlMatrix=[xlMatrix; temp_xlMatrix]; %#ok<AGROW>
                clear temp_xlMatrix tempBehav
            end% neuron
            fprintf('Done\n')
            
        end % session
        [r,~]=size(xlMatrix);
        xlwrite(exptdata.excel_ES_RS,xlMatrix,labelBrain,['A2:I',num2str(r+1)])
        clear r c xlMatrix labelBrain
    end % brain area
    
end % syncDir

%% =====================================================================================================================
if strcmp(syncDir,'READ') || strcmp (syncDir,'BOTH')
    fprintf('\nReading data from Excel Spreadsheet (Manual Classifications and Quality Ratings)...')
    for brainArea=1:2
        if brainArea==1
            labelBrain='V4';
        else
            labelBrain='TE';
        end
        clear temp*
        [tempNum,~,tempRaw]=xlsread(exptdata.excel_ES_RS,labelBrain,'A2:L2000');
        
        if brainArea==1
            ESRS_NeuronData_V4.ClassiferData=tempNum;
            ESRS_NeuronData_V4.RawData=tempRaw;
        else
            ESRS_NeuronData_TE.ClassiferData=tempNum;
            ESRS_NeuronData_TE.RawData=tempRaw;
        end
    end % brain Area
    save([exptdata.projectdir,'es_rs_NeuronalData.mat'],'ESRS_NeuronData_V4','ESRS_NeuronData_TE')
    fprintf('done.\n\n')

    %% ==================================================================================================================================================
    % Paste data into AllSpikes and save again...
    % Load individual monkey data and paste in new data
    fprintf('\nInserting Manual Classifications, Quality Ratings, and Stimulus Identity into Massive Matrices...\n')
    temp1=load([exptdata.projectdir,'Vortex_',exptdata.analysisName,'_megaMatrix.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
    temp2=load([exptdata.projectdir,'Vulcan_',exptdata.analysisName,'_megaMatrix.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs');
    V4_megaMatrix500NoScrubsALL=[temp1.V4_megaMatrix500NoScrubs; temp2.V4_megaMatrix500NoScrubs];
    TE_megaMatrix500NoScrubsALL=[temp1.TE_megaMatrix500NoScrubs; temp2.TE_megaMatrix500NoScrubs];
    clear temp*
    
    % for now, scroll through each session, etc.
    % V4
    sessions=unique(V4_megaMatrix500NoScrubsALL(:,2));
    for sess=1:length(unique(sessions))
        fprintf(['...V4 Session: ',num2str(sessions(sess)),'\n'])
        tempNeuronsSession=unique(V4_megaMatrix500NoScrubsALL(V4_megaMatrix500NoScrubsALL(:,2)==sessions(sess),3)); % select unique numbers PER session
        for nn=1:length(tempNeuronsSession)
            % Paste manual classification and quality
            tempVals=ESRS_NeuronData_V4.ClassiferData(ESRS_NeuronData_V4.ClassiferData(:,2)==sessions(sess) & ESRS_NeuronData_V4.ClassiferData(:,4)==tempNeuronsSession(nn),7:9);
            tempRowCounter=length(V4_megaMatrix500NoScrubsALL(V4_megaMatrix500NoScrubsALL(:,2)==sessions(sess) & V4_megaMatrix500NoScrubsALL(:,3)==tempNeuronsSession(nn),21));
            tempPaste=repmat(tempVals,tempRowCounter,1);
            V4_megaMatrix500NoScrubsALL(V4_megaMatrix500NoScrubsALL(:,2)==sessions(sess) & V4_megaMatrix500NoScrubsALL(:,3)==tempNeuronsSession(nn),20:22)=tempPaste;
            clear tempVals tempRowCounter tempPaste
        end
        tempStim=exptdata.dms500.stimIDs(exptdata.dms500.stimIDs(:,1)==sessions(sess),2:5);
        tempRowCounter=length(V4_megaMatrix500NoScrubsALL(V4_megaMatrix500NoScrubsALL(:,2)==sessions(sess),23));
        tempPaste=repmat(tempStim,tempRowCounter,1);
        V4_megaMatrix500NoScrubsALL(V4_megaMatrix500NoScrubsALL(:,2)==sessions(sess),23:26)=tempPaste;
        clear tempStim tempNeuronsSession tempRowCounter tempPaste
    end
    
    % TE (need to add once V4 is ok)
    sessions=unique(TE_megaMatrix500NoScrubsALL(:,2));
    for sess=1:length(unique(sessions))
        fprintf(['...TE Session: ',num2str(sessions(sess)),'\n'])
        tempNeuronsSession=unique(TE_megaMatrix500NoScrubsALL(TE_megaMatrix500NoScrubsALL(:,2)==sessions(sess),3)); % select unique numbers PER session
        for nn=1:length(tempNeuronsSession)
            % Paste manual classification and quality
            tempVals=ESRS_NeuronData_TE.ClassiferData(ESRS_NeuronData_TE.ClassiferData(:,2)==sessions(sess) & ESRS_NeuronData_TE.ClassiferData(:,4)==tempNeuronsSession(nn),7:9);
            tempRowCounter=length(TE_megaMatrix500NoScrubsALL(TE_megaMatrix500NoScrubsALL(:,2)==sessions(sess) & TE_megaMatrix500NoScrubsALL(:,3)==tempNeuronsSession(nn),21));
            tempPaste=repmat(tempVals,tempRowCounter,1);
            TE_megaMatrix500NoScrubsALL(TE_megaMatrix500NoScrubsALL(:,2)==sessions(sess) & TE_megaMatrix500NoScrubsALL(:,3)==tempNeuronsSession(nn),20:22)=tempPaste;
            clear tempRowPointer tempRowCounter tempPaste
        end
        tempStim=exptdata.dms500.stimIDs(exptdata.dms500.stimIDs(:,1)==sessions(sess),2:5);
        tempRowCounter=length(TE_megaMatrix500NoScrubsALL(TE_megaMatrix500NoScrubsALL(:,2)==sessions(sess),23));
        tempPaste=repmat(tempStim,tempRowCounter,1);
        TE_megaMatrix500NoScrubsALL(TE_megaMatrix500NoScrubsALL(:,2)==sessions(sess),23:26)=tempPaste;
        clear tempStim tempNeuronsSession tempRowCounter tempPaste
    end
    
    fprintf('...done\n\n')
    fprintf('Saving MegaMatrices...')
    % monkey 1
    V4_megaMatrix500NoScrubs=V4_megaMatrix500NoScrubsALL(V4_megaMatrix500NoScrubsALL(:,1)==1,:); % select only monkey
    TE_megaMatrix500NoScrubs=TE_megaMatrix500NoScrubsALL(TE_megaMatrix500NoScrubsALL(:,1)==1,:); %#ok<*NASGU>
    save([exptdata.projectdir,'Vortex_',exptdata.analysisName,'_megaMatrix_synced.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs','-v7.3');
    % monkey 2
    V4_megaMatrix500NoScrubs=V4_megaMatrix500NoScrubsALL(V4_megaMatrix500NoScrubsALL(:,1)==2,:);
    TE_megaMatrix500NoScrubs=TE_megaMatrix500NoScrubsALL(TE_megaMatrix500NoScrubsALL(:,1)==2,:);
    save([exptdata.projectdir,'Vulcan_',exptdata.analysisName,'_megaMatrix_synced.mat'],'V4_megaMatrix500NoScrubs','TE_megaMatrix500NoScrubs','-v7.3');
    fprintf('done.\n\n')
end % syncDir
return