function ndist=ephys_multipleRegression_shuffle(test_evs,test_ev_labels,test_modelName)
%% This version will shuffle the classification indices and run through multiple times
% The ndist will be a null distribution of betaWeights

global timeEpochs500 sampleRate epochNum evdata currentSpikeData currentBrainLabel exptdata


% NEW WAY:
% Contact jules re: run
% instead of running glm on every neuron, every session - select at random
% 100 neurons across the population and run on those only:
% tmpAll_uniqueNeurons=unique([sessionNum neuronNum])
% = datasample(tmpAll,numSamplesPerShuffle)

    



%% Specify number of shuffles
numShuffles=100;
numSamplesPerShuffle=100;


if nargin<4
    regressType='linear';
end

sessions=unique(currentSpikeData(:,2)); % identify individual sessions
warning off

ndist=struct('model',struct('test_evs',test_evs,'test_ev_labels',test_ev_labels,'test_modelName',test_modelName),...
    'sessions',struct('betaWeights',[]));
ndist_thresholds=nan(numShuffles,length(test_ev_labels)+1);


tmp_unique([currentSpikeData(:,2) currentSpikeData(:,3)],'rows')


for shuff=1:numShuffles
    fprintf(['Shuffling Model "',test_modelName,'" (',currentBrainLabel,') SHUFFLE ',num2str(shuff),'of',num2str(numShuffles),' => Session: ']);
    for sess=1:length(sessions)
        clear tmpSess_*
        fprintf([num2str(sessions(sess)),'.'])
        tmpSess_indxB=find(evdata.session == sessions(sess));
        tmpSess_ev=test_evs(tmpSess_indxB,:);
        
        %%
        tmpSess_ev=shuffle(tmpSess_ev);
        %%
        
        tmpSess_numNeuronsSession=max(currentSpikeData(currentSpikeData(:,2)==sessions(sess),10)); % how many neurons for this session?
        for nn=1:tmpSess_numNeuronsSession
            tmpSess_indxN=find(currentSpikeData(:,2) == sessions(sess) & currentSpikeData(:,3) == nn & ismember(currentSpikeData(:,5),evdata.trialNumber(tmpSess_indxB)));
            if ~isempty(tmpSess_indxN)
                tmpSess_newTimeRange=timeEpochs500(epochNum,1):sampleRate:timeEpochs500(epochNum,2);
                % Check to see if length of EV matches data
                if size(tmpSess_ev,1)~=length(tmpSess_indxN)
                    tmpSess_ev=tmpSess_ev(1:length(tmpSess_indxN),:); % fix needed in case behav and neuro trials don't match
                    fprintf('Prob! ')
                end
                % Run actual regression
                for tt=1:length(tmpSess_newTimeRange)
                    if strcmp(regressType,'linear')
                        ndist(shuff).sessions(sess).betaWeights(nn,tt,:)=...
                            glmfit(tmpSess_ev, currentSpikeData(tmpSess_indxN,tmpSess_newTimeRange(tt)));
                        %BetaStuff(brainArea).sessions(sess).epoch(epochNum).betaWeights(nn,tt,:)= glmfit(predictrz, ztransfnan_ahb(AllspikeData(indx,tt)))';
                    elseif strcmp(regressType,'logistic')
                        % NEED TO FIX
                    elseif strcmp(regressType,'stepwise')
                        % NEED TO FIX
                    end
                end
            end
            clear tt 
        end
    end
end

for shuff=1:numShuffles
    %% Calculate statistics for each shuffle
    % Prep Distribution Data for this shuffle
    % ** This treats unit of replication as SESSION **
    
    % maybe instead of having unit of replication be SESSION, instead
    % randomly select 50-100 neurons across all sessions? OR 5 from each
    % session?
    
    tmp_shuffData=[]; tmp_shuff_Stuff=[]; %#ok<NASGU>
    for sess=1:length(sessions)
        tmp_shuff_Stuff=squeeze(mean(ndist(shuff).sessions(sess).betaWeights,1))'; %
        for ev=1:length(test_ev_labels)+1 % scroll through each regressor including the intercept
            tmp_shuffData(sess,ev,:)=tmp_shuff_Stuff(ev,:);  %#ok<AGROW>
        end
    end
    
    % Testing Significance
    for ev=1:length(test_ev_labels)+1
        [~, tmp_shuff_pval]=masst(squeeze(tmp_shuffData(:,ev,:)));
        
        tmpSig=tmp_shuff_pval(:,ev)<0.025 | tmp_shuff_pval(:,ev) > 0.975;
        tmpcounter=0; tmpbig_counter=0; % default to cluster of one
        
        % Chris: this is a dumb way to do it, but can't remember the clever one
        % Andy: seems to work.  Goes through each timepoints and tallys whether it
        % is significant or not. If the next one is also sig, the counter goes up.
        % Counter counts each instance, "bigcounter" tracks the current largest
        % cluster, which will then serve as the threshold.
        % (method of Maris & Oostenveld, 2007 J Neurosci Methods)
        
        for tp=1:size(tmp_shuff_pval,1)
            if tmpSig(tp)==1 && tmpcounter==0
                tmpcounter=1;
            elseif tmpSig(tp)==1 && tmpcounter>0
                tmpcounter=tmpcounter+1;
            elseif tmpSig(tp)==0 && tmpcounter>0
                if tmpcounter>=tmpbig_counter
                    tmpbig_counter=tmpcounter;
                end
                tmpcounter=0;
            elseif tmpSig(tp)==0 && tmpcounter==0
                tmpcounter=0;
            end
        end
        ndist_thresholds(shuff, ev)=tmpbig_counter;
        
    end
    
    fprintf(') => done\n')
end


%% Save Data
save([exptdata.projectdir,'d500_nullDist_model',test_modelName,'.mat'],'ndist','nullDist','ndist_tval','ndist_pval','ndist_threshold')

%% Plot Null Distribution Data
% Figure
figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial')
for sp=1:length(test_ev_labels)+1 % create subplot for each regressor
    if sp==1
        subplot(length(test_ev_labels)+1,1,sp);  hold on
        plot(squeeze(mean(ndist.nullDist(:,1,:),1))','k-','LineWidth',1); % intercept
        plot([750/sampleRate 750/sampleRate],[0 25],'k-','LineWidth',1.5)
        plot([1500/sampleRate 1500/sampleRate],[0 25],'k-','LineWidth',1.5)
        plot([250/sampleRate 250/sampleRate],[0 25],'k:','LineWidth',1.5)
        plot([1000/sampleRate 1000/sampleRate],[0 25],'k:','LineWidth',1.5)
        plot([1750/sampleRate 1750/sampleRate],[0 25],'k:','LineWidth',1.5)
        fill([25 25 50 50],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
        fill([100 100 125 125],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
        fill([175 175 200 200],[0 25 25 0],'k','Edgecolor','none','FaceAlpha',0.3)
        title({['MODEL: ',test_modelName],['Intercept: Area ',currentBrainLabel,' - numSessions: ',num2str(length(sessions))]},'FontSize',16); ylim([0 20])
        xlabel(['Time from Event (downsampled to x/',num2str(sampleRate),' ms, 1st Stim | 2nd Stim | Choice Array )']);
        ylabel('Parameter Estimate')
    else

        subplot(length(test_ev_labels)+1,1,sp);  hold on
        fillsteplot(squeeze(ndist.nullDist(:,sp,:)),2,'-','v1',1:size(squeeze(ndist.nullDist(:,sp,:)),2),{'r'});
        plot([750/sampleRate 750/sampleRate],[-25 25],'k-','LineWidth',1.5)
        plot([1500/sampleRate 1500/sampleRate],[-25 25],'k-','LineWidth',1.5)
        plot([250/sampleRate 250/sampleRate],[-25 25],'k:','LineWidth',1.5)
        plot([1000/sampleRate 1000/sampleRate],[-25 25],'k:','LineWidth',1.5)
        plot([1750/sampleRate 1750/sampleRate],[-25 25],'k:','LineWidth',1.5)
        fill([25 25 50 50],[-5 5 5 -5],'k','Edgecolor','none','FaceAlpha',0.1)
        fill([100 100 125 125],[-5 5 5 -5],'k','Edgecolor','none','FaceAlpha',0.1)
        fill([175 175 200 200],[-5 5 5 -5],'k','Edgecolor','none','FaceAlpha',0.1)
        ylim([ -.25 .25]);
        %legend([h2,ev_labels,'Location','south','NumColumns',numel(ev_labels))
        xlabel(['Time from Event (downsampled to x/',num2str(sampleRate),' ms, 1st Stim | 2nd Stim | Choice Array )']);
        ylabel('Parameter Estimate')
        title(test_ev_labels(sp-1),'FontSize',16)
    end
end
return