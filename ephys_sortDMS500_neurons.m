function [output,newheader]=lsn_ephys_sortDMS500_neurons(data,header,wG,filelabel,generateFig)
% function to classify neurons as visual, etc. and to perform simple
% univariate analysis

global exptdata
numUnits=unique(header(:,4));
newheader=header;
for un=1:length(numUnits)
    % is this neuron visually response?
    output(un).trialBaseline=nanmean(data(header(:,4)==un,wG(2)+abs(min(exptdata.xrange_psths))+exptdata.baselineResponseWindow(1):wG(2)+abs(min(exptdata.xrange_psths))+exptdata.baselineResponseWindow(2)),2); %#ok<*AGROW>
    output(un).trial1stStim=nanmean(data(header(:,4)==un,wG(2)+abs(min(exptdata.xrange_psths))+exptdata.sensoryResponseWindow(1):wG(2)+abs(min(exptdata.xrange_psths))+exptdata.sensoryResponseWindow(2)),2);
    output(un).trial2ndStim=nanmean(data(header(:,4)==un,wG(3)+abs(min(exptdata.xrange_psths))+exptdata.sensoryResponseWindow(1):wG(3)+abs(min(exptdata.xrange_psths))+exptdata.sensoryResponseWindow(2)),2);
        
    newheader(newheader(:,4)==un,20)=output(un).trialBaseline;
    newheader(newheader(:,4)==un,16)=output(un).trial1stStim;
    newheader(newheader(:,4)==un,17)=output(un).trial2ndStim;
    
    baselineActivity=nanmean(data(header(:,4)==un,wG(2)+abs(min(exptdata.xrange_psths))+exptdata.baselineResponseWindow(1):wG(2)+abs(min(exptdata.xrange_psths))+exptdata.baselineResponseWindow(2)),2);
    sensoryResponse_1stPresStim1=nanmean(data(header(:,4)==un & header(:,8)==1,wG(2)+abs(min(exptdata.xrange_psths))+exptdata.sensoryResponseWindow(1):wG(2)+abs(min(exptdata.xrange_psths))+exptdata.sensoryResponseWindow(2)),2);
    sensoryResponse_1stPresStim2=nanmean(data(header(:,4)==un & header(:,8)==2,wG(2)+abs(min(exptdata.xrange_psths))+exptdata.sensoryResponseWindow(1):wG(2)+abs(min(exptdata.xrange_psths))+exptdata.sensoryResponseWindow(2)),2);
    
    output(un).summaryStats(1)=nanmean(output(un).trialBaseline);
    output(un).summaryStats(2)=nanstd(output(un).trialBaseline);
    output(un).summaryStats(3)=nanmean(sensoryResponse_1stPresStim1);
    output(un).summaryStats(4)=nanstd(sensoryResponse_1stPresStim1);
    output(un).summaryStats(5)=nanmean(sensoryResponse_1stPresStim2);
    output(un).summaryStats(6)=nanstd(sensoryResponse_1stPresStim2);
    
    [output(un).ttestBaselineV1stStim1(1),output(un).ttestBaselineV1stStim1(2)]=ttest2(baselineActivity,sensoryResponse_1stPresStim1,'vartype','equal','alpha',0.05,'tail','both');
    [output(un).ttestBaselineV1stStim2(1),output(un).ttestBaselineV1stStim2(2)]=ttest2(baselineActivity,sensoryResponse_1stPresStim2,'vartype','equal','alpha',0.05,'tail','both');
    [output(un).ttest1stStim1vStim2(1),output(un).ttest1stStim1vStim2(2)]=ttest2(sensoryResponse_1stPresStim1,sensoryResponse_1stPresStim2,'vartype','equal','alpha',0.05,'tail','both');
        
    newheader(newheader(:,4)==un,21)=output(un).ttestBaselineV1stStim1(2);
    newheader(newheader(:,4)==un,22)=output(un).ttestBaselineV1stStim2(2);
    newheader(newheader(:,4)==un,23)=output(un).ttest1stStim1vStim2(2);
    
    
    %% Plot Figure
    if generateFig==1
        % First Stimulus = Stim1
        figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial')
        
        nml=max([output(un).summaryStats(3) output(un).summaryStats(5)]); % normalise everything to maximum 1stPres response
        %% StimA
        subplot(2,4,1) % FPonset
        hold on % stimulus onset
        tmp_SpikeData=data(header(:,4)==un,wG(1):wG(2)-1);
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.dms500.stimAstimA_expectRepeat_actualRepeat,:),1)/nml,'r-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.dms500.stimAstimB_expectRepeat_actualAltern,:),1)/nml,'r:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.dms500.stimAstimA_expectAltern_actualRepeat,:),1)/nml,'b:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.dms500.stimAstimB_expectAltern_actualAltern,:),1)/nml,'b-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.dms500.stimAstimA_expectEither_actualRepeat,:),1)/nml,'g-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.dms500.stimAstimB_expectEither_actualAltern,:),1)/nml,'g:','LineWidth',2); catch; end
        grid on; ylim([0 5])
        legend('expR/actR','expR/actA','expA/actR','expR/actA','expE/actR','expE/actA')
        title([filelabel,' Neuron # ',num2str(un),' - Responses to First Stimulus = 1']);
        xlabel('Time From FP Onset (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
        %[t p]=masst(squeeze(tmpx(:,1,:))-squeeze(tmpx(:,2,:)));
        %patchyvect(p>0.999,[6.5 7],{'w','r'});
        
        subplot(2,4,2) % First Stim
        hold on % stimulus onset
        tmp_SpikeData=data(header(:,4)==un,wG(2):wG(3)-1);
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimA_expectRepeat_actualRepeat,:),1)/nml,'r-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimB_expectRepeat_actualAltern,:),1)/nml,'r:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimA_expectAltern_actualRepeat,:),1)/nml,'b:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimB_expectAltern_actualAltern,:),1)/nml,'b-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimA_expectEither_actualRepeat,:),1)/nml,'g-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimB_expectEither_actualAltern,:),1)/nml,'g:','LineWidth',2); catch; end
        grid on;  ylim([0 5])
        legend('expR/actR','expR/actA','expA/actR','expR/actA','expE/actR','expE/actA')
        xlabel('Time From First Stimulus Presentation (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
        %[t p]=masst(squeeze(tmpx(:,1,:))-squeeze(tmpx(:,2,:)));
        %patchyvect(p>0.999,[6.5 7],{'w','r'});
        
        subplot(2,4,3) % Second Stim
        hold on % stimulus onset
        tmp_SpikeData=data(header(:,4)==un,wG(3):wG(4)-1);
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimA_expectRepeat_actualRepeat,:),1)/nml,'r-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimB_expectRepeat_actualAltern,:),1)/nml,'r:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimA_expectAltern_actualRepeat,:),1)/nml,'b:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimB_expectAltern_actualAltern,:),1)/nml,'b-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimA_expectEither_actualRepeat,:),1)/nml,'g-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimB_expectEither_actualAltern,:),1)/nml,'g:','LineWidth',2); catch; end
        grid on;  ylim([0 5])
        legend('expR/actR','expR/actA','expA/actR','expR/actA','expE/actR','expE/actA')
        xlabel('Time From Second Stimulus Presentation (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
        %[t p]=masst(squeeze(tmpx(:,1,:))-squeeze(tmpx(:,2,:)));
        %patchyvect(p>0.999,[6.5 7],{'w','r'});
        
        subplot(2,4,4) % Reward
        hold on % stimulus onset
        tmp_SpikeData=data(header(:,4)==un,wG(4):wG(5)-1);
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimA_expectRepeat_actualRepeat,:),1)/nml,'r-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimB_expectRepeat_actualAltern,:),1)/nml,'r:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimA_expectAltern_actualRepeat,:),1)/nml,'b:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimB_expectAltern_actualAltern,:),1)/nml,'b-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimA_expectEither_actualRepeat,:),1)/nml,'g-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimAstimB_expectEither_actualAltern,:),1)/nml,'g:','LineWidth',2); catch; end
        grid on;  ylim([0 5])
        legend('expR/actR','expR/actA','expA/actR','expR/actA','expE/actR','expE/actA')
        xlabel('Time From Reward Delivery (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
        %[t p]=masst(squeeze(tmpx(:,1,:))-squeeze(tmpx(:,2,:)));
        %patchyvect(p>0.999,[6.5 7],{'w','r'});
        
        %% StimB
        subplot(2,4,5) % FPonset
        hold on % stimulus onset
        tmp_SpikeData=data(header(:,4)==un,wG(1):wG(2)-1);
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.stimBstimB_expectRepeat_actualRepeat,:),1)/nml,'r-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.stimBstimA_expectRepeat_actualAltern,:),1)/nml,'r:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.stimBstimB_expectAltern_actualRepeat,:),1)/nml,'b:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.stimBstimA_expectAltern_actualAltern,:),1)/nml,'b-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.stimBstimB_expectEither_actualRepeat,:),1)/nml,'g-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.stimBstimA_expectEither_actualAltern,:),1)/nml,'g:','LineWidth',2); catch; end
        grid on;  ylim([0 5])
        legend('expR/actR','expR/actA','expA/actR','expR/actA','expE/actR','expE/actA')
        title([filelabel,' Neuron # ',num2str(un),' - Responses to First Stimulus = 2']);
        xlabel('Time From FP Onset (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
        %[t p]=masst(squeeze(tmpx(:,1,:))-squeeze(tmpx(:,2,:)));
        %patchyvect(p>0.999,[6.5 7],{'w','r'});
        
        subplot(2,4,6) % First Stim
        hold on % stimulus onset
        tmp_SpikeData=data(header(:,4)==un,wG(2):wG(3)-1);
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimB_expectRepeat_actualRepeat,:),1)/nml,'r-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimA_expectRepeat_actualAltern,:),1)/nml,'r:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimB_expectAltern_actualRepeat,:),1)/nml,'b:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimA_expectAltern_actualAltern,:),1)/nml,'b-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimB_expectEither_actualRepeat,:),1)/nml,'g-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimA_expectEither_actualAltern,:),1)/nml,'g:','LineWidth',2); catch; end
        grid on;  ylim([0 5])
        legend('expR/actR','expR/actA','expA/actR','expR/actA','expE/actR','expE/actA')
        xlabel('Time From First Stimulus Presentation (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
        %[t p]=masst(squeeze(tmpx(:,1,:))-squeeze(tmpx(:,2,:)));
        %patchyvect(p>0.999,[6.5 7],{'w','r'});
        
        subplot(2,4,7) % Second Stim
        hold on % stimulus onset
        tmp_SpikeData=data(header(:,4)==un,wG(3):wG(4)-1);
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimB_expectRepeat_actualRepeat,:),1)/nml,'r-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimA_expectRepeat_actualAltern,:),1)/nml,'r:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimB_expectAltern_actualRepeat,:),1)/nml,'b:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimA_expectAltern_actualAltern,:),1)/nml,'b-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimB_expectEither_actualRepeat,:),1)/nml,'g-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimA_expectEither_actualAltern,:),1)/nml,'g:','LineWidth',2); catch; end
        grid on;  ylim([0 5])
        legend('expR/actR','expR/actA','expA/actR','expR/actA','expE/actR','expE/actA')
        xlabel('Time From Second Stimulus Presentation (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
        %[t p]=masst(squeeze(tmpx(:,1,:))-squeeze(tmpx(:,2,:)));
        %patchyvect(p>0.999,[6.5 7],{'w','r'});
        
        subplot(2,4,8) % Reward
        hold on % stimulus onset
        tmp_SpikeData=data(header(:,4)==un,wG(4):wG(5)-1);
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimB_expectRepeat_actualRepeat,:),1)/nml,'r-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimA_expectRepeat_actualAltern,:),1)/nml,'r:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimB_expectAltern_actualRepeat,:),1)/nml,'b:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimA_expectAltern_actualAltern,:),1)/nml,'b-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimB_expectEither_actualRepeat,:),1)/nml,'g-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(exptdata.rsvp300.stimBstimA_expectEither_actualAltern,:),1)/nml,'g:','LineWidth',2); catch; end
        grid on;  ylim([0 5])
        legend('expR/actR','expR/actA','expA/actR','expR/actA','expE/actR','expE/actA')
        xlabel('Time From Reward Delivery (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
        %[t p]=masst(squeeze(tmpx(:,1,:))-squeeze(tmpx(:,2,:)));
        %patchyvect(p>0.999,[6.5 7],{'w','r'});
        savefig([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_',filelabel,'Neuron',num2str(un),'.fig'])
        export_fig([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_',filelabel,'Neuron',num2str(un),'.jpg'],'-jpg')
        print(gcf,[exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_',filelabel,'Neuron',num2str(un),'.jpg'],'-djpeg') % generates an JPEG file of the figure
       % export_fig([exptdata.projectdir,exptdata.monkeyname,'_',exptdata.analysisName,'_',filelabel,'Neuron',num2str(un),'.eps'], '-eps', '-rgb', '-transparent');
        close all
    end
end



return