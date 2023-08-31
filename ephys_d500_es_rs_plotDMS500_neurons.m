function ephys_d500_es_rs_plotDMS500_neurons(data,header,windowGap,filelabel)
global exptdata
numUnits=unique(header(:,4));

for un=1:length(numUnits)
    if exist([exptdata.figuredir500,filesep,'neuronPrintouts',filesep,exptdata.monkeyname,'_',exptdata.analysisName,'_',filelabel,'_Neuron',num2str(un),'.jpg'],'file') %&& exptdata.reprocess~=1
        disp(['Printout for ',filelabel,', ','Neuron',num2str(un),' exists.  Skipping...'])
        continue
    else
        %% First Stimulus = Stim1
        ylimit=round(max(nanmean(data(header(:,4)==un,21:end),1))/50)*50; if ylimit==0||isnan(ylimit), ylimit=50; end % round peak firing rate to nearest 50sp/s
        figure; set(gcf,'Units','Normalized');  set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial')
        tmp_SpikeData=data(header(:,4)==un,windowGap(1):end); %#ok<NASGU>
        ptr=find(header(:,4)==un,1); % find first row that features UN as unit number
        sessNum=header(ptr,2);
        sessUnitNum=header(ptr,3); clear ptr

        %% Passive Stim1
        subplot(2,3,1) % Stim1
        hold on % stimulus onset
        tmp_SpikeData=data(header(:,4)==un,windowGap(2):windowGap(3)-1);
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimA_expectRepeat_actualRepeat),:),1),'r-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimB_expectRepeat_actualAltern),:),1),'r:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimA_expectAltern_actualRepeat),:),1),'b:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimB_expectAltern_actualAltern),:),1),'b-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimA_expectEither_actualRepeat),:),1),'g-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimB_expectEither_actualAltern),:),1),'g:','LineWidth',2); catch; end
        grid on; ylim([0 ylimit])
        legend('expR/actR','expR/actA','expA/actR','expR/actA','expE/actR','expE/actA')
        title([filelabel,' Neuron No.',num2str(un),' (Session#',num2str(sessNum),' / Unit#',num2str(sessUnitNum),') - Responses to First Stimulus (A) : PASSIVE']);
        xlabel('Time From Stim 1 Onset (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
     
        subplot(2,3,2) % Stim2
        hold on % stimulus onset
        tmp_SpikeData=data(header(:,4)==un,windowGap(3):windowGap(4)-1);
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimA_expectRepeat_actualRepeat),:),1),'r-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimB_expectRepeat_actualAltern),:),1),'r:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimA_expectAltern_actualRepeat),:),1),'b:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimB_expectAltern_actualAltern),:),1),'b-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimA_expectEither_actualRepeat),:),1),'g-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimB_expectEither_actualAltern),:),1),'g:','LineWidth',2); catch; end
        grid on; ylim([0 ylimit])
        legend('expR/actR','expR/actA','expA/actR','expR/actA','expE/actR','expE/actA') 
        xlabel('Time From Stim 2 Onset (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
        
        subplot(2,3,3) % Stim2
        hold on % stimulus onset
        tmp_SpikeData=data(header(:,4)==un,windowGap(4):windowGap(5)-1);
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimA_expectRepeat_actualRepeat),:),1),'r-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimB_expectRepeat_actualAltern),:),1),'r:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimA_expectAltern_actualRepeat),:),1),'b:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimB_expectAltern_actualAltern),:),1),'b-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimA_expectEither_actualRepeat),:),1),'g-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.passiveCond,exptdata.dms500.stimAstimB_expectEither_actualAltern),:),1),'g:','LineWidth',2); catch; end
        grid on; ylim([0 ylimit])
        legend('expR/actR','expR/actA','expA/actR','expR/actA','expE/actR','expE/actA') 
        xlabel('Time From Choice Array Onset (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
     
        %% Active Stim1
        subplot(2,3,4) % Stim1
        hold on % stimulus onset
        tmp_SpikeData=data(header(:,4)==un,windowGap(2):windowGap(3)-1);
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimA_expectRepeat_actualRepeat),:),1),'r-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimB_expectRepeat_actualAltern),:),1),'r:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimA_expectAltern_actualRepeat),:),1),'b:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimB_expectAltern_actualAltern),:),1),'b-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimA_expectEither_actualRepeat),:),1),'g-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimB_expectEither_actualAltern),:),1),'g:','LineWidth',2); catch; end
        grid on; ylim([0 ylimit])
        legend('expR/actR','expR/actA','expA/actR','expR/actA','expE/actR','expE/actA')
        title([filelabel,' Neuron No.',num2str(un),' (Session#',num2str(sessNum),' / Unit#',num2str(sessUnitNum),') - Responses to First Stimulus (A) : ACTIVE']);
        xlabel('Time From Stim 1 Onset (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
     
        subplot(2,3,5) % Stim2
        hold on % stimulus onset
        tmp_SpikeData=data(header(:,4)==un,windowGap(3):windowGap(4)-1);
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimA_expectRepeat_actualRepeat),:),1),'r-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimB_expectRepeat_actualAltern),:),1),'r:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimA_expectAltern_actualRepeat),:),1),'b:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimB_expectAltern_actualAltern),:),1),'b-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimA_expectEither_actualRepeat),:),1),'g-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimB_expectEither_actualAltern),:),1),'g:','LineWidth',2); catch; end
        grid on; ylim([0 ylimit])
        legend('expR/actR','expR/actA','expA/actR','expR/actA','expE/actR','expE/actA') 
        xlabel('Time From Stim 2 Onset (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
        
        subplot(2,3,6) % Choice Array
        hold on % stimulus onset
        tmp_SpikeData=data(header(:,4)==un,windowGap(4):windowGap(5)-1);
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimA_expectRepeat_actualRepeat),:),1),'r-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimB_expectRepeat_actualAltern),:),1),'r:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimA_expectAltern_actualRepeat),:),1),'b:','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimB_expectAltern_actualAltern),:),1),'b-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimA_expectEither_actualRepeat),:),1),'g-','LineWidth',2); catch; end
        try plot(exptdata.xrange_psths,nanmean(tmp_SpikeData(intersect(exptdata.dms500.activeCond,exptdata.dms500.stimAstimB_expectEither_actualAltern),:),1),'g:','LineWidth',2); catch; end
        grid on; ylim([0 ylimit])
        legend('expR/actR','expR/actA','expA/actR','expR/actA','expE/actR','expE/actA') 
        xlabel('Time From Choice Array Onset (ms)'); ylabel('Avg SPDEN (sp/s)'); xlim([min(exptdata.xrange_psths) max(exptdata.xrange_psths)]);
     
        %savefig([exptdata.figuredir500,filesep,'neuronPrintouts',filesep,exptdata.monkeyname,'_',exptdata.analysisName,'_',filelabel,'_Neuron',num2str(un),'.fig'])
        jpgfigname=[exptdata.figuredir500,filesep,'neuronPrintouts',filesep,exptdata.monkeyname,'_',exptdata.analysisName,'_',filelabel,'_Neuron',num2str(un),'.jpg'];
        print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
        close all
    end
end
return
