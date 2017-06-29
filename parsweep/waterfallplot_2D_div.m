appender = 'divided';
load([ pwd '\parsweep\saves\divsimpar'])
load([ pwd '\parsweep\saves\divparlists'])
cmap = linspecer(4); cmap = cmap(:,1:3);
m = 12;
set(gcf,'Position',[974.5000   98.5000  706 868],'Color',[1 1 1]);
cv = 0;
figure
for n = [1 6 8]
    cv = cv + 1;
    np = 2*n;
    load([pwd '\parsweep\saves\MNxstore_' num2str([m np]) appender])
    
    xstore = tmpst;
    sbi = 0;
    for i = [1 2 5]
        sbi = sbi+1;
        subplot(3,1,sbi)
        if i==1 || i == 3
            plot(tvec,log10(xstore(i,:)/100),'color',cmap(cv,:),'linewidth',2);
        else
            plot(tvec,log10(xstore(i,:)),'color',cmap(cv,:),'linewidth',2);
        end
        title(titname{i},'fontweight','bold','FontSize',12)
        hold on
        
    end
    
    kleg{cv} = ['K_{aa} = 10^{' num2str(log10(K_aa(np))) '} mol.dm^{-3}'];
end
subplot(3,1,1)
ylabel({'log Mean crystal'; 'volume (dm^3)'})
%     ylabel({'$log \hspace{1mm} \bar{V}_{crys}^{cyto}$'; '$(cm^3)$'},'interpreter','latex')
ylim([-21.5 -17.5])
%    set(gca,'Position',[ 974.5000   98.5000  857.5000  868.0000])

subplot(3,1,2)
ylabel({'log [crystal] in '; 'membrane (mol dm^{-3})'})
%     ylabel({'$log \hspace{1mm} [crys]^{mem}$'; '$(mol \hspace{1mm} dm^{-3})$'},'interpreter','latex')
ylim([-9.5 -5])
%     set(gca,'Position',[ 974.5000   98.5000  857.5000  868.0000])

subplot(3,1,3)
ylabel({'log SA of Cytosol'; '(cm^2)'})
%     ylabel({'$log \hspace{1mm} SA_{cyto} \hspace{1mm}$'; '$(cm^2)$'},'interpreter','latex')
xlabel('Time (days)')
%     xlabel('$Time \hspace{1mm}(days)$','interpreter','latex')
ylim([-9.4 -8.9])
%     set(gca,'Position',[   974.5000   98.5000  857.5000  868.0000])

%     h = legend(kleg,'FontSize',10);
%     v = get(h,'title');
%     set(v,'string','K_{aa}','fontweight','bold');
%     set(h,'Position',[  0.6306    0.8100    0.2613    0.0993])
% h.Position = [0.8142    0.7665    0.1589    0.1837];
% h.Title.String = 'K_{aa}';
% h.Title.FontWeight = 'bold';
shg

% saveallfiguresSFLAP([pwd '\parsweep\figures\protocell_divider'],'-tif','-r600'); close all