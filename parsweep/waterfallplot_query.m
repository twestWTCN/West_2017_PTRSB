% close all
load([pwd '\parsweep\saves\simpar'])
load([pwd '\parsweep\saves\parlists'])


% load('C:\Users\Tim\Documents\Work\OriginsofLife\cleanform\FinalGoAt\parsweep\saves\simpar')
figure
set(gcf,'color','w');

% set(gcf,'Position',[  680.0000  280.5000  784.0000  697.5000]);
qpnts = [15 2; 15 15; 2 2; 2 15; 8 8]; 
% qpnts = [1 1; 1 4; 4 1; 2 2; 4 4]; 

cmap = linspecer(size(qpnts,1));
porder = [1 6 3 2 5];
for j = 1:size(qpnts,1)
    load([pwd '\parsweep\saves\MNxstore_' num2str([qpnts(j,1) qpnts(j,2)]) 'query'])
    xstore = tmpst;
    for i = 1:size(porder,2)
        subplot(size(porder,2),1,i)
        if i<7
             plot(tvec,log10(xstore(porder(i),:)),'color',cmap(j,:),'linewidth',2);
        else
            plot(tvec,xstore(porder(i),:),'color',cmap(j,:),'linewidth',2);
        end
        title(titname{porder(i)},'FontSize',12,'fontweight','bold')
        hold on
       xlim([0 100])
    end
  qpvals(j,:) = [log10(R_orgs_cat(qpnts(j,1))) log10(K_aa(qpnts(j,2)))];
  kleg{j} = ['Parameter set ' num2str(j)];
end
subplot(size(porder,2),1,1)
ylabel({'log Mean crystal'; 'volume (dm^3)'})
% ylabel({'$log \bar{V}_{crys}^{cyto}$'; '$(dm^3)$'},'interpreter','latex')
ylim([-19.5 -15.5])
subplot(size(porder,2),1,2)
ylabel({'log [crystal] in'; 'cytosol (mol dm^{-3})'})
% ylabel({'$log [crys]^{cyto}$'; '$(mol \hspace{1mm} dm^{-3})$'},'interpreter','latex')
ylim([-9 -4])
subplot(size(porder,2),1,3)
ylabel({'log [amino acids] in'; 'cytosol (mol dm^{-3})'})
% ylabel({'$log [aa]^{cyto}$'; '$(mol \hspace{1mm} dm^{-3})$'},'interpreter','latex')
ylim([-6.5 -0.5])
subplot(size(porder,2),1,4)
ylabel({'log [crystal] in '; 'membrane (mol dm^{-3})'})
% ylabel({'$log [crys]^{mem}$'; '$(mol \hspace{1mm} dm^{-3})$'},'interpreter','latex')
ylim([-9.5 -4.5])
subplot(size(porder,2),1,5)
ylabel({'log SA';'(cm^2)'})
% ylabel({'$log SA^{cyto}$';'$(cm^2)$'},'interpreter','latex')
ylim([-10 -7])
xlabel('Time (days)')
% xlabel('$Time (days)$','interpreter','latex')
% h = legend(kleg)
% v = get(h,'title');
% set(v,'string','Rcat (M)');
% set(h,'Position',[0.8206    0.0282    0.1594    0.1671])
shg
set(gcf,'Position',[543   101   548   858])
saveallfiguresSFLAP([pwd '\parsweep\figures\protocell_reworked_query_series'],'-tif','-r600'); close all