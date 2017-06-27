clear all; %close all
load([pwd '\parsweep\saves\parlists'])
load([pwd '\parsweep\saves\simpar'])

% load(['C:\Users\Tim\Documents\Work\OriginsofLife\cleanform\FinalGoAt\parsweep\saves\parlists'])

for n =1:Nmn
    for m = 1:Nmn
        load([pwd '\parsweep\saves\MNxstore_' num2str([n m]) 'query'])
        %         load(['C:\Users\Tim\Documents\Work\OriginsofLife\cleanform\FinalGoAt\parsweep\saves\MNxstore_' num2str([n m]) 'query'])
        SAmax(n,m) = max(tmpst(5,:));
        if  max(tmpst(5,:)) >0.9e-9
            L = min(find(tmpst(5,:)> 0.99*max(tmpst(5,:))));
            SAtime(n,m) =  SAmax(n,m)/tvec(L);
        else
            SAtime(n,m) = NaN;
        end
    end
end
figure
colormap jet
set(gcf,'color','w');
colormap jet

% [Xqflat Yqflat Vq] = interpmap(log10(kd_cys),log10(Rcat(1:15)),log10(SAmax(1:15,:)),0.05,'linear')
% imagesc(Xqflat,Yqflat,Vq)
% SAmax(SAmax<1e-7) = NaN;
imagesc2(log10(K_aa),log10(R_orgs_cat(1:Nmn)),log10(SAmax(1:Nmn,:)))
dat = log10(SAmax(1:Nmn,:)); [nr,nc] = size(dat);
%  pcolor([dat nan(nr,1); nan(1,nc+1)]);
% imagesc(1:16,1:16,log10(SAmax(1:Nmn,:)))

set(gca,'XAxisLocation','bottom')
set(gca,'YDir','normal')
xlabel({'log A.A. Binding Constant(mol dm^{-3})'},'FontSize',12)
ylabel('log Catalytic Rate (mol cm^{-2} s^{-1})','FontSize',12)

% xlabel({'$log \hspace{1mm} K_{aa} - (mol \hspace{1mm} dm^{-3})$'},'FontSize',12,'interpreter','latex')
% ylabel('$log \hspace{1mm} R_{cat}^{org} - (mol \hspace{1mm} cm^{-2} \hspace{1mm} s^{-1})$','FontSize',12,'interpreter','latex')
axis square
title('Membrane surface area at equilibrium','FontSize',12,'fontweight','bold')
c = colorbar;
ylabel(c,'log Membrane surface area(cm^2)','FontSize',12)
% ylabel(c,'$log \hspace{1mm} SA^{cyto} \hspace{1mm} (cm^2)$','FontSize',12,'interpreter','latex')

figure
set(gcf,'color','w');
colormap jet
% colormap(flipud(colormap))
SAtime_sc = SAtime;
% SAtime_sc(SAmax<1e-7) = NaN;
% [Xqflat Yqflat Vq] = interpmap(log10(kd_cys),log10(Rcat(1:15)),log10(SAtime_sc(1:15,:)),0.05,'linear')
% imagesc(Xqflat,Yqflat,Vq)
imagesc2(log10(K_aa),log10(R_orgs_cat(1:Nmn)),log10(SAtime_sc(1:Nmn,:)))
% imagesc(1:16,1:16,log10(SAtime_sc(1:Nmn,:)))
set(gca,'XAxisLocation','bottom')
set(gca,'YDir','normal')
xlabel({'log A.A. Binding Constant(mol dm^{-3})'},'FontSize',12)
ylabel('log Catalytic Rate (mol cm^{-2} s^{-1})','FontSize',12)
% xlabel({'$log \hspace{1mm} K_{aa} - (mol \hspace{1mm} dm^{-3})$'},'FontSize',12,'interpreter','latex')
% ylabel('$log \hspace{1mm} R_{cat}^{org} - (mol \hspace{1mm} cm^{-2} \hspace{1mm} s^{-1})$','FontSize',12,'interpreter','latex')
axis square
title('Rate of membrane surface area increase','FontSize',12,'fontweight','bold')

c = colorbar;
ylabel(c,{'log Rate of growth(cm^2 day^-1)'},'FontSize',12)
% ylabel(c,'$log \hspace{1mm} (cm^2 \hspace{1mm} day^{-1})$','FontSize',12,'interpreter','latex')
saveallfiguresSFLAP([pwd '\parsweep\figures\protocell_parspace'],'-tif','-r600'); close all