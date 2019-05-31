clear all
en = readNPY('trop_RHT_7km_pmax_psum_endjf.npy');
ln = readNPY('trop_RHT_7km_pmax_psum_lndjf.npy');

% save the values read in from the numpy files
eps = 0.01802/0.02897;
pmaxEN = en(6,:); 
pmaxLN = ln(6,:);
sstEN = en(3,:); 
sstLN = ln(3,:);
rhEN = en(9,:)./(eps.*SATVPLIQUID(en(8,:))./(40613 - SATVPLIQUID(en(8,:))));
rhLN = ln(9,:)./(eps.*SATVPLIQUID(ln(8,:))./(40613 - SATVPLIQUID(ln(8,:))));

% filter for instances where pmax is non-zero and there is underlying SST
ii = find(~isnan(pmaxEN) & pmaxEN ~= 0 & sstEN > 0);
pmaxEN = pmaxEN(ii);
t2mEN = en(1,ii);
sktEN = en(2,ii);
sstEN = en(3,ii);
dptEN = en(4,ii);
depthEN = sstEN - en(5,ii);

ii = find(~isnan(pmaxLN) & pmaxLN ~= 0 & sstLN > 0);
pmaxLN = pmaxLN(ii);
t2mLN = ln(1,ii);
sktLN = ln(2,ii);
sstLN = ln(3,ii);
dptLN = ln(4,ii);
depthLN = sstLN - ln(5,ii);
clear ii

%%
% filter for d1, d2, d3 systems
pr = struct('ENd1',0,'ENd2',0,'ENd3',0,'LNd1',0,'LNd2',0','LNd3',0,...
    'EN99d1',0,'EN99d2',0,'EN99d3',0,'LN99d1',0,'LN99d2',0,'LN99d3',0);
sst = struct('ENd1',0,'ENd2',0,'ENd3',0,'LNd1',0,'LNd2',0','LNd3',0,...
    'EN99d1',0,'EN99d2',0,'EN99d3',0,'LN99d1',0,'LN99d2',0,'LN99d3',0);
t2m = struct('ENd1',0,'ENd2',0,'ENd3',0,'LNd1',0,'LNd2',0','LNd3',0,...
    'EN99d1',0,'EN99d2',0,'EN99d3',0,'LN99d1',0,'LN99d2',0,'LN99d3',0);
dpt = struct('ENd1',0,'ENd2',0,'ENd3',0,'LNd1',0,'LNd2',0','LNd3',0,...
    'EN99d1',0,'EN99d2',0,'EN99d3',0,'LN99d1',0,'LN99d2',0,'LN99d3',0);
skt = struct('ENd1',0,'ENd2',0,'ENd3',0,'LNd1',0,'LNd2',0','LNd3',0,...
    'EN99d1',0,'EN99d2',0,'EN99d3',0,'LN99d1',0,'LN99d2',0,'LN99d3',0);

d1 = find(depthEN < 65);
d2 = find(depthEN >= 65 & depthEN < 85);
d3 = find(depthEN >= 65);
pr.ENd1 = pmaxEN(d1); pr.ENd2 = pmaxEN(d2); pr.ENd3 = pmaxEN(d3);
t2m.ENd1 = t2mEN(d1); t2m.ENd2 = t2mEN(d2); t2m.ENd3 = t2mEN(d3); 
sst.ENd1 = sstEN(d1); sst.ENd2 = sstEN(d2); sst.ENd3 = sstEN(d3);
dpt.ENd1 = dptEN(d1); dpt.ENd2 = dptEN(d2); dpt.ENd3 = dptEN(d3);
skt.ENd1 = sktEN(d1); skt.ENd2 = sktEN(d2); skt.ENd3 = sktEN(d3);

d1 = find(depthLN < 65);
d2 = find(depthLN >= 65 & depthLN < 85);
d3 = find(depthLN >= 85);
pr.LNd1 = pmaxLN(d1); pr.LNd2 = pmaxLN(d2); pr.LNd3 = pmaxLN(d3);
t2m.LNd1 = t2mLN(d1); t2m.LNd2 = t2mLN(d2); t2m.LNd3 = t2mLN(d3); 
sst.LNd1 = sstLN(d1); sst.LNd2 = sstLN(d2); sst.LNd3 = sstLN(d3);
dpt.LNd1 = dptLN(d1); dpt.LNd2 = dptLN(d2); dpt.LNd3 = dptLN(d3);
skt.LNd1 = sktLN(d1); skt.LNd2 = sktLN(d2); skt.LNd3 = sktLN(d3);

clear d1 d2 d3 pmaxEN pmaxLN t2mEN t2mLN ii depthEN sstEN sstLN depthLN
clear en ln dptEN dptLN sktEN sktLN

%%
% separate values into temperature bins and take the maximum 25 precip 
% values from each
depth = {'d1','d2','d3'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE BASED ON WHICH TEMPERATURE TO USE FOR SCALING
tstruct = sst;
% CHANGE BASED ON HOW MANY POINTS TO TAKE FROM EACH BIN
num = 5;
% CHANGE BASED ON HOW MANY TEMPERATURE BINS TO SET UP
num2 = 30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ll = 1:3
    % EL NINO
    temps = linspace(min(tstruct.(strcat('EN',depth{ll}))),max(tstruct.(strcat('EN',...
        depth{ll}))),num2);
    t99 = []; pr99 = []; 
    for ii = 1:length(temps)-1
        jj = find(tstruct.(strcat('EN',depth{ll})) >= temps(ii) & tstruct.(strcat('EN',...
            depth{ll})) < temps(ii+1));
        subset = pr.(strcat('EN',depth{ll}))(jj);
        group = tstruct.(strcat('EN',depth{ll}))(jj);
        pr99 = [pr99, maxk(subset,num)];
        t99 = [t99, maxk(group,num)];
    end
    pr.(strcat('EN99',depth{ll})) = pr99; 
    tstruct.(strcat('EN99',depth{ll})) = t99; 
    
    % LA NINA
    temps = linspace(min(tstruct.(strcat('LN',depth{ll}))),max(tstruct.(strcat('LN',...
        depth{ll}))),num2);
    t99 = []; pr99 = [];
    for ii = 1:length(temps)-1
        jj = find(tstruct.(strcat('LN',depth{ll})) >= temps(ii) & tstruct.(strcat('LN',...
            depth{ll})) < temps(ii+1));
        subset = pr.(strcat('LN',depth{ll}))(jj);
        group = tstruct.(strcat('LN',depth{ll}))(jj);
        pr99 = [pr99, maxk(subset,num)];
        t99 = [t99, maxk(group,num)];
    end
    pr.(strcat('LN99',depth{ll})) = pr99; 
    tstruct.(strcat('LN99',depth{ll})) = t99; 
end
clear depth group subset ii jj ll pr99 t99 temps num num2

%%
clf
Opt.Interpreter = 'tex';

figure(10)
subplot_tight(1,3,1,[0.1 0.05])
hold on
ll = find(tstruct.EN99d2 <= 303);
scatter(tstruct.EN99d2(ll),pr.EN99d2(ll),25,[0 0.95 0],'filled','markeredgecolor','k')
[fitEN,gofEN,~] = fit(tstruct.EN99d2(ll)',pr.EN99d2(ll)','exp1');
h = plot(fitEN);
set(h,'Color',[0 0.9 0],'LineWidth',1.25);

ll = find(tstruct.EN99d2 > 303);
scatter(tstruct.EN99d2(ll),pr.EN99d2(ll),20,[0 0.95 0],'*')
[fitEN,~,~] = fit(tstruct.EN99d2(ll)',pr.EN99d2(ll)','exp1');
h = plot(fitEN);
set(h,'Color',[0 0.9 0],'LineWidth',1.25,'LineStyle','--');

ll = find(tstruct.LN99d2 <= 300);
scatter(tstruct.LN99d2(ll),pr.LN99d2(ll),25,[0 0.35 0],'filled','d',...
    'markeredgecolor','k')
[fitLN,gofLN,~] = fit(tstruct.LN99d2(ll)',pr.LN99d2(ll)','exp1');
h = plot(fitLN);
set(h,'Color',[0 0.35 0],'LineWidth',1.25);

ll = find(tstruct.LN99d2 > 300);
scatter(tstruct.LN99d2(ll),pr.LN99d2(ll),30,[0 0.35 0],'x')
ll = find(tstruct.LN99d2 > 303);
[fitLN,~,~] = fit(tstruct.LN99d2(ll)',pr.LN99d2(ll)','exp1');
h = plot(fitLN);
set(h,'Color',[0 0.35 0],'LineWidth',1.25,'LineStyle','--');

% Clausius-Clapeyron relation
Rv = 461.52;        % water vapor gas constant [J kg-1 K-1]
Lv = 2458.3*1000;   % heat of enthalpy [J kg-1]
xx = linspace(280,310,50);
plot(xx,exp(-Lv/Rv*(1./xx-1/290)),'color',[0.5 0.5 0.5],...
        'linestyle','--','handlevisibility','off')

text(0.05,0.95,'\bf{a}','fontsize',16,'units','normalized')
text(0.45,0.3,sprintf('r_{EN}^2 = %0.3f',gofEN.rsquare),'color',[0 0.95 0],...
    'fontsize',11,'units','normalized')
text(0.45,0.23,sprintf('RMSE_{EN} = %0.3f',gofEN.rmse),'color',[0 0.95 0],...
    'fontsize',11,'units','normalized')
text(0.45,0.16,sprintf('r_{LN}^2 = %0.3f',gofLN.rsquare),'color',[0 0.35 0],...
    'fontsize',11,'units','normalized')
text(0.45,0.09,sprintf('RMSE_{LN} = %0.3f',gofLN.rmse),'color',[0 0.35 0],...
    'fontsize',11,'units','normalized')

legend('off')
xlabel('Underlying SST [K]')
ylabel('P_{max} [mm h^{-1}]')
set(gca,'yscale','log','fontsize',13)
ylim([0.1 150])

%%
subplot_tight(1,3,2,[0.1 0.05])
hold on
ll = find(tstruct.EN99d1 <= 295);
scatter(tstruct.EN99d1(ll),pr.EN99d1(ll),25,[0 0 1],'filled','markeredgecolor','k')
[fitEN,gofEN,~] = fit(tstruct.EN99d1(ll)',pr.EN99d1(ll)','exp1');
h = plot(fitEN);
set(h,'Color',[0 0 1],'LineWidth',1.25);

ll = find(tstruct.EN99d1 > 295);
scatter(tstruct.EN99d1(ll),pr.EN99d1(ll),25,[0 0 1],'*')
[fitEN,~,~] = fit(tstruct.EN99d1(ll)',pr.EN99d1(ll)','exp1');
h = plot(fitEN);
set(h,'Color',[0 0 1],'LineWidth',1.25,'LineStyle','--');

ll = find(tstruct.LN99d1 <= 295);
scatter(tstruct.LN99d1(ll),pr.LN99d1(ll),25,[0 0 0.5],'filled','d','markeredgecolor','k')
[fitLN,gofLN,~] = fit(tstruct.LN99d1(ll)',pr.LN99d1(ll)','exp1');
h = plot(fitLN);
set(h,'Color',[0 0 0.5],'LineWidth',1.25);

ll = find(tstruct.LN99d1 > 295);
scatter(tstruct.LN99d1(ll),pr.LN99d1(ll),25,[0 0 0.5],'x')
[fitLN,~,~] = fit(tstruct.LN99d1(ll)',pr.LN99d1(ll)','exp1');
h = plot(fitLN);
set(h,'Color',[0 0 0.5],'LineWidth',1.25,'LineStyle','--');

Rv = 461.52;   % water vapor gas constant [J kg-1 K-1]
Lv = 2458.3*1000;   % heat of enthalpy [J kg-1]
xx = linspace(280,310,50);
plot(xx,exp(-Lv/Rv*(1./xx-1/290)),'color',[0.5 0.5 0.5],...
        'linestyle','--','handlevisibility','off')
    
text(0.05,0.95,'\bf{b}','fontsize',16,'units','normalized')
text(0.45,0.3,sprintf('r_{EN}^2 = %0.3f',gofEN.rsquare),'color',[0 0 1],...
    'fontsize',11,'units','normalized')
text(0.45,0.23,sprintf('RMSE_{EN} = %0.3f',gofEN.rmse),'color',[0 0 1],...
    'fontsize',11,'units','normalized')
text(0.45,0.16,sprintf('r_{LN}^2 = %0.3f',gofLN.rsquare),'color',[0 0 0.5],...
    'fontsize',11,'units','normalized')   
text(0.45,0.09,sprintf('RMSE_{LN} = %0.3f',gofLN.rmse),'color',[0 0 0.5],...
    'fontsize',11,'units','normalized')   
legend('off')

xlabel('Underlying SST [K]'); ylabel('');
set(gca,'yscale','log','fontsize',13)
ylim([0.1 150])

%%
subplot_tight(1,3,3,[0.1 0.05])
hold on

ll = find(tstruct.EN99d3 <= 304);
scatter(tstruct.EN99d3(ll),pr.EN99d3(ll),25,[1 0 0],'filled','markeredgecolor','k')
[fitEN,gofEN,~] = fit(tstruct.EN99d3(ll)',pr.EN99d3(ll)','exp1');
h = plot(fitEN);
set(h,'Color',[1 0 0],'LineWidth',1.25);

ll = find(tstruct.EN99d3 > 304);
scatter(tstruct.EN99d3(ll),pr.EN99d3(ll),25,[1 0 0],'*')
[fitEN,~,~] = fit(tstruct.EN99d3(ll)',pr.EN99d3(ll)','exp1');
h = plot(fitEN);
set(h,'Color',[1 0 0],'LineWidth',1.25,'LineStyle','--');

ll = find(tstruct.LN99d3 <= 300);
scatter(tstruct.LN99d3(ll),pr.LN99d3(ll),25,[0.5 0 0],'filled','d','markeredgecolor','k')
[fitLN,gofLN,~] = fit(tstruct.LN99d3(ll)',pr.LN99d3(ll)','exp1');
h = plot(fitLN);
set(h,'Color',[0.5 0 0],'LineWidth',1.25);

ll = find(tstruct.LN99d3 > 300);
scatter(tstruct.LN99d3(ll),pr.LN99d3(ll),30,[0.5 0 0],'x')
[fitLN,~,~] = fit(tstruct.LN99d3(ll)',pr.LN99d3(ll)','exp1');
h = plot(fitLN);
set(h,'Color',[0.5 0 0],'LineWidth',1.25,'LineStyle','--');

Rv = 461.52;   % water vapor gas constant [J kg-1 K-1]
Lv = 2458.3*1000;   % heat of enthalpy [J kg-1]
xx = linspace(280,310,50);
plot(xx,exp(-Lv/Rv*(1./xx-1/290)),'color',[0.5 0.5 0.5],...
        'linestyle','--','handlevisibility','off')

text(0.05,0.95,'\bf{c}','fontsize',16,'units','normalized')
text(0.55,0.3,sprintf('r_{EN}^2 = %0.3f',gofEN.rsquare),'color',[1 0 0],...
    'fontsize',11,'units','normalized')
text(0.55,0.23,sprintf('RMSE_{EN} = %0.3f',gofEN.rmse),'color',[1 0 0],...
    'fontsize',11,'units','normalized')
text(0.55,0.16,sprintf('r_{LN}^2 = %0.3f',gofLN.rsquare),'color',[0.5 0 0],...
    'fontsize',11,'units','normalized')   
text(0.55,0.09,sprintf('RMSE_{LN} = %0.3f',gofLN.rmse),'color',[0.5 0 0],...
    'fontsize',11,'units','normalized')   
legend('off')
    
xlabel('Underlying SST [K]'); ylabel('');
set(gca,'yscale','log','fontsize',13)
ylim([0.1 150])
legend('off')

orient(gcf,'landscape')
set(gcf,'PaperPosition',[0.25,0.1,11,4.5])
% print(gcf,'-dpdf','-r250','CC-scalings-sst-redone2')  %'-bestfit
