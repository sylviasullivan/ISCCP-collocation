clear all
en = readNPY('trop_RHTWCAPE_7km_pmax_psum_endjf_pre.npy');
ln = readNPY('trop_RHTWCAPE_7km_pmax_psum_lndjf_pre.npy');

% save the values read in from the numpy files
eps = 0.01802/0.02897;
pmaxEN = en(6,:); pmaxLN = ln(6,:);
sstEN = en(3,:); sstLN = ln(3,:);
rhEN = en(9,:)./(eps.*SATVPLIQUID(en(8,:))./(40613 - SATVPLIQUID(en(8,:))));
rhLN = ln(9,:)./(eps.*SATVPLIQUID(ln(8,:))./(40613 - SATVPLIQUID(ln(8,:))));
qvEN = en(9,:); qvLN = ln(9,:);
wEN = en(10,:); wLN = ln(10,:);
capeEN = en(11,:); capeLN = ln(11,:);

% filter for instances where pmax is non-zero and there is underlying SST
ii = find(~isnan(pmaxEN) & pmaxEN ~= 0 & sstEN > 0 & rhEN < 1 & capeEN > 1);
pmaxEN = pmaxEN(ii);
t2mEN = en(1,ii);
sktEN = en(2,ii);
sstEN = en(3,ii);
dptEN = en(4,ii);
depthEN = sstEN - en(5,ii);
rhEN = rhEN(ii);
qvEN = qvEN(ii);
wEN = wEN(ii);
capeEN = capeEN(ii);

ii = find(~isnan(pmaxLN) & pmaxLN ~= 0 & sstLN > 0  & rhLN < 1 & capeLN > 1);
pmaxLN = pmaxLN(ii);
t2mLN = ln(1,ii);
sktLN = ln(2,ii);
sstLN = ln(3,ii);
dptLN = ln(4,ii);
depthLN = sstLN - ln(5,ii);
rhLN = rhLN(ii);
qvLN = qvLN(ii);
wLN = wLN(ii);
capeLN = capeLN(ii);
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
rh = struct('ENd1',0,'ENd2',0,'ENd3',0,'LNd1',0,'LNd2',0','LNd3',0,...
    'EN99d1',0,'EN99d2',0,'EN99d3',0,'LN99d1',0,'LN99d2',0,'LN99d3',0);
qv = struct('ENd1',0,'ENd2',0,'ENd3',0,'LNd1',0,'LNd2',0','LNd3',0,...
    'EN99d1',0,'EN99d2',0,'EN99d3',0,'LN99d1',0,'LN99d2',0,'LN99d3',0);
w = struct('ENd1',0,'ENd2',0,'ENd3',0,'LNd1',0,'LNd2',0','LNd3',0,...
    'EN99d1',0,'EN99d2',0,'EN99d3',0,'LN99d1',0,'LN99d2',0,'LN99d3',0);
cape = struct('ENd1',0,'ENd2',0,'ENd3',0,'LNd1',0,'LNd2',0','LNd3',0,...
    'EN99d1',0,'EN99d2',0,'EN99d3',0,'LN99d1',0,'LN99d2',0,'LN99d3',0);

d1 = find(depthEN < 65);
d2 = find(depthEN >= 65 & depthEN < 85);
d3 = find(depthEN >= 65);
pr.ENd1 = pmaxEN(d1); pr.ENd2 = pmaxEN(d2); pr.ENd3 = pmaxEN(d3);
t2m.ENd1 = t2mEN(d1); t2m.ENd2 = t2mEN(d2); t2m.ENd3 = t2mEN(d3); 
sst.ENd1 = sstEN(d1); sst.ENd2 = sstEN(d2); sst.ENd3 = sstEN(d3);
dpt.ENd1 = dptEN(d1); dpt.ENd2 = dptEN(d2); dpt.ENd3 = dptEN(d3);
skt.ENd1 = sktEN(d1); skt.ENd2 = sktEN(d2); skt.ENd3 = sktEN(d3);
rh.ENd1 = rhEN(d1); rh.ENd2 = rhEN(d2); rh.ENd3 = rhEN(d3);
qv.ENd1 = qvEN(d1); qv.ENd2 = qvEN(d2); qv.ENd3 = qvEN(d3);
w.ENd1 = wEN(d1); w.ENd2 = wEN(d2); w.ENd3 = wEN(d3);
cape.ENd1 = capeEN(d1); cape.ENd2 = capeEN(d2); cape.ENd3 = capeEN(d3);

d1 = find(depthLN < 65);
d2 = find(depthLN >= 65 & depthLN < 85);
d3 = find(depthLN >= 85);
pr.LNd1 = pmaxLN(d1); pr.LNd2 = pmaxLN(d2); pr.LNd3 = pmaxLN(d3);
t2m.LNd1 = t2mLN(d1); t2m.LNd2 = t2mLN(d2); t2m.LNd3 = t2mLN(d3); 
sst.LNd1 = sstLN(d1); sst.LNd2 = sstLN(d2); sst.LNd3 = sstLN(d3);
dpt.LNd1 = dptLN(d1); dpt.LNd2 = dptLN(d2); dpt.LNd3 = dptLN(d3);
skt.LNd1 = sktLN(d1); skt.LNd2 = sktLN(d2); skt.LNd3 = sktLN(d3);
rh.LNd1 = rhLN(d1); rh.LNd2 = rhLN(d2); rh.LNd3 = rhLN(d3);
qv.LNd1 = qvLN(d1); qv.LNd2 = qvLN(d2); qv.LNd3 = qvLN(d3);
w.LNd1 = wLN(d1); w.LNd2 = wLN(d2); w.LNd3 = wLN(d3);
cape.LNd1 = capeLN(d1); cape.LNd2 = capeLN(d2); cape.LNd3 = capeLN(d3);

clear d1 d2 d3 pmaxEN pmaxLN t2mEN t2mLN ii depthEN sstEN sstLN depthLN
clear en ln dptEN dptLN sktEN sktLN rhEN rhLN qvEN qvLN capeEN capeLN

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
num2 = 25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ll = 1:3
    % EL NINO
    temps = linspace(min(cape.(strcat('EN',depth{ll}))),max(cape.(strcat('EN',...
        depth{ll}))),num2);
%     temps = linspace(min(tstruct.(strcat('EN',depth{ll}))),max(tstruct.(strcat('EN',...
%         depth{ll}))),num2);
    t99 = []; pr99 = []; rh99 = []; qv99 = []; w99 = []; cape99 = [];
    for ii = 1:length(temps)-1
        jj = find(cape.(strcat('EN',depth{ll})) >= temps(ii) & cape.(strcat('EN',...
            depth{ll})) < temps(ii+1));
%         jj = find(tstruct.(strcat('EN',depth{ll})) >= temps(ii) & tstruct.(strcat('EN',...
%             depth{ll})) < temps(ii+1));
        subset = pr.(strcat('EN',depth{ll}))(jj);
        group = tstruct.(strcat('EN',depth{ll}))(jj);
        guy = rh.(strcat('EN',depth{ll}))(jj);
        gal = qv.(strcat('EN',depth{ll}))(jj);
        kid = w.(strcat('EN',depth{ll}))(jj);
        loser = cape.(strcat('EN',depth{ll}))(jj);
        
        [prpr,jj] = maxk(subset,num);
        pr99 = [pr99, prpr];
        t99 = [t99, group(jj)];
        rh99 = [rh99, guy(jj)];
        qv99 = [qv99, gal(jj)];
        w99 = [w99, kid(jj)];
        cape99 = [cape99, loser(jj)]; %maxk(loser,num)
    end
    pr.(strcat('EN99',depth{ll})) = pr99; 
    tstruct.(strcat('EN99',depth{ll})) = t99; 
    rh.(strcat('EN99',depth{ll})) = rh99;
    qv.(strcat('EN99',depth{ll})) = qv99;
    w.(strcat('EN99',depth{ll})) = w99;
    cape.(strcat('EN99',depth{ll})) = cape99;
    
    % LA NINA
    temps = linspace(min(cape.(strcat('LN',depth{ll}))),max(cape.(strcat('LN',...
        depth{ll}))),num2);
%     temps = linspace(min(tstruct.(strcat('LN',depth{ll}))),max(tstruct.(strcat('LN',...
%         depth{ll}))),num2);
    t99 = []; pr99 = []; rh99 = []; qv99 = []; w99 = []; cape99 = [];
    for ii = 1:length(temps)-1
        jj = find(cape.(strcat('LN',depth{ll})) >= temps(ii) & cape.(strcat('LN',...
            depth{ll})) < temps(ii+1));
%         jj = find(tstruct.(strcat('EN',depth{ll})) >= temps(ii) & tstruct.(strcat('EN',...
%             depth{ll})) < temps(ii+1));
        subset = pr.(strcat('LN',depth{ll}))(jj);
        group = tstruct.(strcat('LN',depth{ll}))(jj);
        guy = rh.(strcat('LN',depth{ll}))(jj);
        gal = qv.(strcat('LN',depth{ll}))(jj);
        kid = w.(strcat('LN',depth{ll}))(jj);
        loser = cape.(strcat('LN',depth{ll}))(jj);
        
        [prpr,jj] = maxk(subset,num);
        pr99 = [pr99, prpr];
        t99 = [t99, group(jj)];
        rh99 = [rh99, guy(jj)];
        qv99 = [qv99, gal(jj)];
        w99 = [w99, kid(jj)];
        cape99 = [cape99, loser(jj)]; %maxk(loser,num)
    end
    pr.(strcat('LN99',depth{ll})) = pr99; 
    tstruct.(strcat('LN99',depth{ll})) = t99; 
    rh.(strcat('LN99',depth{ll})) = rh99;
    qv.(strcat('LN99',depth{ll})) = qv99;
    w.(strcat('LN99',depth{ll})) = w99;
    cape.(strcat('LN99',depth{ll})) = cape99;
end
clear depth group subset ii jj ll pr99 t99 temps num num2 guy qv99 rh99 w99
clear gal kid wEN wLN cape99 loser

%%
clf
Opt.Interpreter = 'tex';

f = figure(10);
subplot_tight(1,3,1,[0.11 0.07])
hold on
scatter(log(cape.EN99d2),log(pr.EN99d2),rh.EN99d2.*30,[0 0.95 0],'filled',...
    'markeredgecolor','k')
P1 = polyfit(log(cape.EN99d2),log(pr.EN99d2),1);
plot(log(cape.EN99d2),P1(1)*log(cape.EN99d2)+P1(2),'Color',[0 0.9 0],'LineWidth',1.25);
r = corrcoef(P1(1)*log(cape.EN99d2)+P1(2),log(pr.EN99d2));
text(0.05,0.87,sprintf('r_{EN}^2 = %0.3f',r(1,2)),'color',[0 0.95 0],...
    'fontsize',11,'units','normalized')

P2 = polyfit(log(cape.LN99d2),log(pr.LN99d2),1);
plot(log(cape.LN99d2),P2(1)*log(cape.LN99d2)+P2(2),'Color',[0 0.35 0],'LineWidth',1.25);
scatter(log(cape.LN99d2),log(pr.LN99d2),rh.LN99d2.*30,[0 0.35 0],'filled',...
    'markeredgecolor','k')
r = corrcoef(P2(1)*log(cape.LN99d2)+P2(2),log(pr.LN99d2));
text(0.05,0.8,sprintf('r_{LN}^2 = %0.3f',r(1,2)),'color',[0 0.35 0],...
    'fontsize',11,'units','normalized')

% Clausius-Clapeyron relation
xx = linspace(100,5000,50);
plot(log(xx),-8+1.5.*(log(xx)),'color',[0.5 0.5 0.5],...
    'linestyle','--','handlevisibility','off')

text(0.05,0.95,'\bf{a}','fontsize',16,'units','normalized')
legend('off')
xlabel('CAPE [J kg^{-1}]')
ylabel('P_{max} [mm h^{-1}]')
ss = num2str([exp(-1); exp(0); exp(1); exp(2); exp(3); exp(4); exp(5)],3);
tt = num2str([exp(6); exp(6.5); exp(7); exp(7.5); exp(8); exp(8.5); exp(9)],4);
set(gca,'fontsize',13,'yticklabel',ss,'xticklabel',tt)
ylim([0 5.5]); xlim([log(500),log(5500)])

%%
subplot_tight(1,3,2,[0.1 0.05])
hold on
scatter(log(cape.EN99d1),log(pr.EN99d1),rh.EN99d1.*30,[0 0 1],'filled',...
    'markeredgecolor','k')
P3 = polyfit(log(cape.EN99d1),log(pr.EN99d1),1);
plot(log(cape.EN99d1),P3(1)*log(cape.EN99d1)+P3(2),'Color',[0 0 1],'LineWidth',1.25);
r = corrcoef(P3(1)*log(cape.EN99d1)+P3(2),log(pr.EN99d1));
text(0.05,0.87,sprintf('r_{EN}^2 = %0.3f',r(1,2)),'color',[0 0 1],...
    'fontsize',11,'units','normalized')

scatter(log(cape.LN99d1),log(pr.LN99d1),rh.LN99d1.*30,[0 0 0.5],'filled',...
    'd','markeredgecolor','k')
P4 = polyfit(log(cape.LN99d1),log(pr.LN99d1),1);
plot(log(cape.LN99d1),P4(1)*log(cape.LN99d1)+P4(2),'Color',[0 0 0.5],'LineWidth',1.25);
r = corrcoef(P4(1)*log(cape.LN99d1)+P4(2),log(pr.LN99d1));
text(0.05,0.8,sprintf('r_{LN}^2 = %0.3f',r(1,2)),'color',[0 0 0.5],...
    'fontsize',11,'units','normalized')

% Clausius-Clapeyron relation
xx = linspace(100,5000,50);
plot(log(xx),-2+0.5.*(log(xx)),'color',[0.5 0.5 0.5],...
    'linestyle','--','handlevisibility','off')
        
text(0.05,0.95,'\bf{b}','fontsize',16,'units','normalized')
legend('off')
xlabel('CAPE [J kg^{-1}]'); ylabel('');
ss = num2str([exp(-1); exp(0); exp(1); exp(2); exp(3); exp(4); exp(5)],3);
tt = num2str([exp(6); exp(6.5); exp(7); exp(7.5); exp(8); exp(8.5); exp(9)],4);
set(gca,'fontsize',13,'yticklabel',ss,'xticklabel',tt)
ylim([0 5.5]); xlim([log(500),log(5500)])

%%
subplot_tight(1,3,3,[0.1 0.05])
hold on
scatter(log(cape.EN99d3),log(pr.EN99d3),rh.EN99d3.*30,[1 0 0],'filled',...
    'markeredgecolor','k')
P5 = polyfit(log(cape.EN99d3),log(pr.EN99d3),1);
plot(log(cape.EN99d3),P5(1)*log(cape.EN99d3)+P5(2),'Color',[1 0 0],'LineWidth',1.25);
r = corrcoef(P5(1)*log(cape.EN99d3)+P5(2),log(pr.EN99d3));
text(0.05,0.87,sprintf('r_{EN}^2 = %0.3f',r(1,2)),'color',[1 0 0],...
    'fontsize',11,'units','normalized')

scatter(log(cape.LN99d3),log(pr.LN99d3),rh.LN99d3.*30,[0.5 0 0],'filled',...
    'markeredgecolor','k')
P6 = polyfit(log(cape.LN99d3),log(pr.LN99d3),1);
plot(log(cape.LN99d3),P6(1)*log(cape.LN99d3)+P6(2),'Color',[0.5 0 0],'LineWidth',1.25);
r = corrcoef(P6(1)*log(cape.LN99d3)+P6(2),log(pr.LN99d3));
text(0.05,0.8,sprintf('r_{LN}^2 = %0.3f',r(1,2)),'color',[0.5 0 0],...
    'fontsize',11,'units','normalized')

% Clausius-Clapeyron relation
xx = linspace(100,5000,50);
plot(log(xx),-2+0.5.*(log(xx)),'color',[0.5 0.5 0.5],...
    'linestyle','--','handlevisibility','off')

text(0.05,0.95,'\bf{c}','fontsize',16,'units','normalized')
legend('off') 
xlabel('CAPE [J kg^{-1}]'); ylabel('');
ss = num2str([exp(-1); exp(0); exp(1); exp(2); exp(3); exp(4); exp(5)],3);
tt = num2str([exp(6); exp(6.5); exp(7); exp(7.5); exp(8); exp(8.5); exp(9)],4);
set(gca,'fontsize',13,'yticklabel',ss,'xticklabel',tt)
ylim([0 5.5]); xlim([log(500),log(5500)])

orient(gcf,'landscape')
f.Units = 'inches';
set(gcf,'PaperPosition',[0.25,0.25,11,4.5])
% print(gcf,'-dpdf','-r250','CC-scalings-logCAPE')  %'-bestfit
