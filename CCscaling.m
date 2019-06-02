clear all
en = readNPY('trop_RHTW_5km_pmax_psum_endjf.npy');
ln = readNPY('trop_RHTW_5km_pmax_psum_lndjf.npy');

% save the values read in from the numpy files
eps = 0.01802/0.02897;
pmaxEN = en(6,:); 
pmaxLN = ln(6,:);
sstEN = en(3,:); 
sstLN = ln(3,:);
rhEN = en(9,:)./(eps.*SATVPLIQUID(en(8,:))./(40613 - SATVPLIQUID(en(8,:))));
rhLN = ln(9,:)./(eps.*SATVPLIQUID(ln(8,:))./(40613 - SATVPLIQUID(ln(8,:))));
qvEN = en(9,:);
qvLN = ln(9,:);
% wEN = en(10,:);
% wLN = ln(10,:);

% filter for instances where pmax is non-zero and there is underlying SST
ii = find(~isnan(pmaxEN) & pmaxEN ~= 0 & sstEN > 0 & rhEN < 1);
pmaxEN = pmaxEN(ii);
t2mEN = en(1,ii);
sktEN = en(2,ii);
sstEN = en(3,ii);
dptEN = en(4,ii);
depthEN = sstEN - en(5,ii);
rhEN = rhEN(ii);
qvEN = qvEN(ii);
% wEN = wEN(ii);

ii = find(~isnan(pmaxLN) & pmaxLN ~= 0 & sstLN > 0  & rhLN < 1);
pmaxLN = pmaxLN(ii);
t2mLN = ln(1,ii);
sktLN = ln(2,ii);
sstLN = ln(3,ii);
dptLN = ln(4,ii);
depthLN = sstLN - ln(5,ii);
rhLN = rhLN(ii);
qvLN = qvLN(ii);
% wLN = wLN(ii);
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
% w = struct('ENd1',0,'ENd2',0,'ENd3',0,'LNd1',0,'LNd2',0','LNd3',0,...
%     'EN99d1',0,'EN99d2',0,'EN99d3',0,'LN99d1',0,'LN99d2',0,'LN99d3',0);

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
% w.ENd1 = wEN(d1); w.ENd2 = wEN(d2); w.ENd3 = wEN(d3);

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
% w.LNd1 = wLN(d1); w.LNd2 = wLN(d2); w.LNd3 = wLN(d3);

clear d1 d2 d3 pmaxEN pmaxLN t2mEN t2mLN ii depthEN sstEN sstLN depthLN
clear en ln dptEN dptLN sktEN sktLN rhEN rhLN qvEN qvLN

%%
% separate values into temperature bins and take the maximum 25 precip 
% values from each
depth = {'d1','d2','d3'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WHICH TEMPERATURE TO USE FOR SCALING
tstruct = sst;
% HOW MANY TEMPERATURE BINS
num = 15; num2 = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ll = 1:3
    % EL NINO
    temps = linspace(285,307,num);
    prQ = []; tQ = []; rhQ = []; qvQ = [];
    pr99 = []; t99 = []; rh99 = []; qv99 = [];
    for ii = 1:length(temps)-1
        jj = find(tstruct.(strcat('EN',depth{ll})) >= temps(ii) & tstruct.(strcat('EN',...
            depth{ll})) < temps(ii+1));
        
        % read in the relevant values
        precip = pr.(strcat('EN',depth{ll}))(jj);
        temperature = tstruct.(strcat('EN',depth{ll}))(jj);
        relhumidity = rh.(strcat('EN',depth{ll}))(jj);
        spechumidity = qv.(strcat('EN',depth{ll}))(jj);
        
        % calculate the quantiles of the subsets in this bin
        qpr = quantile(precip,[0.2 0.5 0.75 0.9]);
        qt = quantile(temperature,[0.2 0.5 0.75 0.9]);
        qrh = quantile(relhumidity,[0.2 0.5 0.75 0.9]);
        qqv = quantile(spechumidity,[0.2 0.5 0.75 0.9]);
        
        % append the quantiles to their respective lists
        prQ = [prQ; qpr]; 
        tQ = [tQ; qt]; 
        rhQ = [rhQ; qrh];
        qvQ = [qvQ; qqv];
        
        % extract the five largest precip values in each bin and their
        % correspondings conditions
        [prpr,i] = maxk(precip,num2);
        pr99 = [pr99, prpr];
        t99 = [t99, temperature(i)];
        rh99 = [rh99, relhumidity(i)];
        qv99 = [qv99, spechumidity(i)];   
    end
    pr.(strcat('ENq',depth{ll})) = prQ; 
    tstruct.(strcat('ENq',depth{ll})) = tQ; 
    rh.(strcat('ENq',depth{ll})) = rhQ;
    qv.(strcat('ENq',depth{ll})) = qvQ;
    
    pr.(strcat('EN99',depth{ll})) = pr99; 
    tstruct.(strcat('EN99',depth{ll})) = t99; 
    rh.(strcat('EN99',depth{ll})) = rh99;
    qv.(strcat('EN99',depth{ll})) = qv99;
    
    % LA NINA
    temps = linspace(285,307,num);
    prQ = []; tQ = []; rhQ = []; qvQ = [];
    pr99 = []; t99 = []; rh99 = []; qv99 = [];
    for ii = 1:length(temps)-1
        jj = find(tstruct.(strcat('LN',depth{ll})) >= temps(ii) & tstruct.(strcat('LN',...
            depth{ll})) < temps(ii+1));
        
        % read in the relevant values
        precip = pr.(strcat('LN',depth{ll}))(jj);
        temperature = tstruct.(strcat('LN',depth{ll}))(jj);
        relhumidity = rh.(strcat('LN',depth{ll}))(jj);
        spechumidity = qv.(strcat('LN',depth{ll}))(jj);
        
        % calculate the quantiles of the subsets in this bin
        qpr = quantile(precip,[0.2 0.5 0.75 0.9]);
        qt = quantile(temperature,[0.2 0.5 0.75 0.9]);
        qrh = quantile(relhumidity,[0.2 0.5 0.75 0.9]);
        qqv = quantile(spechumidity,[0.2 0.5 0.75 0.9]);
        
        % append the quantiles to their respective lists
        prQ = [prQ; qpr]; 
        tQ = [tQ; qt]; 
        rhQ = [rhQ; qrh];
        qvQ = [qvQ; qqv];
        
        % extract the five largest precip values in each bin and their
        % correspondings conditions
        [prpr,i] = maxk(precip,num2);
        pr99 = [pr99, prpr];
        t99 = [t99, temperature(i)];
        rh99 = [rh99, relhumidity(i)];
        qv99 = [qv99, spechumidity(i)];
    end
    pr.(strcat('LNq',depth{ll})) = prQ; 
    tstruct.(strcat('LNq',depth{ll})) = tQ; 
    rh.(strcat('LNq',depth{ll})) = rhQ;
    qv.(strcat('LNq',depth{ll})) = qvQ;
    
    pr.(strcat('LN99',depth{ll})) = pr99; 
    tstruct.(strcat('LN99',depth{ll})) = t99; 
    rh.(strcat('LN99',depth{ll})) = rh99;
    qv.(strcat('LN99',depth{ll})) = qv99;
end
clear depth ii jj pr99 t99 qv99 rh99 prQ tQ rhQ qvQ i ll num num2 precip
clear qpr qt qrh qqv  prpr relhumidity spechumidity temperature temps

%%
clf
Opt.Interpreter = 'tex';

f = figure(10);
subplot_tight(1,3,1,[0.11 0.07])
hold on
scatter(tstruct.EN99d2,log(pr.EN99d2),rh.EN99d2.*30,[0 0.95 0],'filled',...
    'markeredgecolor','k')
[dataout, ~, ~, ~] = lowess([tstruct.EN99d2; log(pr.EN99d2)]',0.5,1);
[peak, peakloc] = max(dataout(:,3));
plot(dataout(:,1),dataout(:,3),'Color',[0 0.9 0],'LineWidth',1.25);
plot([dataout(peakloc,1),dataout(peakloc,1)],[log(0.1),peak],'color',[0 0.95 0],...
    'linestyle','--')
r = corrcoef(dataout(:,2),dataout(:,3));
text(0.05,0.87,sprintf('r_{EN}^2 = %0.3f',r(1,2)),'color',[0 0.95 0],...
    'fontsize',11,'units','normalized')

scatter(tstruct.LN99d2,log(pr.LN99d2),rh.LN99d2.*30,[0 0.35 0],'filled',...
    'markeredgecolor','k')
[dataout, ~, ~, ~] = lowess([tstruct.LN99d2; log(pr.LN99d2)]',0.5,1);
[peak, peakloc] = max(dataout(:,3));
plot(dataout(:,1),dataout(:,3),'Color',[0 0.35 0],'LineWidth',1.25);
plot([dataout(peakloc,1),dataout(peakloc,1)],[log(0.1),peak],'color',[0 0.35 0],...
    'linestyle','--')
r = corrcoef(dataout(:,2),dataout(:,3));
text(0.05,0.8,sprintf('r_{LN}^2 = %0.3f',r(1,2)),'color',[0 0.35 0],...
    'fontsize',11,'units','normalized')

% Clausius-Clapeyron relation
% xx = linspace(260,310,50);
% for b = 260:10:320
%     plot(xx,log(8.73)+0.08.*(xx - b),'color',[0.5 0.5 0.5],...
%         'linestyle','--','handlevisibility','off')
% end

text(0.05,0.95,'\bf{a}','fontsize',16,'units','normalized')
legend('off')
xlabel('Underlying SST [K]')
ylabel('P_{max} [mm h^{-1}]')
ss = num2str([exp(-1); exp(0); exp(1); exp(2); exp(3); exp(4); exp(5)],3);
set(gca,'fontsize',13,'yticklabel',ss)
ylim([-1 5.5])

%     ,m/.%%
subplot_tight(1,3,2,[0.1 0.05])
hold on
scatter(tstruct.EN99d1,log(pr.EN99d1),rh.EN99d1.*30,[0 0 1],'filled',...
    'markeredgecolor','k')
[dataout, ~, ~, ~] = lowess([tstruct.EN99d1; log(pr.EN99d1)]',0.5,1);
[peak, peakloc] = max(dataout(:,3));
plot(dataout(:,1),dataout(:,3),'Color',[0 0 1],'LineWidth',1.25);
plot([dataout(peakloc,1),dataout(peakloc,1)],[log(0.1),peak],'Color',[0 0 1],...
    'linestyle','--')
r = corrcoef(dataout(:,2),dataout(:,3));
text(0.05,0.87,sprintf('r_{EN}^2 = %0.3f',r(1,2)),'color',[0 0 1],...
    'fontsize',11,'units','normalized')

scatter(tstruct.LN99d1,log(pr.LN99d1),rh.LN99d1.*30,[0 0 0.5],'filled',...
    'd','markeredgecolor','k')
[dataout, ~, ~, ~] = lowess([tstruct.LN99d1; log(pr.LN99d1)]',0.5,1);
[peak, peakloc] = max(dataout(:,3));
plot(dataout(:,1),dataout(:,3),'Color',[0 0 0.5],'LineWidth',1.25);
plot([dataout(peakloc,1),dataout(peakloc,1)],[log(0.1),peak],'Color',[0 0 0.5],...
    'linestyle','--')
r = corrcoef(dataout(:,2),dataout(:,3));
text(0.05,0.8,sprintf('r_{LN}^2 = %0.3f',r(1,2)),'color',[0 0 0.5],...
    'fontsize',11,'units','normalized')

% Clausius-Clapeyron relation
% xx = linspace(260,310,50);
% for b = 260:10:320
%     plot(xx,log(8.73)+0.07.*(xx - b),'color',[0.5 0.5 0.5],...
%         'linestyle','--','handlevisibility','off')
% end
        
text(0.05,0.95,'\bf{b}','fontsize',16,'units','normalized')
legend('off')
xlabel('Underlying SST [K]'); ylabel('');
ss = num2str([exp(-1); exp(0); exp(1); exp(2); exp(3); exp(4); exp(5)],3);
set(gca,'fontsize',13,'yticklabel',ss)
ylim([-1 5.5])

%%
subplot_tight(1,3,3,[0.1 0.05])
hold on
scatter(tstruct.EN99d3,log(pr.EN99d3),rh.EN99d3.*30,[1 0 0],'filled',...
    'markeredgecolor','k')
[dataout, ~, ~, ~] = lowess([tstruct.EN99d3; log(pr.EN99d3)]',0.5,1);
[peak, peakloc] = max(dataout(:,3));
plot(dataout(:,1),dataout(:,3),'Color',[1 0 0],'LineWidth',1.25);
plot([dataout(peakloc,1),dataout(peakloc,1)],[log(0.1),peak],'Color',[1 0 0],...
    'linestyle','--')
r = corrcoef(dataout(:,2),dataout(:,3));
text(0.05,0.87,sprintf('r_{EN}^2 = %0.3f',r(1,2)),'color',[1 0 0],...
    'fontsize',11,'units','normalized')

scatter(tstruct.LN99d3,log(pr.LN99d3),rh.LN99d3.*30,[0.5 0 0],'filled',...
    'markeredgecolor','k')
[dataout, ~, ~, ~] = lowess([tstruct.LN99d3; log(pr.LN99d3)]',0.5,1);
[peak, peakloc] = max(dataout(:,3));
plot(dataout(:,1),dataout(:,3),'Color',[0.5 0 0],'LineWidth',1.25);
plot([dataout(peakloc,1),dataout(peakloc,1)],[log(0.1),peak],'Color',[0.5 0 0],...
    'linestyle','--')
text(0.05,0.8,sprintf('r_{LN}^2 = %0.3f',r(1,2)),'color',[0.5 0 0],...
    'fontsize',11,'units','normalized')

% Clausius-Clapeyron relation
Rv = 461.52;   % water vapor gas constant [J kg-1 K-1]
% Lv = 2458.3*1000;   % heat of enthalpy [J kg-1]
% xx = linspace(280,310,50);
% plot(xx,-Lv/Rv*(1./xx-1/290),'color',[0.5 0.5 0.5],...
%         'linestyle','--','handlevisibility','off')
% xx = linspace(260,310,50);
% for b = 260:10:320
%     plot(xx,log(8.73)+0.07.*(xx - b),'color',[0.5 0.5 0.5],...
%         'linestyle','--','handlevisibility','off')
% end
text(0.05,0.95,'\bf{c}','fontsize',16,'units','normalized')
legend('off') 
xlabel('Underlying SST [K]'); ylabel('');
ss = num2str([exp(-1); exp(0); exp(1); exp(2); exp(3); exp(4); exp(5)],3);
set(gca,'fontsize',13,'yticklabel',ss)
ylim([-1 5.5])
legend('off')

orient(gcf,'landscape')
f.Units = 'inches';
set(gcf,'PaperPosition',[0.25,0.1,11,4.5])
% print(gcf,'-dpdf','-r250','CC-scalings-sst-redone4')  %'-bestfit
