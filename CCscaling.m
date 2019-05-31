clear all
endjf = readNPY('trop_SST_pavg_pmax_endjf_local_cores.npy');
lndjf = readNPY('trop_SST_pavg_pmax_lndjf_local_cores.npy');

fields = {'depth','pavg','pmax','CTT','SST'};
valsEN = [endjf(2,:) - endjf(1,:); endjf(3,:); endjf(4,:); endjf(1,:); ...
    endjf(2,:)];
valsLN = [lndjf(2,:) - lndjf(1,:); lndjf(3,:); lndjf(4,:); lndjf(1,:); ...
    lndjf(2,:)];
Systems = struct(fields{1},0,fields{2},0,fields{3},0,fields{4},0,fields{5},0);

for ii = 1:length(fields)
    if ii == 1
       Systems.(fields{ii}) = struct('EN',valsEN(ii,:),'LN',valsLN(ii,:));
    else if ii == 2 || ii == 3
       Systems.(fields{ii}) = struct('EN',0,'LN',0);
       choo = valsEN(ii,Systems.(fields{1}).EN < 65); 
       chee = valsEN(ii,Systems.(fields{1}).EN >= 65 & Systems.(fields{1}).EN < 85); 
       chii = valsEN(ii,Systems.(fields{1}).EN >= 85);
       i1 = find(choo ~= 0); i2 = find(chee ~= 0); i3 = find(chii ~= 0);
       Systems.(fields{ii}).EN = struct('shallow',choo(i1),'middle',chee(i2),...
           'deep',chii(i3));
       
       choo = valsLN(ii,Systems.(fields{1}).LN < 65); 
       chee = valsLN(ii,Systems.(fields{1}).LN >= 65 & Systems.(fields{1}).LN < 85);
       chii = valsLN(ii,Systems.(fields{1}).LN >= 85);
       i4 = find(choo ~= 0); i5 = find(chee ~= 0); i6 = find(chii ~= 0);
       Systems.(fields{ii}).LN = struct('shallow',choo(i4),'middle',chee(i5),...
           'deep',chii(i6));
        else
            Systems.(fields{ii}) = struct('EN',0,'LN',0);
            choo = valsEN(ii,Systems.(fields{1}).EN < 65); 
            choo = choo(i1);
            chee = valsEN(ii,Systems.(fields{1}).EN >= 65 & Systems.(fields{1}).EN < 85); 
            chee = chee(i2);
            chii = valsEN(ii,Systems.(fields{1}).EN >= 85);
            chii = chii(i3);
            Systems.(fields{ii}).EN = struct('shallow',choo,'middle',chee,...
                'deep',chii);
            
            choo = valsLN(ii,Systems.(fields{1}).LN < 65);
            choo = choo(i4);
            chee = valsLN(ii,Systems.(fields{1}).LN >= 65 & Systems.(fields{1}).LN < 85); 
            chee = chee(i5);
            chii = valsLN(ii,Systems.(fields{1}).LN >= 85);
            chii = chii(i6);
            Systems.(fields{ii}).LN = struct('shallow',choo,'middle',chee,...
                'deep',chii);
        end
    end
end
clear endjf lndjf fields valsEN valsLN chee choo chii 
clear i1 i2 i3 i4 i5 i6 ii jj

%%
clf
fig = figure(5);
ii = 25; kk = 250; ll = 100;
fs = 11.5;

subplot_tight(2,2,1,[0.09,0.09])
hold on
xx = linspace(285,308,50);
guy = logspace(-1,1.4,8);
for jj = 1:length(guy)
    plot(xx,guy(jj)*exp(-7000.*(1./xx-1/xx(1))),'color',[0.5 0.5 0.5],'linestyle','--')
end
scatter(Systems.SST.EN.shallow(1:ii:end),Systems.pmax.EN.shallow(1:ii:end),...
    20,'blue','filled','markeredgecolor','black')
scatter(Systems.SST.EN.middle(1:kk:end),Systems.pmax.EN.middle(1:kk:end),...
    20,'green','filled','markeredgecolor','black','markerfacealpha',0.5)
scatter(Systems.SST.EN.deep(1:ll:end),Systems.pmax.EN.deep(1:ll:end),...
    20,'red','filled','markeredgecolor','black','markerfacealpha',0.5)
xlim([286 304]); ylim([0.1 15])
set(gca,'fontsize',fs); set(gca,'yscale','log');
ylabel('Maximum precipitation [mm h^{-1}]','fontsize',fs-1)
text(0.05,0.95,'\bf{(a)}','fontsize',fs+2,'units','normalized')

subplot_tight(2,2,2,[0.09,0.09])
hold on
for jj = 1:length(guy)
    plot(xx,guy(jj)*exp(-7000.*(1./xx-1/xx(1))),'color',[0.5 0.5 0.5],'linestyle','--')
end
scatter(Systems.SST.EN.shallow(1:ii:end),Systems.pavg.EN.shallow(1:ii:end),...
    20,'blue','filled','markeredgecolor','black')
scatter(Systems.SST.EN.middle(1:kk:end),Systems.pavg.EN.middle(1:kk:end),...
    20,'green','filled','markeredgecolor','black','markerfacealpha',0.5)
scatter(Systems.SST.EN.deep(1:ll:end),Systems.pavg.EN.deep(1:ll:end),...
    20,'red','filled','markeredgecolor','black','markerfacealpha',0.5)
xlim([286 304]); ylim([0.1 5])
set(gca,'fontsize',fs); set(gca,'yscale','log');
ylabel('Mean precipitation [mm h^{-1}]','fontsize',fs-1)
text(0.05,0.95,'\bf{(b)}','fontsize',fs+2,'units','normalized')

subplot_tight(2,2,3,[0.09,0.09])
hold on
for jj = 1:length(guy)
    plot(xx,guy(jj)*exp(-7000.*(1./xx-1/xx(1))),'color',[0.5 0.5 0.5],'linestyle','--')
end
scatter(Systems.SST.LN.shallow(1:ii:end),Systems.pmax.LN.shallow(1:ii:end),...
    20,'blue','filled','markeredgecolor','black')
scatter(Systems.SST.LN.middle(1:kk:end),Systems.pmax.LN.middle(1:kk:end),...
    20,'green','filled','markeredgecolor','black','markerfacealpha',0.5)
scatter(Systems.SST.LN.deep(1:ll:end),Systems.pmax.LN.deep(1:ll:end),...
    20,'red','filled','markeredgecolor','black','markerfacealpha',0.5)
xlim([286 304]); ylim([0.1 15])
set(gca,'fontsize',fs); set(gca,'yscale','log');
ylabel('Maximum precipitation [mm h^{-1}]','fontsize',fs-1)
text(0.05,0.95,'\bf{(c)}','fontsize',fs+2,'units','normalized')

subplot_tight(2,2,4,[0.07,0.07])
hold on
for jj = 1:length(guy)
    plot(xx,guy(jj)*exp(-7000.*(1./xx-1/xx(1))),'color',[0.5 0.5 0.5],'linestyle','--')
end
scatter(Systems.SST.LN.shallow(1:ii:end),Systems.pavg.LN.shallow(1:ii:end),...
    20,'blue','filled','markeredgecolor','black')
scatter(Systems.SST.LN.middle(1:kk:end),Systems.pavg.LN.middle(1:kk:end),...
    20,'green','filled','markeredgecolor','black','markerfacealpha',0.5)
scatter(Systems.SST.LN.deep(1:ll:end),Systems.pavg.LN.deep(1:ll:end),...
    20,'red','filled','markeredgecolor','black','markerfacealpha',0.5)
xlim([286 304]); ylim([0.1 5])
set(gca,'fontsize',fs); set(gca,'yscale','log');
ylabel('Mean precipitation [mm h^{-1}]','fontsize',fs-1)
text(0.05,0.95,'\bf{(d)}','fontsize',fs+2,'units','normalized')
