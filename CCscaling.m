clear all
endjf = readNPY('trop_SST_pcore_endjf_local_cores.npy'); 
lndjf = readNPY('trop_SST_pcore_lndjf_local_cores.npy');

fields = {'depth','pmax','SST'};
valsEN = [endjf(2,:) - endjf(1,:); endjf(3,:); endjf(2,:)];
valsLN = [lndjf(2,:) - lndjf(1,:); lndjf(3,:); lndjf(2,:)];
Systems = struct(fields{1},0,fields{2},0,fields{3},0);

for ii = 1:length(fields)
    if ii == 1
       Systems.(fields{ii}) = struct('EN',valsEN(ii,:),'LN',valsLN(ii,:));
    else if ii == 2
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
       i4 = find(choo ~= 0 & ~isnan(choo)); i5 = find(chee ~= 0 & ~isnan(chee)); 
       i6 = find(chii ~= 0 & ~isnan(chii));
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
set(fig,'PaperOrientation','landscape')
ii = 10; kk = 100; ll = 50;
fs = 13.5;

subplot_tight(1,2,1,[0.12,0.09])
hold on
xx = linspace(290,305,50);
guy = logspace(0.2,1.1,6);
Rv = 461.52;   % water vapor gas constant [J kg-1 K-1]
Lv = 2458.3*1000;   % heat of enthalpy [J kg-1]
for jj = 1:length(guy)
    plot(xx,guy(jj)*exp(-Lv/Rv*(1./xx-1/290)),'color',[0.5 0.5 0.5],...
        'linestyle','--','handlevisibility','off')
end
pmaxEN1 = Systems.pmax.EN.shallow > prctile(Systems.pmax.EN.shallow,98);
pmaxEN2 = Systems.pmax.EN.middle > prctile(Systems.pmax.EN.middle,99.9);
pmaxEN3 = Systems.pmax.EN.deep > prctile(Systems.pmax.EN.deep,99.5);
scatter(Systems.SST.EN.shallow(pmaxEN1),Systems.pmax.EN.shallow(pmaxEN1),...
    25,'blue','filled','markeredgecolor','black','displayname','least deep');
scatter(Systems.SST.EN.middle(pmaxEN2),Systems.pmax.EN.middle(pmaxEN2),...
    25,'green','filled','markeredgecolor','black','markerfacealpha',0.5,...
    'displayname','intermediate');
scatter(Systems.SST.EN.deep(pmaxEN3),Systems.pmax.EN.deep(pmaxEN3),...
    25,'red','filled','markeredgecolor','black','markerfacealpha',0.5,...
    'displayname','deepest');
% xlim([290 305]); ylim([5 20])
set(gca,'fontsize',fs); set(gca,'yscale','log');
ylabel('Maximum precipitation [mm h^{-1}]','fontsize',fs-1)
text(0.05,0.95,'\bf{(a)} - El Niño','fontsize',fs+2,'units','normalized')
xlabel('Sea surface temperature [K]','fontsize',fs-1)
ll = legend;
legend('boxoff')
set(ll,'location','southeast')

subplot_tight(1,2,2,[0.12,0.09])
hold on
for jj = 1:length(guy)
    plot(xx,guy(jj)*exp(-Lv/Rv*(1./xx-1/xx(1))),'color',[0.5 0.5 0.5],'linestyle','--')
end
pmaxLN1 = Systems.pmax.LN.shallow > prctile(Systems.pmax.LN.shallow,98);
pmaxLN2 = Systems.pmax.LN.middle > prctile(Systems.pmax.LN.middle,99.9);
pmaxLN3 = Systems.pmax.LN.deep > prctile(Systems.pmax.LN.deep,99.5);
scatter(Systems.SST.LN.shallow(pmaxLN1),Systems.pmax.LN.shallow(pmaxLN1),...
    25,'blue','filled','markeredgecolor','black')
scatter(Systems.SST.LN.middle(pmaxLN2),Systems.pmax.LN.middle(pmaxLN2),...
    25,'green','filled','markeredgecolor','black','markerfacealpha',0.5)
scatter(Systems.SST.LN.deep(pmaxLN3),Systems.pmax.LN.deep(pmaxLN3),...
    25,'red','filled','markeredgecolor','black','markerfacealpha',0.5)
% xlim([290 305]); ylim([5 20])
set(gca,'fontsize',fs); set(gca,'yscale','log');
text(0.05,0.95,'\bf{(b)} - La Niña','fontsize',fs+2,'units','normalized')
xlabel('Sea surface temperature [K]','fontsize',fs-1)

print(fig,'-dpdf','-r200','CC-scaling-extremes')
