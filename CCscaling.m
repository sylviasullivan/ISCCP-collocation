clear all
en = readNPY('trop_SST_CWVC_pmax_psum_endjf.npy');
ln = readNPY('trop_SST_CWVC_pmax_psum_lndjf.npy');

% save the values read in from the numpy files
sstEN = en(1,:); 
sstLN = ln(1,:);
pmaxEN = en(3,:); 
pmaxLN = ln(3,:);
cwvcEN = en(5,:);
cwvcLN = ln(5,:);

% filter for instances where pmax is non-zero and there is underlying SST
ii = find(~isnan(pmaxEN) & pmaxEN ~= 0 & sstEN > 0);
pmaxEN = pmaxEN(ii);
sstEN = en(1,ii);
depthEN = sstEN - en(2,ii);
cwvcEN = cwvcEN(ii);

ii = find(~isnan(pmaxLN) & pmaxLN ~= 0 & sstLN > 0);
pmaxLN = pmaxLN(ii);
sstLN = ln(1,ii);
depthLN = sstLN - ln(2,ii);
cwvcLN = cwvcLN(ii);
clear ii

%%
% filter for d1, d2, d3 systems
pr1 = struct('ENd1',0,'ENd2',0,'ENd3',0,'LNd1',0,'LNd2',0','LNd3',0,...
    'EN99d1',0,'EN99d2',0,'EN99d3',0,'LN99d1',0,'LN99d2',0,'LN99d3',0);
pr2 = struct('ENd1',0,'ENd2',0,'ENd3',0,'LNd1',0,'LNd2',0','LNd3',0,...
    'EN99d1',0,'EN99d2',0,'EN99d3',0,'LN99d1',0,'LN99d2',0,'LN99d3',0);
sst = struct('ENd1',0,'ENd2',0,'ENd3',0,'LNd1',0,'LNd2',0','LNd3',0,...
    'EN99d1',0,'EN99d2',0,'EN99d3',0,'LN99d1',0,'LN99d2',0,'LN99d3',0);
cwvc = struct('ENd1',0,'ENd2',0,'ENd3',0,'LNd1',0,'LNd2',0','LNd3',0,...
    'EN99d1',0,'EN99d2',0,'EN99d3',0,'LN99d1',0,'LN99d2',0,'LN99d3',0);

d1 = find(depthEN < 65);
d2 = find(depthEN >= 65 & depthEN < 85);
d3 = find(depthEN >= 65);
pr1.ENd1 = pmaxEN(d1); pr1.ENd2 = pmaxEN(d2); pr1.ENd3 = pmaxEN(d3);
pr2.ENd1 = pmaxEN(d1); pr2.ENd2 = pmaxEN(d2); pr2.ENd3 = pmaxEN(d3);
sst.ENd1 = sstEN(d1); sst.ENd2 = sstEN(d2); sst.ENd3 = sstEN(d3);
cwvc.ENd1 = cwvcEN(d1); cwvc.ENd2 = cwvcEN(d2); cwvc.ENd3 = cwvcEN(d3);

d1 = find(depthLN < 65);
d2 = find(depthLN >= 65 & depthLN < 85);
d3 = find(depthLN >= 85);
pr1.LNd1 = pmaxLN(d1); pr1.LNd2 = pmaxLN(d2); pr1.LNd3 = pmaxLN(d3);
pr2.LNd1 = pmaxLN(d1); pr2.LNd2 = pmaxLN(d2); pr2.LNd3 = pmaxLN(d3);
sst.LNd1 = sstLN(d1); sst.LNd2 = sstLN(d2); sst.LNd3 = sstLN(d3);
cwvc.LNd1 = cwvcLN(d1); cwvc.LNd2 = cwvcLN(d2); cwvc.LNd3 = cwvcLN(d3);

clear d1 d2 d3 pmaxEN pmaxLN ii depthEN sstEN sstLN depthLN en ln cwvcEN cwvcLN

%%
% SST binning and extraction of the maximum 5 precip vals from each
depth = {'d1','d2','d3'};
% how many bins?
num = 15; num2 = 5;
for ll = 1:3
    % EL NINO
    temps = linspace(285,307,num);
    prQ = []; tQ = []; pr99 = []; t99 = [];
    for ii = 1:length(temps)-1
        jj = find(sst.(strcat('EN',depth{ll})) >= temps(ii) & ...
            sst.(strcat('EN',depth{ll})) < temps(ii+1));
        
        % read in the relevant values
        precip = pr1.(strcat('EN',depth{ll}))(jj);
        temperature = sst.(strcat('EN',depth{ll}))(jj);
        
        % calculate the quantiles of the subsets in this bin
        qpr = quantile(precip,[0.2 0.5 0.75 0.9]);
        qt = quantile(temperature,[0.2 0.5 0.75 0.9]);

        % append the quantiles to their respective lists
        prQ = [prQ; qpr]; 
        tQ = [tQ; qt]; 
        
        % extract the five largest precip values in each bin and their
        % correspondings conditions
        [prpr,i] = maxk(precip,num2);
        pr99 = [pr99, prpr];
        t99 = [t99, temperature(i)];
    end
    pr1.(strcat('ENq',depth{ll})) = prQ; 
    sst.(strcat('ENq',depth{ll})) = tQ; 
    
    pr1.(strcat('EN99',depth{ll})) = pr99; 
    sst.(strcat('EN99',depth{ll})) = t99; 
    
    % LA NINA
    temps = linspace(285,307,num);
    prQ = []; tQ = []; pr99 = []; t99 = []; 
    for ii = 1:length(temps)-1
        jj = find(sst.(strcat('LN',depth{ll})) >= temps(ii) & ...
            sst.(strcat('LN',depth{ll})) < temps(ii+1));
        
        % read in the relevant values
        precip = pr1.(strcat('LN',depth{ll}))(jj);
        temperature = sst.(strcat('LN',depth{ll}))(jj);
        
        % calculate the quantiles of the subsets in this bin
        qpr = quantile(precip,[0.2 0.5 0.75 0.9]);
        qt = quantile(temperature,[0.2 0.5 0.75 0.9]);
        
        % append the quantiles to their respective lists
        prQ = [prQ; qpr]; 
        tQ = [tQ; qt]; 
        
        % extract the five largest precip values in each bin and their
        % correspondings conditions
        [prpr,i] = maxk(precip,num2);
        pr99 = [pr99, prpr];
        t99 = [t99, temperature(i)];
    end
    pr1.(strcat('LNq',depth{ll})) = prQ; 
    sst.(strcat('LNq',depth{ll})) = tQ; 
    
    pr1.(strcat('LN99',depth{ll})) = pr99; 
    sst.(strcat('LN99',depth{ll})) = t99; 
end

%%
for ll = 1:3
    % EL NINO
    temps = linspace(5/100.,65/100.,num);
    prQ = []; cwvcQ = []; pr99 = []; cwvc99 = [];
    for ii = 1:length(temps)-1
        jj = find(cwvc.(strcat('EN',depth{ll})) >= temps(ii) & ...
            cwvc.(strcat('EN',depth{ll})) < temps(ii+1));
        
        % read in the relevant values
        precip = pr2.(strcat('EN',depth{ll}))(jj);
        colwater = cwvc.(strcat('EN',depth{ll}))(jj);
        
        % calculate the quantiles of the subsets in this bin
        qpr = quantile(precip,[0.2 0.5 0.75 0.9]);
        qcwvc = quantile(colwater,[0.2 0.5 0.75 0.9]);

        % append the quantiles to their respective lists
        prQ = [prQ; qpr]; 
        cwvcQ = [cwvcQ; qcwvc];
        
        % extract the five largest precip values in each bin and their
        % correspondings conditions
        [prpr,i] = maxk(precip,num2);
        pr99 = [pr99, prpr];
        cwvc99 = [cwvc99, colwater(i)];  
    end
    pr2.(strcat('ENq',depth{ll})) = prQ; 
    cwvc.(strcat('ENq',depth{ll})) = cwvcQ;
    
    pr2.(strcat('EN99',depth{ll})) = pr99; 
    cwvc.(strcat('EN99',depth{ll})) = cwvc99;
    
    % LA NINA
    temps = linspace(5/100.,65/100.,num);
    prQ = []; cwvcQ = []; pr99 = []; cwvc99 = [];
    for ii = 1:length(temps)-1
        jj = find(cwvc.(strcat('LN',depth{ll})) >= temps(ii) & ...
            cwvc.(strcat('LN',depth{ll})) < temps(ii+1));
        
        % read in the relevant values
        precip = pr2.(strcat('LN',depth{ll}))(jj);
        colwater = cwvc.(strcat('LN',depth{ll}))(jj);
        
        % calculate the quantiles of the subsets in this bin
        qpr = quantile(precip,[0.2 0.5 0.75 0.9]);
        qcwvc = quantile(colwater,[0.2 0.5 0.75 0.9]);
        
        % append the quantiles to their respective lists
        prQ = [prQ; qpr]; 
        cwvcQ = [cwvcQ; qcwvc];
        
        % extract the five largest precip values in each bin and their
        % correspondings conditions
        [prpr,i] = maxk(precip,num2);
        pr99 = [pr99, prpr];
        cwvc99 = [cwvc99, colwater(i)];
    end
    pr2.(strcat('LNq',depth{ll})) = prQ; 
    cwvc.(strcat('LNq',depth{ll})) = cwvcQ;
    
    pr2.(strcat('LN99',depth{ll})) = pr99; 
    cwvc.(strcat('LN99',depth{ll})) = cwvc99;
end

clear depth ii jj pr99 t99 cwvc99 prQ tQ cwvcQ i ll num num2 precip
clear qpr qt qcwvc prpr colwater temperature temps
