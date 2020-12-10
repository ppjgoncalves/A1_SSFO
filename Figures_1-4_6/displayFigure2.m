%% code and data for Figure 2
function [n_spikes,stat_troughs,stat_widths,s_width,n] = displayFigure2
% load:
%   - spikeforms: normalized waveforms of all spikes
%   - zeroIdx: index for single-units only (no 0-cluster)
%   - peakIdx: correct position of waveform's peak (appropriate spikewaves)
%   - widht: width of waveform in ms
%   - trough: normalized distance from baseline to trough
%   - dataArray: data of all recorded data
%   - audIdx: data for analysis only
%   - newPosIdx: index for sign. enhanced units
%   - newNegIdx: index for sign. reduced units

% output:
%   - n_spikes: number of enhanced and reduced (2.row) units with spike
%   waveforms
%   - stat_troughs: mean, std p trough values (rows) for enhanced and reduced units (col.)
%   - stat_widths: mean, std p width values (rows) for enhanced and reduced units (col.)
%   - s_width: mean, std, n, p of effect values (rows) for units with a
%   small and large spike widths (col.)
%   - n: number of units in scatter plot for enhanced and reduced units

load('spikeforms.mat','spikeforms','zeroIdx','peakIdx','width','trough')  
load('sigDataArray.mat','dataArray','audIdx')
load('newNegPosIdcs.mat','newPosIdx','newNegIdx') 
percCh = [dataArray{:,10}]';   %#ok<IDISVAR,USENS>
stIdx = [newPosIdx,newNegIdx];

% figure A 
n_spikes = local_showSpikes(spikeforms,zeroIdx,peakIdx,audIdx,stIdx);

% figure B
[stat_troughs,stat_widths] = local_showRatioHistogram(trough,width,...
    zeroIdx,peakIdx,audIdx,stIdx);

% figure C
s_width = local_showDeltaFRHistogram(percCh,width,zeroIdx,peakIdx,audIdx,stIdx);

% figure D
n = local_plotWidthTrough(width,trough,audIdx,zeroIdx,peakIdx,stIdx);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
% figure A
function n_spikes = local_showSpikes(data,zeroIdx,peakIdx,audIdx,stIdx)
Col = [0.2 0.6 0;0.6 0.4 0.8];
Fs = 30e3;

figure
h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape')
% setting
nBins = size(data,2);
x = 1/Fs*1e3 * (1:nBins); % time line in ms

% sum up total number of units each
n = NaN(1,2);

for i = 1:2
    theseUnits = stIdx(:,i);

    tmp_idx = [theseUnits,audIdx,~zeroIdx,peakIdx];
    tmp_idx = logical(sum(tmp_idx,2) == size(tmp_idx,2)); % all indices musst be true
    
    n(i) = sum(tmp_idx);
    
    thisData = data(tmp_idx,:);
    y = mean(thisData);
    SEM = std(thisData) / sqrt(size(thisData,1)); % standard error of the mean
    
    % plot + - SEM first
    low = y - SEM;
    high = y + SEM;
    h =  plot(x,high,'Color',Col(i,:),'LineWidth',0.2);
    h.Color(4) = 0.4;
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on
    h = plot(x,low,'Color',Col(i,:),'LineWidth',0.2);
    h.Color(4) = 0.4;
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    % fill area in between
    x2 = [x,fliplr(x)];
    inBetween = [high,fliplr(low)];
    h = fill(x2,inBetween,Col(i,:),'LineStyle','none');
    set(h,'facealpha',.4)
    
    % insert mean
    h = plot(x,y,'LineWidth',1,'Color',Col(i,:));
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
end

hold off
axis tight
set(gca,'FontSize',12)
set(gca,'TickDir','out')
ylabel('voltage [normalized]','FontSize',14)
xlabel('time [ms]','FontSize',14)
box off

legend({'sign.enhanced','sign.reduced'},'location','northeast','FontSize',14)
legend boxoff

n_spikes(1) = n(1);
n_spikes(2) = n(2);
end

% figure B
function [stat_troughs,stat_widths] = local_showRatioHistogram(trough,width,...
    zeroIdx,peakIdx,audIdx,stIdx)

figure
h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape')
Col = [0.2 0.6 0;0.6 0.4 0.8];
n =  zeros(1,2);

% find the right egdes
tmp_idx = [stIdx,audIdx,~zeroIdx,peakIdx];
tmp_idx = logical(sum(tmp_idx,2) == 4); % all indices musst be true

trough_edges = min(-trough(tmp_idx)) : 0.1 : max(-trough(tmp_idx)); 
width_edges = 0 : 0.02 : max(width(tmp_idx));

allTrough = []; %statistics
allWidth = []; %statistics
groups = []; %statistics

for i = 1:2
    theseUnits = stIdx(:,i);
    tmp_idx = [theseUnits,audIdx,~zeroIdx,peakIdx];
    tmp_idx = logical(sum(tmp_idx,2) == size(tmp_idx,2)); % all indices musst be true
    n(i) = sum(tmp_idx);
    allTrough = [allTrough;-trough(tmp_idx)];
    allWidth = [allWidth;width(tmp_idx);];
    groups = [groups;i*ones(n(i),1)]; %statistics
    
    % take troughs
    theseTroughs = -trough(tmp_idx);
    N = histcounts(theseTroughs,trough_edges,'Normalization', 'probability');
    
    subplot(2,1,1)
    plot(trough_edges(2:end),N,'LineWidth',1.5,'Color',Col(i,:))
    hold on
    if i == 3
        hold off
        axis tight
        set(gca,'FontSize',12)
        ylabel({'number of units';'[probability]'},'FontSize',14)
        xlabel('trough amplitude [normalized]','FontSize',14)
        box off
    end
    
    % take widths
    theseWidths = width(tmp_idx);
    N = histcounts(theseWidths,width_edges,'Normalization', 'probability');
    
    subplot(2,1,2)
    plot(width_edges(2:end),N,'LineWidth',1.5,'Color',Col(i,:))
    hold on
    if i == 3
        hold off
        axis tight
        set(gca,'FontSize',12)
        ylabel({'number of units';'[probability]'},'FontSize',14)
        xlabel('spike width [ms]','FontSize',14)
        box off
    end
end

subplot(2,1,1)
set(gca,'XLim',[-0.3627 1.2])
subplot(2,1,2)
legend({'sign.enhanced','sign.reduced'},'location','northeast','FontSize',14)
legend boxoff

% mean, std p values (rows) for enhanced and reduced units (col.)
% trough values
stat_troughs = NaN(3,2);
stat_troughs(1,1) = mean(allTrough(groups == 1));
stat_troughs(1,2) = mean(allTrough(groups == 2));
stat_troughs(2,1) = std(allTrough(groups == 1));
stat_troughs(2,2) = std(allTrough(groups == 2));
[~,~,STATS] = kruskalwallis(allTrough, groups);
[c,~] = multcompare(STATS); % Tukey Test
stat_troughs(3,1) = c(6);

% trough values
stat_widths = NaN(3,2);
stat_widths(1,1) = mean(allWidth(groups == 1));
stat_widths(1,2) = mean(allWidth(groups == 2));
stat_widths(2,1) = std(allWidth(groups == 1));
stat_widths(2,2) = std(allWidth(groups == 2));
[~,~,STATS] = kruskalwallis(allWidth, groups);
[c,~] = multcompare(STATS); % Tukey Test
stat_widths(3,1) = c(6);


end

% figure C
function s = local_showDeltaFRHistogram(effect,width, zeroIdx,peakIdx,audIdx,stIdx)
% histogram of firing rates for each group, distinguished bei 0.14octave
% width

% borderwidth to split the data set
border = 0.12; 

figure('units', 'centimeters', 'position', [12 8 14 8]);
% get(gca,'colororder')

% split units - index vector
tmp_idx = [stIdx,audIdx,~zeroIdx,peakIdx];
tmp_idx = logical(sum(tmp_idx,2) == 4); % all indices musst be true

small = width <= border;
huge = ~small;
small = small & tmp_idx;
huge = huge & tmp_idx;

% insert stuff
% median(effect(small))
% mad(effect(small))
% mean(effect(small))
% std(effect(small))
% median(effect(huge))
% mad(effect(huge))
% mean(effect(huge))
% std(effect(huge))

% effect: everything higher than 100 should scaled down to 100
effect(effect>100) = 100;

% histogram
edges = -100 : 5 : 100;
histogram(effect(small),edges,'Normalization','probability','FaceColor',[.1,.1,.1])
hold on
histogram(effect(huge),edges,'Normalization','probability','FaceColor',[.9,.9,.9]);
hold off
set(gca,'FontSize',12)
ylabel({'number of units';'[probability]'},'FontSize',14)
xlabel('firing rate change after activation [%]','FontSize',14)
legend({'width <= 0.12ms','width > 0.12ms'},'location','northeast','FontSize',14)
legend boxoff

% statistics: mean, std, n, p (rows) for small and large widths (col.)
s =  zeros(4,2);
% small width mean, std, n
s(1,1) = mean(effect(small));
s(2,1) = std(effect(small));
s(3,1) = sum(small);
% large width mean, std, n
s(1,2) = mean(effect(huge));
s(2,2) = std(effect(huge));
s(3,2) = sum(huge);
% Mann-Whitney U Test
[p,~,~] = ranksum(effect(small),effect(huge));
s(4,1) = p;
end

% figure D
function n = local_plotWidthTrough(width,trough,audIdx,zeroIdx,peakIdx,stIdx)

figure
h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape')
Col = [0.6 0.4 0.8;0.2 0.6 0];
n =  zeros(1,2);

for i = 1:2
    if i == 1
        iG = 2;
    else
        iG = 1;
    end
    
    theseUnits = stIdx(:,iG);
    
    tmp_idx = [theseUnits,audIdx,~zeroIdx,peakIdx];
    tmp_idx = logical(sum(tmp_idx,2) == size(tmp_idx,2)); % all indices musst be true
    n(i) = sum(tmp_idx);
    
    % take widths
    theseWidths = width(tmp_idx);
    ratio =  -trough(tmp_idx);
    
    if i == 1
        scatter(theseWidths,ratio,60,Col(i,:)) % (x,y,MarkerSize)
        hold on
    else
        scatter(theseWidths,ratio,60,Col(i,:),'filled')
        hold off   
    end
end

%set(gca,'YLim',[0.5 1])
set(gca,'FontSize',12)
ylabel('trough amplitude [normalized]','FontSize',14)
xlabel('spike width [ms]','FontSize',14)

legend({'sign.reduced','sign.enhanced'},'location','northeast','FontSize',14)
legend boxoff
end