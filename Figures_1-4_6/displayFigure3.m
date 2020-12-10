%% code and data for Figure 3
function [m,n,p] = displayFigure3
% load:
% data - all recorded data
% audIdx - data for analysis only
% tuningcurves - tuningcurves of all units
% tuningproperties - best frequency and bandwidth of each tuning curve
% output:
%   - m: median change changes in octave for BF and BW (rows) for enhanced
%   and reduced units (columns)
%   - p: p-value wilcoxon ranksum test for BF and BW (rows) for enhanced
%   and reduced units (columns)
%   - n: number of used units for BF and BW (rows) for enhanced
%   and reduced units (columns)

load('sigDataArray.mat','dataArray','audIdx')
load('TuningProperties.mat')     %#ok<LOAD> % best frequency and bandwidth of all units
load('tuningcurves.mat','tuningcurve_onset','f')         % tuning curves of all units
load('newNegPosIdcs.mat','newPosIdx','newNegIdx') % load significant groups

% only take cells that are actually active
posIdx = newPosIdx & audIdx;
negIdx = newNegIdx & audIdx;

tc = tuningcurve_onset; % abbrev.
bf = BF{1}(:,1);  %#ok<IDISVAR,USENS>
nU = length(bf);
nF = length(f);


%% figure 2A 
% tuningcurve difference- overview
allDiff = NaN(nU,2*nF -1);

% organize the overview plot
for iU = 1:nU 
    
    thisTune(1,:) = tc{iU,1}(1:34);
    thisTune(2,:) = tc{iU,2}(1:34);
    
    % normalize
    thisTune = thisTune/max(thisTune(1,:));
    
    % find BF-index
    [~,bfIdx] = min(abs(f-bf(iU)));
    idx = nF-bfIdx + (1:nF);
    
    % tuning difference
    tuneDiff = diff(thisTune);%./thisTune(1,:);
    
    % now put it in the right place in the matrix
    allDiff(iU,idx) = tuneDiff;
end

% sort by the change at BF
[~,sIdx] = sort(allDiff(:,nF));
allDiff = allDiff(sIdx,:);
allDiff(allDiff<-1) = -1;
allDiff(allDiff>1) = 1;


cmap = [1 1 1; cool];
octTicks = 2:8:(2*nF-1); % tick positions for octaves
octLbls = -4:4;

figure
colormap(cmap)
subplot(1,7,1:3)
imagesc(allDiff(posIdx(sIdx),:));
hold on
plot(nF,36,'bo',nF,99,'bo') % example cells
set(gca,'CLim',[-1.1 1.1],'XTick',octTicks,'XTickLabel',octLbls,'box', 'off','tickdir','out')
xlabel('frequency rel. to BF')
ylabel('units sorted by the change at BF')
title('enhanced units')

subplot(1,7,4:6)
imagesc(allDiff(negIdx(sIdx),:));
hold on
plot(nF,88,'bo',nF,928,'bo') % example cells
set(gca,'CLim',[-1.1 1.1],'XTick',octTicks,'XTickLabel',octLbls,'box', 'off','tickdir','out')
xlabel('frequency rel. to BF')
title('reduced units')

subplot(1,7,7)
c = colorbar;
set(gca,'CLim',[-1.1 1.1])
c.Label.String = 'rate change after SSFO activation (a.u.)';
axis off


%% Figure 2B and 2C
% BW and BF change

bw = BW{1}; %#ok<IDISVAR,USENS> % onset only
bf = BF{1}; % onset only

f2 = figure;
subplot(2,2,1)
scatter(1e-3*bf(posIdx,1),1e-3*bf(posIdx,2),'k')
set(gca,'xscale','log','yscale','log')
set(gca,'XLim',[4 100],'YLim',[4 100])
axis square
box on
local_matchAxis
xlabel('best frequency [kHz], control')
ylabel('best frequency [kHz], PV-SSFO')
set(gca,'TickDir','out')

subplot(2,2,3)
scatter(bw(posIdx,1),bw(posIdx,2),'k')
ylim([0 4])
xlim([0 4])
axis square
box on
local_matchAxis
xlabel('bandwidth [oct], control')
ylabel('bandwidth [oct], PV-SSFO')
set(gca,'TickDir','out')


edgs = logspace(0,2); % logarthmic -> 50 bins from 1 to 100 Hz

subplot(2,2,2)
histogram2(1e-3*bf(negIdx,1),1e-3*bf(negIdx,2),edgs,edgs,'DisplayStyle','tile','EdgeColor','none');
set(gca,'xscale','log','yscale','log')
set(gca,'XLim',[4 100],'YLim',[4 100])
set(gca,'clim',[0 15])
axis square
box on
local_matchAxis
xlabel('best frequency [kHz], control')
set(gca,'TickDir','out')

edgs = 0:0.1:4;
subplot(2,2,4)
histogram2(bw(negIdx,1),bw(negIdx,2),edgs,edgs,'DisplayStyle','tile','EdgeColor','none');
set(gca,'clim',[0 15])
ylim([0 4])
xlim([0 4])
axis square
box on
local_matchAxis
xlabel('bandwidth [oct], control')
set(gca,'TickDir','out')

colorbar

f2.Renderer = 'Painters';


%% Figure 2 examples
exUnits = [36, 99, 88, 928]; 
isNeg = [0 0 1 1];

tC_on = tc(sIdx,2); % must have run % tuningCurveDifferences to have sIdx in the workspace
tC_off = tc(sIdx,1);


for iExpl = 1:4

    neg = isNeg(iExpl);
    exNo = exUnits(iExpl);
    if neg
        negT_on = tC_on(negIdx(sIdx,1)); % this is on the original full data set
        negT_off = tC_off(negIdx(sIdx,1)); % this is on the original full data set
    else
        negT_on = tC_on(posIdx(sIdx,1)); % this is on the original full data set
        negT_off = tC_off(posIdx(sIdx,1)); % this is on the original full data set
    end
    
    figure
    
    ex1 = negT_off{exNo};
    ex1(2,:) = negT_on{exNo};
    ex1 = ex1/max(ex1(1,:));
    
    nf = length(ex1);
    
    % tuning curve plot
    subplot(4,2,[1 3 5])
    
    semilogx(1e-3*f(1:nf),ex1')
    axis tight
    box off
    set(gca,'XTick', [10 20 40], 'tickdir','out')
    
    l=legend('contr.','PV act.');
    l.Box='off';
    ylabel('onset response (a.u.)')
    xlabel('frequency [kHz]')

    % difference
    subplot(4,2,7)
    dxpos = diff(ex1);
    dxpos(dxpos<0) = 0;
    dxneg = diff(ex1);
    dxneg(dxpos>0) = 0;
    
    area(dxpos,'FaceColor','g','EdgeColor','none'), hold on
    area(dxneg,'FaceColor','r','EdgeColor','none')
    axis tight
    box off
    set(gca,'XTickLabel',round(1e-3*f(5:5:30)),'tickdir','out','Color','none','XColor','none')
    ylim(max(abs(ylim))*[-1 1])
    ylabel('rel. change')

    % fits
    ex1(:,any(ex1<0.05)) = []; % omit values less than 10% of maximum rate
    subplot(122)
    scatter(ex1(1,:),ex1(2,:),'k')
    xlim([0 max(xlim)]) % make sure it starts at 0
    axis square
    local_matchAxis
    % add: points which are used for the fit (colored)
    invalid = ex1(2,:) < 0.1 | ex1(1,:) < 0.1;
    ex1(:,invalid) = [];
    hold on
    scatter(ex1(1,:),ex1(2,:),'k','filled')
    hold off
    xlim([0 max(xlim)]) % make sure it starts at 0
    local_matchAxis
    
    [m,b,~]=local_lsqfitma(ex1(1,:),ex1(2,:));
    
    line(xlim,m*xlim+b,'linestyle','--');
    xlabel('spike rate control')
    ylabel('spike rate PV activated')
    set(gca,'tickdir','out')
    t = get(gca,'XTick');
    set(gca,'YTick',t)
    box on
end


%% statistics
% bandwidth 
diffBW = bw(:,2)-bw(:,1);
bw_red = diffBW(negIdx); % reduced units only
bw_red = bw_red(~isnan(bw_red));
bw_enh = diffBW(posIdx); % enhanced units only
bw_enh = bw_enh(~isnan(bw_enh));

m = zeros(2,2); % like figure
m(2,1) = median(bw_enh);
m(2,2) = median(bw_red);

n = zeros(2,2); % like figure
n(2,1) = length(bw_enh);
n(2,2) = length(bw_red);

% wilcoxon ranksum test 
% reduced units
bw_red = bw(negIdx,:);
bw_red = bw_red(~isnan(bw_red(:,1)),:);
[p_red,~,~] = signrank(bw_red(:,1),bw_red(:,2));
% enhanced units
bw_enh = bw(posIdx,:);
bw_enh = bw_enh(~isnan(bw_enh(:,1)),:);
[p_enh,~,~] = signrank(bw_enh(:,1),bw_enh(:,2));

p = zeros(2,2);
p(2,1) = p_enh;
p(2,2) = p_red;

% best frequency
% calculate difference of BF in octave
diffBF = log2(bf(:,2)./bf(:,1)); 
bf_red = diffBF(negIdx,:); % reduced units only
bf_red = bf_red(~isnan(bf_red(:,1)),:);
bf_enh = diffBF(posIdx,:); % enhanced units only
bf_enh = bf_enh(~isnan(bf_enh(:,1)),:);

m(1,1) = median(bf_enh);
m(1,2) = median(bf_red);

% wilcoxon ranksum test 
% reduced units
bf_red = bf(negIdx,:);
bf_red = bf_red(~isnan(bf_red(:,1)),:);
[p_red,~,~] = signrank(bf_red(:,1),bf_red(:,2));
% enhanced units
bf_enh = bf(posIdx,:);
bf_enh = bf_enh(~isnan(bf_enh(:,1)),:);
[p_enh,~,~] = signrank(bf_enh(:,1),bf_enh(:,2));
p(1,1) = p_enh;
p(1,2) = p_red;

n(1,1) = length(bf_enh);
n(1,2) = length(bf_red);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
function local_matchAxis(ax)

if nargin == 0
    ax = gca;
end

axis(ax);
%axis tight
if strcmp(get(ax,'xscale'), 'log') || strcmp(get(ax,'yscale'), 'log')
    xlim([max(min(axis),0) max(axis)]);
else
    xlim([min(axis) max(axis)]);
end

ylim(xlim)

line(xlim,ylim,'color','r')
end

function [m,b,r,sm,sb]=local_lsqfitma(X,Y)
n = length(X);
 % Calculate sums and other re-used expressions
 
Sx = sum(X);
Sy = sum(Y);
xbar = Sx/n;
ybar = Sy/n;
U = X - xbar;
V = Y - ybar;
 
Suv = sum(U .* V);
Su2 = sum(U .^2);
Sv2 = sum(V .^2);
 
sigx = sqrt(Su2/(n-1));
sigy = sqrt(Sv2/(n-1));
 
% Calculate m, b, r, sm, and sb
 
m = (Sv2 - Su2 + sqrt(((Sv2 - Su2)^2) + (4 * Suv^2)))/(2 * Suv);
b = ybar - m * xbar;
r = Suv / sqrt(Su2 * Sv2);
 
sm = (m/r) * sqrt((1 - r^2)/n);
sb1 = (sigy - sigx * m)^2;
sb2 = (2 * sigx * sigy) + ((xbar^2 * m * (1 + r))/r^2);
sb = sqrt((sb1 + ((1 - r) * m * sb2))/n);
end
