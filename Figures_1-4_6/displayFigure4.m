%% code and data for Figure 4
function [n_histo,m_histo,iqr_histo,fraction_data,n_position,n_single] = displayFigure4
% load:
% data - all recorded data
% audIdx - data for analysis only
% tuningcurves - tuningcurves of all units
% tuningproperties - best frequency and bandwidth of each tuning curve
% output:
%   - n_histo: number of enhanced (1.row) and reduced (2.row) units in
%   total (left) and with r>0.5 (right)
%   - m_histo: median values for slope (left) and intercept (right) for 
%       enhanced (1.row) and reduced (2.row) units 
%   - iqr_histo: inter quartile range for slope (left) and intercept (right) for 
%       enhanced (1.row) and reduced (2.row) units 
%   - fraction_data: mean (1.row), std (2.row) and SEM (3.row) for 
%       divisive (1.col), subtractive (2.col), both (3.col) or none (4.col)
%       fraction
%   - n_position: number of recording positions
%   - n_single: number of units

load('sigDataArray.mat','dataArray','audIdx')
load('TuningProperties.mat')     %#ok<LOAD> % best frequency and bandwidth of all units
load('tuningcurves.mat','tuningcurve_onset','f')         % tuning curves of all units

% effect change
percCh = [dataArray{:,10}]';  %#ok<*IDISVAR,NODEF>

% load new sign Subgroups and replace the old ones
load('newNegPosIdcs.mat','newPosIdx','newNegIdx')

% unit IDs
imun = dataArray(:,1);
imTbl = local_imun2table(imun);

nU = length(percCh);

idcs = {newPosIdx & audIdx,newNegIdx& audIdx};
lbls = {'enhanced units', 'reduced units'};

% pre-allocate
pVal_slope = nan(nU,1);
pVal_interc = nan(nU,1);


%% collect slope and intercept 

for iU = 1:nU
    thisOn = tuningcurve_onset{iU,2}; %#ok<USENS>
    thisOff = tuningcurve_onset{iU,1};
    m = max(thisOff);
    
    % normalize
    thisOn = thisOn/m;
    thisOff = thisOff/m;
    invalid = thisOn < 0.1 | thisOff < 0.1;
    thisOn(invalid) = [];
    thisOff(invalid) = [];
    nSamples = sum(~invalid);
    if nSamples>3 % on
        [r_tuning(iU),p_sign(iU)] = corr(thisOff',thisOn');
        [m,b,~,sm,sb]=local_lsqfitma(thisOff,thisOn);
        p_tuning(iU,:) = [m,b];
        sd_tuning(iU,:) = [sm,sb];
        
        if ~isnan(m) && ~isnan(sm) && sm>0
            pVal_slope(iU) = local_testMaRegress(m,sm,1,nSamples,'right');
        end
        
        if ~isnan(b) && ~isnan(sb) && sb>0
            pVal_interc(iU) = local_testMaRegress(b,sb,0,nSamples,'right');
        end
    end
end

%% plot figure 3A
fHist = figure;
col = [0.2 0.6 0;0.6 0.4 0.8];
gr = [0.4, 0.4, 0.4];

val = abs(r_tuning') > 0.5; % exclude low correlations from the plots
n_histo = zeros(2,2);
m_histo = zeros(2,2);
iqr_histo = zeros(2,2);

nCat = numel(idcs);
for iI = 1:nCat
    figure(fHist)
    thisIdcs = idcs{iI};
    n_histo(iI,1) = sum(thisIdcs);
    n_histo(iI,2) = sum(thisIdcs&val);
    
    subplot(nCat,3,(iI-1)*3+1);
    histogram(r_tuning(thisIdcs),-0.5:0.1:1,'FaceColor',col(iI,:))
    ylabel(lbls{iI});
    if iI == 1
        title('correlations coeff.')
    elseif iI == 3
        xlabel('r')
    end
    xlim([-0.5 1])
    set(gca,'tickdir','out','box','off')
    
    subplot(nCat,3,(iI-1)*3+2)
    if iI == 1
        histogram(p_tuning(thisIdcs&val,1),-0.5:0.15:2.5,'FaceColor',col(iI,:))
    else
        histogram(p_tuning(thisIdcs&val,1),-0.5:0.1:2.5,'FaceColor',col(iI,:))
    end
    if iI == 1
        title('slope')
    elseif iI == 3
        xlabel('slope [Hz/Hz]')
    end
    xlim([-0.5 2.5])
    set(gca,'tickdir','out','box','off')
    m_histo(iI,1) = median(p_tuning(thisIdcs&val,1));
    iqr_histo(iI,1) = iqr(p_tuning(thisIdcs&val,1));
    
    subplot(nCat,3,(iI-1)*3+3)
    if iI == 1
        histogram(p_tuning(thisIdcs&val,2),-0.5:0.075:2.5,'FaceColor',col(iI,:))
    else
        histogram(p_tuning(thisIdcs&val,2),-0.5:0.075:2.5,'FaceColor',col(iI,:))
    end
    if iI == 1
        title('intercept, rel. to max')
    elseif iI == 3
        xlabel('norm. intercept')
    end
    xlim([-0.5 1])
    set(gca,'tickdir','out','box','off')
    m_histo(iI,2) = median(p_tuning(thisIdcs&val,2));
    iqr_histo(iI,2) = iqr(p_tuning(thisIdcs&val,2));    
    
end

% adapt yLim for each subplot
yEnh = [0,20];
yRed = [0,200];
subplot(2,3,1)
set(gca,'YLim',yEnh)
line([0.5 0.5],yEnh,'color','k','LineStyle','--')
subplot(2,3,2)
set(gca,'YLim',yEnh)
line([1 1],yEnh,'color','r')
subplot(2,3,3)
set(gca,'YLim',yEnh)
line([0 0],yEnh,'color','r')
subplot(2,3,4)
set(gca,'YLim',yRed)
line([0.5 0.5],yRed,'color','k','LineStyle','--')
subplot(2,3,5)
set(gca,'YLim',yRed)
line([1 1],yRed,'color','r')
subplot(2,3,6)
set(gca,'YLim',yRed)
line([0 0],yRed,'color','r')


%% look at the significance levels in each animal and make the summary plot of fraction of cells
% figure 3B
pIdx = idcs{2}; % take only reduced units
 a = findgroups(imTbl.animal(pIdx));

sigData =  [];
sigData(:,1) = splitapply(@(x,y) mean(x>=0.05&y<0.05),pVal_interc(pIdx),pVal_slope(pIdx),a);
sigData(:,2) = splitapply(@(x,y) mean(x<0.05&y>=0.05),pVal_interc(pIdx),pVal_slope(pIdx),a);
sigData(:,3) = splitapply(@(x,y) mean(x<0.05&y<0.05),pVal_interc(pIdx),pVal_slope(pIdx),a);
sigData(:,4) = splitapply(@(x,y) mean(x>=0.05&y>=0.05),pVal_interc(pIdx),pVal_slope(pIdx),a);

fraction_data(1,:) = mean(sigData);
fraction_data(2,:) = std(sigData);
fraction_data(3,:) = std(sigData)./sqrt(length(sigData)); %SEM

figure
Colors=[gr;gr;gr;gr];
UnivarScatter(sigData,'Label',{'Divisive','Subtractive','Both','None'},...
    'MarkerFaceColor',Colors,'Compression',3, 'Width', 0.3,'SEMColor',Colors*1.5,'StdColor',Colors*2);
box off
set(gca,'tickdir','out','fontsize',14)
ylabel('Fraction')
ylim([0 1])
xlim([0.6 4.3])
title('Fraction divisive/subtractive')
set(findall(gcf,'-property','FontSize'),'FontSize',16)

%%  check change in slope and intercept for each position, plot against rate reduction in that position
p = findgroups(imTbl.animal,imTbl.position); % each position in an animal defines a group
mSlope = splitapply(@(x) median(x,'omitnan'),p_tuning(:,1),p);
mInterc = splitapply(@(x) median(x,'omitnan'),p_tuning(:,2),p);
mChange = splitapply(@(x) median(x,'omitnan'),percCh,p);
n_position = length(mChange);

h =@(a,b,c,d) median(((a-1)./b)-(c./d),'omitnan'); 
mDivisiveness  = splitapply(h,p_tuning(:,1),sd_tuning(:,1),p_tuning(:,2),sd_tuning(:,2),p);
nDiv = numel(mDivisiveness);
n = numel(mChange);

figure
lw = 1; % linewidth

subplot(2,1,1)
plot(mChange,mSlope,'o','Color',gr,'LineWidth',lw)
set(gca,'YLim',[0 2.3])
ylabel('slope coeeficient [Hz/Hz]')
xlim([-70 40])
[m,b,r] = local_lsqfitma(mChange,mSlope);
line(xlim,m*xlim+b,'linestyle','--');
v = n-1;
tval = -r*sqrt((n-2)/(1-r^2));
tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution
p = 1-tdist2T(tval,v);

tx = -65;
ty = max(ylim)-0.1*diff(ylim);

text(tx,ty,sprintf('r=%0.2f, p=%0.2d',r,p),'fontsize',12)
set(gca,'tickdir','out','fontsize',14,'box','off')
ylim([0 2.1])
title('slope')

subplot(2,1,2)
plot(mChange,mInterc,'o','Color',gr,'LineWidth',lw)
ylabel(' y-intercept (norm.)')
xlim([-70 40])
[~,~,r] = local_lsqfitma(mChange,mInterc);
v = n-1;
tval = -r*sqrt((n-2)/(1-r^2));
tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution
p = 1-tdist2T(tval,v);

tx = -15;
ty = min(ylim)+0.1*diff(ylim);
text(tx,ty,sprintf('r=%0.2f, p=%0.3f',r,p),'fontsize',12)

set(gca,'YLim',[-0.23 0.28])
set(gca,'tickdir','out','fontsize',14,'box','off')
ylim([-0.51 0.51])
title('y-intercept')

set(findall(gcf,'-property','FontSize'),'FontSize',16)

%% plot slope and intercept for each unit as a function of its % change
figure
sz = 20; %markersize
c = gr; % color
mea = 0.6; % marker edge alpha transpancy

ax(1) = subplot(221);
scatter(percCh(val&percCh<0),p_tuning(val&percCh<0,1),sz,c,'MarkerEdgeAlpha',mea);
ylabel('slope')
line(xlim,[1 1],'Color','r')
[~,p(1,1)] = corr(percCh(val&percCh<0),p_tuning(val&percCh<0,1));
[m,b,r(1,1),sm,~]=local_lsqfitma(percCh(val&percCh<0),p_tuning(val&percCh<0,1));
line(xlim,b+m*xlim,'color','k')
title('negative change, slope')
thisP = local_testMaRegress(m,sm,1,sum(val&percCh<0)); %#ok<*NASGU>
text(-90,1.8,sprintf('r = %0.2f, p=%0.2d',r(1,1),p(1,1)));
n_single(1,1) = length(percCh(val&percCh<0));

ax(2) = subplot(222);
scatter(percCh(val&percCh>0),p_tuning(val&percCh>0,1),sz,c,'MarkerEdgeAlpha',mea);
ylabel('slope')
line(xlim,[1 1],'Color','r')
[~,p(1,2)] = corr(percCh(val&percCh>0),p_tuning(val&percCh>0,1));
[m,b,r(1,2),sm,~]=local_lsqfitma(percCh(val&percCh>0),p_tuning(val&percCh>0,1));
line(xlim,b+m*xlim,'color','k')
title('positve change, slope')
thisP = local_testMaRegress(m,sm,1,sum(val&percCh>0));
text(20,0.3,sprintf('r = %0.2f, p=%0.3f',r(1,2),p(1,2)));
n_single(1,2) = length(percCh(val&percCh>0));

ax(3) = subplot(223);
scatter(percCh(val&percCh<0),p_tuning(val&percCh<0,2),sz,c,'MarkerEdgeAlpha',mea);
ylabel('intercept')
line(xlim,[0 0],'Color','r')
[~,p(2,1)] = corr(percCh(val&percCh<0),p_tuning(val&percCh<0,2));
[m,b,r(2,1),~,~]=local_lsqfitma(percCh(val&percCh<0),p_tuning(val&percCh<0,2));
line(xlim,b+m*xlim,'color','k')
xlabel('modulation [%]')
title('negative change, intercept')
text(-90,0.4,sprintf('r = %0.2f, p=%0.3f',r(2,1),p(2,1)));

ax(4) = subplot(224);
scatter(percCh(val&percCh>0),p_tuning(val&percCh>0,2),sz,c,'MarkerEdgeAlpha',mea);
ylabel('intercept')
line(xlim,[0 0],'Color','r')
[~,p(2,2)] = corr(percCh(val&percCh>0),p_tuning(val&percCh>0,2));
[m,b,r(2,2),~,~]=local_lsqfitma(percCh(val&percCh>0),p_tuning(val&percCh>0,2));
line(xlim,b+m*xlim,'color','k')
xlabel('modulation [%]')
title('positve change, intercept')
text(20,-0.3,sprintf('r = %0.2f, p=%0.3f',r(2,2),p(2,2)));

subplot(221), ylim([0 2.1])
subplot(222), axis([0 100 0 2.1]) 

subplot(223), ylim([-0.51 0.51])
subplot(224), axis([0 100 -0.51 0.51]) 

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
function [m,b,r,sm,sb]=local_lsqfitma(X,Y)

% Determine the size of the vector
 
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


function local_matchAxis(ax) %#ok<DEFNU>
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

line(xlim,ylim,'color','k')
end

function imTbl = local_imun2table(imun)
if ~iscell(imun)
    imun = {imun};
end

nI = length(imun);
imunStr(nI).animal = [];
for iI=1:nI
    imunStr(iI).animal = ['h0' imun{iI}(1:3)];
    imunStr(iI).position = - 1e-4*eval(imun{iI}(4:8));
    imunStr(iI).sorting = eval(imun{iI}(9));
    imunStr(iI).probName = imun{iI}(10:17);
    imunStr(iI).unitNo = eval(imun{iI}(end-1:end));
end

imTbl = struct2table(imunStr);
end

function p = local_testMaRegress(m,sd,refVal,N,tail)
%TESTMAREGRESS test the results of the Major axis regression (t-test)
%   m  -  mean of the regression parameter (slope or intercept)
%   sd - standard deviation of the regression paramter
%   N  - number of points that went into the regression
%   ref 


if nargin <= 4
    tail = 'both';
end

v = N-1;
tval = (m-refVal) / sqrt(sd^2/N);       % Calculate T-Statistic
tdist2 = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution



switch tail
    case 'left'
        tdist1 = @(t,v) 1-(1-tdist2(t,v))/2; 
        p = abs((1-tval<0) - tdist1(tval,v));
    case 'right'
        tdist1 = @(t,v) 1-(1-tdist2(t,v))/2;              % 1-tailed t-distribution
        p = abs((1-tval>0) - tdist1(tval,v));
    case 'both'
        
        p = 1-tdist2(tval,v);
end

end