%% code and data for Figure 1
function [m,s,r,pcoff,p,n_auditory,n] = displayFigure1
% load:
% data - all recorded data
% audIdx - data for analysis only
% output:
%   - m: mean percentage effect changes over all auditory units for ST and DRC
%   - s: standard deviation of percentage effect changes over all auditory units for ST and DRC
%   - r: correlation coefficient between both stimulus paradigms
%   - p: p-value, Wilcoxon Ranksum test of effect change between both stimulus paradigms
%   - n_auditory: number of all auditory units
%   - n: % number of auditory units with a significant decrease and increase of
            % firing rate after PV action and with no significant change in single
            % tones stimulus paradigm

load('sigDataArray.mat','dataArray','audIdx')
load('newNegPosIdcs.mat') %#ok<LOAD> % load significant groups

%% figure 1c
% examples PSTHs of 2 units for both stimulus paradigms
% unit IDs:
example_1 = '1533887501718192003';
example_2 = '1534262500910111202';

% find example units in data array
idx_1 = strcmp(dataArray(:,1),example_1);  %#ok<*IDISVAR,NODEF>
idx_2 = strcmp(dataArray(:,1),example_2);

local_displayExamplePSTHS(dataArray,idx_1,idx_2)

%% figure 1d
% joint histogram
n_auditory = local_displayJointHistogram(dataArray,audIdx);

%% figure 1e
% distribution figure for ST only for significant effect changes 
local_displayStDistribution(dataArray,audIdx,newNegIdx,newPosIdx);
[m,s,r,pcoff,p,n] = local_information(dataArray,audIdx,newNegIdx,newPosIdx);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1c

function local_displayExamplePSTHS(dataArray,idx_1,idx_2)
% example 1
figure('units', 'centimeters', 'position', [12 8 20 6]);

ST_off = mean(dataArray{idx_1,2},2);
ST_on = mean(dataArray{idx_1,3},2);
t_ST = 0.01:0.01:1;

DRC_off = dataArray{idx_1,4};
DRC_on = dataArray{idx_1,5};
t_DRC = 0.01:0.01:1;

start = 400;

subplot(1,2,1)
plot(t_ST,ST_off,'b','LineWidth',1.5)
hold on
plot(t_ST,ST_on,'r','LineWidth',1.5)
hold off
set(gca,'FontSize',12)
ylabel('firing rate [Hz]','FontSize',14)
xlabel('time [s]','FontSize',14)
title('example 1','FontSize',14)

subplot(1,2,2)
plot(t_DRC,DRC_off(start:start+99),'b','LineWidth',1.5)
hold on
plot(t_DRC,DRC_on(start:start+99),'r','LineWidth',1.5)
set(gca,'XLim',[0 1])
set(gca,'FontSize',12)
ylabel('firing rate [Hz]','FontSize',14)
xlabel('time [s]','FontSize',14)

% example 2
figure('units', 'centimeters', 'position', [12 8 20 6]);

ST_off = mean(dataArray{idx_2,2},2);
ST_on = mean(dataArray{idx_2,3},2);
t_ST = 0.01:0.01:1;

DRC_off = dataArray{idx_2,4};
DRC_on = dataArray{idx_2,5};
t_DRC = 0.01:0.01:1;

start = 400;

subplot(1,2,1)
plot(t_ST,ST_off,'b','LineWidth',1.5)
hold on
plot(t_ST,ST_on,'r','LineWidth',1.5)
hold off
set(gca,'FontSize',12)
ylabel('firing rate [Hz]','FontSize',14)
xlabel('time [s]','FontSize',14)
legend({'control','PV+ activation'},'location','northeast','FontSize',12)
legend boxoff
title('example 2','FontSize',14)

subplot(1,2,2)
plot(t_DRC,DRC_off(start:start+99),'b','LineWidth',1.5)
hold on
plot(t_DRC,DRC_on(start:start+99),'r','LineWidth',1.5)
set(gca,'XLim',[0 1])
set(gca,'FontSize',12)
ylabel('firing rate [Hz]','FontSize',14)
xlabel('time [s]','FontSize',14)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1d
function n_audIdx = local_displayJointHistogram(dataArray,audIdx)
eff = cell2mat(dataArray(:,10:11));

% discard nonsign
eff = eff(audIdx,:);
figure('units', 'centimeters', 'position', [10 8 14 11]);

joint_edges = -100:5:50;
histogram2(eff(:,1),eff(:,2),joint_edges,joint_edges,'DisplayStyle','tile','ShowEmptyBins','off');
axis square
hold on
minL = min(joint_edges);
maxL = max(joint_edges);
line([minL maxL],[minL maxL],'Color','r') % bisecting
hold off
set(gca,'FontSize',12)
xlabel({'ST';'Firing rate change';'after activation [%]'},'FontSize',12)
ylabel({'Firing rate change';'after activation [%]';'DRC'},'FontSize',12)
 
% colorbar
h = colorbar;
caxis([0 10])
thisLabel = get(h,'YTickLabel');
thisLabel{end} = ['> ',thisLabel{end}];
set(h,'YTickLabel',thisLabel)
ylabel(h, 'Counts','FontSize',12)

n_audIdx = sum(audIdx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1e

function local_displayStDistribution(dataArray,audIdx,negIdx,posIdx)
% extract already calculated percentage effect changes between SSFO on and
% control condition for single Tones and DRC for each unit
eff = cell2mat(dataArray(:,10));

% make already subroups: 
signR = eff(negIdx & audIdx);
signE = eff(posIdx & audIdx);
others = ~negIdx & ~posIdx & audIdx;
nsign = eff(others);

histo_edges = [-100:5:100,max(max(eff(audIdx)))];

% barplot
histo_R = histcounts(signR,histo_edges);   
histo_E = histcounts(signE,histo_edges);
histo_N = histcounts(nsign,histo_edges);

width = 14; 
height = 10;
figure('units', 'centimeters', 'position', [5 10 width height]);
stacked = [histo_R;histo_E;histo_N]';
H = bar(histo_edges(1:end-1),stacked,'stacked');
% change colors
col = [0.6 0.4 0.8;0.2 0.6 0;0.8 0.8 0.8];
colorSet = [];
for i = 1:3
    myColors = col(i,:);
    colorSet = [colorSet myColors];
    H(i).FaceColor = 'flat';
    H(i).CData = myColors;
  %  if i == 3
  %      H(i).EdgeColor = 'none';
  %  end
end
legend({'sign.reduced','sign.enhanced','non-sign.'})
legend boxoff
box off
set(gca,'YLim',[0 120])
set(gca,'FontSize',14)
ylabel('Counts','FontSize',12)
xlabel({'Firing rate change';'after activation [%]'},'FontSize',12)

legend boxoff

end

function [m,s,r,pcoff,p,n] = local_information(dataArray,audIdx,negIdx,posIdx)
% fig D
% mean and standard deviation for percentage effect changes between SSFO on and
% control condition of both stimulus paradigms 
eff = cell2mat(dataArray(:,10:11));
% discard nonsign
%signIdx = negIdx & audIdx | posIdx & audIdx;
eff = eff(audIdx,:);

m = mean(eff);
s = std(eff);
% correlation coefficient and Wilcoxon Ranksum test of effect change between both stimulus paradigms
[r,pcoff] = corrcoef(eff(:,1),eff(:,2));
pcoff = pcoff(2,1);
r = r(2,1);
[p,~,~] = signrank(eff(:,1),eff(:,2));

% fig E
% number of auditory units with a significant decrease and increase of
% firing rate after PV action and with no significant change in single
% tones stimulus paradigm
n = NaN(1,3); % decrease, increase, no change
n(1) = sum(negIdx & audIdx);
n(2) = sum(posIdx & audIdx);
n(3) = sum(~negIdx & ~posIdx & audIdx);

end