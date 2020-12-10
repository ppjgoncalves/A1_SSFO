%% code and data for Figure 6
function displayFigure6
% load:
% compTable: now we have the variable compTbl in our workspace


load('compTable.mat')  %#ok<LOAD>

paradigms = {'Tones','DRC','Voc'};
parms = {'slope','interc'};
pName= {'pSlope','pInterc'};

idxPl = reshape(1:6,3,2)';

f=figure;
for iPrm = 1:2
    x = compTbl.([parms{iPrm} paradigms{1}]); % alway compare to tones
    xOutl = isoutlier(x);
    xSig = compTbl.([pName{iPrm} paradigms{1}])<0.05;
    for iPrd = 2:3
        y = compTbl.([parms{iPrm} paradigms{iPrd}]);
        sig = compTbl.([pName{iPrm} paradigms{iPrd}])<0.05 & xSig;
        outl = isoutlier(y) | xOutl;
        
        subplot(2,3,idxPl(iPrm,iPrd-1)) %  put it in the right place
        
        scatter(x(~outl & ~sig),y(~outl & ~sig),4,0.5*[1 1 1]), hold on
        scatter(x(~outl & sig),y(~outl & sig),4,'k','markerfacecolor','k')
        
        axis tight
        local_matchAxis
        axis square
        set(gca,'tickdir','out')
        
        if iPrm == 1
            title(paradigms{iPrd})
        else
            xlabel(paradigms{1})
        end
        if iPrd ==2
            ylabel(parms{iPrm})
        
        end
    end
end


%% histoggrams
outl = any(isoutlier(compTbl(:,[6 8 13])),2);
bins = linspace(0,2,20);
hT = histcounts(compTbl.slopeTones(~outl),bins);
hD = histcounts(compTbl.slopeDRC(~outl),bins);
hV = histcounts(compTbl.slopeVoc(~outl),bins);


subplot(233)
plot(bins(1:end-1),[hT;hD;hV],'LineWidth',2)
xlabel('slope')
ylabel('# of units')
set(gca,'tickdir','out','yaxisloc','right')
box off
axis tight
line([1 1],ylim,'color','k','linestyle',':')

outl = any(isoutlier(compTbl(:,[7 9 14])),2);
bins = linspace(-0.35,0.35,20);
hT = histcounts(compTbl.intercTones(~outl),bins);
hD = histcounts(compTbl.intercDRC(~outl),bins);
hV = histcounts(compTbl.intercVoc(~outl),bins);

l = legend('tones','DRC','Voc');
l.Box = 'off';
l.ItemTokenSize(1) = 5;
l.Position = [0.8120 0.7612 0.0925 0.1610];


subplot(236)
plot(bins(1:end-1),[hT;hD;hV],'LineWidth',2)

set(gca,'tickdir','out','yaxisloc','right')
box off
axis tight
line([0 0],ylim,'color','k','linestyle',':')
xlabel('intercept')
ylabel('# of units')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local function

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

line(xlim,ylim,'color','k')
end