clear;
showGraphs = true;
editTable = false;

channelCount = 32;
prePhase = zeros(channelCount,channelCount);
postPhase = zeros(channelCount,channelCount);
delPhase = zeros(channelCount,channelCount);

load('/Users/Geoff/Documents/MATLAB/EEG/phaseSlopeInfo.mat');
delSlope = table;
i = 1;
j = 2;
for pairCounter = 1:length(summaries(1).channelPairs)
    
    preSlope(pairCounter) = summaries(1).channelPairs(pairCounter).slope;
    postSlope(pairCounter) = summaries(2).channelPairs(pairCounter).slope;
    label = sprintf('%s-%s',summaries(1).channelPairs(pairCounter).channel1,summaries(1).channelPairs(pairCounter).channel2);
    if(editTable)
        delSlope{pairCounter,1} = postSlope(pairCounter) - preSlope(pairCounter);
        delSlope{pairCounter,2} = {label};
        delSlope{pairCounter,3} = preSlope(pairCounter);
        delSlope{pairCounter,4} = postSlope(pairCounter);
    end
    prePhase(i,j) = preSlope(pairCounter);
    postPhase(i,j) = postSlope(pairCounter);
    delPhase(i,j) = postSlope(pairCounter) - preSlope(pairCounter);
    prePhase(j,i) = preSlope(pairCounter);
    postPhase(j,i) = postSlope(pairCounter);
    delPhase(j,i) = postSlope(pairCounter) - preSlope(pairCounter);
    if(showGraphs)
        %debug
        close all;
        figure;
        hold on;
        plot(summaries(1).channelPairs(pairCounter).phaseAngle);
        x = 20:100;
        y = x .* summaries(1).channelPairs(pairCounter).poly(1) + summaries(1).channelPairs(pairCounter).poly(2);
        plot(x, y, 'r');
        
        title(sprintf('%s pre', label));
        %end debug
    end
    j = j + 1;
    if(j > channelCount)
        i = i + 1;
        j = i + 1;        
    end
end
if(editTable)
    delSlope.Properties.VariableNames{1} = 'changeInSlope';
    delSlope.Properties.VariableNames{2} = 'channelPair';
    delSlope.Properties.VariableNames{3} = 'preSlope';
    delSlope.Properties.VariableNames{4} = 'postSlope';
    sortedTable = sortrows(delSlope,1);
end
preMin = min(min(prePhase));
postMin = min(min(postPhase));
delMin = min(min(delPhase));
preMax = max(max(prePhase));
postMax = max(max(postPhase));
delMax = max(max(delPhase));
minmin = min([preMin,postMin,delMin]);
maxmax = max([preMax, postMax, delMax]);
close all;

[labels] = antChannelLocs;

figure;
imagesc(prePhase);
title('phase slope index, pre-neurofeedback');
cbar = colorbar();
set(gca, 'CLim', [minmin, maxmax]);
set(gca, 'xtick', 1:32);
set(gca, 'xticklabel', labels(1:32));
set(gca, 'ytick', 1:32);
set(gca, 'yticklabel', labels(1:32));

figure;
imagesc(postPhase);
title('phase slope index, post-neurofeedback');
cbar = colorbar();
set(gca, 'CLim', [minmin, maxmax]);
set(gca, 'xtick', 1:32);
set(gca, 'xticklabel', labels(1:32));
set(gca, 'ytick', 1:32);
set(gca, 'yticklabel', labels(1:32));

figure;
imagesc(delPhase);
title('change in phase slope index');
cbar = colorbar();
set(gca, 'CLim', [minmin, maxmax]);
set(gca, 'xtick', 1:32);
set(gca, 'xticklabel', labels(1:32));
set(gca, 'ytick', 1:32);
set(gca, 'yticklabel', labels(1:32));

