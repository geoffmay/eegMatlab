function task = guessWahbehTaskFromAlpha(surfCoh)
% 
% inputFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\Oregon\artifactRemovedAnalytics';
% 
% inputFile = 'PM101Surface.mat';
% 
% 
% if(exist('surfCoh', 'var'))
%     keepVar = 'surfCoh';
%     allVar = who;
%     allVar{end+1} = 'allVar';
%     allVar(strcmp(allVar, keepVar)) = [];
%     clear(allVar{:});
% else
%     
%     load(fullfile(inputFolder, inputFile));
% end
ind = find(strcmp(surfCoh.labels, 'abs Oz 9Hz-12Hz'));
doPlot = false;

OzAlpha = surfCoh.matrix(:, ind);
filterSize = 2500;
a = 1;
b = repmat(1/filterSize, [1, filterSize]);
filtAlpha1 = filtfilt(b, a, OzAlpha);
alphaMean = mean(filtAlpha1);
alphaStd = std(filtAlpha1);
if(doPlot)
    figure;
    hold on;
%     plot(filtAlpha)
    plot(filtAlpha1)
    plot(OzAlpha)
    zoom xon
end

highInd = find(filtAlpha1 > (alphaMean + 1 .* alphaStd));
lowInd = find(filtAlpha1 < (alphaMean + 0 .* alphaStd));

%signal starts low (at zero)
low = lowInd(1);
high = highInd(1);
counter = 0;
keepGoing = true;

while(keepGoing)
    counter = counter + 1;
    nextLowIndirect = min(find(lowInd > high));
    if(length(nextLowIndirect) > 0)
        nextLow = lowInd(nextLowIndirect);
        falling(counter) = nextLow;
        low = nextLow;
    else
        keepGoing = false;
    end
    nextHighIndirect = min(find(highInd > nextLow));
    nextHigh = highInd(nextHighIndirect);
    rising(counter) = high;
    high = nextHigh;
end
task = zeros(1, length(filtAlpha1));
wiggleRoom = 1000;
for i = 1:length(falling)
    %shift = + filterSize / 2;
    %shift = - filterSize / 2;
    shift = 0;
    start1 = rising(i) + shift;
    start2 = rising(i+1) + shift;
    end1 = falling(i) + shift;
    task(start1+wiggleRoom:end1-wiggleRoom) = 1;
    task(end1+wiggleRoom:start2-wiggleRoom) = -1;
end
task(start2+wiggleRoom:end) = 1;

if(doPlot)
    plot(task);
end