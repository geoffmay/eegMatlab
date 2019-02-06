inputFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\Oregon\artifactRemovedAnalytics';
oldFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\Oregon\artifactRemovedReliability';

inputFile = 'PM101Surface.mat';


if(exist('surfCoh', 'var'))
    keepVar = 'surfCoh';
    allVar = who;
    allVar{end+1} = 'allVar';
    allVar(strcmp(allVar, keepVar)) = [];
    clear(allVar{:});
else
    
    load(fullfile(inputFolder, inputFile));
end
doPlot = true;

task = guessWahbehTaskFromAlpha(surfCoh);

% if(doPlot)
%     plot(simp);
% end

% measureCounter = 1;
% for targetMeasure = targetMeasures
%     plotLabel{measureCounter} = sprintf('eyes closed %s', targetMeasure{1});
%     measureCounter = measureCounter + 1;
%     plotLabel{measureCounter} = sprintf('eyes open %s', targetMeasure{1});
%     measureCounter = measureCounter + 1;
%     plotLabel{measureCounter} = sprintf('combined %s', targetMeasure{1});
%     measureCounter = measureCounter + 1;
%     
% end


if(fewMeasures)
    targetMeasures = {'abs Oz 9Hz-12Hz', 'abs Fp1 5Hz-8Hz'};
    measureCounter = 1;
    for targetMeasure = targetMeasures
        
        ind = find(strcmp(surfCoh.labels, targetMeasure));
        measureData = surfCoh.matrix(:, ind);
        eyesClosed = measureData(task == 1);
        eyesOpen = measureData(task == -1);
        %     if(doPlot)
        %         figure;
        %         hold on;
        %
        %         plot((1:length(eyesClosed))./128, eyesClosed);
        %         plot((1:length(eyesOpen))./128, eyesOpen);
        %     end
        
        %compute reliability for eyes open and eyes closed states
        minDeltaFraction = 1e-4;
        i = 1;
        frameStep = 5;
        frameCount = 128*frameStep*i;
        while(frameCount < (length(eyesOpen)*.4) || frameCount < (length(eyesClosed)*.4))
            reliabilities.frameCounts(i) = frameCount;
            if( frameCount < (length(eyesClosed)*.4))
                reliabilities.eyesClosed(i) = asymptoteSignalReliability(eyesClosed, frameCount);
                fprintf(' (eyes closed)')
            end
            if( frameCount < (length(eyesOpen)*.4))
                reliabilities.eyesOpen(i) = asymptoteSignalReliability(eyesOpen, frameCount);
                fprintf(' (eyes open)')
            end
            i = i + 1;
            frameCount = 128*5*i;
        end
        
        %load old combined reliability data
        oldData = load(fullfile(oldFolder, inputFile));
        oldInd = find(strcmp(oldData.summary.surfLabels, targetMeasure));
        stdCoeff = norminv(.95);
        for i = 1:length(oldData.summary.surfResample)
            a = oldData.summary.surfResample(i);
            oldPlot(i) = a.averageDifference(oldInd) + stdCoeff * a.stddevDifference(oldInd);
            xOld(i) = a.frameCount / 128;
        end
        
        a = [reliabilities.eyesOpen];
        a = [a.percentiles];
        for i = 1:length(a)
            eyesOpen95(i) = a(i).percentileValues(3);
        end
        xOpen = (1:length(eyesOpen95)) .* frameStep;
        
        a = [reliabilities.eyesClosed];
        a = [a.percentiles];
        for i = 1:length(a)
            eyesClosed95(i) = a(i).percentileValues(3);
        end
        xClosed = (1:length(eyesClosed95)) .* frameStep;
        
        plotX{measureCounter} = xClosed;
        plotY{measureCounter} = eyesClosed95;
        plotLabel{measureCounter} = sprintf('eyes closed %s', targetMeasure{1});
        measureCounter = measureCounter + 1;
        
        plotX{measureCounter} = xOpen;
        plotY{measureCounter} = eyesOpen95;
        plotLabel{measureCounter} = sprintf('eyes open %s', targetMeasure{1});
        measureCounter = measureCounter + 1;
        
        plotX{measureCounter} = xOld;
        plotY{measureCounter} = oldPlot;
        plotLabel{measureCounter} = sprintf('combined %s', targetMeasure{1});
        measureCounter = measureCounter + 1;
        
    end
else
%     ind = find(strcmp(surfCoh.labels, targetMeasure));
    eyesClosed = surfCoh.matrix(task == 1, :);
    eyesOpen = surfCoh.matrix(task == -1, :);
    
    %compute reliability for eyes open and eyes closed states
    minDeltaFraction = 1e-4;
    i = 1;
    frameStep = 5;
    frameCount = 128*frameStep*i;
    while(frameCount < (length(eyesOpen)*.4) || frameCount < (length(eyesClosed)*.4))
        reliabilities.frameCounts(i) = frameCount;
        if( frameCount < (length(eyesClosed)*.4))
            reliabilities.eyesClosed(i) = asymptoteSignalReliability(eyesClosed, frameCount);
            fprintf(' (eyes closed)')
        end
        if( frameCount < (length(eyesOpen)*.4))
            reliabilities.eyesOpen(i) = asymptoteSignalReliability(eyesOpen, frameCount);
            fprintf(' (eyes open)')
        end
        i = i + 1;
        frameCount = 128*5*i;
    end
    
    %load old combined reliability data
    oldData = load(fullfile(oldFolder, inputFile));
    oldInd = find(strcmp(oldData.summary.surfLabels, targetMeasure));
    stdCoeff = norminv(.95);
    for i = 1:length(oldData.summary.surfResample)
        a = oldData.summary.surfResample(i);
        oldPlot(i) = a.averageDifference(oldInd) + stdCoeff * a.stddevDifference(oldInd);
        xOld(i) = a.frameCount / 128;
    end
    
    a = [reliabilities.eyesOpen];
    a = [a.percentiles];
    for i = 1:length(a)
        eyesOpen95(i) = a(i).percentileValues(3);
    end
    xOpen = (1:length(eyesOpen95)) .* frameStep;
    
    a = [reliabilities.eyesClosed];
    a = [a.percentiles];
    for i = 1:length(a)
        eyesClosed95(i) = a(i).percentileValues(3);
    end
    xClosed = (1:length(eyesClosed95)) .* frameStep;
end

%plot reliabilities
if(doPlot)
    figure;
    hold on;
    for i = 1:length(plotX)
        plot(plotX{i}, plotY{i});
    end
    legend(plotLabel);
    xlabel('subsample duration (seconds)');
    ylabel('power difference (log(uv^2))');
    %     plot(xClosed, eyesClosed95);
    %     plot(xOpen, eyesOpen95);
    %     plot(xOld, oldPlot);
    %     legend({'alpha closed','alpha open','alpha combined'});
end


save('C:\Users\Neuro\Documents\MATLAB\processed\Oregon\ptsdTaskSliceDemo.mat');
% 
% monotonic = 1:1000;
% stepfunction = zeros(1, 1000);
% stepfunction(501:end) = 1;
% squarewave = zeros(1, 1000);
% squareCount = 10;
% for i = 1:length(squareCount)
%     squareStart = 1 + (i-1) * (length(squarewave) / squareCount);
%     squareEnd = squareStart + (length(squarewave) / squareCount / 2);
%     squarewave(squareStart:squareEnd) = 1;
% end
% i = 1;
% frameCount = 5*i;
% while(frameCount < (length(monotonic)*.5))
%     dummyReliability.frameCounts(i) = frameCount;
%     dummyReliability.monotonic(i) = asymptoteSignalReliability(squarewave, frameCount);
%     fprintf(' (monotonic)')
%     i = i + 1;
%     frameCount = 5*i;
% end
% 
% a=[dummyReliability.monotonic.percentiles]
% for i = 1:length(a)
%     b = a(i).percentileValues;
%     toPlot(i) = b(4);
% end
% close all;
% figure;
% plot(toPlot);