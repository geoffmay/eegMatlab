resampleFolder = '/Users/Geoff/Documents/reliability5';
valueFolder = '/Users/Geoff/Documents/avgValues';
surface = true;
copyOnly = false;
doPlot = false;

files = dir(resampleFolder);
files([files.isdir]) = [];
valueFiles = dir(valueFolder);
valueFiles ([valueFiles .isdir]) = [];

fileCounter = 1;
completeFileCounter = 1;
clear summaries;
mat = NaN(50, 2480, length(files));
while(fileCounter <= length(files))
    if(copyOnly)
        fprintf('\n file %d of %d, measure %d of %d', fileCounter, length(files), measureCounter, length(data.surfLabels));
    end
    valueIndex = find(strcmp({valueFiles.name}, files(fileCounter).name));
    if(length(valueIndex) > 0)
        filename = files(fileCounter).name;
        data = load(fullfile(resampleFolder, filename));
        if(isfield(data, 'summary'))
            data = data.summary;
            values = load(fullfile(valueFolder, valueFiles(valueIndex).name));
            values = values.summary;
            cohInd = cellfun(@length, strfind(data.surfLabels, 'coh')) > 0;
            powInd = cellfun(@length, strfind(data.surfLabels, 'abs')) > 0;
            cohAvgValue = mean(values.surfSummary.meanValue(cohInd));
            powAvgValue = mean(values.surfSummary.meanValue(powInd));
            clear plot95Raw;
            clear plot95Norm;
            
            measureCounter = 1;
            while(measureCounter <= length(data.surfLabels))
                checkThis = true;
                if(~copyOnly)
                    fprintf('\n file %d of %d, measure %d of %d', fileCounter, length(files), measureCounter, length(data.surfLabels));
                end
                if(surface)
                    thisData = data.surfResample;
                    if(isfield(data, 'surfPermResample'))
                        thisDataNoise = data.surfNoiseResample;
                        thisDataPerm = data.surfPermResample;
                    end
                    thisValue = values.surfSummary;
                    %label = thisValue.labels{measureCounter};
                    valueIndex = measureCounter;
                    label = data.surfLabels{measureCounter};
                    %                     %if(~strcmp(data.surfLabels{measureCounter}, label))
                    %                     if(~strcmp(values.surfSummary.labels{measureCounter}, label))
                    %                         %dataIndex = find(cellfun(@length, strfind(data.surfLabels, label)));
                    valueIndex = find(cellfun(@length, strfind(thisValue.labels, label)));
                    if(length(valueIndex) == 0)
                        checkThis = false;
                    end
                    %                     end
                else
                    thisData = data.icaResample;
                    thisDataNoise = data.icaNoiseResample;
                    thisDataPerm = data.icaPermResample;
                    thisValue = values.icaSummary;
                    label = data.icaLabels{measureCounter};
                    %                     valueIndex = measureCounter;
                    %                     if(~strcmp(thisValue.labels{measureCounter}, label))
                    valueIndex = find(cellfun(@length, strfind(thisValue.labels, label)));
                    if(length(valueIndex) == 0)
                        checkThis = false;
                    end
                    %                     end
                end
                if(checkThis)
                    if(copyOnly)
                        for sampleCounter = 1:length(thisData)
                            row = thisData(sampleCounter).averageDifference(measureCounter);
                            rowPerm = thisDataPerm(sampleCounter).averageDifference(measureCounter);
                            rowNoise = thisDataNoise(sampleCounter).averageDifference(measureCounter);
                            mat(sampleCounter, measureCounter, fileCounter) = row;
                            matPerm = thisDataPerm(sampleCounter).averageDifference(measureCounter);
                            matNoise = thisDataNoise(sampleCounter).averageDifference(measureCounter);
                        end
                    else
                        secondsDuration = NaN(1, length(thisData));
                        avgDif = NaN(1, length(thisData));
                        avgDifCheck = NaN(1, length(thisData));
                        stdDevDif = NaN(1, length(thisData));
                        finalDelDel = NaN(1, length(thisData));
                        for i = 1:length(thisData)
                            secondsDuration(i) = thisData(i).frameCount / 128;
                            avgDif(i) = thisData(i).averageDifference(measureCounter);
                            if(exist('thisDataPerm', 'var'))
                                avgDifPerm(i) = thisDataPerm(i).averageDifference(measureCounter);
                                avgDifNoise(i) = thisDataNoise(i).averageDifference(measureCounter);
                                stdDevDifPerm(i) = thisDataPerm(i).stddevDifference(measureCounter);
                                stdDevDifNoise(i) = thisDataNoise(i).stddevDifference(measureCounter);
                            end
                            avgDifCheck(i) = thisData(i).averageDifferenceCheck(measureCounter);
                            stdDevDif(i) = thisData(i).stddevDifference(measureCounter);
                            finalDelDel(i) = thisData(i).finalDelDel(measureCounter);
                        end
                        ninetyFifthZ = norminv(.95);
                        conf = stdDevDif.*ninetyFifthZ + avgDif;
                        if(exist('thisDataPerm', 'var'))
                            confPerm = stdDevDifPerm.*ninetyFifthZ + avgDifPerm;
                            confNoise = stdDevDifNoise.*ninetyFifthZ + avgDifNoise;
                        end
                        if(strfind(data.surfLabels{measureCounter}, 'coh'))
                            avgValue = cohAvgValue;
                        else
                            avgValue = powAvgValue;
                        end
                        %                     normalized = conf ./ thisValue.meanValue(measureCounter);
                        normalized = conf ./ avgValue;
                        
                        %         [F1, G] = fit(secondsDuration', normalized', 'exp2');
                        %         coeffs = coeffvalues(F1);
                        %                     smoothAvg = slidingAverage(avgDif, 5)
                        %                     smoothConf = slidingAverage(conf, 5)
                        %                     cohDiff = diff(smoothConf)
                        plot95Raw{measureCounter} = conf;
                        plot95Norm{measureCounter} = normalized;
                        
                        if(doPlot)
                            plotTitle = 'coh Fp1-Fp2 5Hz-8Hz';
                            smoothConf = slidingAverage(conf, 5);
                            plotInd = find(strcmp(data.surfLabels, plotTitle));
                            if(measureCounter == plotInd)
                                fontSize = 20;
                                close all;
                                x = secondsDuration'
                                y1 = avgDif'.*100;
                                y2 = conf'.*100;
                                if(exist('thisDataPerm', 'var'))
                                    y1a = avgDifPerm'.*100;
                                    y1b = avgDifNoise'.*100;
                                    y2a = confPerm'.*100;
                                    y2b = confNoise'.*100;
                                end
                                y3 = stdDevDif.*ninetyFifthZ + avgDif
                                meanFit = fit(x,y1,'exp2');
                                confFit = fit(x,y2,'exp2');
                                figure;
                                hold on;
                                scatter(x, y1);
                                scatter(x, y2);
                                if(exist('y1a','var'))
                                    scatter(x, y1a);
                                    scatter(x, y1b);
                                end
                                
                                legend('mean', 'average', 'permuted', 'noise');
                                plot(meanFit);
                                plot(confFit);
                                %                             plot(x, [y1, y2].*100, 'linewidth', 5);
                                legend('average', '95 %ile');
                                ylabel('difference between sample means (% coherence)', 'fontsize', fontSize)
                                xlabel('sample duration (seconds)', 'fontsize', fontSize);
                                title('Difference between subsamples of theta coherence at Fp1-Fp2', 'fontsize', fontSize);
                                set(gca, 'fontsize', fontSize);
                            end
                            
                            plotTitle = 'abs Fp1 5Hz-8Hz';
                            plotInd = find(strcmp(data.surfLabels, plotTitle));
                            if(measureCounter == plotInd)
                                fontSize = 20;
                                close all;
                                x = secondsDuration'
                                y1 = avgDif';
                                y2 = conf';
                                if(exist('thisDataPerm', 'var'))
                                    y1a = avgDifPerm';
                                    y1b = avgDifNoise';
                                end
                                y3 = stdDevDif.*ninetyFifthZ + avgDif
                                figure;
                                if(exist('thisDataPerm', 'var'))
                                    plot(x, [y1, y2, y1a, y1b], 'linewidth', 5);
                                else
                                    plot(x, [y1, y2], 'linewidth', 5);
                                end
                                legend('average', '95 %ile', 'permuted', 'noise');
                                ylabel('difference between sample means (microvolts^2)', 'fontsize', fontSize)
                                xlabel('sample duration (seconds)', 'fontsize', fontSize);
                                title('Difference between subsamples of theta power at Fp1', 'fontsize', fontSize);
                                set(gca, 'fontsize', fontSize);
                            end
                        end
                        normAvg = avgDif ./ thisValue.meanValue(valueIndex)';
                        neighbors = 2;
                        
                        avgAround5Min = avgDif(20-neighbors:20+neighbors);
                        poly = polyfit(1:length(avgAround5Min), avgAround5Min,1);
                        meanSummary.slopeAt5Minutes = poly(1) * 4;
                        meanSummary.valueAt5Minutes = mean(avgAround5Min);
                        normAvgAround5Min = normAvg(20-neighbors:20+neighbors);
                        poly = polyfit(1:length(normAvgAround5Min), normAvgAround5Min,1);
                        meanSummary.normSlopeAt5Minutes = poly(1) * 4;
                        meanSummary.normValueAt5Minutes = mean(normAvgAround5Min);
                        
                        avgAround10Min = avgDif(40-neighbors:40+neighbors);
                        poly = polyfit(1:length(avgAround10Min), avgAround10Min,1);
                        meanSummary.slopeAt10Minutes = poly(1) * 4;
                        meanSummary.valueAt10Minutes = mean(avgAround10Min);
                        normAvgAround10Min = normAvg(40-neighbors:40+neighbors);
                        poly = polyfit(1:length(normAvgAround10Min), normAvgAround10Min,1);
                        meanSummary.normSlopeAt10Minutes = poly(1) * 4;
                        meanSummary.normValueAt10Minutes = mean(normAvgAround10Min);
                        
                        avgAround5Min = conf(20-neighbors:20+neighbors);
                        poly = polyfit(1:length(avgAround5Min), avgAround5Min,1);
                        conf95Summary.slopeAt5Minutes = poly(1) * 4;
                        conf95Summary.valueAt5Minutes = mean(avgAround5Min);
                        normAvgAround5Min = normAvg(20-neighbors:20+neighbors);
                        poly = polyfit(1:length(normAvgAround5Min), normAvgAround5Min,1);
                        conf95Summary.normSlopeAt5Minutes = poly(1) * 4;
                        conf95Summary.normValueAt5Minutes = mean(normAvgAround5Min);
                        
                        avgAround10Min = normalized(40-neighbors:40+neighbors);
                        poly = polyfit(1:length(avgAround10Min), avgAround10Min,1);
                        conf95Summary.slopeAt10Minutes = poly(1) * 4;
                        conf95Summary.valueAt10Minutes = mean(avgAround10Min);
                        normAvgAround10Min = normAvg(40-neighbors:40+neighbors);
                        poly = polyfit(1:length(normAvgAround10Min), normAvgAround10Min,1);
                        conf95Summary.normSlopeAt10Minutes = poly(1) * 4;
                        conf95Summary.normValueAt10Minutes = mean(normAvgAround10Min);
                        
                        summary.conf95Summary = conf95Summary;
                        summary.avgSummary = meanSummary;
                        summary.label = label;
                        summaries(measureCounter) = summary;
                        completeFileCounter = completeFileCounter + 1;
                        clear summary
                    end
                end
                measureCounter = measureCounter + 1;
            end
        else
            fprintf('\nfile not complete, skipping...');
        end
    else
        fprintf('\nfile not complete, skipping...');
    end
    if(exist('summaries', 'var'))
        digestFolder = fullfile(resampleFolder, 'digest 3');
        if(~exist(digestFolder, 'dir'))
            mkdir(digestFolder);
        end
        save(fullfile(digestFolder, filename), 'summaries', 'plot95Raw', 'plot95Norm', 'filename');
    end
    fileCounter = fileCounter + 1;
    clear summaries;
    clear plot95Raw;
    clear plot95Norm;
end


if(doPlot)
    for i = 1:size(summaries,1)
        for j = 1:size(summaries,2)
            data = summaries(i,j).conf95Summary;
            value5(i,j) = data.valueAt5Minutes;
            value10(i,j) = data.valueAt10Minutes;
            nValue5(i,j) = data.normValueAt5Minutes;
            nValue10(i,j) = data.normValueAt10Minutes;
            slope5(i,j) = data.slopeAt5Minutes;
            slope10(i,j) = data.slopeAt10Minutes;
            nSlope5(i,j) = data.normSlopeAt5Minutes;
            nSlope10(i,j) = data.normSlopeAt10Minutes;
            
            
        end
    end
    
    close all;
    figure;
    imagesc([value5, zeros(2480, 5), value10])
    colorbar;
    title('values');
    figure;
    imagesc([nValue5, zeros(2480, 5), nValue10])
    title('norm values');
    colorbar;
    
    figure;
    imagesc([slope5, zeros(2480, 5), slope10])
    colorbar;
    title('slopes');
    figure;
    imagesc([nSlope5, zeros(2480, 5), nSlope10])
    title('norm slopes');
    colorbar;
    tilefigs
    
    
    meanN10 = mean(nValue10, 2);
    meanN5 = mean(nValue5, 2);
    labels = {summaries(1:2480).label};
    %plotChannelPairs({summaries(1:2480).label}, nValue5(:,1))
    plotCoherencePca(meanN5, labels)
    tilefigs;
    fprintf('\n');
    dummy = 1;
end