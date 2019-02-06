resampleFolder = '/home/data/EEG/processed/Oregon/reliability5';
valueFolder = '/home/data/EEG/processed/Oregon/avgValues';
surface = false;

files = dir(resampleFolder);
files([files.isdir]) = [];
valueFiles = dir(valueFolder);
valueFiles ([valueFiles .isdir]) = [];
doPlot = true;

fileCounter = 1;
completeFileCounter = 1;
while(fileCounter <= length(files))
    valueIndex = find(strcmp({valueFiles.name}, files(fileCounter).name));
    if(length(valueIndex) > 0)
        data = load(fullfile(resampleFolder, files(fileCounter).name));
        if(isfield(data, 'summary'))
        values = load(fullfile(valueFolder, valueFiles(valueIndex).name));
        values = values.summary;
        data = data.summary;
        if(doPlot)
          cohTitle = 'coh Fp1-Fp2 5Hz-8Hz';
          cohInd = find(strcmp(data.surfLabels, cohTitle));
        end
        measureCounter = 1;
        while(measureCounter <= length(data.surfLabels))
            fprintf('\n file %d of %d, measure %d of %d', fileCounter, length(files), measureCounter, length(data.surfLabels));
            if(surface)
                thisData = data.surfResample;
                thisValue = values.surfSummary;
                label = data.surfLabels{measureCounter};
                if(~strcmp(thisValue.labels{measureCounter}, label))
                    error('label mismatch');
                end
            else
                thisData = data.icaResample;
                thisValue = values.icaSummary;
                label = data.icaLabels{measureCounter};
                if(~strcmp(thisValue.labels{measureCounter}, label))
                    error('label mismatch');
                end
            end
            secondsDuration = NaN(1, length(thisData));
            avgDif = NaN(1, length(thisData));
            avgDifCheck = NaN(1, length(thisData));
            stdDevDif = NaN(1, length(thisData));
            finalDelDel = NaN(1, length(thisData));
            for i = 1:length(thisData)
                secondsDuration(i) = thisData(i).frameCount / 128;
                avgDif(i) = thisData(i).averageDifference(measureCounter);
                avgDifCheck(i) = thisData(i).averageDifferenceCheck(measureCounter);
                stdDevDif(i) = thisData(i).stddevDifference(measureCounter);
                finalDelDel(i) = thisData(i).finalDelDel(measureCounter);
            end
            ninetyFifthZ = norminv(.95);
            conf = stdDevDif.*ninetyFifthZ + avgDif;
            normalized = conf ./ thisValue.meanValue(measureCounter);
            %         [F1, G] = fit(secondsDuration', normalized', 'exp2');
            %         coeffs = coeffvalues(F1);
            if(doPlot)
              if(measureCounter == cohInd)
                close all;
                figure;
                x = secondsDuration';
                y1 = avgDif' .* 100;
                y2 = conf' .* 100;
                plot(x,[y1, y2], 'linewidth', 5);
                set(gca, 'fontsize', 20);
                legend('average', '95th %ile');         
                xlabel('duration of sample (seconds)');
                ylabel('difference between sample means (%)');
                %plot(F1, secondsDuration', [normalized', normAvg']);
                %title(sprintf('%s difference', cohTitle));
                title('Agreement between subsamples of theta coherence between Fp1-Fp2');
              end
            end
            normAvg = avgDif ./ thisValue.meanValue(measureCounter)';
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
            
            summaries(measureCounter, completeFileCounter) = summary;
            completeFileCounter = completeFileCounter + 1;
            clear summary
            
            measureCounter = measureCounter + 1;
        end
    else
        fprintf('\nfile not complete, skipping...');
    end
    else
        fprintf('\nfile not complete, skipping...');
    end
    fileCounter = fileCounter + 1;
end



for i = 1:size(summaries,1)
    for j = 1:size(summaries,2)
        data = summaries(i,j).conf95Summary;
        value5(i,j) = data.valueAt5Minutes;
        value10(i,j) = data.valueAt10Minutes;
        nValue5(i,j) = data.normValueAt5Minutes;
        nValue10(i,j) = data.normValueAt10Minutes;
        
        
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

