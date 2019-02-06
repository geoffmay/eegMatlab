onlyPresent = 0;

%artifactRemovedFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\Oregon\artifactRemovedAnalytics';
artifactRemovedFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\Oregon\artifactRemovedReliability';
artifactPresentFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\Oregon\artifactPresentReliability';
artifactPresentValueFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\Oregon\reliability9';
outputFilename = 'pairs';
artifactRemovedFiles = dir(artifactRemovedFolder);
artifactRemovedFiles([artifactRemovedFiles.isdir]) = [];
artifactPresentFiles = dir(artifactPresentFolder);
artifactPresentValueFiles = dir(artifactPresentValueFolder);
artifactPresentValueFiles([artifactPresentValueFiles.isdir]) = [];
artifactPresentFiles([artifactPresentFiles.isdir]) = [];
stdDev95 = norminv(0.95);

cohFitOptions = fitoptions('exp2');
cohFitOptions.StartPoint = [0.0931   -0.0266    0.0948   -0.0009];
cohFitOptions.Lower = [0   -1    0   -1];
cohFitOptions.MaxIter = 2000;
cohFitOptions.Upper = [3   0    3   0];
powFitOptions = cohFitOptions;
powFitOptions.Upper = [10   0    10   0];
asyFitOptions = powFitOptions;

if(~onlyPresent)
    for fileCounter = 1:length(artifactRemovedFiles);
        removedFilename = artifactRemovedFiles(fileCounter).name;
        removedReliability = loadReliabilityMatrix(fullfile(artifactRemovedFolder, removedFilename));
        presentReliability = loadReliabilityMatrix(fullfile(artifactPresentFolder, removedFilename), removedReliability.measureLabels);
        isCoh = cellfun(@length, strfind(removedReliability.measureLabels, 'coh')) > 0;
        isAsy = cellfun(@length, strfind(removedReliability.measureLabels, 'asy')) > 0;
        
        for measureIndex = 1:size(removedReliability.averageDifference, 1)
            fprintf('\nfile %d of %d, measure %d of %d', fileCounter, length(artifactRemovedFiles), measureIndex, size(removedReliability.averageDifference, 1));
            x = removedReliability.sampleMinuteDurations';
            y = [presentReliability.percentile95Difference(measureIndex,:); removedReliability.percentile95Difference(measureIndex,:)]';
            if(isCoh(measureIndex))
                fitOptions = cohFitOptions;
            elseif(isAsy(measureIndex))
                fitOptions = asyFitOptions;
            else
                fitOptions = powFitOptions;
            end
            [iteration.presentCurve, iteration.presentGoodness] = fit(x, y(:,1), 'exp2', fitOptions);
            [iteration.removedCurve, iteration.removedGoodness] = fit(x, y(:,2), 'exp2', fitOptions);
            iteration.measureLabel = removedReliability.measureLabels{measureIndex};
            iteration.sourceFile = removedFilename;
            if(fileCounter == 1 && measureIndex == 1)
                fits = repmat(iteration, [length(artifactRemovedFiles), size(removedReliability.averageDifference, 1)]);
            end
            fits(fileCounter, measureIndex) = iteration;
        end
    end
    outputFilename = 'C:\Users\Neuro\Documents\MATLAB\processed\Oregon\fits\pairs.mat';
    save(outputFilename, 'fits', '-v7.3');
    fprintf('\n');
else
    doPlot= false;
    if(doPlot)
        close all;
        figure;
        x = removedReliability.sampleMinuteDurations;
        measureLabel = 'abs Pz 1Hz-4Hz';
        measureIndex = find(strcmp(removedReliability.measureLabels, measureLabel));
        y = [presentReliability.percentile95Difference(measureIndex,:); removedReliability.percentile95Difference(measureIndex,:)];
        plot(x,y);
        title(sprintf('95%% %s',removedReliability.measureLabels{measureIndex}));
        legend({'present', 'absent'});
    end
    
    
    
    for fileCounter = 1:length(artifactPresentFiles);
        presentFilename = artifactRemovedFiles(fileCounter).name;
        %     removedFilename = artifactRemovedFiles(fileCounter).name;
        %     removedReliability = loadReliabilityMatrix(fullfile(artifactRemovedFolder, removedFilename));
        presentReliability = loadReliabilityMatrix(fullfile(artifactPresentFolder, presentFilename));
        isCoh = cellfun(@length, strfind(presentReliability.measureLabels, 'coh')) > 0;
        isAsy = cellfun(@length, strfind(presentReliability.measureLabels, 'asy')) > 0;
        for measureIndex = 1:size(presentReliability.averageDifference, 1)
            fprintf('\nfile %d of %d, measure %d of %d', fileCounter, length(artifactRemovedFiles), measureIndex, size(presentReliability.averageDifference, 1));
            x = presentReliability.sampleMinuteDurations';
            %         y = [presentReliability.percentile95Difference(measureIndex,:); presentReliability.percentile95Difference(measureIndex,:)]';
            y = presentReliability.percentile95Difference(measureIndex,:)';
            if(~any(isnan(y)))
                if(isCoh(measureIndex))
                    fitOptions = cohFitOptions;
                elseif(isAsy(measureIndex))
                    fitOptions = asyFitOptions;
                else
                    fitOptions = powFitOptions;
                end
                [iteration.presentCurve, iteration.presentGoodness] = fit(x, y(:,1), 'exp2', fitOptions);
                %         [iteration.presentCurve, iteration.presentGoodness] = fit(x, y(:,2), 'exp2', fitOptions);
                iteration.measureLabel = presentReliability.measureLabels{measureIndex};
                iteration.sourceFile = presentFilename;
                if(fileCounter == 1 && measureIndex == 1)
                    fits = repmat(iteration, [length(artifactPresentFiles), size(presentReliability.averageDifference, 1)]);
                    remove = zeros(size(fits));
                end
                fits(fileCounter, measureIndex) = iteration;
            else
                remove(fileCounter, measureIndex) = 1;
            end
        end
    end
    outputFilename = 'C:\Users\Neuro\Documents\MATLAB\processed\Oregon\fits\artifactPresent.mat';
    save(outputFilename, 'fits', '-v7.3');
end

