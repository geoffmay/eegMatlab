load('/Users/Geoff/Downloads/zScoreAnalysis.mat')
plotImages = false;
plotElectrodes = true;

verum = analysis.verum.meanImprove(analysis.verum.keep);
sham = analysis.sham.meanImprove(analysis.sham.keep);
difference = analysis.meanImprovementDiff;
labels = analysis.labels;

for i = 1:length(labels)
    label = labels{i};
    dashes = strfind(label, '-');
    if(length(dashes) == 2)
        newSpaces = dashes;
    elseif(length(dashes) == 3)
        newSpaces = dashes([1, 3]);
    else
        newSpaces = [];
    end
    if(length(newSpaces > 1))
        label(newSpaces) = ' ';
        labels{i} = label;
    end
end

%plotCoherencePca(difference, labels);
%tilefigs;

if(plotImages | plotElectrodes)
    imagesc(analysis.verum.improvementFactors(:, analysis.verum.keep));
    imagesc(analysis.verum.improvementFactors(:, analysis.verum.keep), [0 1]);

    allFiles = [analysis.verum.filenames, analysis.sham.filenames];
    allImprovement = [analysis.verum.improvementFactors; analysis.sham.improvementFactors];
    isVerum = [ones(1,length(analysis.verum.filenames)), zeros(1,length(analysis.sham.filenames))];
    close all;
    for i = 1:25
        searchString = sprintf('_%03d-', i)
        %        ptInd = find(cellfun(@length, strfind(analysis.verum.filenames, searchString)));
        ptInd = find(cellfun(@length, strfind(allFiles, searchString)));
        if(length(ptInd) > 0)
            isThisVerum = sum(isVerum(ptInd)) / length(ptInd);
            fprintf('\n%s: %d', searchString);
            for j = 1:length(ptInd)
                fprintf(' %d', ptInd(j));
            end
            if(plotImages)
                figure;
                %im = imagesc(analysis.verum.improvementFactors(ptInd, analysis.verum.keep));
                im = imagesc(allImprovement(ptInd, analysis.verum.keep), [-1 1]);
                colormap('parula');
                title(sprintf('ROBI_%s (%f)', searchString, isThisVerum));
                %             colorbar;
            end
            if(plotElectrodes)
                close all;
                toPlot = mean(allImprovement(ptInd, analysis.verum.keep));
                plotCoherencePca(toPlot, labels);
                tilefigs;
                fprintf('\n%s', searchString);
            end
        end
    end
    tilefigs;
end
