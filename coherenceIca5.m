%based on coherencePca5, this compares components from each of 3 references
%and finds the timecourse that is maximally correlated for each component;
%this is a way of quantifying the biological "truth" of the component;
%components with low correlations with other references are either
%artifacts of the refrence montage or obscured by the other reference
%montage, with higher likelihood assumed to be the former; mapping the
%component in question will likely help sort this out.

clear;

debug = 0;
downSample = 1;
downSampleFactor = 16;

if(~exist('versions', 'var'))
    inputFolder = '/home/data/EEG/processed/Robi/coherenceIca/old';
    files = dir(inputFolder);
    files([files.isdir]) = [];
    filter = 'Fast';
    files(cellfun(@length, strfind({files.name}, filter)) == 0) = [];
    
    for i = 1:length(files)
        versions{i} = load(fullfile(inputFolder, files(i).name));
        if(i == 1)
            remove = isnan(versions{i}.cohIca.out1(1,:));
        else
            remove = remove | isnan(versions{i}.cohIca.out1(1,:));
        end
    end
    componentCount = size(versions{1}.cohIca.out1,1);
    referenceCount = length(versions);
    
    if(downSample)
        rowCount = ceil(sum(~remove) / downSampleFactor);
        shortSamples = cell(1, length(versions));
        for referenceNumber = 1:referenceCount
            shortMatrix = NaN(rowCount, componentCount);
            for componentNumber = 1:componentCount
                column = downsample(versions{referenceNumber}.cohIca.out1(componentNumber, ~remove), downSampleFactor)';
                shortMatrix(:, componentNumber) = column;
            end
            shortSamples{referenceNumber} = shortMatrix;
        end
    end
end

if(debug)
    componentCount = 100;
end


rho = NaN(componentCount,referenceCount,componentCount,referenceCount);
p = NaN(componentCount,referenceCount,componentCount,referenceCount);
matchRho = NaN(componentCount, referenceCount);
matchRhoIndex = cell(componentCount, referenceCount);
referencePairCounter = 1;
for reference1 = 1:referenceCount
    if(~downSample)
        matrix1 = versions{reference1}.cohIca.out1(:, ~remove)';        
    else
        matrix1 = shortSamples{reference1};
    end
    for reference2 = 1:referenceCount
        if(reference1 ~= reference2)
            if(~downSample)
                matrix2 = versions{reference2}.cohIca.out1(:, ~remove)';
            else
                matrix2 = shortSamples{reference2};
            end
            for component1 = 1:componentCount
                fprintf('\n%s reference pair %d, component %d of %d', char(datetime), referencePairCounter, component1, componentCount);
                for component2 = component1+1:componentCount
                    if(reference1 ~= reference2)
                        [thisRho, thisP] = corr(matrix1(:,component1), matrix2(:,component2));
                        rho(component1, reference1, component2, reference2) = thisRho;
                        rho(component2, reference1, component1, reference2) = thisRho;
                        rho(component1, reference2, component2, reference1) = thisRho;
                        rho(component2, reference2, component1, reference1) = thisRho;
                        
                        p(component1, reference1, component2, reference2) = thisP;
                        p(component2, reference1, component1, reference2) = thisP;
                        p(component1, reference2, component2, reference1) = thisP;
                        p(component2, reference2, component1, reference1) = thisP;
                    end
                end
                
            end
            referencePairCounter = referencePairCounter + 1;
        end
    end
end
referencePairCounter = 1;
for reference1 = 1:referenceCount
    for reference2 = 1:referenceCount
        if(reference1 ~= reference2)
            for component1 = 1:componentCount
                row = rho(component1, reference1, :, reference2);
                row(isnan(row)) = [];
                row = abs(row);
                [match.maxRho(component1, referencePairCounter), match.maxRhoIndex(component1, referencePairCounter)] = max(row);
            end
            match.referencePairLabels{referencePairCounter} = sprintf('%s vs. %s', files(reference1).name, files(reference2).name);
            referencePairCounter = referencePairCounter + 1;
        end
    end
end
figure;
imagesc(match.maxRho');
title('Pairwise correlation of ICA signal timecourses');
ylabel('reference pair')
set(gca, 'yticklabel', {'Avg-Cpz','Avg-Mast','Cpz-Avg','Cpz-Mast','Mast-Avg','Mast-Cpz'});
set(gca, 'xtick', []);
xlabel('maximum correlation for each component');
cbar = colorbar;
ylabel(cbar, 'R^2');

save('robi003IcaReferenceCrossCorrelation.mat', 'match', 'rho', 'p', '-v7.3');

if(false)
    componentNumber = 1;
    for i = 1:length(versions)
        myTitle = strrep(files(i).name, '_', ' ');
        plotCoherencePca(versions{i}.cohPca, versions{i}.cohPca.labelChannelFreq, componentNumber, myTitle);
    end
end