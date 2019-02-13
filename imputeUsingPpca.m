function [output, ppcaResult] = imputeUsingPpca(input, maxMemory)

verbose = 0;

%error checking
if(size(input, 1) < size(input, 2))
    error('number of observations (rows) should be larger than the number of variables (columns)');
end


%figure out max components based on available memory
if(~exist('maxMemory', 'var'))
    [memInfo, memSys] = memory();
    %maxMemory = memInfo.MaxPossibleArrayBytes;
    maxMemory = memSys.PhysicalMemory.Available;
end
%maxComponents = floor(sqrt(maxMemory / size(input,1) / 8));

componentCount = size(input, 2);
maxSamples = floor(maxMemory / (componentCount * componentCount * 8 * 2));
minSamples = size(input, 2);
if(maxSamples < minSamples)
    maxSamples = minSamples;
end
downsampleRate = ceil(size(input, 1) / maxSamples);

%if(verbose)
% if(maxComponents >= componentCount)
%     fprintf('imputing with all %d available components', componentCount);
% else
%     componentCount = maxComponents;
%     fprintf('trimming to %d components to conserve memory', componentCount);
% end
%end

if(downsampleRate == 1)
    if(verbose)
        fprintf('imputing with all %d observations\n', size(input, 1));
    end
else
    oldSize = size(input,1);
    dsInput = downsample(input, downsampleRate);
    if(verbose)
        fprintf('downsapmling to %d observations (%0.1f%%) based on available physical memory (%0.2f GB)\n', size(input, 1), size(input,1) / oldSize * 100, maxMemory / 1024  / 1024  / 1024);
    end
end


%impute using mean values to produce a starting point for nan
meanedInput = dsInput;
for i = 1:size(meanedInput, 2)
    nans = isnan(meanedInput(:,i));
    meanValue = mean(meanedInput(~nans, i));
    meanedInput(nans, i) = meanValue;
end
[pcaResult0.COEFF,pcaResult0.SCORE,pcaResult0.LATENT,pcaResult0.MU,pcaResult0.V,pcaResult0.S] = pca(meanedInput, 'NumComponents', componentCount);


ppcaOptions.Display = 'iter';
ppcaOptions.MaxIter = 500;
ppcaOptions.TolFun = 1e-6;
ppcaOptions.TolX = 1e-6;

tic;
[ppcaResult.COEFF,ppcaResult.SCORE,ppcaResult.LATENT,ppcaResult.MU,ppcaResult.V,ppcaResult.S] = ppcaQuicker(dsInput, componentCount, 'W0', pcaResult0.COEFF, 'Options', ppcaOptions);
ppcaResult.elapsedSeconds = toc;
ppcaResult.maxComponennts = componentCount;

%check for nans in score
if(false)
    scoreNanInd = find(isnan(pcaResult.SCORE));
    scoreNanCoord(:,1) = ceil(scoreNanInd ./ size(pcaResult.SCORE,1))
    scoreNanCoord(:,2) = ceil(mod(scoreNanInd, size(pcaResult.SCORE,1)))
    [scoreNanRow, scoreNanCol] = ind2sub(size(pcaResult.SCORE), scoreNanInd);
end

%apply coefficients to entire input span
if(downsampleRate > 1)
    if(verbose)
        fprintf('%s: using pca coefficients from downsampled data to full sample\n', char(datetime));
    end
    ppcaResult.downsampleRate = downsampleRate;
    ppcaResult.samplesUsed = size(dsInput, 1);
    tic;
    [output, recomputeResult] = imputeUsingFixedComponents(input, ppcaResult);
    ppcaResult.extrapolationSeconds = toc;
else
    output = ppcaResult.SCORE * ppcaResult.COEFF';
end

for i = 1:size(input, 2)
    %     meanValue = nanmean(input(:, i));
    %     reconstructed(:, i) = reconstructed(:, i) + meanValue;
    realMeasure = ~isnan(input(:, i));
    [rhos(i), ps(i)] = corr(input(realMeasure, i), output(realMeasure, i));
end
ppcaResult.correlations = rhos;


if(false)
    %reconstruct missing dataset and compare to known values
    reconstructed = ppcaResult.SCORE * ppcaResult.COEFF';
    output = input;
    for i = 1:size(input, 2)
        meanValue = nanmean(input(:, i));
        reconstructed(:, i) = reconstructed(:, i) + meanValue;
        realMeasure = ~isnan(input(:, i));
        [rhos(i), ps(i)] = corr(input(realMeasure, i), reconstructed(realMeasure, i));
        output(~realMeasure, i) = reconstructed(~realMeasure, i);
    end
    ppcaResult.correlations = rhos;
end

% meanRho = mean(abs(rhos));

if(false)
    %plot a characteristically missing variable
    missingCount = sum(isnan(input));
    ind = find(missingCount == max(missingCount))
    indMin = find(missingCount == min(missingCount))
    column = ind(1);
    toPlot = [input(:, column), output(:, column)];
    figure;
    plot(toPlot);
    legend('actual', 'predicted');
end

% save('C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\imputedByPpca.mat', 'pcaResult', '-v7.3');
%
% save('C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\imputed.mat', '-v7.3');

if(false)
    rng(1);
    comps = rand(100, 2);
    coeffs = [10 1 2 4; 3 5 1 2];
    sigs = comps * coeffs;
    missingFraction = 0.05;
    makeNan = rand(size(sigs)) < missingFraction;
    keepRow = mod(1:size(sigs,1), 2) == 1;
    sigDrop = sigs;
    sigDrop(makeNan) = NaN;
    subSample = sigDrop(keepRow, :);
    
    [pcam.COEFF,pcam.SCORE,pcam.LATENT,pcam.MU,pcam.V,pcam.S] = ppca(subSample, 2);
    
    [pca0.COEFF,pca0.SCORE,pca0.LATENT,pca0.MU,pca0.V,pca0.S] = pca(sigs);
    reconstructed = pca0.SCORE*pca0.COEFF';
    sig2comp = inv(pca0.COEFF);
    %     reconstructedScore = sigs * sig2comp';
    [reconstructed, stats] = imputeUsingFixedComponents(sigDrop, pcam);
    
    
    for i = 1:size(sigs,2)
        meanValue = mean(sigs(:,i));
        reconstructed(:,i) = reconstructed(:,i) + meanValue;
    end
    [rho, p] = corr(sigs(:), reconstructed(:));
    
    if(false)
        epsilon = 1e-20;
        emptyColumn = pca0.LATENT < epsilon;
        %     emptyColumn = all(abs(reconstructedScore) < epsilon);
        reconstructedScore(:, emptyColumn) = [];
        
        %     ratio = reconstructedScore ./ comps;
        
        for i = 1:size(reconstructedScore,2)
            scaledScore(:,i) = reconstructedScore(:, i) .* pca0.LATENT(i);
        end
        
        recMean = mean(reconstructedScore);
        compMean = mean(comps);
        
        recStd = std(reconstructedScore);
        compStd = std(comps);
        
        
        normComp = comps(:,1);
        normComp = normComp - mean(normComp(:,1));
        normComp = normComp ./ std(normComp(:,1));
        normRec = reconstructedScore(:,1);
        normRec = normRec - mean(normRec(:,1));
        normRec = normRec ./ std(normRec(:,1));
        
        
        %     normRec = reconstructedScore(:,1) - mean(reconstructedScore(:,1)) ./ std(reconstructedScore(:,1));
        figure;
        plot([normComp, normRec]);
        legend({'normComp', 'normRec'});
        
        [rho2, p2] = corr(normComp, normRec);
        
        close all;
        scatter(comps(:,1), reconstructedScore(:,1));
        
        plot([comps(:,1), scaledScore(:,1), reconstructedScore(:,1)])
        legend({'original', 'scaled', 'reconstructed'});
    end
end

