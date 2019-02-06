

missing = fullData.eeg.signals;
missingCount = sum(sum(isnan(missing)));
alphaLabels = fullData.eeg.labels;
alphaInd = cellfun(@length, strfind(alphaLabels, '9-12Hz')) > 0;
alphaLabels = alphaLabels(alphaInd);
missing = missing(:, alphaInd);

tic;
[pcaResult.COEFF,pcaResult.SCORE,pcaResult.LATENT,pcaResult.MU,pcaResult.V,pcaResult.S] = ppca(missing, 20);
pcaResult.timeTaken = toc;