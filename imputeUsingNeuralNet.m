function input = imputeUsingNeuralNet(input, hiddenLayerSize)

if(~exist('hiddenLayerSize', 'var'))
    hiddenLayerSize = 4;
end

if(size(input, 1) > size(input, 2))
    error('number of observations (columns) should be larger than the number of variables (rows)');
end
rhos = NaN(size(input, 1), 1);
missingCount = sum(isnan(input), 2);

%impute data using available variables with no missing data
%this could be modified to use partial observations; basically whatever
%supplies the most data; the current method will work when many
%observatiosn are complete.
completeVariable = false(size(missingCount));
completeVariable(missingCount == 0) = true;
originalMissing = sum(missingCount == 0);

%iterate through variables missing from the fewest observations first
sortCount = unique(missingCount);
sortCount(sortCount == 0) = [];
trainingCounter = 1;

for countIndex = 1:length(sortCount)
    count = sortCount(countIndex);
    missingVars = find(missingCount == count);
    net = feedforwardnet(hiddenLayerSize);
    net.trainParam.showWindow = false;
    for varIndex = 1:length(missingVars)
        outputFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\neuralnet\imputed\hiddenLayer4\';
        outputFilename = sprintf('%s%d.mat', outputFolder, trainingCounter);
        if(~exist(outputFilename, 'file'))
            %train the network using the observations that are complete for the
            %variable in question
            var = missingVars(varIndex);
            fprintf('%s: imputing variable number %d of %d\n', char(datetime), trainingCounter, originalMissing);
            
            fullObservations = ~isnan(input(var,:));
            completeData = input(completeVariable, fullObservations);
            missingData = input(var, fullObservations);
            [r, net, trainResultStats] = neuralNetworkTweakExisting(net, completeData, missingData);
            
            %impute missing data using trained network
            imputed = net(input(completeVariable, ~fullObservations));
            input(var, ~fullObservations) = imputed;
            
            %catalogue results
            result.rho = r;
            result.net = net;
            result.trainingStats = trainResultStats;
            result.variableIndex = var;
            %             results(trainingCounter) = result;
            save(outputFilename, 'result');
        end
        
        trainingCounter = trainingCounter  + 1;
        save('C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\imputed.mat', '-v7.3');
    end
    %only update this at the end of the loop because we want to reuse the
    %network for variables that are missing the same number of observations
    completeVariable(sum(isnan(input), 2) == 0) = true;
    
end


folder = 'C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\neuralnet\imputed\hiddenLayer4';
files = dir(folder);
files([files.isdir]) = [];
fullData = load('combinedEegMri.mat');
alphaLabels = fullData.eeg.labels;
alphaLabels = alphaLabels(cellfun(@length, strfind(alphaLabels, '9-12Hz')) > 0);

for i = 1:length(files)
    data = load(fullfile(folder, files(i).name));
    training(i) = data.result.rho.training;
    test(i) = data.result.rho.test;
    validation(i) = data.result.rho.validation;
    imputedLabels(i) = alphaLabels(data.result.variableIndex);
end

alphaPowerInd = find(cellfun(@length, strfind(imputedLabels, 'absolute')) > 0);
alphaPowerLabels = imputedLabels(alphaPowerInd);
alphaPowerRho = test(alphaPowerInd);
[sortRho, sortInd] = sort(test);
sortLabels = imputedLabels(sortInd)

alphaPowerRho = test(alphaPowerInd);

if(false)
    
    for i = 1:length(alphaPowerLabels)
        label = alphaPowerLabels{i};
        items = strsplit(label);
        chanlocs(i).labels = items{2};
    end
    dummyEeg.chanlocs = chanlocs;
    dummyEeg = setStandardLocations(dummyEeg);
    
    figure;
    plot(sortRho);
    ylabel('Pearson''s R, predicted vs. actual test data');
    
    figure;
    topoplot(alphaPowerRho, dummyEeg.chanlocs, 'maplimits', 'maxmin', 'conv', 'on');
    colorbar;
    title('Pearson''s R, predicted vs. actual test data');
end


%comparison to existing (pca-based) imputation

if(false)
    fullData = load('C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\combinedEegMri.mat');
    alphaInd = find(cellfun(@length, strfind(fullData.eeg.labels, '9-12Hz')) > 0);
    alphaSignal = fullData.eeg.signals(:, alphaInd);
    alphaLabls = fullData.eeg.labels(alphaInd);
    
    componentCount = 20;
    %impute using mean values to produce a starting point for nan
    alphaMeanedSignal = alphaSignal;
    for i = 1:size(alphaMeanedSignal, 2)
        nans = isnan(alphaMeanedSignal(:,i));
        meanValue = mean(alphaMeanedSignal(~nans, i));
        alphaMeanedSignal(nans, i) = meanValue;
    end
    [pcaResult0.COEFF,pcaResult0.SCORE,pcaResult0.LATENT,pcaResult0.MU,pcaResult0.V,pcaResult0.S] = pca(alphaMeanedSignal, 'NumComponents', componentCount);
    
    %     [pcaResult.COEFF,pcaResult.SCORE,pcaResult.LATENT,pcaResult.MU,pcaResult.V,pcaResult.S] = ppcaQuicker(alphaSignal, componentCount, 'W0', pcaResult0.COEFF);
    ppcaOptions.Display = 'iter';
    ppcaOptions.MaxIter = 500;
    ppcaOptions.TolFun = 1e-6;
    ppcaOptions.TolX = 1e-6;
    
    tic;
    [pcaResult.COEFF,pcaResult.SCORE,pcaResult.LATENT,pcaResult.MU,pcaResult.V,pcaResult.S] = ppcaQuicker(alphaSignal, componentCount, 'W0', pcaResult0.COEFF, 'Options', ppcaOptions);
    pcaResult.elapsedSeconds = toc;
    pcaResult.maxComponennts = componentCount;
    
    
    %reconstruct missing dataset and compare to known values
    reconstructed = pcaResult.SCORE * pcaResult.COEFF';
    for i = 1:size(alphaSignal, 2)
        meanValue = nanmean(alphaSignal(:, i));
        reconstructed(:, i) = reconstructed(:, i) + meanValue;
    end
    for i = 1:size(alphaSignal, 2)
        realMeasure = ~isnan(alphaSignal(:, i));
        [rhos(i), ps(i)] = corr(alphaSignal(realMeasure, i), reconstructed(realMeasure, i));
    end
    meanRho = mean(abs(rhos))
    
    %plot a characteristically missing variable
    missingCount = sum(isnan(alphaSignal));
    ind = find(missingCount == max(missingCount))
    indMin = find(missingCount == min(missingCount))
    column = ind(1);
    toPlot = [alphaMeanedSignal(:, column), reconstructed(:, column)];
    figure;
    plot(toPlot);
    legend('actual', 'predicted');
    
    
    save('C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\imputedByPpca.mat', 'pcaResult', '-v7.3');
end

if(false)
    %test reconstruction on meaned data (with none missing)
    %newComponentCount = size(alphaMeanedSignal, 2);
    newComponentCount = 20;
    [pcaResult0.COEFF,pcaResult0.SCORE,pcaResult0.LATENT,pcaResult0.MU,pcaResult0.V,pcaResult0.S] = pca(alphaMeanedSignal, 'NumComponents', newComponentCount);
    reconstructed0 = pcaResult0.SCORE * pcaResult0.COEFF';
    for i = 1:size(alphaMeanedSignal, 2)
        meanValue = mean(alphaMeanedSignal(:, i));
        reconstructed0(:, i) = reconstructed0(:, i) + meanValue;
    end
    bigDiff0 = reconstructed0 - alphaMeanedSignal;
    for i = 1:size(alphaMeanedSignal, 2)
        [rhos0(i), ps0(i)] = corr(alphaMeanedSignal(:, i), reconstructed0(:, i));
    end
    meanRho = mean(abs(rhos0))

    
    
    reconstructed = pcaResult.SCORE * pcaResult.COEFF';
    bigDiff = reconstructed
    
    a = rand(100,10);
    [b.COEFF,b.SCORE,b.LATENT,b.MU,b.V,b.S] = pca(a);
    
end

save('C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\imputed.mat', '-v7.3');
