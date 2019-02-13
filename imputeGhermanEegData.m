rootFolder = 'C:\Users\Neuro\Documents\MATLAB\data\GhermanPhiliastides\1.0.1\';
outputFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\imputed';

subFolders = dir(rootFolder);
subFolders(cellfun(@length, strfind({subFolders.name}, 'sub-')) == 0) = [];
if(~exist(outputFolder, 'file'))
    mkdir(outputFolder);
end


chanlocs = readlocs(fullfile(fileparts(which('imputeGhermanEegData')), 'GhermanEfMRIcap64.elp'));
for i = 1:length(chanlocs)
    chanlabels{i} = chanlocs(i).labels;
end

%load the data
fileCounter = 1;
sampleCounter = 0;
for (folderCounter = 1:length(subFolders))
    subFolder = fullfile(rootFolder, subFolders(folderCounter).name, 'EEG');
    files = dir(subFolder);
    files(cellfun(@length, strfind({files.name}, 'EEG_data')) == 0) = [];
    for fileCounter = 1:length(files)
        fprintf('loading folder %d file %d', folderCounter, fileCounter);
        filename = fullfile(subFolder, files(fileCounter).name);
        %         outputFilename = fullfile(outputFolder, files(fileCounter).name);
        %         if(~exist(outputFilename, 'file'))
        %placeholder = sprintf('started on %s', char(datetime));
        %         save(outputFilename, 'placeholder');
        if(length(strfind(filename, '.html')) == 0)
            eeg = loadGhermanEeg(filename);
            data = eeg.data';
            data = downsample(data, 100);
            firstFrame = sampleCounter+1;
            lastFrame = sampleCounter + size(data,1);
            bigEeg(firstFrame:lastFrame, :) = NaN;
            for i = 1:eeg.nbchan
                ind = find(strcmp(chanlabels, eeg.chanlocs(i).labels));
                bigEeg(firstFrame:lastFrame, ind) = data(:, i); %eeg.data(i, :)';
            end
            sampleCounter = lastFrame;
            lastFrames(fileCounter) = lastFrame;
            fileCounter = fileCounter + 1;
        end
        %             tic;
        %             phaseSlopeTopographies = phaseSlopeTopography(eeg);
        %             elapsed = toc;
        %             save(outputFilename, 'phaseSlopeTopographies', 'elapsed', '-v7.3');
    end
end

save('C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\eegSrate100.mat', 'bigEeg', '-v7.3');


[mem, sys] = memory();
physMem = sys.PhysicalMemory.Available;

dsEeg = downsample(bigEeg, 4);
memCap = 2048;
while(memCap < physMem)
    [output, pcaResult] = imputeUsingPpca(dsEeg, memCap);
    outputFilename = sprintf('%s%d.mat', 'C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\imputed\ppcaImputed', memCap);
    save(outputFilename, 'output', 'pcaResult', '-v7.3');
    memCap = memCap * 2;
end

outputFiles = dir(outputFolder);
outputFiles([outputFiles.isdir]) = [];
if(false)
    for i = 1:length(outputFiles)
        fprintf('loading file %d of %d\n', i, length(outputFiles));
        a = load(fullfile(outputFolder, outputFiles(i).name));
        rhos(i) = mean(a.pcaResult.correlations);
    end
    figure;
    plot(rhos);
    
    good = find(rhos > .9);
    clear coeffs;
    for i = 1:length(good)
        a = load(fullfile(outputFolder, outputFiles(i).name));
        rhos(i) = mean(a.pcaResult.correlations);
        coeffs(:,:,i) = a.pcaResult.COEFF;
        latent(:,i) = a.pcaResult.LATENT;
        chan1(:,i) = a.output(:,1);
    end
    
    
    figure;
    counter = 1;
    side = ceil(sqrt(size(coeffs,3)));
    for i = 1:side
        for j = 1:side
            if(counter < size(coeffs,3))
                subplot(side, side, counter);
                imagesc(coeffs(:,:,counter));
                colorbar;
                title(sprintf('%d', counter));
                counter = counter + 1;
            end
        end
    end
    
    
    coeffFlat = reshape(coeffs, [size(coeffs,1)*size(coeffs,2), size(coeffs,3)]);
    for i = 1:size(coeffFlat, 2)
        for j = i+1:size(coeffFlat, 2)
            [coeffRho(i,j), p] = corr(coeffFlat(:,i), coeffFlat(:, j));
            [latentRho(i,j), p] = corr(latent(:,i), latent(:, j));
            [chan1Rho(i,j), p] = corr(latent(:,i), latent(:, j));
        end
    end
    
    
end