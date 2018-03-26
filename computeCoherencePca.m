%computeCoherencePca( inputFilename, outputFilename )
%COMPUTECOHERENCEPCA Summary of this function goes here
%   Detailed explanation goes here

inputFolder = '/media/Seagate Backup Plus Drive/pca';
allFiles = dir(inputFolder);
useBinary = true;


if(~useBinary)
    a = cellfun(@length, strfind({allFiles.name}, '.eegData'));
    inputFiles = allFiles(find(a));
    for fileCounter = 1:length(inputFiles)
        try
            inputFilename = inputFiles(fileCounter).name;% '/media/Seagate Backup Plus Drive/pca/PM102.bdf.coherence.bin';
            inputPath = fullfile(inputFolder, inputFilename);
            intermPath1 = strrep(inputPath, '.eegData', '.coherence.mat');
            outputFilename = strrep(inputPath, '.eegData', '.coherencePca.mat'); %/media/Seagate Backup Plus Drive/pca/PM102.pca.mat';
            if(~exist(intermPath1, 'file'))
                EEG.data = loadRobiDataFile(inputPath)';
                EEG.srate = 1024;
                EEG.nbchan = 34;
                [labels, locs] = wahbehChannelLocs;
                EEG.chanlocs = locs(1:34);
                [ coh.channelPairs, coh.x, coh.channels ] = allChannelCoherence(EEG);
                save(intermPath1, 'coh', '-v7.3');
            else
                load(intermPath1);
            end
            
            if(false)
                data = NaN(sampleCount, sampleSize);
                index = 2;
                counter = 1;
                while(index < length(contents))
                    data(counter, :) = contents(index:index+sampleSize-1);
                    index = index + sampleSize + 1;
                    counter = counter + 1;
                end
                invalidColumn = any(isnan(data),1);
                data(:,invalidColumn) = [];
                
                [p.COEFF, p.SCORE, p.LATENT, p.TSQUARED, p.EXPLAINED, p.MU] = pca(data);
                fprintf('\n%s) %s', char(datetime), outputFilename);
                verboseFilename = strrep(outputFilename, '.mat', '.verbose.mat');
                save(verboseFilename, 'p', '-v7.3');
                p.SCORE = [];
                p.TSQUARED = [];
                save(outputFilename, 'p', '-v7.3');
            end
        catch ex
        end
    end
else
    
    a = cellfun(@length, strfind({allFiles.name}, '.bdf.coherence.bin'));
    
    inputFiles = allFiles(find(a));
    
    
    for fileCounter = 2:length(inputFiles)
        
        try
            
            inputFilename = inputFiles(fileCounter).name;% '/media/Seagate Backup Plus Drive/pca/PM102.bdf.coherence.bin';
            inputPath = fullfile(inputFolder, inputFilename);
            outputFilename = strrep(inputPath, '.bdf.coherence.bin', '.pca.mat'); %/media/Seagate Backup Plus Drive/pca/PM102.pca.mat';
            if(~exist(outputFilename, 'file'))
                fileId = fopen(inputPath);
                fileInfo = dir(inputPath);
                contents = fread(fileId, fileInfo.bytes / 8, 'double');
                fclose(fileId);
                
                chanCount = 32;
                sampleSize = chanCount * (chanCount-1) / 2 * 5;
                
                sampleCount = length(contents) / (sampleSize + 1);
                fprintf('\n%f\n', sampleCount)
                
                data = NaN(sampleCount, sampleSize);
                
                index = 2;
                counter = 1;
                while(index < length(contents))
                    data(counter, :) = contents(index:index+sampleSize-1);
                    index = index + sampleSize + 1;
                    counter = counter + 1;
                end
                invalidColumn = any(isnan(data),1);
                data(:,invalidColumn) = [];
                
                [p.COEFF, p.SCORE, p.LATENT, p.TSQUARED, p.EXPLAINED, p.MU] = pca(allData);
                fprintf('\n%s) %s', char(datetime), outputFilename);
                verboseFilename = strrep(outputFilename, '.mat', '.verbose.mat');
                save(verboseFilename, 'p', '-v7.3');
                p.SCORE = [];
                p.TSQUARED = [];
                save(outputFilename, 'p', '-v7.3');
            end
        catch ex
            save(outputFilename, 'ex');
        end
    end
    
end
%end


