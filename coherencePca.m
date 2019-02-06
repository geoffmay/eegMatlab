inputFolder = '/media/eegDrive/';
inputFiles = dir(inputFolder);
inputFiles([inputFiles.isdir]) = [];

if(false)
%filter
outcome = cellfun(@length, strfind({inputFiles.name},'outcome'));
baseline = cellfun(@length, strfind({inputFiles.name},'baseline'));
remove = find(~baseline & ~outcome);
inputFiles(remove) = [];
end

for fileNumber = 1:length(inputFiles)
    inputFilename = fullfile(inputFolder, inputFiles(fileNumber).name);
    fprintf('\n%s file# %d of %d', char(datetime), fileNumber, length(inputFiles));    
    coherenceData = load(inputFilename);
    badReferenceFrame = checkForBadReference(coherenceData.channelPairs(31).coherence(:,1));
    if(badReferenceFrame < size(coherenceData.channelPairs(31).coherence, 1))
        for i = 1:length(coherenceData.channelPairs)
            coherenceData.channelPairs(i).coherence(badReferenceFrame:end,:) = [];
        end 
    end
    if(isfield(coherenceData, 'channelPairs'))
        
        %determine output filename
        [folder, file, ext] = fileparts(coherenceData.filename);
        file = strrep(folder, '/home/data/EEG/data/ROBI/', '');
        file = strrep(file, '/', ' ');
        outputFolder = '/home/data/EEG/processed/Robi/coherencePca/';
        outputFile = sprintf('%s%s.mat', outputFolder, file);
        
        
        cohSize = size(coherenceData.channelPairs(1).coherence);
        cohSize(2) = cohSize(2) * length(coherenceData.channelPairs);
        
        allCoh = NaN(cohSize);
        freqLabels = {'delta','theta','alpha','beta','hibeta'};
        
        start = 1;
        finish = start + 4;
        for i = 1:length(coherenceData.channelPairs)
            allCoh(:, start:finish) = coherenceData.channelPairs(i).coherence;
            for(j = 1:length(freqLabels))
                allLabels{start + j - 1} = sprintf('%s %s', coherenceData.channelPairs(i).label, freqLabels{j});
            end
            
            start = start + 5;
            finish = start + 4;
            
        end
        
        %remove mastoids
        m1 = cellfun(@length, strfind(allLabels, 'M1'));
        m2 = cellfun(@length, strfind(allLabels, 'M2'));
        remove = find(m1 | m2);
        allCoh(:,remove) = [];
        allLabels(:,remove) = [];
        
        %perform pca
        [cohPca.COEFF, cohPca.SCORE, cohPca.LATENT, cohPca.TSQUARED, cohPca.EXPLAINED, cohPca.MU] = pca(allCoh);
        close all;
        
        %get channel lables only (no frequency)
        onlyLabels = allLabels;
        for i = 1:length(onlyLabels)
            label = onlyLabels{i};
            space = strfind(label, ' ');
            if(length(space) > 0)
                label(space(1):end) = [];
            end
            onlyLabels{i} = label;
        end
        cohPca.labelChannelOnly = onlyLabels;
        cohPca.labelChannelFreq = allLabels;
        cohPca.filename = coherenceData.filename;
        
        %display image of pair coefficients
        if(false)
            imagesc(cohPca.COEFF(:, 1:10));
            colorbar;
        end
        
        %check top values
        doPlot = 0;
        doTopo = 0;
        maxPair = length(cohPca.labelChannelOnly);
        maxCoefficient = 20;
        
        if(doTopo)
            [chanLabs, chanlocs] = antChannelLocs;
        else
            [chanLabs] = antChannelLocs;
        end
        chanCoefficients = zeros(length(chanLabs), maxCoefficient);
        
        M1 = strcmp(chanLabs, 'M1');
        M2 = strcmp(chanLabs, 'M2');
        CPz = strcmp(chanLabs, 'CPz');
        remove = M1 | M2 | CPz;
        
        for i = 1:maxCoefficient
            a = cohPca.COEFF(:, i);
            [a1, ind] = sort(a, 1, 'descend');
            if(doPlot)
                if(~doTopo)
                    figure;
                    plot(a1);
                end
            end
            sortLabels = cohPca.labelChannelOnly(ind);
            maxWeight = max(cohPca.COEFF(:, i));
            alpha = cohPca.COEFF(:, i) ./ maxWeight;
            i
            if(true)
                %plot sums from top channels
                for chanPairCounter = 1:maxPair
                    chans = strsplit(sortLabels{chanPairCounter}, '-');
                    ind1 = find(strcmp(chanLabs, chans{1}));
                    ind2 = find(strcmp(chanLabs, chans{2}));
                    
                    chanCoefficients(ind1, i) = chanCoefficients(ind1, i) + cohPca.COEFF(chanPairCounter, i);
                    chanCoefficients(ind2, i) = chanCoefficients(ind2, i) + cohPca.COEFF(chanPairCounter, i);
                end
                if(doPlot)
                    figure;
                    topoplot(chanCoefficients(~remove,i), chanlocs(~remove), 'maplimits', [-3 3]);
                    colorbar;
                end
            end
            sortLabels(1:20)
            if(doPlot)
                if(doTopo)
                    plotChannelPairs(sortLabels(1:100));
                    title(sprintf('component %d %s', i, cohPca.filename));
                end
            end
        end
        %reduce size to save disk space
        cohPca.COEFF(:, maxCoefficient+1:end) = [];
        cohPca.SCORE(:, maxCoefficient+1:end) = [];
        cohPca.EXPLAINED(:, maxCoefficient+1:end) = [];
        
        save(outputFile, 'cohPca', 'chanCoefficients');
    end
end

