%compute the timecourse for each component during neurofeedback, using the
%separating matrix derived from the first resting state session.

outputFolder = '/media/eegDrive/003IcaTimeCourse';
cohFolder = '/media/eegDrive';
icaFolder = '/home/data/EEG/processed/Robi/coherenceIca/old2';
icaFile = 'ROBI_003_baseline eyes open-MastoidsFast.mat';
cohFiles = dir(cohFolder);
cohFiles([cohFiles.isdir]) = [];
fileFilter = '_003';
remove = find(~cellfun(@length, strfind({cohFiles.name}, fileFilter)));
cohFiles(remove) = [];

load(fullfile(icaFolder,icaFile));
sep = cohIca.separatingMatrix;
freqLabels = {' delta',' theta',' alpha',' beta',' hibeta'};
chanLabels = antChannelLocs;
m1 = find(strcmp(chanLabels,'M1'));
m2 = find(strcmp(chanLabels,'M2'));
cpz = find(strcmp(chanLabels,'CPz'));
chanLabels([m1 m2 cpz:end]) = [];
for fileCounter = 1:length(cohFiles)
    outputPath = fullfile(outputFolder, cohFiles(fileCounter).name);
    if(~exist(outputPath, 'file'))
        startDate = char(datetime);
        save(outputPath, 'startDate');
        data = load(fullfile(cohFolder, cohFiles(fileCounter).name));
        sampleCount = size(data.channelPairs(1).coherence,1);
        measureCount = length(cohIca.labelChannelFreq);
        timeCourse = NaN(measureCount, sampleCount);
        powerChannelCounter = 1;
        powerFreqCounter = 1;
        for measureCounter = 1:measureCount
            label = cohIca.labelChannelFreq{measureCounter};
            freqIndex = -1;
            for freqCounter = 1:length(freqLabels)
                if(strfind(label, freqLabels{freqCounter}))
                    freqIndex = freqCounter;
                end
            end
            if(length(strfind(label, 'abs')) == 0) %coherence
                space = strfind(label, ' ');
                thisLabel = label(1:space(1)-1);
                pairIndex = find(strcmp({data.channelPairs.label}, thisLabel));
                timeCourse(measureCounter,:) = data.channelPairs(pairIndex).coherence(:,freqIndex)';
            else %power
                timeCourse(measureCounter,:) = data.channels(powerChannelCounter).absolutePower(:,powerFreqCounter)';
                powerFreqCounter = powerFreqCounter + 1;
                if(powerFreqCounter > 5)
                    powerFreqCounter = 1;
                    powerChannelCounter = powerChannelCounter + 1;
                end
            end
        end
        componentTimecourses = cohIca.separatingMatrix * timeCourse;
        save(outputPath, 'componentTimecourses', '-v7.3');
    end
    if(length(strfind(outputPath, '_tx ')))
        eventFilename = cohFiles(fileCounter).name;
        eventFilename(9) = '/';
        slash = strfind(eventFilename, '_63');
        eventFilename(slash) = '/';
        eventFilename = strrep(eventFilename, 'coherenceStats.mat', 'Events.txt');
        eventFilename = ['/home/data/EEG/data/ROBI/', eventFilename];
        interpolated = loadInterpolatedZScore(eventFilename);
        clear componentTimecourses;
        load(outputPath); %componentTimecourse (matrix: componentCount x sampleCount)
        components = componentTimecourses';
        smoothZ = interpolated.SmoothZ';
        minLength = min(size(components,1), size(smoothZ,1));
        for i = 1:size(components,2)
            fprintf('\n%d of %d',i,size(components,2));
            [rho(i), p(i)] = corr(smoothZ(1:minLength, 1), components(1:minLength, i));
        end
        a = find(abs(rho') > .05);
        mat = cohIca.separatingMatrix(a,:);
        for i = 1:length(a)
            mat(i,:) = mat(i,:) .* rho(a(i));
        end
        Z = linkage(mat');
        perm = clusterPermutation(Z);       
        b = mat(:,perm);
        close all;
        lim = prctile(reshape(b,[1, numel(b)]), 90);
        imagesc(b,[-lim,lim]);
        colorbar;
        c = mean(mat,1);
        plotCoherencePca(c);

        a = sort(rho);
        close all;
        plot(a);
        dummy = 1;
        
    end
end
