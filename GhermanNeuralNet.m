

%if false, adds channels one at a time, and includes all coherence pairs
individualCoherence = false;


%pick up where we left off

dataFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides';
outFolder = fullfile(dataFolder, 'neuralnet');
combinedFilename = fullfile(dataFolder, 'combinedEegMri.mat');

if(~exist(combinedFilename, 'file'))
    eegFolder = '/home/data/EEG/processed/GhermanPhilastides/convolved';
    mriFolder = '/home/data/EEG/processed/GhermanPhilastides/networkTimecourses';
    eegFiles = dir(eegFolder);
    eegFiles = eegFiles(cellfun(@length, strfind({eegFiles.name}, '.mat')) > 0);
    
    %load mri data
    volumeCounter = 0;
    validFileCounter = 0;
    for i = 1:length(eegFiles)
        numberString = strrep(eegFiles(i).name, 'EEG_data_sub-', '');
        if(strcmp(numberString, eegFiles(i).name))
            error('unexpected eeg filename: %s', eegFiles(i).name);
        end
        subjectString = numberString(1:2);
        runString = numberString(8:9);
        if(strcmp(runString, '01'))
            validFileCounter = validFileCounter + 1;
            %             mriFilename = sprintf('/home/data/subjects/sub-GhermanPhiliastides%s/infomap/main/groupnetworktimecourses.txt', subjectString);
            %             mriTimeCourse = textread(mriFilename);
            %             networkNumbers = mriTimeCourse(1,:);
            %             mriTimeCourse = mriTimeCourse(2:end,:);
            
            mriFilename = fullfile(mriFolder, sprintf('sub-GhermanPhiliastides%s.mat', subjectString));
            mriTimeCourse = load(mriFilename);
            
            %             nonNanVolumes(validFileCounter) = sum(~all(isnan(mriTimeCourse), 2));
            %             volumeCount = size(mriTimeCourse, 1);
            volumeCount = size(mriTimeCourse.networks.timecourse, 1);
            mri.networkNumbers = mriTimeCourse.networks.numbers;
            mri.subjects(volumeCounter+1:volumeCounter+volumeCount) = str2num(subjectString);
            %     mri.runs(volumeCounter+1:volumeCounter+volumeCount) = str2num(runString);
            mri.networkActivity(volumeCounter+1:volumeCounter+volumeCount, :) = mriTimeCourse.networks.timecourse;
            volumeCounter = volumeCounter + volumeCount;
        end
    end
    
    %get eeg labels
    labels = cell(0);
    for i = 1:length(eegFiles)
        fprintf('loading labels file %d of %d\n', i, length(eegFiles));
        eegFilename = fullfile(eegFolder, eegFiles(i).name);
        eeg = load(eegFilename);
        labels = unique([labels, eeg.convolvedData.labels]);
    end
    
    sampleEeg = load('/home/data/subjects/sub-GhermanPhiliastides02/EEG/EEG_events_sub-02_run-02.mat');
    
    %load eeg data
    volumeCounter = 0;
    clear eeg;
    eeg.labels = labels;
    eeg.signals = NaN(size(mri.networkActivity, 1), length(labels));
    for i = 1:length(eegFiles)
        fprintf('organizing eeg data, file %d of %d\n', i, length(eegFiles));
        numberString = strrep(eegFiles(i).name, 'EEG_data_sub-', '');
        if(strcmp(numberString, eegFiles(i).name))
            error('unexpected eeg filename: %s', eegFiles(i).name);
        end
        subjectString = numberString(1:2);
        runString = numberString(8:9);
        eegFilename = fullfile(eegFolder, eegFiles(i).name);
        thisEeg = load(eegFilename);
        volumeCount = size(thisEeg.convolvedData.signals, 2);
        for j = 1:size(thisEeg.convolvedData.signals, 1)
            destInd = find(strcmp(eeg.labels, thisEeg.convolvedData.labels{j}));
            eeg.signals(volumeCounter+1:volumeCounter + volumeCount, destInd) = thisEeg.convolvedData.signals(j, :)';
        end
        eeg.subjects(volumeCounter+1:volumeCounter+volumeCount) = str2num(subjectString);
        eeg.runs(volumeCounter+1:volumeCounter+volumeCount) = str2num(runString);
        eeg.times(volumeCounter+1:volumeCounter+volumeCount) = thisEeg.convolvedData.times;
        volumeCounter = volumeCounter + volumeCount;
    end
    
    save(combinedFilename, 'eeg', 'mri', '-v7.3');
    
    %debug
    subjs = unique(eeg.subjects);
    for i = 1:length(subjs)
        fprintf('%d of %d\n', i, length(subjs));
        subCount(i, 1) = sum(eeg.subjects == subjs(i) & eeg.runs == 1);
        subCount(i, 2) = sum(eeg.subjects == subjs(i) & eeg.runs == 2);
        mri1 = sprintf('/home/data/subjects/sub-GhermanPhiliastides%02d/ses-01/func/sub-%02d_task-main_run-01_bold.nii.gz', subjs(i), subjs(i));
        mri2 = sprintf('/home/data/subjects/sub-GhermanPhiliastides%02d/ses-01/func/sub-%02d_task-main_run-02_bold.nii.gz', subjs(i), subjs(i));
        nif1 = niftiread(mri1);
        nif2 = niftiread(mri2);
        subCount(i, 3) = size(nif1, 4);
        subCount(i, 4) = size(nif2, 4);
        subCount(i, 5) = subCount(i,1)-subCount(i,3);
        subCount(i, 6) = subCount(i,2)-subCount(i,4);
        subCount(i, 7) = i;
    end
    for i = 1:length(subjs)
        run1StartTime(i, 1) = min(eeg.times(eeg.subjects == subjs(i) & eeg.runs == 1));
        run2StartTime(i, 1) = min(eeg.times(eeg.subjects == subjs(i) & eeg.runs == 2));
        run1EndTime(i, 1) = max(eeg.times(eeg.subjects == subjs(i) & eeg.runs == 1));
        run2EndTime(i, 1) = max(eeg.times(eeg.subjects == subjs(i) & eeg.runs == 2));
        run1MaxInterval(i, 1) = max(diff(eeg.times(eeg.subjects == subjs(i) & eeg.runs == 1)));
        run2MaxInterval(i, 1) = max(diff(eeg.times(eeg.subjects == subjs(i) & eeg.runs == 2)));
    end
    subCount
    
    subjectNumber = subCount(:,7);
    run1EegVolumes = subCount(:,1);
    run2EegVolumes = subCount(:,2);
    run1MriVolumes = subCount(:,3);
    run2MriVolumes = subCount(:,4);
    run1Difference = subCount(:,5);
    run2Difference = subCount(:,6);
    metadata.volumeCounts = table(subjectNumber, run1EegVolumes, run1MriVolumes, run1Difference, run2EegVolumes, run2MriVolumes, run2Difference);
    metadata.timingInfo = table(subjectNumber, run1StartTime, run1MaxInterval, run1EndTime, run2StartTime, run2EndTime, run2MaxInterval);
    
    save(combinedFilename, 'eeg', 'mri', 'metadata', '-v7.3');
    
    if(false)
        for i = 1:length(subjs)
            for j = 1:2
                counts(i,j) = sum(eeg.subjects == i & eeg.runs == j);
                counts(i,j+2) = sum(mri.subjects == i);
            end
        end
        
        subsample = eeg.times(eeg.subjects == 15 & eeg.runs == 2);
        subsample1 = eeg.times(eeg.subjects == 15 & eeg.runs == 1);
        d = diff(subsample);
        tabulate(d);
        d1 = diff(subsample1);
        tabulate(d1);
        fprintf('%f\n',max(subsample1));
        fprintf('%f\n',max(subsample));
        
        
        eegnans = sum(isnan(eeg.signals));
        tab = tabulate(eegnans);
        tab(tab(:,2)==0,:) = [];
        
        volumenans = sum(isnan(eeg.signals), 2);
        tabV = tabulate(volumenans);
        tabV(tabV(:,2)==0,:) = [];
        
        sampleFile = '/home/data/subjects/sub-GhermanPhiliastides01/ses-01/func/sub-01_task-main_run-01_bold.nii.gz';
        sampleFile2 = '/home/data/subjects/sub-GhermanPhiliastides01/ses-01/func/sub-01_task-main_run-02_bold.nii.gz';
        
        subjs = unique(mri.subjects);
        for s = 1:length(subjs)
            fprintf('%d: %d\n', subjs(s), sum(mri.subjects == subjs(s)));
        end
        %cifti = ft_read_cifti(sampleFile);
        a = niftiread(sampleFile);
        b = niftiread(sampleFile2);
        %     sampleFile2 = '/home/data/subjects/sub-GhermanPhiliastides23/ses-01/func/sub-23_task-main_run-01_bold.nii.gz';
        %cifti = ft_read_cifti(sampleFile);
        %     b = niftiread(sampleFile2);
        %end debug
    end
else
    if(~exist('eeg', 'var'))
        load(combinedFilename);
    end
end

subjectsWithExplainedProblems = [5, 14, 15];
%subject 5: mising 10 mri volumes from run 2
%subject 14: mising 5 mri volumes from each run
%subject 15: mising 5 mri volumes from run 2


subjectsWithUnexplainedProblems = [1, 7, 12, 23];
%subject 1: mising 1 eeg volume marker from run 1, 3 from run 2
%subject 7: mising 8 mri volumes from each run
%subject 12: mising 1 eeg volume marker from run 1
%subject 24: timing information is compressed to 1/5 time

perfectSubjects = 1:24;
perfectSubjects = setdiff(perfectSubjects, subjectsWithExplainedProblems);
perfectSubjects = setdiff(perfectSubjects, subjectsWithUnexplainedProblems);



%drop extra eeg volume markers to synchronize with MRI
discardLastVolumeMarkers = true;
onlyAlpha = true;

dropEeg = false(size(eeg.signals, 1), 1);
remove1 = [];
remove2 = [];
for i = 1:size(metadata.volumeCounts)
    subj = metadata.volumeCounts{i, 'subjectNumber'};
    run1Ind = find(eeg.subjects == subj & eeg.runs == 1);
    run2Ind = find(eeg.subjects == subj & eeg.runs == 2);
    run1VolumeCount = metadata.volumeCounts{i, 'run1MriVolumes'};
    run2VolumeCount = metadata.volumeCounts{i, 'run2MriVolumes'};
    if(discardLastVolumeMarkers)
        remove1 = run1Ind(run1VolumeCount+1:end);
        remove2 = run2Ind(run2VolumeCount+1:end);
    end
    dropEeg(remove1) = true;
    dropEeg(remove2) = true;
end
eeg.signals(dropEeg, :) = [];
eeg.subjects(dropEeg) = [];
eeg.runs(dropEeg) = [];
eeg.times(dropEeg) = [];


if(onlyAlpha)
    keepInd = cellfun(@length, strfind(eeg.labels, '9-12Hz')) > 0;
    eeg.labels = eeg.labels(keepInd);
    eeg.signals = eeg.signals(:, keepInd);
end

%drop motion frames
eegInput = eeg.signals';
mriOutput = mri.networkActivity';
dropMri = all(isnan(mriOutput), 1);
mriOutput(:, dropMri) = [];
eegInput(:, dropMri) = [];

%help neural net handle nans
eegNan = isnan(eegInput);

anyNans = any(eegNan, 2);

fastCheck = true;
if(fastCheck)
    neuralNet.cleanInput = eegInput(~anyNans,:);
    neuralNet.bigNetworks = mriOutput(1:10,:);
    neuralNet.net = feedforwardnet(30);
    neuralNet.eegLabels = eeg.labels(~anyNans);
    [neuralNet.r, neuralNet.net, neuralNet.trainResultStats] = neuralNetworkTweakExisting(neuralNet.net, neuralNet.cleanInput, neuralNet.bigNetworks);
    
        networks.colors = (...
        [.5, .5, .5; ... %0: unassigned = gray
        1, 0, 0; ... %1: default = red
        0, 0, .5; ... %2: secondary visual = blue
        1, 1, 0; ... %3: frontoparietal = yellow
        227/255, 172/255, 117/255; ... %4: primary visual = peach
        0, 1, 0; ... %5: dorsal attention = green
        230/255, 156/255, 231/255; ... %6: sensory(?) = lavender
        0, .5, .5; ... %7: ventral attention = teal
        0, 0, 0; ... %8: salience = black
        .5, 0, .5; ... %9: ciguloopercular = violet
        0, 1, 1; ... %10: motor hand = light blue
        227/255, 133/255, 45/255; ... %11: motor mouth = orange
        150/255, 79/255, 219/255; ... %12: opercular = lighter purple
        41/255, 82/255, 123/255; ... %13: ? = pale blue
        77/255, 236/255, 78/255; ... %14: ? = lime green
        42/255, 42/255, 252/255; ... %15: sky blue = parietal memory (?)
        1, 1, 1; ... %16: ? = white
        41/255, 129/255, 42/255 % 1,1,1 ... %17: motor foot = dark green
        ]);
    
    networks.labels = {'unassigned', 'default', 'secondary visual', 'frontoparietal',...
        'primary visual', 'dorsal attention', 'sensory', 'ventral attention', ...
        'salience', 'cingulooperular', 'motor hand', 'motor mouth', 'opercular', ...
        'pale blue', 'lime green', 'parietal memory', 'white', 'motor foot'}';
    neuralNet.networks = networks;
    
    save('C:\Users\Neuro\Documents\MATLAB\processed\GhermanPhilastides\neuralnet\cleanonly.mat', 'neuralNet', '-v7.3');
    
end


eegInput(eegNan) = 0;
eegNanBinary = ones(size(eegNan));
eegNanBinary(eegNan) = -1;
eegInput = [eegInput; eegNanBinary];
sampleOutput = mriOutput([2, 4],:);

if(false)
    sampleOutput = mriOutput([2, 4],:);
    sampleInput = eegInput(1:10, :);
    [regressed, neuralNet] = neuralNetworkAdvancedGeneration(sampleInput, sampleOutput, 4, true);
end

%do the fit
[regressed, neuralNet] = neuralNetworkAdvancedGeneration(eegInput, sampleOutput, 10, true);

%save results
outputFilename = fullfile(outFolder, 'firstPass.mat');
save(outputFilename, 'regressed', 'neuralNet', '-v7.3');



if(individualCoherence)
    % keep = cellfun(@length, strfind(trainingData.eegColumnLabels, 'ower')) > 0;
    oldKeep = false(size(trainingData.eegColumnLabels));
    
    outputData = trainingData.mriOutput';
    
    %add one eeg parameter at a time, the one with the highest test R value
    loopCounter = 1;
    while(loopCounter < length(oldKeep))
        maxR = 0;
        maxI = -1;
        clear maxTrial;
        for i = 1:length(oldKeep)
            %         fprintf('%s: loop %d, training %d of %d (maxI = %d, maxR = %f)\n', char(datetime), loopCounter, i, length(oldKeep), maxI, maxR);
            if(~oldKeep(i))
                keep = oldKeep;
                keep(i) = true;
                inputData = trainingData.eegInput(:, keep)';
                trial.r = neuralNetworkAdvancedGeneration(inputData, outputData);
                score = trial.r.test + trial.r.validation + 0.1 * trial.r.training;
                if(score > maxR)
                    trial.keep = find(keep);
                    trial.newIndex = i;
                    maxR = score;
                    maxI = i;
                    maxTrial = trial;
                end
            end
        end
        trials(loopCounter) = maxTrial;
        oldKeep = false(size(trainingData.eegColumnLabels));
        oldKeep(maxTrial.keep) = true;
        outFile = fullfile(outFolder, sprintf('sparseEegMriNeuralNet%d.mat', loopCounter));
        save(outFile, 'trials', '-v7.3');
        fprintf('%s: loop %d, maxR = %f, measure = %s\n', char(datetime), loopCounter, maxR, trainingData.eegColumnLabels{maxTrial.newIndex});
        loopCounter = loopCounter + 1;
    end
    
    inputData = trainingData.eegInput(:, keep)';
    outputData = trainingData.mriOutput';
    regressed = neuralNetworkAdvancedGeneration(inputData, outputData);
    
    rs = [trials.r];
    plot([[rs.training]', [rs.test]', [rs.validation]']);
    legend({'training', 'test', 'validation'});
else % individualChannels == false
    distributionSize = 21;
    outputData = trainingData.mriOutput';
    
    %get unique channel names
    channelItemIndex = 1;
    measureItemIndex = 2;
    allChans = cell(0);
    for i = 1:length(trainingData.eegColumnLabels)
        label = trainingData.eegColumnLabels{i};
        items = strsplit(label, ' ');
        if(i == 1)
            if(strcmp(items{1}, 'coherence') || (length(strfind(items{1}, 'ower')) > 0))
                channelItemIndex = 2;
                measureItemIndex = 1;
            elseif(strcmp(items{2}, 'coherence') || (length(strfind(items{2}, 'ower')) > 0))
                channelItemIndex = 1;
                measureItemIndex = 2;
            else
                error('unhandled eeg label format: %s', label);
            end
        end
        
        %if false, adds channels one at a time, and includes all coherence pairs
        individualCoherence = true;
        
        
        %pick up where we left off
        % outFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\MrEegLink';
        outFolder = '/home/data/EEG/processed/GhermanPhilastides/neuralnet';
        combinedFilename = '/home/data/EEG/processed/GhermanPhilastides/combinedEegMri.mat';
        
        if(exist(combinedFilename, 'file'))
            eegFolder = '/home/data/EEG/processed/GhermanPhilastides/convolved';
            mriFolder = '/home/data/EEG/processed/GhermanPhilastides/convolved';
            eegFiles = dir(eegFolder);
            eegFiles = eegFiles(cellfun(@length, strfind({eegFiles.name}, '.mat')) > 0);
            
            %load mri data
            volumeCounter = 0;
            for i = 1:length(eegFiles)
                numberString = strrep(eegFiles(i).name, 'EEG_data_sub-', '');
                if(strcmp(numberString, eegFiles(i).name))
                    error('unexpected eeg filename: %s', eegFiles(i).name);
                end
                subjectString = numberString(1:2);
                runString = numberString(8:9);
                mriFilename = sprintf('/home/data/subjects/sub-GhermanPhiliastides%s/infomap/main/groupnetworktimecourses.txt', subjectString);
                mriTimeCourse = textread(mriFilename);
                networkNumbers = mriTimeCourse(1,:);
                mriTimeCourse = mriTimeCourse(2:end,:);
                volumeCount = size(mriTimeCourse, 1);
                if(~exist('mri', 'var'))
                    mri.networkNumbers = networkNumbers;
                else
                    if(any(mri.networkNumbers ~= networkNumbers))
                        error('inconsistent network numbers across files at: %s', eegFilename);
                    end
                end
                mri.subjects(volumeCounter+1:volumeCounter+volumeCount) = str2num(subjectString);
                %     mri.runs(volumeCounter+1:volumeCounter+volumeCount) = str2num(runString);
                mri.networkActivity(volumeCounter+1:volumeCounter+volumeCount, :) = mriTimeCourse;
                volumeCounter = volumeCounter + volumeCount;
            end
            
            %get eeg labels
            labels = cell(0);
            for i = 1:length(eegFiles)
                fprintf('loading labels file %d of %d\n', i, length(eegFiles));
                eegFilename = fullfile(eegFolder, eegFiles(i).name);
                eeg = load(eegFilename);
                labels = unique([labels, eeg.convolvedData.labels]);
            end
            
            %load eeg data
            volumeCounter = 0;
            clear eeg;
            eeg.labels = labels;
            eeg.signals = NaN(size(mri.networkActivity, 1), length(labels));
            for i = 1:length(eegFiles)
                fprintf('organizing eeg data, file %d of %d\n', i, length(eegFiles));
                numberString = strrep(eegFiles(i).name, 'EEG_data_sub-', '');
                if(strcmp(numberString, eegFiles(i).name))
                    error('unexpected eeg filename: %s', eegFiles(i).name);
                end
                subjectString = numberString(1:2);
                runString = numberString(8:9);
                eegFilename = fullfile(eegFolder, eegFiles(i).name);
                thisEeg = load(eegFilename);
                volumeCount = size(thisEeg.convolvedData.signals, 2);
                for j = 1:size(thisEeg.convolvedData.signals, 1)
                    destInd = find(strcmp(eeg.labels, thisEeg.convolvedData.labels{j}));
                    eeg.signals(volumeCounter+1:volumeCounter + volumeCount, destInd) = thisEeg.convolvedData.signals(j, :)';
                end
                eeg.subjects(volumeCounter+1:volumeCounter+volumeCount) = str2num(subjectString);
                eeg.runs(volumeCounter+1:volumeCounter+volumeCount) = str2num(runString);
                eeg.times(volumeCounter+1:volumeCounter+volumeCount) = thisEeg.convolvedData.times;
                volumeCounter = volumeCounter + volumeCount;
            end
            
            save(combinedFilename, 'eeg', 'mri', '-v7.3');
            
            %debug
            eegnans = sum(isnan(eeg.signals));
            tab = tabulate(eegnans);
            tab(tab(:,2)==0,:) = [];
            
            volumenans = sum(isnan(eeg.signals), 2);
            tabV = tabulate(volumenans);
            tabV(tabV(:,2)==0,:) = [];
            %end debug
        end
        
        
        if(individualCoherence)
            % keep = cellfun(@length, strfind(trainingData.eegColumnLabels, 'ower')) > 0;
            oldKeep = false(size(trainingData.eegColumnLabels));
            
            outputData = trainingData.mriOutput';
            
            %add one eeg parameter at a time, the one with the highest test R value
            loopCounter = 1;
            while(loopCounter < length(oldKeep))
                maxR = 0;
                maxI = -1;
                clear maxTrial;
                for i = 1:length(oldKeep)
                    %         fprintf('%s: loop %d, training %d of %d (maxI = %d, maxR = %f)\n', char(datetime), loopCounter, i, length(oldKeep), maxI, maxR);
                    if(~oldKeep(i))
                        keep = oldKeep;
                        keep(i) = true;
                        inputData = trainingData.eegInput(:, keep)';
                        trial.r = neuralNetworkAdvancedGeneration(inputData, outputData);
                        score = trial.r.test + trial.r.validation + 0.1 * trial.r.training;
                        if(score > maxR)
                            trial.keep = find(keep);
                            trial.newIndex = i;
                            maxR = score;
                            maxI = i;
                            maxTrial = trial;
                        end
                    end
                end
                trials(loopCounter) = maxTrial;
                oldKeep = false(size(trainingData.eegColumnLabels));
                oldKeep(maxTrial.keep) = true;
                outFile = fullfile(outFolder, sprintf('sparseEegMriNeuralNet%d.mat', loopCounter));
                save(outFile, 'trials', '-v7.3');
                fprintf('%s: loop %d, maxR = %f, measure = %s\n', char(datetime), loopCounter, maxR, trainingData.eegColumnLabels{maxTrial.newIndex});
                loopCounter = loopCounter + 1;
            end
            
            inputData = trainingData.eegInput(:, keep)';
            outputData = trainingData.mriOutput';
            regressed = neuralNetworkAdvancedGeneration(inputData, outputData);
            
            rs = [trials.r];
            plot([[rs.training]', [rs.test]', [rs.validation]']);
            legend({'training', 'test', 'validation'});
        else % individualChannels == false
            distributionSize = 21;
            outputData = trainingData.mriOutput';
            
            %get unique channel names
            channelItemIndex = 1;
            measureItemIndex = 2;
            allChans = cell(0);
            for i = 1:length(trainingData.eegColumnLabels)
                label = trainingData.eegColumnLabels{i};
                items = strsplit(label, ' ');
                if(i == 1)
                    if(strcmp(items{1}, 'coherence') || (length(strfind(items{1}, 'ower')) > 0))
                        channelItemIndex = 2;
                        measureItemIndex = 1;
                    elseif(strcmp(items{2}, 'coherence') || (length(strfind(items{2}, 'ower')) > 0))
                        channelItemIndex = 1;
                        measureItemIndex = 2;
                    else
                        error('unhandled eeg label format: %s', label);
                    end
                end
                chans = strsplit(items{channelItemIndex}, '-');
                allChans = unique([allChans, chans]);
            end
            
            oldKeep = [];
            if(exist('savedData', 'var'))
                
                
                %if false, adds channels one at a time, and includes all coherence pairs
                individualCoherence = true;
                
                
                %pick up where we left off
                % outFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\MrEegLink';
                outFolder = '/home/data/EEG/processed/GhermanPhilastides/neuralnet';
                combinedFilename = '/home/data/EEG/processed/GhermanPhilastides/combinedEegMri.mat';
                
                if(exist(combinedFilename, 'file'))
                    eegFolder = '/home/data/EEG/processed/GhermanPhilastides/convolved';
                    mriFolder = '/home/data/EEG/processed/GhermanPhilastides/convolved';
                    eegFiles = dir(eegFolder);
                    eegFiles = eegFiles(cellfun(@length, strfind({eegFiles.name}, '.mat')) > 0);
                    
                    %load mri data
                    volumeCounter = 0;
                    for i = 1:length(eegFiles)
                        numberString = strrep(eegFiles(i).name, 'EEG_data_sub-', '');
                        if(strcmp(numberString, eegFiles(i).name))
                            error('unexpected eeg filename: %s', eegFiles(i).name);
                        end
                        subjectString = numberString(1:2);
                        runString = numberString(8:9);
                        mriFilename = sprintf('/home/data/subjects/sub-GhermanPhiliastides%s/infomap/main/groupnetworktimecourses.txt', subjectString);
                        mriTimeCourse = textread(mriFilename);
                        networkNumbers = mriTimeCourse(1,:);
                        mriTimeCourse = mriTimeCourse(2:end,:);
                        volumeCount = size(mriTimeCourse, 1);
                        if(~exist('mri', 'var'))
                            mri.networkNumbers = networkNumbers;
                        else
                            if(any(mri.networkNumbers ~= networkNumbers))
                                error('inconsistent network numbers across files at: %s', eegFilename);
                            end
                        end
                        mri.subjects(volumeCounter+1:volumeCounter+volumeCount) = str2num(subjectString);
                        %     mri.runs(volumeCounter+1:volumeCounter+volumeCount) = str2num(runString);
                        mri.networkActivity(volumeCounter+1:volumeCounter+volumeCount, :) = mriTimeCourse;
                        volumeCounter = volumeCounter + volumeCount;
                    end
                    
                    %get eeg labels
                    labels = cell(0);
                    for i = 1:length(eegFiles)
                        fprintf('loading labels file %d of %d\n', i, length(eegFiles));
                        eegFilename = fullfile(eegFolder, eegFiles(i).name);
                        eeg = load(eegFilename);
                        labels = unique([labels, eeg.convolvedData.labels]);
                    end
                    
                    %load eeg data
                    volumeCounter = 0;
                    clear eeg;
                    eeg.labels = labels;
                    eeg.signals = NaN(size(mri.networkActivity, 1), length(labels));
                    for i = 1:length(eegFiles)
                        fprintf('organizing eeg data, file %d of %d\n', i, length(eegFiles));
                        numberString = strrep(eegFiles(i).name, 'EEG_data_sub-', '');
                        if(strcmp(numberString, eegFiles(i).name))
                            error('unexpected eeg filename: %s', eegFiles(i).name);
                        end
                        subjectString = numberString(1:2);
                        runString = numberString(8:9);
                        eegFilename = fullfile(eegFolder, eegFiles(i).name);
                        thisEeg = load(eegFilename);
                        volumeCount = size(thisEeg.convolvedData.signals, 2);
                        for j = 1:size(thisEeg.convolvedData.signals, 1)
                            destInd = find(strcmp(eeg.labels, thisEeg.convolvedData.labels{j}));
                            eeg.signals(volumeCounter+1:volumeCounter + volumeCount, destInd) = thisEeg.convolvedData.signals(j, :)';
                        end
                        eeg.subjects(volumeCounter+1:volumeCounter+volumeCount) = str2num(subjectString);
                        eeg.runs(volumeCounter+1:volumeCounter+volumeCount) = str2num(runString);
                        eeg.times(volumeCounter+1:volumeCounter+volumeCount) = thisEeg.convolvedData.times;
                        volumeCounter = volumeCounter + volumeCount;
                    end
                    
                    save(combinedFilename, 'eeg', 'mri', '-v7.3');
                    
                    %debug
                    eegnans = sum(isnan(eeg.signals));
                    tab = tabulate(eegnans);
                    tab(tab(:,2)==0,:) = [];
                    
                    volumenans = sum(isnan(eeg.signals), 2);
                    tabV = tabulate(volumenans);
                    tabV(tabV(:,2)==0,:) = [];
                    %end debug
                end
                
                
                if(individualCoherence)
                    % keep = cellfun(@length, strfind(trainingData.eegColumnLabels, 'ower')) > 0;
                    oldKeep = false(size(trainingData.eegColumnLabels));
                    
                    outputData = trainingData.mriOutput';
                    
                    %add one eeg parameter at a time, the one with the highest test R value
                    loopCounter = 1;
                    while(loopCounter < length(oldKeep))
                        maxR = 0;
                        maxI = -1;
                        clear maxTrial;
                        for i = 1:length(oldKeep)
                            %         fprintf('%s: loop %d, training %d of %d (maxI = %d, maxR = %f)\n', char(datetime), loopCounter, i, length(oldKeep), maxI, maxR);
                            if(~oldKeep(i))
                                keep = oldKeep;
                                keep(i) = true;
                                inputData = trainingData.eegInput(:, keep)';
                                trial.r = neuralNetworkAdvancedGeneration(inputData, outputData);
                                score = trial.r.test + trial.r.validation + 0.1 * trial.r.training;
                                if(score > maxR)
                                    trial.keep = find(keep);
                                    trial.newIndex = i;
                                    maxR = score;
                                    maxI = i;
                                    maxTrial = trial;
                                end
                            end
                        end
                        trials(loopCounter) = maxTrial;
                        oldKeep = false(size(trainingData.eegColumnLabels));
                        oldKeep(maxTrial.keep) = true;
                        outFile = fullfile(outFolder, sprintf('sparseEegMriNeuralNet%d.mat', loopCounter));
                        save(outFile, 'trials', '-v7.3');
                        fprintf('%s: loop %d, maxR = %f, measure = %s\n', char(datetime), loopCounter, maxR, trainingData.eegColumnLabels{maxTrial.newIndex});
                        loopCounter = loopCounter + 1;
                    end
                    
                    inputData = trainingData.eegInput(:, keep)';
                    outputData = trainingData.mriOutput';
                    regressed = neuralNetworkAdvancedGeneration(inputData, outputData);
                    
                    rs = [trials.r];
                    plot([[rs.training]', [rs.test]', [rs.validation]']);
                    legend({'training', 'test', 'validation'});
                else % individualChannels == false
                    distributionSize = 21;
                    outputData = trainingData.mriOutput';
                    
                    %get unique channel names
                    channelItemIndex = 1;
                    measureItemIndex = 2;
                    allChans = cell(0);
                    for i = 1:length(trainingData.eegColumnLabels)
                        label = trainingData.eegColumnLabels{i};
                        items = strsplit(label, ' ');
                        if(i == 1)
                            if(strcmp(items{1}, 'coherence') || (length(strfind(items{1}, 'ower')) > 0))
                                channelItemIndex = 2;
                                measureItemIndex = 1;
                            elseif(strcmp(items{2}, 'coherence') || (length(strfind(items{2}, 'ower')) > 0))
                                channelItemIndex = 1;
                                measureItemIndex = 2;
                            else
                                error('unhandled eeg label format: %s', label);
                            end
                        end
                        chans = strsplit(items{channelItemIndex}, '-');
                        allChans = unique([allChans, chans]);
                    end
                    
                    %if false, adds channels one at a time, and includes all coherence pairs
                    individualCoherence = true;
                    
                    
                    %pick up where we left off
                    % outFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\MrEegLink';
                    outFolder = '/home/data/EEG/processed/GhermanPhilastides/neuralnet';
                    combinedFilename = '/home/data/EEG/processed/GhermanPhilastides/combinedEegMri.mat';
                    
                    if(exist(combinedFilename, 'file'))
                        eegFolder = '/home/data/EEG/processed/GhermanPhilastides/convolved';
                        mriFolder = '/home/data/EEG/processed/GhermanPhilastides/convolved';
                        eegFiles = dir(eegFolder);
                        eegFiles = eegFiles(cellfun(@length, strfind({eegFiles.name}, '.mat')) > 0);
                        
                        %load mri data
                        volumeCounter = 0;
                        for i = 1:length(eegFiles)
                            numberString = strrep(eegFiles(i).name, 'EEG_data_sub-', '');
                            if(strcmp(numberString, eegFiles(i).name))
                                error('unexpected eeg filename: %s', eegFiles(i).name);
                            end
                            subjectString = numberString(1:2);
                            runString = numberString(8:9);
                            mriFilename = sprintf('/home/data/subjects/sub-GhermanPhiliastides%s/infomap/main/groupnetworktimecourses.txt', subjectString);
                            mriTimeCourse = textread(mriFilename);
                            networkNumbers = mriTimeCourse(1,:);
                            mriTimeCourse = mriTimeCourse(2:end,:);
                            volumeCount = size(mriTimeCourse, 1);
                            if(~exist('mri', 'var'))
                                mri.networkNumbers = networkNumbers;
                            else
                                if(any(mri.networkNumbers ~= networkNumbers))
                                    error('inconsistent network numbers across files at: %s', eegFilename);
                                end
                            end
                            mri.subjects(volumeCounter+1:volumeCounter+volumeCount) = str2num(subjectString);
                            %     mri.runs(volumeCounter+1:volumeCounter+volumeCount) = str2num(runString);
                            mri.networkActivity(volumeCounter+1:volumeCounter+volumeCount, :) = mriTimeCourse;
                            volumeCounter = volumeCounter + volumeCount;
                        end
                        
                        %get eeg labels
                        labels = cell(0);
                        for i = 1:length(eegFiles)
                            fprintf('loading labels file %d of %d\n', i, length(eegFiles));
                            eegFilename = fullfile(eegFolder, eegFiles(i).name);
                            eeg = load(eegFilename);
                            labels = unique([labels, eeg.convolvedData.labels]);
                        end
                        
                        %load eeg data
                        volumeCounter = 0;
                        clear eeg;
                        eeg.labels = labels;
                        eeg.signals = NaN(size(mri.networkActivity, 1), length(labels));
                        for i = 1:length(eegFiles)
                            fprintf('organizing eeg data, file %d of %d\n', i, length(eegFiles));
                            numberString = strrep(eegFiles(i).name, 'EEG_data_sub-', '');
                            if(strcmp(numberString, eegFiles(i).name))
                                error('unexpected eeg filename: %s', eegFiles(i).name);
                            end
                            subjectString = numberString(1:2);
                            runString = numberString(8:9);
                            eegFilename = fullfile(eegFolder, eegFiles(i).name);
                            thisEeg = load(eegFilename);
                            volumeCount = size(thisEeg.convolvedData.signals, 2);
                            for j = 1:size(thisEeg.convolvedData.signals, 1)
                                destInd = find(strcmp(eeg.labels, thisEeg.convolvedData.labels{j}));
                                eeg.signals(volumeCounter+1:volumeCounter + volumeCount, destInd) = thisEeg.convolvedData.signals(j, :)';
                            end
                            eeg.subjects(volumeCounter+1:volumeCounter+volumeCount) = str2num(subjectString);
                            eeg.runs(volumeCounter+1:volumeCounter+volumeCount) = str2num(runString);
                            eeg.times(volumeCounter+1:volumeCounter+volumeCount) = thisEeg.convolvedData.times;
                            volumeCounter = volumeCounter + volumeCount;
                        end
                        
                        save(combinedFilename, 'eeg', 'mri', '-v7.3');
                        
                        %debug
                        eegnans = sum(isnan(eeg.signals));
                        tab = tabulate(eegnans);
                        tab(tab(:,2)==0,:) = [];
                        
                        volumenans = sum(isnan(eeg.signals), 2);
                        tabV = tabulate(volumenans);
                        tabV(tabV(:,2)==0,:) = [];
                        %end debug
                    end
                    
                    
                    if(individualCoherence)
                        % keep = cellfun(@length, strfind(trainingData.eegColumnLabels, 'ower')) > 0;
                        oldKeep = false(size(trainingData.eegColumnLabels));
                        
                        outputData = trainingData.mriOutput';
                        
                        %add one eeg parameter at a time, the one with the highest test R value
                        loopCounter = 1;
                        while(loopCounter < length(oldKeep))
                            maxR = 0;
                            maxI = -1;
                            clear maxTrial;
                            for i = 1:length(oldKeep)
                                %         fprintf('%s: loop %d, training %d of %d (maxI = %d, maxR = %f)\n', char(datetime), loopCounter, i, length(oldKeep), maxI, maxR);
                                if(~oldKeep(i))
                                    keep = oldKeep;
                                    keep(i) = true;
                                    inputData = trainingData.eegInput(:, keep)';
                                    trial.r = neuralNetworkAdvancedGeneration(inputData, outputData);
                                    score = trial.r.test + trial.r.validation + 0.1 * trial.r.training;
                                    if(score > maxR)
                                        trial.keep = find(keep);
                                        trial.newIndex = i;
                                        maxR = score;
                                        maxI = i;
                                        maxTrial = trial;
                                    end
                                end
                            end
                            trials(loopCounter) = maxTrial;
                            oldKeep = false(size(trainingData.eegColumnLabels));
                            oldKeep(maxTrial.keep) = true;
                            outFile = fullfile(outFolder, sprintf('sparseEegMriNeuralNet%d.mat', loopCounter));
                            save(outFile, 'trials', '-v7.3');
                            fprintf('%s: loop %d, maxR = %f, measure = %s\n', char(datetime), loopCounter, maxR, trainingData.eegColumnLabels{maxTrial.newIndex});
                            loopCounter = loopCounter + 1;
                        end
                        
                        inputData = trainingData.eegInput(:, keep)';
                        outputData = trainingData.mriOutput';
                        regressed = neuralNetworkAdvancedGeneration(inputData, outputData);
                        
                        rs = [trials.r];
                        plot([[rs.training]', [rs.test]', [rs.validation]']);
                        legend({'training', 'test', 'validation'});
                    else % individualChannels == false
                        distributionSize = 21;
                        outputData = trainingData.mriOutput';
                        
                        %get unique channel names
                        channelItemIndex = 1;
                        measureItemIndex = 2;
                        allChans = cell(0);
                        for i = 1:length(trainingData.eegColumnLabels)
                            label = trainingData.eegColumnLabels{i};
                            items = strsplit(label, ' ');
                            if(i == 1)
                                
                                
                                %if false, adds channels one at a time, and includes all coherence pairs
                                individualCoherence = true;
                                
                                
                                %pick up where we left off
                                % outFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\MrEegLink';
                                outFolder = '/home/data/EEG/processed/GhermanPhilastides/neuralnet';
                                combinedFilename = '/home/data/EEG/processed/GhermanPhilastides/combinedEegMri.mat';
                                
                                if(exist(combinedFilename, 'file'))
                                    eegFolder = '/home/data/EEG/processed/GhermanPhilastides/convolved';
                                    mriFolder = '/home/data/EEG/processed/GhermanPhilastides/convolved';
                                    eegFiles = dir(eegFolder);
                                    eegFiles = eegFiles(cellfun(@length, strfind({eegFiles.name}, '.mat')) > 0);
                                    
                                    %load mri data
                                    volumeCounter = 0;
                                    for i = 1:length(eegFiles)
                                        numberString = strrep(eegFiles(i).name, 'EEG_data_sub-', '');
                                        if(strcmp(numberString, eegFiles(i).name))
                                            error('unexpected eeg filename: %s', eegFiles(i).name);
                                        end
                                        subjectString = numberString(1:2);
                                        runString = numberString(8:9);
                                        mriFilename = sprintf('/home/data/subjects/sub-GhermanPhiliastides%s/infomap/main/groupnetworktimecourses.txt', subjectString);
                                        mriTimeCourse = textread(mriFilename);
                                        networkNumbers = mriTimeCourse(1,:);
                                        mriTimeCourse = mriTimeCourse(2:end,:);
                                        volumeCount = size(mriTimeCourse, 1);
                                        if(~exist('mri', 'var'))
                                            mri.networkNumbers = networkNumbers;
                                        else
                                            if(any(mri.networkNumbers ~= networkNumbers))
                                                error('inconsistent network numbers across files at: %s', eegFilename);
                                            end
                                        end
                                        mri.subjects(volumeCounter+1:volumeCounter+volumeCount) = str2num(subjectString);
                                        %     mri.runs(volumeCounter+1:volumeCounter+volumeCount) = str2num(runString);
                                        mri.networkActivity(volumeCounter+1:volumeCounter+volumeCount, :) = mriTimeCourse;
                                        volumeCounter = volumeCounter + volumeCount;
                                    end
                                    
                                    %get eeg labels
                                    labels = cell(0);
                                    for i = 1:length(eegFiles)
                                        fprintf('loading labels file %d of %d\n', i, length(eegFiles));
                                        eegFilename = fullfile(eegFolder, eegFiles(i).name);
                                        eeg = load(eegFilename);
                                        labels = unique([labels, eeg.convolvedData.labels]);
                                    end
                                    
                                    %load eeg data
                                    volumeCounter = 0;
                                    clear eeg;
                                    eeg.labels = labels;
                                    eeg.signals = NaN(size(mri.networkActivity, 1), length(labels));
                                    for i = 1:length(eegFiles)
                                        fprintf('organizing eeg data, file %d of %d\n', i, length(eegFiles));
                                        numberString = strrep(eegFiles(i).name, 'EEG_data_sub-', '');
                                        if(strcmp(numberString, eegFiles(i).name))
                                            error('unexpected eeg filename: %s', eegFiles(i).name);
                                        end
                                        subjectString = numberString(1:2);
                                        runString = numberString(8:9);
                                        eegFilename = fullfile(eegFolder, eegFiles(i).name);
                                        thisEeg = load(eegFilename);
                                        volumeCount = size(thisEeg.convolvedData.signals, 2);
                                        for j = 1:size(thisEeg.convolvedData.signals, 1)
                                            destInd = find(strcmp(eeg.labels, thisEeg.convolvedData.labels{j}));
                                            eeg.signals(volumeCounter+1:volumeCounter + volumeCount, destInd) = thisEeg.convolvedData.signals(j, :)';
                                        end
                                        eeg.subjects(volumeCounter+1:volumeCounter+volumeCount) = str2num(subjectString);
                                        eeg.runs(volumeCounter+1:volumeCounter+volumeCount) = str2num(runString);
                                        eeg.times(volumeCounter+1:volumeCounter+volumeCount) = thisEeg.convolvedData.times;
                                        volumeCounter = volumeCounter + volumeCount;
                                    end
                                    
                                    save(combinedFilename, 'eeg', 'mri', '-v7.3');
                                    
                                    %debug
                                    eegnans = sum(isnan(eeg.signals));
                                    tab = tabulate(eegnans);
                                    tab(tab(:,2)==0,:) = [];
                                    
                                    volumenans = sum(isnan(eeg.signals), 2);
                                    tabV = tabulate(volumenans);
                                    tabV(tabV(:,2)==0,:) = [];
                                    %end debug
                                end
                                
                                
                                if(individualCoherence)
                                    % keep = cellfun(@length, strfind(trainingData.eegColumnLabels, 'ower')) > 0;
                                    oldKeep = false(size(trainingData.eegColumnLabels));
                                    
                                    outputData = trainingData.mriOutput';
                                    
                                    %add one eeg parameter at a time, the one with the highest test R value
                                    loopCounter = 1;
                                    while(loopCounter < length(oldKeep))
                                        maxR = 0;
                                        maxI = -1;
                                        clear maxTrial;
                                        for i = 1:length(oldKeep)
                                            %         fprintf('%s: loop %d, training %d of %d (maxI = %d, maxR = %f)\n', char(datetime), loopCounter, i, length(oldKeep), maxI, maxR);
                                            if(~oldKeep(i))
                                                keep = oldKeep;
                                                keep(i) = true;
                                                inputData = trainingData.eegInput(:, keep)';
                                                trial.r = neuralNetworkAdvancedGeneration(inputData, outputData);
                                                score = trial.r.test + trial.r.validation + 0.1 * trial.r.training;
                                                if(score > maxR)
                                                    trial.keep = find(keep);
                                                    trial.newIndex = i;
                                                    maxR = score;
                                                    maxI = i;
                                                    maxTrial = trial;
                                                end
                                            end
                                        end
                                        trials(loopCounter) = maxTrial;
                                        oldKeep = false(size(trainingData.eegColumnLabels));
                                        oldKeep(maxTrial.keep) = true;
                                        outFile = fullfile(outFolder, sprintf('sparseEegMriNeuralNet%d.mat', loopCounter));
                                        save(outFile, 'trials', '-v7.3');
                                        fprintf('%s: loop %d, maxR = %f, measure = %s\n', char(datetime), loopCounter, maxR, trainingData.eegColumnLabels{maxTrial.newIndex});
                                        loopCounter = loopCounter + 1;
                                    end
                                    
                                    inputData = trainingData.eegInput(:, keep)';
                                    outputData = trainingData.mriOutput';
                                    regressed = neuralNetworkAdvancedGeneration(inputData, outputData);
                                    
                                    rs = [trials.r];
                                    plot([[rs.training]', [rs.test]', [rs.validation]']);
                                    legend({'training', 'test', 'validation'});
                                else % individualChannels == false
                                    distributionSize = 21;
                                    outputData = trainingData.mriOutput';
                                    
                                    %get unique channel names
                                    channelItemIndex = 1;
                                    measureItemIndex = 2;
                                    allChans = cell(0);
                                    for i = 1:length(trainingData.eegColumnLabels)
                                        label = trainingData.eegColumnLabels{i};
                                        items = strsplit(label, ' ');
                                        if(i == 1)
                                            if(strcmp(items{1}, 'coherence') || (length(strfind(items{1}, 'ower')) > 0))
                                                channelItemIndex = 2;
                                                measureItemIndex = 1;
                                            elseif(strcmp(items{2}, 'coherence') || (length(strfind(items{2}, 'ower')) > 0))
                                                channelItemIndex = 1;
                                                measureItemIndex = 2;
                                            else
                                                error('unhandled eeg label format: %s', label);
                                            end
                                        end
                                        chans = strsplit(items{channelItemIndex}, '-');
                                        allChans = unique([allChans, chans]);
                                    end
                                    
                                    oldKeep = [];
                                    if(exist('savedData', 'var'))
                                        oldKeep = savedData.bestChan;
                                    end
                                    powerIndex = cellfun(@length, strfind(trainingData.eegColumnLabels, 'power')) > 0;
                                    coherenceIndex = cellfun(@length, strfind(trainingData.eegColumnLabels, 'coherence')) > 0;
                                    
                                    for h = (length(oldKeep)+1):length(allChans)
                                        bestScore = 0;
                                        clear best loops
                                        loopCounter = 1;
                                        for i = 1:length(allChans)
                                            keepChan = oldKeep;
                                            if(~any(keepChan == i))
                                                
                                                %get input eeg data that includes exclusively channels of
                                                %interest
                                                keepChan(end + 1) = i;
                                                keepMeasure = false(size(trainingData.eegColumnLabels));
                                                for j = 1:length(keepChan)
                                                    jLabel = allChans{keepChan(j)};
                                                    hasJ = cellfun(@length, strfind(trainingData.eegColumnLabels, jLabel)) > 0;
                                                    keepMeasure(powerIndex & hasJ) = true;
                                                    for k = (j+1):length(keepChan)
                                                        kLabel = allChans{keepChan(k)};
                                                        hasK = cellfun(@length, strfind(trainingData.eegColumnLabels, kLabel)) > 0;
                                                        keepMeasure(coherenceIndex & hasJ & hasK) = true;
                                                    end
                                                end
                                                inputData = trainingData.eegInput(:, keepMeasure)';
                                                
                                                %generate a distribution of trials
                                                clear regressed;
                                                for j = 1:distributionSize
                                                    fprintf('%s: input size %d fo %d, testing channel %d of %d, iteration %d of %d', char(datetime), h, length(allChans), i, length(allChans), j, distributionSize);
                                                    regressed(j) = neuralNetworkAdvancedGeneration(inputData, outputData);
                                                    quickScores(j) = regressed(j).training * .1 + regressed(j).test * .45 + regressed(j).validation * .45;
                                                    fprintf(' score: %0.3f\n', quickScores(j));
                                                end
                                                medScore = median(quickScores);
                                                if(medScore > bestScore)
                                                    bestChan = keepChan;
                                                    bestScore = medScore;
                                                    best.Score = bestScore;
                                                    best.Inputs = allChans(keepChan);
                                                    best.Measures = trainingData.eegColumnLabels(keepMeasure);
                                                    best.Regressed = regressed;
                                                end
                                                loop.score = medScore;
                                                loop.channel = allChans{i};
                                                loops(loopCounter) = loop;
                                                loopCounter = loopCounter + 1;
                                            end %if new channel not in oldKeep
                                        end %i loop
                                        oldKeep = bestChan;
                                        best.Inputs
                                        outputFilename = fullfile(outFolder, sprintf('networkDisribution%d.mat', h));
                                        save(outputFilename, 'best', 'loops', 'bestChan', '-v7.3');
                                    end %h loop
                                end
                                
                                if(false)
                                    %make a topoplot of the electrodes in order
                                    %topoplot(
                                end
                                
                                
                                
                                
                                if(strcmp(items{1}, 'coherence') || (length(strfind(items{1}, 'ower')) > 0))
                                    channelItemIndex = 2;
                                    measureItemIndex = 1;
                                elseif(strcmp(items{2}, 'coherence') || (length(strfind(items{2}, 'ower')) > 0))
                                    channelItemIndex = 1;
                                    measureItemIndex = 2;
                                else
                                    error('unhandled eeg label format: %s', label);
                                end
                            end
                            chans = strsplit(items{channelItemIndex}, '-');
                            allChans = unique([allChans, chans]);
                        end
                        
                        oldKeep = [];
                        if(exist('savedData', 'var'))
                            oldKeep = savedData.bestChan;
                        end
                        powerIndex = cellfun(@length, strfind(trainingData.eegColumnLabels, 'power')) > 0;
                        coherenceIndex = cellfun(@length, strfind(trainingData.eegColumnLabels, 'coherence')) > 0;
                        
                        for h = (length(oldKeep)+1):length(allChans)
                            bestScore = 0;
                            clear best loops
                            loopCounter = 1;
                            for i = 1:length(allChans)
                                keepChan = oldKeep;
                                if(~any(keepChan == i))
                                    
                                    %get input eeg data that includes exclusively channels of
                                    %interest
                                    keepChan(end + 1) = i;
                                    keepMeasure = false(size(trainingData.eegColumnLabels));
                                    for j = 1:length(keepChan)
                                        jLabel = allChans{keepChan(j)};
                                        hasJ = cellfun(@length, strfind(trainingData.eegColumnLabels, jLabel)) > 0;
                                        keepMeasure(powerIndex & hasJ) = true;
                                        for k = (j+1):length(keepChan)
                                            kLabel = allChans{keepChan(k)};
                                            hasK = cellfun(@length, strfind(trainingData.eegColumnLabels, kLabel)) > 0;
                                            keepMeasure(coherenceIndex & hasJ & hasK) = true;
                                        end
                                    end
                                    inputData = trainingData.eegInput(:, keepMeasure)';
                                    
                                    %generate a distribution of trials
                                    clear regressed;
                                    for j = 1:distributionSize
                                        fprintf('%s: input size %d fo %d, testing channel %d of %d, iteration %d of %d', char(datetime), h, length(allChans), i, length(allChans), j, distributionSize);
                                        regressed(j) = neuralNetworkAdvancedGeneration(inputData, outputData);
                                        quickScores(j) = regressed(j).training * .1 + regressed(j).test * .45 + regressed(j).validation * .45;
                                        fprintf(' score: %0.3f\n', quickScores(j));
                                    end
                                    medScore = median(quickScores);
                                    if(medScore > bestScore)
                                        bestChan = keepChan;
                                        bestScore = medScore;
                                        best.Score = bestScore;
                                        best.Inputs = allChans(keepChan);
                                        best.Measures = trainingData.eegColumnLabels(keepMeasure);
                                        best.Regressed = regressed;
                                    end
                                    loop.score = medScore;
                                    loop.channel = allChans{i};
                                    loops(loopCounter) = loop;
                                    loopCounter = loopCounter + 1;
                                end %if new channel not in oldKeep
                            end %i loop
                            oldKeep = bestChan;
                            best.Inputs
                            outputFilename = fullfile(outFolder, sprintf('networkDisribution%d.mat', h));
                            save(outputFilename, 'best', 'loops', 'bestChan', '-v7.3');
                        end %h loop
                    end
                    
                    if(false)
                        %make a topoplot of the electrodes in order
                        %topoplot(
                    end
                    
                    
                    
                    
                    
                    
                    oldKeep = [];
                    if(exist('savedData', 'var'))
                        oldKeep = savedData.bestChan;
                    end
                    powerIndex = cellfun(@length, strfind(trainingData.eegColumnLabels, 'power')) > 0;
                    coherenceIndex = cellfun(@length, strfind(trainingData.eegColumnLabels, 'coherence')) > 0;
                    
                    for h = (length(oldKeep)+1):length(allChans)
                        bestScore = 0;
                        clear best loops
                        loopCounter = 1;
                        for i = 1:length(allChans)
                            keepChan = oldKeep;
                            if(~any(keepChan == i))
                                
                                %get input eeg data that includes exclusively channels of
                                %interest
                                keepChan(end + 1) = i;
                                keepMeasure = false(size(trainingData.eegColumnLabels));
                                for j = 1:length(keepChan)
                                    jLabel = allChans{keepChan(j)};
                                    hasJ = cellfun(@length, strfind(trainingData.eegColumnLabels, jLabel)) > 0;
                                    keepMeasure(powerIndex & hasJ) = true;
                                    for k = (j+1):length(keepChan)
                                        kLabel = allChans{keepChan(k)};
                                        hasK = cellfun(@length, strfind(trainingData.eegColumnLabels, kLabel)) > 0;
                                        keepMeasure(coherenceIndex & hasJ & hasK) = true;
                                    end
                                end
                                inputData = trainingData.eegInput(:, keepMeasure)';
                                
                                %generate a distribution of trials
                                clear regressed;
                                for j = 1:distributionSize
                                    fprintf('%s: input size %d fo %d, testing channel %d of %d, iteration %d of %d', char(datetime), h, length(allChans), i, length(allChans), j, distributionSize);
                                    regressed(j) = neuralNetworkAdvancedGeneration(inputData, outputData);
                                    quickScores(j) = regressed(j).training * .1 + regressed(j).test * .45 + regressed(j).validation * .45;
                                    fprintf(' score: %0.3f\n', quickScores(j));
                                end
                                medScore = median(quickScores);
                                if(medScore > bestScore)
                                    bestChan = keepChan;
                                    bestScore = medScore;
                                    best.Score = bestScore;
                                    best.Inputs = allChans(keepChan);
                                    best.Measures = trainingData.eegColumnLabels(keepMeasure);
                                    best.Regressed = regressed;
                                end
                                loop.score = medScore;
                                loop.channel = allChans{i};
                                loops(loopCounter) = loop;
                                loopCounter = loopCounter + 1;
                            end %if new channel not in oldKeep
                        end %i loop
                        oldKeep = bestChan;
                        best.Inputs
                        outputFilename = fullfile(outFolder, sprintf('networkDisribution%d.mat', h));
                        save(outputFilename, 'best', 'loops', 'bestChan', '-v7.3');
                    end %h loop
                end
                
                if(false)
                    %make a topoplot of the electrodes in order
                    %topoplot(
                end
                
                
                
                
                ep = savedData.bestChan;
            end
            powerIndex = cellfun(@length, strfind(trainingData.eegColumnLabels, 'power')) > 0;
            coherenceIndex = cellfun(@length, strfind(trainingData.eegColumnLabels, 'coherence')) > 0;
            
            for h = (length(oldKeep)+1):length(allChans)
                bestScore = 0;
                clear best loops
                loopCounter = 1;
                for i = 1:length(allChans)
                    keepChan = oldKeep;
                    if(~any(keepChan == i))
                        
                        %get input eeg data that includes exclusively channels of
                        %interest
                        keepChan(end + 1) = i;
                        keepMeasure = false(size(trainingData.eegColumnLabels));
                        for j = 1:length(keepChan)
                            jLabel = allChans{keepChan(j)};
                            hasJ = cellfun(@length, strfind(trainingData.eegColumnLabels, jLabel)) > 0;
                            keepMeasure(powerIndex & hasJ) = true;
                            for k = (j+1):length(keepChan)
                                kLabel = allChans{keepChan(k)};
                                hasK = cellfun(@length, strfind(trainingData.eegColumnLabels, kLabel)) > 0;
                                keepMeasure(coherenceIndex & hasJ & hasK) = true;
                            end
                        end
                        inputData = trainingData.eegInput(:, keepMeasure)';
                        
                        %generate a distribution of trials
                        clear regressed;
                        for j = 1:distributionSize
                            fprintf('%s: input size %d fo %d, testing channel %d of %d, iteration %d of %d', char(datetime), h, length(allChans), i, length(allChans), j, distributionSize);
                            regressed(j) = neuralNetworkAdvancedGeneration(inputData, outputData);
                            quickScores(j) = regressed(j).training * .1 + regressed(j).test * .45 + regressed(j).validation * .45;
                            fprintf(' score: %0.3f\n', quickScores(j));
                        end
                        medScore = median(quickScores);
                        if(medScore > bestScore)
                            bestChan = keepChan;
                            bestScore = medScore;
                            best.Score = bestScore;
                            best.Inputs = allChans(keepChan);
                            best.Measures = trainingData.eegColumnLabels(keepMeasure);
                            best.Regressed = regressed;
                        end
                        loop.score = medScore;
                        loop.channel = allChans{i};
                        loops(loopCounter) = loop;
                        loopCounter = loopCounter + 1;
                    end %if new channel not in oldKeep
                end %i loop
                oldKeep = bestChan;
                best.Inputs
                outputFilename = fullfile(outFolder, sprintf('networkDisribution%d.mat', h));
                save(outputFilename, 'best', 'loops', 'bestChan', '-v7.3');
            end %h loop
        end
        
        if(false)
            %make a topoplot of the electrodes in order
            %topoplot(
        end
        
        
        
        
        
        chans = strsplit(items{channelItemIndex}, '-');
        allChans = unique([allChans, chans]);
    end
    
    oldKeep = [];
    if(exist('savedData', 'var'))
        oldKeep = savedData.bestChan;
    end
    powerIndex = cellfun(@length, strfind(trainingData.eegColumnLabels, 'power')) > 0;
    coherenceIndex = cellfun(@length, strfind(trainingData.eegColumnLabels, 'coherence')) > 0;
    
    for h = (length(oldKeep)+1):length(allChans)
        bestScore = 0;
        clear best loops
        loopCounter = 1;
        for i = 1:length(allChans)
            keepChan = oldKeep;
            if(~any(keepChan == i))
                
                %get input eeg data that includes exclusively channels of
                %interest
                keepChan(end + 1) = i;
                keepMeasure = false(size(trainingData.eegColumnLabels));
                for j = 1:length(keepChan)
                    jLabel = allChans{keepChan(j)};
                    hasJ = cellfun(@length, strfind(trainingData.eegColumnLabels, jLabel)) > 0;
                    keepMeasure(powerIndex & hasJ) = true;
                    for k = (j+1):length(keepChan)
                        kLabel = allChans{keepChan(k)};
                        hasK = cellfun(@length, strfind(trainingData.eegColumnLabels, kLabel)) > 0;
                        keepMeasure(coherenceIndex & hasJ & hasK) = true;
                    end
                end
                inputData = trainingData.eegInput(:, keepMeasure)';
                
                %generate a distribution of trials
                clear regressed;
                for j = 1:distributionSize
                    fprintf('%s: input size %d fo %d, testing channel %d of %d, iteration %d of %d', char(datetime), h, length(allChans), i, length(allChans), j, distributionSize);
                    regressed(j) = neuralNetworkAdvancedGeneration(inputData, outputData);
                    quickScores(j) = regressed(j).training * .1 + regressed(j).test * .45 + regressed(j).validation * .45;
                    fprintf(' score: %0.3f\n', quickScores(j));
                end
                medScore = median(quickScores);
                if(medScore > bestScore)
                    bestChan = keepChan;
                    bestScore = medScore;
                    best.Score = bestScore;
                    best.Inputs = allChans(keepChan);
                    best.Measures = trainingData.eegColumnLabels(keepMeasure);
                    best.Regressed = regressed;
                end
                loop.score = medScore;
                loop.channel = allChans{i};
                loops(loopCounter) = loop;
                loopCounter = loopCounter + 1;
            end %if new channel not in oldKeep
        end %i loop
        oldKeep = bestChan;
        best.Inputs
        outputFilename = fullfile(outFolder, sprintf('networkDisribution%d.mat', h));
        save(outputFilename, 'best', 'loops', 'bestChan', '-v7.3');
    end %h loop
end

if(false)
    %make a topoplot of the electrodes in order
    %topoplot(
end





