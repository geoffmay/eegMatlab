if(false)
    [pairs, inputFolder] = getNativeZScoreCompleters;
    outputFolder = '/home/psy/Documents/MATLAB/data/nativeZScore';
else
    inputFiles = dir('/media/psy/NIC_DRIVE01/zScoresComplete');    
    inputFiles([inputFiles.isdir]) = [];
    outputFolder = '/home/psy/Documents/MATLAB/data/nativeZScore2';
    pairs = {inputFiles.name};
    hasnotes = cellfun(@length, strfind(pairs, 'zScorenotes.txt')) > 0;
    pairs(hasnotes) = [];
end
 
for timepointCounter = 1:size(pairs,2)
    for subjectCounter = 1:size(pairs,1)
        fprintf('\ntime %d, subject %d', timepointCounter, subjectCounter);
        filePath = fullfile(inputFolder, pairs{subjectCounter, timepointCounter});
        zScore = loadNativeZScores(filePath);
        summary.filename = filePath;
        summary.labels = zScore.labels;
        summary.meanZScores = mean(zScore.data, 2);
        summary.stdZScores = std(zScore.data, 0, 2);
        summary.skewZScores = skewness(zScore.data, [], 2);
        summary.kurtZScores = kurtosis(zScore.data, [], 2);
        binSize = 128*15;
        binCounter = 1;
        for i = 1:binSize:size(zScore.data,2)
            piece.meanZScores = mean(zScore.data(i:i+binSize-1), 2);
            piece.stdZScores = std(zScore.data(i:i+binSize-1), 0, 2);
            piece.startSeconds = i / 128;
            piece.endSeconds = (i+binSize-1)/128;
            summary.pieces(binCounter) = piece;
            binCounter = binCounter + 1;
        end
        summaries(subjectCounter, timepointCounter) = summary;
        save(fullfile(outputFolder, 'summaries.mat'), 'summaries');
    end
end
 
 
labels =summaries(1).labels;
isAbs = cellfun(@length, strfind(labels, 'abs')) > 0;
isRat = cellfun(@length, strfind(labels, 'rat')) > 0;
 
 
for subjectCounter = 1:size(pairs,1)
    pre = summaries(subjectCounter, 1);
    post = summaries(subjectCounter, 2);
    eegMeasure = pre.labels;
    preZ = pre.meanZScores;
    postZ = post.meanZScores;
    preZ(isAbs | isRat) = [];
    postZ(isAbs | isRat) = [];
    abnormal = abs(preZ) > 2;
    improved = abs(preZ) > abs(postZ);
    worsened = abs(preZ) < abs(postZ);
    changeMagnitude = postZ - preZ ./ pre.stdZScores;
    remove = isnan(changeMagnitude);
    tab = table(eegMeasure, preZ, postZ, improved, worsened, changeMagnitude);
    tab(remove,:) = [];
    tables{subjectCounter} = tab;
end
 
save(fullfile(outputFolder,'tables.mat'), 'tables');