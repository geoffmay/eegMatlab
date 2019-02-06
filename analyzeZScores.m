minPieceCount = 25;
piecesToSample = 4;

filename = '/Users/Geoff/Downloads/allSessionsZScores.mat';
filename = 'J:\zScoresCompletesummary\summaries.mat';
filename = 'J:\zScoresCompletesummary\summaries.mat';
if(~exist('data', 'var'))
    data = load(filename);
end

shams = {'ROBI_008', 'ROBI_009', 'ROBI_010', 'ROBI_011', 'ROBI_012', 'ROBI_016'};
isSham = false(1,length(data.summaries));
for i = 1:length(isSham)
    for j = 1:length(shams)
        if(length(strfind(data.summaries(i).filename, shams{j})) > 0)
            isSham(i) = true;
        end
    end
end

isTx = cellfun(@length, strfind({data.summaries.filename}, 'tx')) > 0;
isRobi = cellfun(@length, strfind({data.summaries.filename}, 'ROBI')) > 0;
labels = data.summaries(1).labels;
isCoherence = cellfun(@length, strfind(labels, 'coh-')) > 0;
isRelPow = cellfun(@length, strfind(labels, 'rel-')) > 0;


isGood = isTx & isRobi;
txInd = find(isTx & isRobi);

lastRealPiece = NaN(1, length(txInd));
hasEnoughPieces = false(1, length(txInd));
pieceCounts = NaN(1, length(txInd));
txFilenames = cell(1, length(txInd));

pieceChanges = NaN(length(txInd), 5700);
startPieces = NaN(length(txInd), 5700);
endPieces = NaN(length(txInd), 5700);
improvementFactors = NaN(length(txInd), 5700);

fileCounter = 1;
txCounter = 1;
while(fileCounter < length(data.summaries))
    summary = data.summaries(fileCounter);
    
    if(isGood(fileCounter))
        %get the index of the last "good" piece
        if(txCounter == 1)
            lastRealPiece(txCounter) = length(summary.pieces);
        else
            prevSummary = data.summaries(txInd(txCounter-1));
            if(txInd(txCounter) ~= fileCounter)
                error('mismatched counter');
            end
            j = length(summary.pieces);
            match = true;
            while(match)
                if(length(prevSummary.pieces) >= j)
                    prevPiece = [prevSummary.pieces(j).meanZScores];
                    thisPiece = [summary.pieces(j).meanZScores];
                    difference = thisPiece -  prevPiece;
                    if(any(difference))
                        match = false;
                        lastRealPiece(txCounter) = j;
                    end
                end
                j = j - 1;
            end
            
        end
        
        %get pieces from the start and end of the session
        if(lastRealPiece(txCounter) > minPieceCount)
            hasEnoughPieces(txCounter) = true;
            pieceCount = length(summary.pieces);
            pieceCounts(txCounter) = pieceCount;
            if(pieceCount > 10)
                startBits = summary.pieces(2:piecesToSample+1);
                endBits = summary.pieces(end-piecesToSample:end-1);
                startSum = zeros(5700,1);
                endSum = zeros(5700,1);
                for i = 1:piecesToSample
                    startSum = startSum + startBits(i).meanZScores(1:end-6);
                    endSum = endSum + endBits(i).meanZScores(1:end-6);
                end
                %                 startPiece = summary.pieces(3);
                %                 endPiece = summary.pieces(lastRealPiece(txCounter)-2);
                %                 startMean = [startPiece.meanZScores(1:end-6)];
                %                 endMean = [endPiece.meanZScores(1:end-6)];
                startMean = startSum ./ piecesToSample;
                endMean = endSum ./ piecesToSample;
                pieceChange = endMean - startMean;
                improvementFactor = NaN(size(pieceChange));
                for i = 1:length(improvementFactor)
                    improvementFactor(i) = -1 * startMean(i) * pieceChange(i);
                end
                improvementFactors(txCounter, :) = improvementFactor';
                pieceChanges(txCounter, :) = pieceChange';
                startPieces(txCounter, :) = startMean';
                endPieces(txCounter, :) = endMean';
                if(false)
                    pieceMean = [startPiece.meanZScores];
                    little = pieceMean < -10;
                end
                dummy = 1;
            end
        end
        txFilenames{txCounter} = summary.filename;
        txSham(txCounter) = isSham(fileCounter);
        txCounter = txCounter + 1;
    end
    fileCounter = fileCounter + 1;
end
if(false)
    imagesc(pieceChanges);
    colorbar;
end

pieceChanges(~hasEnoughPieces, :) = [];
startPieces(~hasEnoughPieces, :) = [];
endPieces(~hasEnoughPieces, :) = [];
verum.improvementFactors = improvementFactors(hasEnoughPieces & ~txSham, :);
sham.improvementFactors = improvementFactors(hasEnoughPieces & txSham, :);
verum.filenames = txFilenames(hasEnoughPieces & ~txSham);
sham.filenames = txFilenames(hasEnoughPieces & txSham);
verum.meanImprove = mean(verum.improvementFactors);
sham.meanImprove = mean(sham.improvementFactors);

improvementFactors(~hasEnoughPieces, :) = [];


% 
% meanChange = mean(pieceChanges);
meanImprovement = mean(improvementFactors);
meanVerumImprove = mean(verumImprove);
meanShamImprove = mean(shamImprove);

% figure;
% plot(pieceChanges(isCoherence));
% title('coherence');
% figure;
% plot(pieceChanges(isRelPow));
% title('relative power');
% figure;
% plot(meanImprovement(isCoherence));
% title('coherence improvement');
% figure;
% plot(meanImprovement(isRelPow));
% title('relative power improvement');

isInteresting = isCoherence | isRelPow;
% 
% figure;
% toPlot = improvementFactors(:,isInteresting);
% imagesc(toPlot, [-5, 5]);
% colorbar;
% title('improvementFactors');
% 
% figure;
% toPlot1 = startPieces(:,isInteresting);
% imagesc(toPlot1, [-5, 5]);
% colorbar;
% title('start Z score');
% 
% figure;
% toPlot2 = pieceChanges(:,isInteresting);
% imagesc(toPlot2, [-5, 5]);
% colorbar;
% title('Z score changes');
% 
% figure;
% toPlot3 = endPieces(:,isInteresting);
% imagesc(toPlot3, [-5, 5]);
% colorbar;
% title('end Z scores');
if(false)
    figure;
    plot(meanImprovement(isInteresting));
    title('improvement factors');
end
% 
% figure;
% meanChanges = mean(pieceChanges);
% startMeans = mean(startPieces);
% plot([meanChanges(isInteresting)', startMeans(isInteresting)']);
% legend({'change', 'start'});
% title('change vs start');


starts = verum.meanImprove(isInteresting);
% changes = meanChanges(isInteresting);
isVeryInteresting = starts ~= 0;

keepIndices = find(isInteresting);
keepIndices = keepIndices(isVeryInteresting);

verum.keep = keepIndices; 
sham.keep = keepIndices; 
difference = verum.meanImprove(verum.keep) - sham.meanImprove(sham.keep);

analysis.meanImprovementDiff = difference;
analysis.verum = verum;
analysis.sham = sham;
analysis.labels = labels(keepIndices);

figure;
plot(difference);
title('difference');

figure;
hold on;
plot(verum.meanImprove(verum.keep));
plot(sham.meanImprove(sham.keep));
legend('verum', 'sham');


iStarts = starts(isVeryInteresting);
iChanges = changes(isVeryInteresting);


meanVerumImprove1 = meanVerumImprove(isInteresting);
meanShamImprove1 = meanShamImprove(isInteresting);
meanVerumImprove2 = meanVerumImprove1(isVeryInteresting);
meanShamImprove2 = meanShamImprove1(isVeryInteresting);

figure;
plot(meanVerumImprove2);
title('verum');
figure;
plot(meanShamImprove2);
title('sham');

look1 = meanVerumImprove2 - meanShamImprove2;
figure;
plot(look1);
title('look1');

figure;
imagesc(meanVerumImprove2);
colorbar;
title('verum');

figure;
imagesc(meanShamImprove2);
colorbar;
title('sham');





[sortedStart, startSortInd] = sort(iStarts);
sortedChange = iChanges(startSortInd);
axis = zeros(size(sortedStart));
figure;
plot([sortedStart', sortedChange', axis']);
legend({'sortedStart', 'sortedChange'});


figure;
intImp = meanImprovement(isInteresting);
verIntImp = intImp(isVeryInteresting);
plot(verIntImp);
title('improvement factors');
