baselineFilename = '/home/data/EEG/data/RobiPilot/RA/baseline eyes open/63055198303000001.eegData';
treatment1Filename = '/home/data/EEG/data/RobiPilot/RA/tx 1/630564112186073937.eegData';
txNotesFilename = '/home/data/EEG/data/RobiPilot/RA/tx 1/630564112186073937Events.txt';

notes = loadInterpolatedZScore(txNotesFilename);

icaCoh = deriveIcaCoherenceMatrix(baselineFilename);
surfCoh = deriveRobiCoherenceMatrix(baselineFilename);


[icaCohPca.COEFF, icaCohPca.SCORE, icaCohPca.LATENT, icaCohPca.TSQUARED, icaCohPca.EXPLAINED, icaCohPca.MU] = pca(icaCoh.matrix);
[surfCohPca.COEFF, surfCohPca.SCORE, surfCohPca.LATENT, surfCohPca.TSQUARED, surfCohPca.EXPLAINED, surfCohPca.MU] = pca(icaCoh.matrix);


tic;
tx1 = loadRobiDataFile(treatment1Filename);
toc;

cutoff = 30 * 2048 * 60;
tx1(cutoff:end, :) = [];
[~,chanlocs] = antChannelLocs;
removeLocs = [13 19 34];
chanlocs(removeLocs) = [];
labels = {chanlocs.labels};

if(false)
  indices = 1:33;
  x = (1:size(tx1,1)) ./ (2048 * 60);
  plot(x, tx1(:,indices));
  legend(labels(indices));
end
eegTranspose = tx1';
eegTranspose(removeLocs,:) = [];
icaWeights = icaCoh.icaInfo.weights{:,:};

tx1Ic = icaWeights* eegTranspose;
EEG.data = tx1Ic;
EEG.srate = 2048;
EEG.chanlocs =chanlocs;
EEG.nbchan = 31;
[ channelPairs, x, channels, freqInfo ]  = allChannelCoherence(EEG);
[cohMatrix, cohLabels] = convertCoherenceStructToMatrix(channelPairs, freqInfo, channels);

cohIndices = 1:(31*30 / 2 * 5);
cohMatrix = cohMatrix(:, cohIndices)';
pcaWeights = icaCohPca.COEFF(cohIndices, 1)';
pcaLabels = icaCoh.labels(cohIndices);

%todo: plot coherence maps using ica locations


%tx pca score has bumps....
txPcaScore = pcaWeights * cohMatrix * 1000000;

if(false)
  order = 3;
  framelen = 2047;
  
  
  sgf = sgolayfilt(txPcaScore,order,framelen);
  
  x = (1:length(txPcaScore)) ./ (128 * 60);
  
  plot(x, [sgf']);
  scoreFft = abs(fft(txPcaScore));
  scoreFft = scoreFft(2:(length(scoreFft)/2 + 1));
  x2 = (1:length(scoreFft)) .* (64 / length(scoreFft));
  scoreFft = scoreFft .* x2; %still quite a bit of noise from 2 Hz and harmonics...
  scoreFft = scoreFft .* x2; %still quite a bit of noise from 2 Hz and harmonics...
  plot(x2, scoreFft)
  
end

%what if we recalculate from new ica components?


EEG1.data = eegTranspose;
EEG1.srate = 2048;
EEG1.chanlocs =chanlocs;
EEG1.nbchan = 31;
[ channelPairs1, x1, channels1, freqInfo1 ]  = allChannelCoherence(EEG1);
[cohMatrix1, cohLabels1] = convertCoherenceStructToMatrix(channelPairs1, freqInfo1, channels1);

cohMatrix1 = cohMatrix1(:, cohIndices)';
[icaCohPca1.COEFF, icaCohPca1.SCORE, icaCohPca1.LATENT, icaCohPca1.TSQUARED, icaCohPca1.EXPLAINED, icaCohPca1.MU] = pca(cohMatrix1);
pcaWeights1 = icaCohPca1.COEFF(cohIndices, 1)';
txPcaScore1 = pcaWeights1 * cohMatrix1 * 1000000;
scoreFft1 = abs(fft(txPcaScore1));
scoreFft1 = scoreFft1(2:(length(scoreFft1)/2 + 1));
x21 = (1:length(scoreFft1)) .* (64 / length(scoreFft1));
scoreFft1 = scoreFft1 .* x2; %still quite a bit of noise from 2 Hz and harmonics...
scoreFft1 = scoreFft1 .* x2; %still quite a bit of noise from 2 Hz and harmonics...

figure;
plot(x2, scoreFft1)


