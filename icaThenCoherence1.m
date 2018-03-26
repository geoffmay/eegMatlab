
fullFile = 1;

filename = 'E:\eyetracker\AS\baseline eyes open\630491559160578485.eegData';

tic;

data1 = loadRobiDataFile(filename);
[~, chanlocs] = antChannelLocs;
dropChannels = [13 19];
chanlocs(dropChannels)=[];
dropChannels = [13 19 34];
data1(:,dropChannels) = [];
data = data1';
fprintf('\n%s: start ica', char(datetime));
tic;
data = data .* 1000000;

[cohIca.weights,cohIca.sphere,cohIca.compvars,cohIca.bias,cohIca.signs,cohIca.lrates,cohIca.activations] = runica(data, 'maxSteps', 1024 * 16);
toc;
fprintf('\n%s: finish ica', char(datetime));

save('AS runIca.mat', 'cohIca', '-v7.3');
clear EEG;
EEG.data = data;
EEG.srate = 2048;
EEG.nbchan = size(data,1);
EEG.chanlocs = chanlocs;
[ orig.coh, orig.x, orig.pow, orig.freqInfo ] = allChannelCoherence(EEG);

[orig.cohMat, orig.labels] = convertCoherenceStructToMatrix(orig.coh, orig.freqInfo);

[cohPca.coeff, cohPca.score, cohPca.latent, cohPca.tsquared, cohPca.explained] = pca(orig.cohMat);

clear EEG;
EEG.data = cohIca.activations;
EEG.srate = 2048;
EEG.nbchan = size(cohIca.weights,1);
for i = 1:EEG.nbchan
    EEG.chanlocs(i).labels = sprintf('icaComp%d', i);
end
[ coh, x, pow, freqInfo ] = allChannelCoherence(EEG);
[cohMat, labels] = convertCoherenceStructToMatrix(coh, freqInfo);

% allComps = cell(1, length(labels) * 2);
% counter = 1;
% for i = 1:length(labels)
%     items = strsplit(labels{i}, ' ');
%     comps = strsplit(items{2}, '-');
%     allComps(counter:counter+1) = comps;
%     counter = counter + 2;
% end
% comps =unique(allComps);

flipIcaMat = cohIca.weights;
for i = 1:size(flipIcaMat,2)
%     comp = flipIcaMat(:,i);
    comp = flipIcaMat(i,:);
    avg = mean(comp);
    if(avg < 0)
        comp = comp .* -1;
%         flipIcaMat(:,i) = comp;
        flipIcaMat(i,:) = comp;
    end
end

xs = [chanlocs.X];
xs = xs(1:31)';
ys = [chanlocs.Y];
ys = ys(1:31)';
zs = [chanlocs.Z];
zs = zs(1:31)';
for i = 1:size(flipIcaMat,2)
    comp = flipIcaMat(:,i);    
    xSum = 0;
    ySum = 0;
    zSum = 0;
    weightSum = 0;
    for j = 1:length(comp)
        xSum = xSum + xs(j) * comp(j);
        ySum = ySum + ys(j) * comp(j);
        zSum = zSum + zs(j) * comp(j);
        weightSum = weightSum + comp(j) * comp(j);
    end
    weightAvg = sqrt(weightSum);
    compLoc.X = xSum / weightAvg;
    compLoc.Y = ySum / weightAvg;
    compLoc.Z = zSum / weightAvg;
    compLocs(i) = compLoc;
end
figure;
scatter3([compLocs.X], [compLocs.Y], [compLocs.Z]);
xlabel('ant-post');
ylabel('med-lat');
zlabel('sup-inf');

close all;
for i = 1:length(compLocs)
    % i = 0;
    % i = i + 1;
    clust = ones(1, length(compLocs));
    clust(i) = 2;
    % gscatter([compLocs.X], [compLocs.Y], clust);
    figure;
    topoplot(flipIcaMat(:,i), chanlocs(1:31));
    fprintf('\ni = %d', i);
    title(sprintf('comp %d', i));
end
tilefigs

[icaCohPca.coeff, icaCohPca.score, icaCohPca.latent, icaCohPca.tsquared, icaCohPca.explained] = pca(cohMat);

[rho, p] = corr(icaCohPca.score(:,1), cohPca.score(:,1));

timePointCount = length(icaCohPca.score(:,1));
timePoints = (1:timePointCount) ./ 128 ./ 60;

icaFft = fft(icaCohPca.score(:,1));
surfFft = fft(cohPca.score(:,1));

% figure;
% plot([abs(icaFft), abs(surfFft)]);
cohLength = size(cohPca.score,1);
x = (0:(cohLength / 2)) ./ cohLength .* freqInfo.coherenceSampleRateHz;
toPlot = 2:(cohLength/2);

x = x(toPlot);
iF = abs(icaFft(toPlot));
sF = abs(surfFft(toPlot));

iF1 = NaN(1, length(iF));
sF1 = NaN(1, length(sF));
for i= 1:length(iF1)
%     iF1(i) = iF(i) * x(i) * x(i);
%     sF1(i) = sF(i) * x(i) * x(i);
    iF1(i) = iF(i) * x(i);
    sF1(i) = sF(i) * x(i);
end

figure;
hold on;
plot(x, sF1);
plot(x, iF1);
legend({'surface', 'ica'});

figure;
plot(timePoints, icaCohPca.score(:,1));

cohCorrR = ones(size(cohMat, 2));
cohCorrP = ones(size(cohMat, 2));
for i = 1:size(cohMat, 2)
    fprintf('\n%s: row %d of %d', char(datetime), i, size(cohMat, 2));
    for j = (i+1):size(cohMat, 2)
        [rho, p] = corr(cohMat(:, i), cohMat(:,j));
        cohCorrR(i,j) = rho;
        cohCorrR(j,i) = rho;
        cohCorrP(i,j) = p;
        cohCorrP(j,i) = p;        
    end
end

tic;
icaCohCorr = paircorr_mod(cohMat);
[sortMat, perm] = clusterMatrix(icaCohCorr);
permLab = labels(perm);
for i = 1:length(permLab)
  line = permLab{i};
  items = strsplit(line, ' ');
  freqs{i} = items{end};
end
uniqueFreqs = unique(freqs);
for i = 1:length(freqs)
  freqInd(i) = find(strcmp(uniqueFreqs, freqs(i)));
end
plot(freqInd);
set(gca, 'yticklabel', uniqueFreqs);
set(gca, 'ytick', 1:5);
ylim([0.5, 5.5]);
figure;



permLab(200:300)
toc;

toc;

beep
