function output = checkReferencePhaseSlope(filename, doPlots)

if(~exist('doPlots','var'))
  doPlots = false;
end


channelCount = 34;

file = dir(filename);
fileLength = file.bytes / 8;
summary.filename = filename;

fileId = fopen(filename, 'r');
sampleCount = fileLength / channelCount;
if(sampleCount ~= floor(sampleCount))
  sampleCount = floor(sampleCount);
  fileLength = sampleCount * channelCount;
end
contents = fread(fileId, fileLength, 'double');
fclose(fileId);
data = reshape(contents, channelCount, fileLength / channelCount)';
reref = data;
avgref = data;
cpzEEG.nchan = 34;
cpzEEG.srate = 2048;
[labels, chanlocs] = antChannelLocs;
cpzEEG.chanlocs = chanlocs;
mastoidEEG = cpzEEG;

cpzIndex = find(strcmp(labels,'CPz'));
m1Index = find(strcmp(labels,'M1'));
m2Index = find(strcmp(labels,'M2'));
channelCounter = 1;
for i=1:size(data,1)
    sample = data(i,:);
    linkedMastoids = (sample(m1Index) + sample(m2Index)) / 2;
    avg = mean(sample(1:33));
    newSample = sample - linkedMastoids;
    avgSample = sample - avg;
    reref(i,1:33) = newSample(1:33);
    avgref(i,1:33) = avgSample(1:33);
end
% cpzEEG.data = data';
% mastoidEEG.data = reref';
% clear data reref;

mastoidSlope = zeros(33,33);
cpzSlope = zeros(33,33);
avgSlope = zeros(33,33);
mastoidPhasePlot = cell(0);
cpzPhasePlot = cell(0);
avgPhasePlot = cell(0);
pairCounter = 1;
for chan1 = 1:33
    for chan2 = chan1+1:33
        fprintf('.');
        if(mod(pairCounter,100) == 0)
            fprintf('\n%d',pairCounter);
        end
        if(pairCounter == 1)
            tic;
        end
        %         [cpzPsi, cpzPhasePlot{pairCounter}] = phaseSlopeIndex(cpzEEG, chan1, chan2, 1,50);
        [cpzPsi, cpzPhasePlot{pairCounter}] = phaseSlopeIndexOld(data, chan1, chan2);
        cpzSlope(chan1,chan2) = cpzPsi;
        cpzSlope(chan2,chan1) = -cpzPsi;
        %         [mastoidPsi, mastoidPhasePlot{pairCounter}] = phaseSlopeIndex(mastoidEEG, chan1, chan2, 1,50);
        [mastoidPsi, mastoidPhasePlot{pairCounter}] = phaseSlopeIndexOld(reref, chan1, chan2);
        mastoidSlope(chan1,chan2) = mastoidPsi;
        mastoidSlope(chan2,chan1) = -mastoidPsi;

        [avgPsi, avgPhasePlot{pairCounter}] = phaseSlopeIndexOld(data, chan1, chan2);
        avgSlope(chan1,chan2) = avgPsi;
        avgSlope(chan2,chan1) = -avgPsi;

        if(pairCounter == 1)
            toc;
        end
        
        pairCounter = pairCounter + 1;
    end
end
%debug
%end debug
mastoidSlope = mastoidSlope(1:33,1:33);
cpzSlope = cpzSlope(1:33,1:33);
avgSlope = avgSlope(1:33,1:33);
globalMin = min([min(min(mastoidSlope)),min(min(cpzSlope)),min(min(avgSlope))]);
globalMax = max([max(max(mastoidSlope)),max(max(cpzSlope)),max(max(avgSlope))]);
if(doPlots)
  figure;
  imagesc(mastoidSlope,[globalMin,globalMax]);
  title('mastoid');
  figure;
  imagesc(cpzSlope,[globalMin,globalMax]);
  title('cpz');
  figure;
  imagesc(avgSlope,[globalMin,globalMax]);
  title('avg');

  figure;
  difference = mastoidSlope - cpzSlope;
  imagesc(difference,[globalMin,globalMax]);
  title('cpz-mastoid difference');
  
end

cpzTopo = phaseSlopeTopography(cpzSlope);
avgTopo = phaseSlopeTopography(avgSlope);
mastoidTopo = phaseSlopeTopography(mastoidSlope);
[labels,chanlocs] = antChannelLocs;

if(doPlots)
  figure;
  topoplot(cpzTopo,chanlocs(1:32));
  title('cpz');
  
  figure;
  topoplot(mastoidTopo,chanlocs(1:32));
  title('mastoid');
  
  figure;
  topoplot(avgTopo,chanlocs(1:32));
  title('avg');
end

output.cpzSlope = cpzSlope;
output.mastoidSlope = mastoidSlope;
output.avgSlope = avgSlope;
output.cpzPhasePlot = cpzPhasePlot;
output.mastoidPhasePlot = mastoidPhasePlot;
output.avgPhasePlot = avgPhasePlot;


