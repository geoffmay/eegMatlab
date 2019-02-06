function [ icaCoherence, summary ] = deriveIcaCoherenceMatrix( eegDataFilename, downsampleRate , icaWeights )
%DERIVEICACOHERENCEMATRIX describes coherence between EEG ica components
%   Computes ica components for an EEG signal, then computes a
%   timepoint*coherence matrix, downsampled to 128Hz.  The frequency bands
%   used in coherence computation are defined in the 'allChannelCoherence'
%   function.
%
%   INPUTS:
%
%   eegDataFilename can be a char array that leads to a ROBI-style eeg file,
%   a cell array containing multiple char arrays pointing to files (in
%   which case ica will be run on the concatenated files), or an
%   eeglab-style struct that contains data in a channel*timepoint matrix,
%   srate, etc.  All ROBI-style files are a binary formatted list of
%   doubles, 34 at a time, with channel labels defined by the
%   'antChannelLocs' function, sampled at 2048 samples per second.  This
%   function rereferences the data to linked mastoids, then the mastoids 
%   and "counter" channel are dropped.  
%
%   downsampleRate, if defined, will resample the data to the new rate in
%   Hz.  This is useful for very large files or lists of files that add up
%   to a large sum.  The resulting (weight * sphere) matrix is multiplied
%   by the original data to obtain "high fidelity" ica component
%   activations.  Entering 0 for this value will skip the downSample step
%   and allow the user to pass the third variable.
%
%   icaWeights, if given, bypasses the ica computation step.  This is also
%   useful for using component weights from a 'baseline' eeg to process
%   subsequent files if desired.
%
%   OUTPUT: a struct with the following fields:
%
%   

doCpp = false;

if(~exist('runica', 'file'))
  eeglab;
end
if(exist('downsampleRate', 'var'))
  if(downsampleRate == 0)
    clear downsampleRate;
  end
end
if(~exist('eegDataFilename', 'var'))
  error('must supply eegDataFilename');
  eegDataFilename = '/home/data/EEG/data/ROBI/ROBI_003/baseline eyes open/630158995692243270.eegData';
  eegDataFilename = '/home/data/EEG/data/RobiPilot/RA/baseline eyes open/630551983032598584.eegData';
end

if(isstruct(eegDataFilename))
  keep = 1:31;
  data = eegDataFilename.data(keep,:);
  chanlocs = eegDataFilename.chanlocs(keep);
  sourceFileInfo.filename = eegDataFilename.filename;
  sourceFileInfo.dataLength = eegDataFilename.pnts;
else
  
  if(ischar(eegDataFilename))
    eegDataFilename = {eegDataFilename};
  end
  
  %load EEG channels, drop mastoids and counter
  for fileCounter = 1:length(eegDataFilename)
    data1 = loadRobiDataFile(eegDataFilename{fileCounter});
    sourceFileInfo.filename = eegDataFilename{fileCounter};
    sourceFileInfo.dataLength = size(data1,1);
    [~, chanlocs] = antChannelLocs;
    dropChannels = [13 19 34];
    chanlocs(dropChannels)=[];
    data1(:,dropChannels) = [];
    if(fileCounter == 1)
      data = data1';
    else
      data(:, end+1:end+size(data1,1)) = data1';
    end
  end
  
  %compute ica
  data = data .* 1000000;  %convert volts to microvolts; helps prevent early ica convergence
end

srate = 2048;
if(isstruct(eegDataFilename))
  srate = eegDataFilename.srate;
end


if(exist('downsampleRate', 'var'))  
  fprintf('\n(%s) downsampling from %d to %d Hz for ica', char(datetime), srate, downsampleRate);
  for i = 1:size(data,1)
    upSample = double(data(i, :));
    downSample = resample(upSample, downsampleRate, srate);
    data2(i,:) = downSample;
  end
  originalData = data;
  data = data2;
end


if(~exist('icaWeights', 'var'))
  fprintf('\n%s: start ica', char(datetime));
  [ica.weights,ica.sphere,ica.compvars,ica.bias,ica.signs,ica.lrates,ica.activations] = runica(data, 'maxSteps', 1024 * 16);
  sphWeights = ica.weights * ica.sphere;
  if(exist('downsampleRate', 'var'))
    newActivations1 = sphWeights * originalData;
    if(true)
      %verify
      x = double(newActivations1(1,:)');
      y = ica.activations(1,:)';
      xLen = length(x);
      yLen = length(y);
      while(xLen * yLen >= (2^31))
        xLen = floor(xLen/2);
        yLen = floor(yLen/2);
      end
      x1 = resample(x, yLen, xLen);
      if(length(x1) > length(y))
        x1(length(y)+1:end) = [];
      elseif(length(y) > length(x1))
        y(length(x1)+1:end) = [];
      end
      [rho, p] = corr(y, x1);
      if(rho < 0.99)
        icaCoherence.note = sprintf('rebuilt ica activation matrix does not match given activation matrix; rho = %f', rho);
      else
        icaCoherence.note = sprintf('');
      end
    end    
    ica.activations = newActivations1;
  end
  %copy to return structure
  for i = 1:size(ica.weights, 2)
    icaInfo.weights(:,i) = table(ica.weights(:,i));
    icaInfo.weights.Properties.VariableNames{i} = chanlocs(i).labels;
    icaInfo.spheredWeights(:,i) = table(sphWeights(:,i));
    icaInfo.spheredWeights.Properties.VariableNames{i} = chanlocs(i).labels;
  end
else
  %use weights (sphere * weights) if they are provided.
  ica.activations = icaWeights{:,:} * data;
  icaInfo.spheredWeights = icaWeights;
  sphWeights = icaWeights;
end
icaInfo.sourceFileInfo = sourceFileInfo;




%ensure all components are net positive (LORETA does this).
% flipCoefficientsFlag = ones(1,size(ica.weights,2));
flipActivations = ica.activations;
flipCoefficients = sphWeights;
flipSphereCoefficientsFlag = ones(1,size(sphWeights,2));
for i = 1:size(flipCoefficients,2)
  %   comp = flipCoefficients(i,:);
  %   avg = mean(comp);
  sphAvg = mean(sphWeights(i,:));
  %   if(avg < 0)
  %     flipCoefficientsFlag(i) = -1;
  %     comp = comp .* -1;
  %     ica.weights(i,:) = comp;
  %     %     ica.activations(i,:) = ica.activations(i,:) .* -1;
  %   end
  if(sphAvg < 0)
    flipSphereCoefficientsFlag(i) = -1;
    flipCoefficients(i,:) = flipCoefficients(i,:) .* -1;
    flipActivations(i,:) = flipActivations(i,:) .* -1;
  end
end

%compute coherence of ica components
clear EEG;
EEG.srate = srate;
EEG.nbchan = size(data,1);
EEG.data = ica.activations;
EEG.srate = srate;
EEG.nbchan = size(ica.weights,1);
for i = 1:EEG.nbchan
  EEG.chanlocs(i).labels = sprintf('icaComp%d', i);
end

if(doCpp)
  output = eegCppAnalytics(EEG);
  icaCoherence.matrix = output.data;
  icaCoherence.labels = output.measures;
else
  [ coh, x, pow, freqInfo ] = allChannelCoherence(EEG);
  [mat, labels] = convertCoherenceStructToMatrix(coh, freqInfo, pow);
  %   [mat, labels] = allChannelCoherence(EEG);
  icaCoherence.matrix = mat;
  icaCoherence.labels = labels;
%   icaCoherence.timePoints = timePoints;
end





%calculate approximate component locations
xs = [chanlocs.X];
xs = xs(1:31)';
ys = [chanlocs.Y];
ys = ys(1:31)';
zs = [chanlocs.Z];
zs = zs(1:31)';
for i = 1:size(ica.weights,2)
  comp = ica.weights(:,i);
  sphComp = sphWeights(:,i);
  xSum = 0;
  ySum = 0;
  zSum = 0;
  sphXSum = 0;
  sphYSum = 0;
  sphZSum = 0;
  weightSum = 0;
  sphereSum = 0;
  for j = 1:length(comp)
    xSum = xSum + xs(j) * comp(j);
    ySum = ySum + ys(j) * comp(j);
    zSum = zSum + zs(j) * comp(j);
    weightSum = weightSum + comp(j) * comp(j);
    sphXSum = sphXSum + xs(j) * sphComp(j);
    sphYSum = sphYSum + ys(j) * sphComp(j);
    sphZSum = sphZSum + zs(j) * sphComp(j);
    sphereSum = sphereSum + sphComp(j) * sphComp(j);
  end
  weightAvg = sqrt(weightSum);
  sphereAvg = sqrt(sphereSum);
  compLoc.X = xSum / weightAvg;
  compLoc.Y = ySum / weightAvg;
  compLoc.Z = zSum / weightAvg;
  sphereLoc.X = sphXSum / sphereAvg;
  sphereLoc.Y = sphYSum / sphereAvg;
  sphereLoc.Z = sphZSum / sphereAvg;
  compLocs(i) = compLoc;
  compSphLocs(i) = sphereLoc;
end

icaInfo.locations = compLocs;
icaInfo.spheredLocations = compSphLocs;

%compute the 'x' axis (minutes) for coherence timecourse
% timePointCount = size(mat, 1);
% timePoints = (1:timePointCount) ./ 128 ./ 60;

icaCoherence.icaInfo = icaInfo;

if(nargout > 1)  
  summary = summarizeCoherenceMatrix(icaCoherence);
end

%todo: make sure you can reconstruct (code is the start of the false block
%below)
dummy = 1;


if(false)
  
  %check that we can "rebuild" the activations from the data and weights
  rebuilt1 = sphWeights * data;
  given1 = ica.activations;
  rebuilt2 = flipCoefficients * data;
  given2 = flipActivations;
  clear rs ps rs1 ps1 rs2 ps2
  for i = 1:size(data,1)
    fprintf('.');
    for j = 1:size(data,1);
      %       [rs1(i,j), ps1(i,j)] = corr(rebuilt1(j,:)', given1(i,:)');
      [rs2(i,j), ps2(i,j)] = corr(rebuilt2(j,:)', given2(i,:)');
    end
  end
  imagesc(rs2); %this should be 1 on the diagonal and close to zero elsewhere
  colorbar;
  
  %check "manual matrix multiplication"
  rebuiltData = zeros(size(dataIca));
  errIca = NaN(size(dataIca));
  for sampleNumber = 1:size(dataIca, 1)
    fprintf('\n%f', sampleNumber / size(dataIca,1));
    for compNum = 1:size(sphWeights,1)
      for chanNum = 1:size(sphWeights, 2)
        rebuiltData(sampleNumber,compNum) = rebuiltData(sampleNumber,compNum) + data1(sampleNumber,chanNum) * sphWeights(compNum,chanNum);
      end
    end
    %debug
    errIca(sampleNumber,:) = rebuiltData(sampleNumber,:) .* 1000000 - rebuilt1(:,sampleNumber)';
    %end debug
  end
  
  
  %pairwise correlation of component correlations
  [cohCorrR, cohCorrP] = paircorr_mod(cohMat);
  
  
  
  
  
  %mark pairs that have more than two components
  rSquared = cohCorrR .* cohCorrR;
  for i = 1:size(rSquared,1)
    fprintf('\n%s: %d of %d', char(datetime), i, size(rSquared,1));
    pair1 = labels{i};
    items = strsplit(pair1, ' ');
    comps1 = strsplit(items{2}, '-');
    for j = 1:size(rSquared,2)
      pair2 = labels{j};
      items = strsplit(pair2, ' ');
      comps1(3:4) = strsplit(items{2}, '-');
      ucomp = unique(comps1);
      if(length(ucomp) > 2)
        moreThanTwoComponents(i,j) = 1;
      end
    end
  end
  
  %find components
  threeCompRSquared = cohCorrR .*cohCorrR;
  threeCompRSquared(moreThanTwoComponents == 0) = 0;
  [sort3CompRS, perm] = sortMatrix(threeCompRSquared);
  count = 1000;
  x = NaN(count, 1);
  y = NaN(count, 1);
  for i = 1:count
    fprintf('\n%s: %d of %d', char(datetime), i, count);
    threshold = i / count;
    x(i) = threshold;
    y(i) = sum(sum(threeCompRSquared > threshold));
  end
  
  threshold = 0.7;
  supra = (rSquared > threshold) & moreThanTwoComponents;
  supLab = cell(1,sum(sum(supra)));
  counter = 1;
  for i = 1:size(supra,1)
    for j = (i+1):size(supra,2)
      if(supra(i,j))
        supLab{counter} = sprintf('%s vs %s', labels{i}, labels{j});
        counter = counter + 1;
      end
    end
  end
  supLab(counter:end) = [];
  
  clear rs ps
  testCount = 3;
  a = rand(1000,testCount);
  [bica.weights,bica.sphere,bica.compvars,bica.bias,bica.signs,bica.lrates,bica.activations] = runica(a');
  
  c{1} = bica.activations' * bica.weights * bica.sphere;
  c{2} = (bica.weights * bica.sphere * bica.activations)';
  
  d{1} = a * bica.weights;
  d{2} = a * bica.weights * bica.sphere;
  d{3} = bica.activations';
  dbig = [a d{1} d{2} d{3}];
  
  %   plot([d{2}(:,1) d{3}(:,1)]);
  legend({'recalc', 'given'});
  clear rs ps;
  for i= 1:testCount
    for j = 1:testCount
      rs(i,j) = corr(d{2}(:,i), d{3}(:,j));
    end
  end
  rs
  
end

%  ****
%
%
% voltages = loadRobiDataFile(eegDataFilename);
% [ica.weights,ica.sphere,ica.compvars,ica.bias,ica.signs,ica.lrates,ica.activations] = runica(voltages');
% icaComps = ica.weights;
%
% EEG.data = loadRobiDataFile(eegDataFilename);
% EEG.data(:, 32:end) = [];
% EEG.data = EEG.data';
% EEG.srate = 2048;
% [~, EEG.chanlocs] = antChannelLocs;
% EEG.chanlocs(32:end) = [];
% EEG.nbchan = 31;
% coh = allChannelCoherence(EEG);
% mat = NaN(size(coh(1).coherence, 1), size(coh(1).coherence, 2) * length(coh));
% counter = 1;
% for i = 1:length(coh)
%     mat(:, counter:counter+size(coh(1).coherence, 2)-1) = coh(i).coherence;
%     counter = counter + size(coh(1).coherence, 2);
% end