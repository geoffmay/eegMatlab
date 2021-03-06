function [cohInfo, summary] = deriveRobiCoherenceMatrix(eegDataFilename)

doCpp = false;

if(isstruct(eegDataFilename))
  EEG = eegDataFilename;
  %   discardStart = 32;
  %   EEG.data(discardStart:end, :) = [];
  %   EEG.chanlocs(discardStart:end, :) = [];
  %   EEG.nbchan = discardStart - 1;
else
  EEG.data = loadRobiDataFile(eegDataFilename);
  EEG.data(:, 32:end) = [];
  EEG.data = EEG.data';
  EEG.srate = 2048;
  [~, EEG.chanlocs] = antChannelLocs;
  EEG.chanlocs(32:end) = [];
  EEG.nbchan = 31;
end
if(doCpp)
  output = eegCppAnalytics(EEG);
  cohInfo.matrix = output.data;
  cohInfo.labels = output.measures;
else
  [ coh, x, channels, freqInfo ] = allChannelCoherence2(EEG);
  [mat, labels] = convertCoherenceStructToMatrix(coh, freqInfo, channels);
%   [mat, labels] = allChannelCoherence3(EEG);
  cohInfo.timePoints = x;
  timePointCount = size(mat, 1);
  timePoints = (1:timePointCount) ./ 128 ./ 60;
  cohInfo.timePoints = timePoints;
  cohInfo.matrix = mat;
  cohInfo.labels = labels;
end

if(nargout > 1)
  summary = summarizeCoherenceMatrix(cohInfo);
end
