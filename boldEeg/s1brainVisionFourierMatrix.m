function outputFolder = s1brainVisionFourierMatrix(edfFilename, fourierFolder)

% 
% locations = fileLocations;
% ind = strfind(edfFilename, 'EEG_data_sub-');
% if(length(ind) > 0)
%     numText = edfFilename((ind(end) + length('EEG_data_sub-')):end);
%     underscore = strfind(numText, '_');
%     numText(underscore:end) = [];
%     folder = fullfile(locations.ghermanFolder, ['sub-GhermanPhiliastides', numText], 'EEG');
%     eeg = loadGhermanEeg(fullfile(folder, edfFilename));
% else
%     locations = fileLocations();
%     folder = locations.brainVision;
%     eeg = loadBrainvisionEdf(fullfile(folder, edfFilename));
% end
% 
% 
% networkFilename = strrep(eeg.filename, '-edf.edf', '_network_timecourses.txt');
% networkActivity = textread(['C:\Vision\Raw Files\Geoff EEG test\mri\' networkFilename]);
% networkLabels = networkActivity(1, :);
% networkActivity = networkActivity(2:end, :);
% 
% targetNetworks = networkActivity(:, [1 3]);
[~, file] = fileparts(edfFilename);
%outputFolder = fullfile('/home/data/EEG/processed/boldEeg/fourier/', file);
outputFolder = fullfile(fourierFolder, file);



if(~exist(outputFolder, 'file'))
  eeg = s2loadBrainvisionEdf(edfFilename);
  if(false)
      maxIndex = 5000;
      eeg.times((maxIndex + 1):end) = [];
      eeg.pnts = maxIndex;
      eeg.data(:, (maxIndex + 1:end)) = [];
  end
  s3eegFourier(eeg, outputFolder);
%   fourierEeg.filename = edfFilename;
%   save(outputFolder, 'fourierEeg', '-v7.3');
end



if(false)
  [mat, labels] = s4convertCoherenceStructToMatrix(coh, freqInfo, channels);
  
  %   [mat, labels] = allChannelCoherence3(EEG);
  fourierEeg = eeg_emptyset;
  fourierEeg.times = x;
  % cohInfo.timePoints = x;
  % timePointCount = size(mat, 1);
  % timePoints = (1:timePointCount) ./ 128 ./ 60;
  % cohInfo.timePoints = timePoints;
  
  fourierEeg.data = mat;
  % cohInfo.matrix = mat;
  for i = 1:length(labels)
    fourierEeg.chanlocs(i).labels = labels{i};
  end
  % cohInfo.labels = labels;
  
  % s5saveSplitCohFile(cohInfo);
  
  
end