function [eeg] = loadRobiEeg(filename)

% filename = '/Volumes/NIC_DRIVE01/WD Backup.swstor/Neurofeedback/ZWFjZDNiYzU2NWI4NGU1Zm/Volume{550af545-a501-11e5-8846-806e6f6e6963}/Users/Neurofeedback/Documents/EEG/ROBI_013/outcome eyes open/630378445770810582.eegData';
% filename = '/Volumes/NIC_DRIVE01/WD Backup.swstor/Neurofeedback/ZWFjZDNiYzU2NWI4NGU1Zm/Volume{550af545-a501-11e5-8846-806e6f6e6963}/Users/Neurofeedback/Documents/EEG/ROBI_013/outcome impedance/630378441018603929.eegData';
% filename = '/Volumes/NIC_DRIVE01/WD Backup.swstor/Neurofeedback/ZWFjZDNiYzU2NWI4NGU1Zm/Volume{550af545-a501-11e5-8846-806e6f6e6963}/Users/Neurofeedback/Documents/EEG/ROBI_013/outcome eyes closed/630378452749420357.eegData';

fileStarts = 1;
if(iscell(filename))
  counter = 1;
  for i = 1:length(filename)
    fileStarts(i) = counter;
    thisEeg = loadRobiDataFileInternal(filename{i});
    if(i == 1)
      eeg = thisEeg;
    end
    eeg.data(counter:counter + size(thisEeg.data,1) - 1, :) = thisEeg.data;    
    counter = counter + size(thisEeg.data,1);
  end
else
  eeg = loadRobiDataFileInternal(filename);
end


function eeg = loadRobiDataFileInternal(filename)

rereference = true;

watershed = 630438838651680050;
multiplier = 1000 * 1000 / 1024 / 1024;

channelCount = 34;
file = dir(filename);
fileLength = file.bytes / 8;
sampleCount = fileLength / channelCount;
if(floor(sampleCount) ~= sampleCount)
  sampleCount = floor(sampleCount);
  fileLength = sampleCount * channelCount;
end
[folder, file, ext] = fileparts(filename);
fileTimestamp = str2num(file);
if(fileTimestamp > watershed)
  multiplier = 1;
end
fileId = fopen(filename);
contents = fread(fileId, fileLength, 'double');
fclose(fileId);

%truncate if necessary
sampleCount = fileLength / channelCount;
if(sampleCount ~= floor(sampleCount))
  sampleCount = floor(sampleCount);
  fileLength = sampleCount * channelCount;
  contents = contents(1:fileLength);
end
labels = antChannelLocs;

data = reshape(contents, channelCount, fileLength / channelCount)';
clear contents;

if(rereference)
  cpzIndex = find(strcmp(labels,'CPz'));
  m1Index = find(strcmp(labels,'M1'));
  m2Index = find(strcmp(labels,'M2'));
  
  for i=1:size(data,1)
    sample = data(i,:);
    sample = sample .* multiplier;
    linkedMastoids = (sample(m1Index) + sample(m2Index)) / 2;
    %     avg = mean(sample(1:33));
    newSample = sample - linkedMastoids;
    %     avgSample = sample - avg;
    data(i,1:33) = newSample(1:33);
  end
  if(~exist('eeg_emptyset', 'file'))
    eeglab;
  end
  


end

data(:,34) = [];
[chanlabels, chanlocs] = antChannelLocs;
chanlocs(34) = [];

eeg = eeg_emptyset;
eeg.data = data';
eeg.srate = 2048;
eeg.chanlocs = chanlocs;
eeg.ref = '(M1 + M2) * 0.5';
eeg.pnts = size(data, 1);
eeg.times = (1:size(data, 1)) .* (1/eeg.srate);
eeg.filename = filename;
eeg.nbchan = size(data,2);





