function [data, fileStarts] = loadBrainvisionDataFile(filename)

% filename = '/Volumes/NIC_DRIVE01/WD Backup.swstor/Neurofeedback/ZWFjZDNiYzU2NWI4NGU1Zm/Volume{550af545-a501-11e5-8846-806e6f6e6963}/Users/Neurofeedback/Documents/EEG/ROBI_013/outcome eyes open/630378445770810582.eegData';
% filename = '/Volumes/NIC_DRIVE01/WD Backup.swstor/Neurofeedback/ZWFjZDNiYzU2NWI4NGU1Zm/Volume{550af545-a501-11e5-8846-806e6f6e6963}/Users/Neurofeedback/Documents/EEG/ROBI_013/outcome impedance/630378441018603929.eegData';
% filename = '/Volumes/NIC_DRIVE01/WD Backup.swstor/Neurofeedback/ZWFjZDNiYzU2NWI4NGU1Zm/Volume{550af545-a501-11e5-8846-806e6f6e6963}/Users/Neurofeedback/Documents/EEG/ROBI_013/outcome eyes closed/630378452749420357.eegData';

fileStarts = 1;
if(iscell(filename))
  counter = 1;
  for i = 1:length(filename)
    fileStarts(i) = counter;
    thisData = loadRobiDataFileInternal(filename{i});
    data(counter:counter + size(thisData,1) - 1, :) = thisData;    
    counter = counter + size(thisData,1);
  end
else
  data = loadBrainvisionDataFileInternal(filename);
end


function data = loadBrainvisionDataFileInternal(filename)

rereference = true;
channelCount = 64;
file = dir(filename);
fileLength = file.bytes / 8;
sampleCount = fileLength / channelCount;
if(floor(sampleCount) ~= sampleCount)
  sampleCount = floor(sampleCount);
  fileLength = sampleCount * channelCount;
end
[folder, file, ext] = fileparts(filename);
fileTimestamp = str2num(file);
multiplier = 1;
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

data = reshape(contents, fileLength / channelCount, channelCount);
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
end

