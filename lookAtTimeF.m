function [ output_args ] = lookAtTimeF( input_args )
%LOOKATTIMEF Summary of this function goes here
%   Detailed explanation goes here

channelLabel = 'Fz';
ycutoff = 25;

% root = '/home/data/EEG/processed/Oregon/ERP/wavelet/';
% folders = dir(root);
% for i = length(folders):-1:1
%   if(folders(i).name(1) == '.')
%     folders(i) = [];
%   end
% end

eventOfInterest = 209;

root = '';
folders.name = sprintf('/home/data/EEG/processed/Oregon/ERP/wavelet/timeFEvent%d', eventOfInterest);

destination = sprintf('/home/data/EEG/processed/Oregon/ERP/event%dwavelet%s.mat', eventOfInterest, channelLabel);
if(~exist(destination))
  counter = 1;
  for folderIndex = 1:length(folders)
    folder = fullfile(root, folders(folderIndex).name);
    files = dir(folder);
    files([files.isdir])=[];
    for i = 1:length(files)
      load(fullfile(folder,files(i).name));
      channelIndex = find(strcmp({data.channels.label}, channelLabel));
      timef = data.channels(channelIndex).fullTimeF;
      ersp.freqs = timef.freqs;
      ersp.times = timef.times;
      ersp.ersp = timef.ersp;
      ersp.event = data.event;
      ersp.channel = channelLabel;
      ersp.sourceFile = data.file;
      ersp.eventCount = length(timef.allErsp);
      %     includeFreqs = timef.freqs <= ycutoff;
      %   figure;
      %   imagesc(timef.times, timef.freqs(includeFreqs), timef.ersp(includeFreqs,:));
      %   set(gca,'ydir','normal');
      allErsp(counter) = ersp;
      counter = counter + 1;
      fprintf('%d of %d folders, %d of %d files, %d total\n', folderIndex, length(folders), i, length(files), counter);
    end
    tic;save(destination, 'allErsp');toc;
  end
else
  load(destination);
end

bigMatrix = NaN(length(allErsp), numel(allErsp(1).ersp));
for i = 1:length(allErsp)
  bigMatrix(i,:) = reshape(allErsp(i).ersp, [1, numel(allErsp(i).ersp)]);
end

link = linkage(bigMatrix,'average','euclidean');
tic; Y = pdist(bigMatrix); toc;
c = cophenet(link, Y);
dendrogram(link, 0);

maxes = NaN(1,length(allErsp));
counts = NaN(1,length(allErsp));
for i = length(allErsp):-1:1
  newMax = max(max(allErsp(i).ersp));
  if(newMax < 1e10)
    maxes(i) = newMax;
    counts(i) = allErsp(i).eventCount;
  else
    maxes(i) = [];
    counts(i) = [];
  end
end
scatter(counts,log(maxes));

psych = load('/home/data/EEG/data/Oregon/wahbehVariables.mat');

ycutoff = 25;


flatErsp = zeros(size(psych.vetmindData,1), numel(allErsp(1).ersp));

for fileNumber = 1:length(allErsp)
    [~,sourceFile] = fileparts(allErsp(fileNumber).sourceFile);
    subjectId = sourceFile(1:5);
    psychIndex = find(strcmp(psych.vetmindData{:,1}, upper(subjectId)));
    if(length(psychIndex) > 0)
        matrix = allErsp(fileNumber).ersp;
        sigma = 5;
        size1 = 30;
        x = linspace(-size1/2,size1/2,size1);
        gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
        gaussFilter = gaussFilter / sum(gaussFilter);
        %debug
        for i1 = 1:size(matrix,1)
          matrix1(i1,:) = filtfilt(gaussFilter, 1, matrix(i1,:));
        end
        
%         figure;
%         imagesc(matrix);
%         figure;
%         imagesc(matrix1);
        %end debug
        flat = reshape(matrix, 1, numel(matrix));
        flatErsp(psychIndex,:) = flatErsp(psychIndex) + flat;
    end
end

remove = all(flatErsp == 0, 2);
flatErsp(remove,:) = [];
flatLogErsp = log(flatErsp);
psych.vetmindData(remove,:) = [];
caps = psych.vetmindData{:,'CAPS'};

[rho, p] = corr(flatErsp, caps);
[logRho, logP] = corr(flatLogErsp, caps);

minP = min(min(p));
minLogP = min(min(logP));

end
