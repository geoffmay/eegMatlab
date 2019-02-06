eegFolder = '/home/data/EEG/processed/Robi/coherence2';
mriFilename = '/home/data/EEG/processed/Robi/segregation/output.mat';

if(~exist('mriData','var'))
  mriData = load(mriFilename);
end
eegFiles = dir(eegFolder);

baseIndex = cellfun(@length, strfind({eegFiles.name}, 'baseline eyes open')) > 0;
exitIndex = cellfun(@length, strfind({eegFiles.name}, 'outcome eyes open')) > 0;
for i = 1:length(mriData.subjects)
  subjectName = mriData.subjects{i};
  subjectNumber = subjectName(5:end);
  subjectMatch = cellfun(@length, strfind({eegFiles.name}, ['ROBI_' subjectNumber])) > 0;
  
  exitEegFilename{i} = eegFiles(find(exitIndex & subjectMatch)).name;
  baseEegFilename{i} = eegFiles(find(baseIndex & subjectMatch)).name;
  mriData.exitEeg(i) = load(fullfile(eegFolder, exitEegFilename{i}));
  mriData.baseEeg(i) = load(fullfile(eegFolder, baseEegFilename{i}));
  
end

timepoints = {'pre-NF', 'post-NF', 'follow-up'};
for i = 1:length(mriData.subjects)
  fig = figure;
  set(gcf, 'Color', [1 1 1]);
  bar(mriData.FC(1:2,:,i));
  set(gca, 'xticklabel', (mriData.connections_names));
  legend(timepoints);
  title(mriData.subjects{i});
  filename = sprintf('/home/data/EEG/processed/Robi/segregation/plots/%d - %s', i, mriData.subjects{i});
  set(gca, 'fontsize', 24);
%   export_fig(filename, '-png');
  

end
tilefigs;

concave = [4 5 9 14];
convex = [3 8 18];


%plot power
powerInd = cellfun(@length, strfind(mriData.eegLabels, 'abs')) > 0;
freqLabels = {'1Hz-4Hz', '5Hz-8Hz', '9Hz-12Hz', '13Hz-24Hz', '25Hz-30Hz'};
[labs, locs] = antChannelLocs;
dropLocs = [find(cellfun(@length, strfind({locs.labels}, 'M')) > 0), 34];
locs(dropLocs) = [];
for subjNumber = 1:length(mriData.subjects)
  for freqNumber = 1:length(freqLabels)
    freqInd = cellfun(@length, strfind(mriData.eegLabels, freqLabels{freqNumber}));
    ind = find(freqInd & powerInd);
    %     chanLabels = mriData.eegLabels(ind);
    %     for i = 1:length(chanLabels)
    %       spaces = strfind(chanLabels{i}, ' ');
    %       chanLabs{i} = chanLabels{i}(spaces(1)+1:spaces(2)-1);
    %     end
    base = mriData.baseEeg(subjNumber).thisFile.means;
    exit = mriData.exitEeg(subjNumber).thisFile.means;
    values = exit(ind) - base(ind);
    topoplot(values, locs, 'maplimits', [-1.5 1.5]);
    set(gcf, 'Color', [1 1 1]);
    colorbar;
    filename = [mriData.subjects{subjNumber} ' ' freqLabels{freqNumber}];
    export_fig(filename, '-png');
    close all;
  end
end

i = 1;
while(i <=length(mriData.subjects))
  subject = mriData.subjects{i}
  diffEeg(i,:) = mriData.exitEeg(i).thisFile.means - mriData.baseEeg(i).thisFile.means;
%   plotCoherencePca(diffEeg,mriData.eegLabels);
%   title(subject);
   i = i + 1;
  
end

minDiff = min(diffEeg,[], 1)';
maxDiff = max(diffEeg, [],1)';
figure;
plot([minDiff, maxDiff]);

fileId = fopen('diffEeg.csv', 'w');
for i = 1:length(mriData.eegLabels)
  fprintf(fileId, '%s', mriData.eegLabels{i});
  if(i < length(mriData.eegLabels))
    fprintf(fileId, ',');
  end
end
for j = 1:size(diffEeg, 1)
  fprintf(fileId, '\r\n');
  for i = 1:size(diffEeg, 2)
    fprintf(fileId, '%f', diffEeg(j,i));
    if(i < size(diffEeg,2))
      fprintf(fileId, ',');
    end
  end
end
fclose(fileId);

