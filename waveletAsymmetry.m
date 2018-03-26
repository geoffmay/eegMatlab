function [ output_args ] = waveletAsymmetry( input_args )
%WAVELETASYMMETRY Summary of this function goes here
%   Detailed explanation goes here

% folder = '/home/data/EEG/processed/Oregon/ERP/wavelet/timeFMixedEvents';
% files =dir(folder);

files=[];
fileLocs = cell(0);

folders = [{'/home/data/EEG/processed/Oregon/ERP/wavelet/timeFEvent4'}...
  {'/home/data/EEG/processed/Oregon/ERP/wavelet/timeFEvent129'}...
  {'/home/data/EEG/processed/Oregon/ERP/wavelet/timeFEvent130'}...
  {'/home/data/EEG/processed/Oregon/ERP/wavelet/timeFEvent131'}...
  {'/home/data/EEG/processed/Oregon/ERP/wavelet/timeFEvent132'}];

for i = 1:length(folders)
  files = [files; dir(folders{i})];
  files([files.isdir]) = [];
  fileLocs(end+1:length(files), 1) = folders(i)';
end


channel1 = [{'Fz'}];
channel2 = [{'F4'}, {'AF4'}, {'Fp1'}, {'Fp2'}, {'FC6'}, {'T8'}];
allChannels = [channel1 channel2];
summariesFilename = '/home/data/EEG/processed/waveletSummaries.mat';
psych = load('/home/data/EEG/data/Oregon/wahbehVariables.mat');

 if(~exist(summariesFilename))
%if(true)
  task = cell(1,length(files));
  for i = 1:length(files)
    if(length(strfind(files(i).name, 'ANT'))>0)
      task{i} = 'flanker';
    elseif(length(strfind(files(i).name, 'PM')) > 0)
      task{i} = 'PTSD';
    else
      task{i} = 'tones';
    end
  end
  flankerFiles = {files(find(strcmp(task, 'flanker'))).name};
  event44 = cellfun(@length, strfind(flankerFiles, 'Event44'));
  event209 = cellfun(@length, strfind(flankerFiles, 'Event209'));
  flanker44 = flankerFiles(find(event44));
  flanker209 = flankerFiles(find(event209));
  
  caps = [];
  flatAsym=[];
  flatPower=[];
  tic;
  for i = 1:length(files)
    if(strcmp(task{i}, 'flanker'))
      %a = load(fullfile(folder,files(i).name));
      a = load(fullfile(fileLocs{i},files(i).name));
      a = a.data;
      slashes = strfind(a.file, '/');
      target = a.file(slashes(end)+1:slashes(end)+5);
      channels = {a.channels.label};
      psychIndex = find(strcmp(psych.vetmindData{:,1}, target));
      cap = psych.vetmindData{psychIndex, 'CAPS'};
      ind1 = find(strcmp({a.channels.label}, 'Fz'));
      ind2 = find(strcmp({a.channels.label}, 'T8'));
      if(length(ind1) > 0 & length(ind2) > 0)
        timeF1 = a.channels(ind1).fullTimeF;
        timeF2 = a.channels(ind2).fullTimeF;
        if(length(timeF1.allErsp) > 5)
          
          %       asym = a.channels(ind2).timeF.ersp - a.channels(ind1).timeF.ersp;
          %       power = a.channels(ind2).timeF.ersp;
          asym = timeF2.ersp - timeF1.ersp;
          power = timeF2.ersp;
          if(length(cap) > 0)
            caps(end+1,:) = cap;
            flatAsym(end+1,:) = reshape(asym, 1, size(asym,1)*size(asym,2));
            flatPower(end+1,:) = reshape(power, 1, size(power,1)*size(power,2));
          end
        end
      end
    end
    fprintf('%d of %d\n',i, length(files));
  end
  save(summariesFilename, 'caps', 'flatAsym', 'flatPower');
  toc;
else
  load(summariesFilename);
end

[rho, p] = corr(caps, flatAsym);
[rrho, rp] = corr(caps, flatPower);
  keep = p < 0.05;
  rhomask = rho .* keep;
  pLog = -log10(p);
  rhoGraph = reshape(rhomask, 93, 200);
  pGraph = reshape(pLog, 93, 200);
  figure;
  heatmap(rhoGraph);
  title(sprintf('rho %s', channel2Title));
  figure;
  heatmap(pGraph);
  title(sprintf('p %s', channel2Title));



e1 = strcmp({summaries.electrode1}, 'Fz');
for comparator = 1:length(channel2);
  channel2Title = channel2{comparator};
  if(strcmp(channel2Title,'T4'))
    channel2Title = 'T8';
  end
  e2 = strcmp({summaries.electrode2}, channel2Title);
  flanker = strcmp(task, 'flanker');
  comparisons = find(e1 & e2 & flanker);
  caps = [];
  flatAsym = [];
  for i = 1:length(comparisons)
    summary = summaries(comparisons(i));
    targetFile = upper(summary.file1);
    psychIndex = find(strcmp(psych.vetmindData{:,1}, targetFile(1:5)));
    cap = psych.vetmindData{psychIndex, 'CAPS'};
    disp(i);
    if(i == 14)
      dummy = 0;
    end
    if(length(cap) > 0)
      caps(end+1,:) = cap;
      flatAsym(end+1,:) = reshape(summary.asym, 1, size(summary.asym,1)*size(summary.asym,2));
    end
  end
  [rho, p] = corr(caps, flatAsym);
  keep = p < 0.05;
  rhomask = rho .* keep;
  pLog = -log10(p);
  rhoGraph = reshape(rhomask, 24, 200);
  pGraph = reshape(pLog, 24, 200);
  figure;
  heatmap(rhoGraph);
  title(sprintf('rho %s', channel2Title));
  figure;
  heatmap(pGraph);
  title(sprintf('p %s', channel2Title));
end
tilefigs;
end


%     'FFT Amplitude AsymmetryFP1FzAlpha 1'
%     'FFT Amplitude AsymmetryT4FzAlpha'
%     'FFT Amplitude AsymmetryT4FzAlpha 2'
%     'FFT Amplitude AsymmetryFC6FzAlpha'
%     'FFT Amplitude AsymmetryFC6FzAlpha 1'
%     'FFT Amplitude AsymmetryF4FzAlpha'
%     'FFT Amplitude AsymmetryF4FzAlpha 1'
%     'Z Scored FFT Amplitude AsymmetryFP1FzAlpha 1'
%     'Z Scored FFT Amplitude AsymmetryT4FzAlpha'
%     'Z Scored FFT Amplitude AsymmetryF4FzAlpha'
%     'Z Scored FFT Amplitude AsymmetryF4FzAlpha 1'
%     'Z Scored FFT Amplitude AsymmetryFP2FzAlpha 1'
