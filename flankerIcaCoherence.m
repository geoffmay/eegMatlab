preFolder = '/home/data/EEG/data/Oregon/Flanker1';
postFolder = '/home/data/EEG/data/Oregon/Flanker2';

doAllFiles = false;
preFiles = dir(preFolder);
postFiles = dir(postFolder);
preFiles([preFiles.isdir]) = [];
postFiles([postFiles.isdir]) = [];

plotImage = false;
plotMeans = true;

fileCounter = 1;
while(fileCounter < length(preFiles))
  loadFiles = true;
  if(~doAllFiles)
    if(exist('coh', 'var'))
      loadFiles = false;
    end
  end
  %compute ica
  if(loadFiles)
    preFile = preFiles(fileCounter).name;
    preNumber = preFile(3:5);
    postIndex = find(cellfun(@length, strfind({postFiles.name}, preNumber)));
    if(length(postIndex) > 0)
      postFile = postFiles(postIndex).name;
    else
      fprintf('\nno match for %s', preFile);
    end
    
    pre = loadBdf(fullfile(preFolder, preFile));
    preEdit = pre;
    remove = 32:43;
    preEdit.data(remove,:) = [];
    preEdit.chanlocs(remove) = [];
    preEdit.nbchan = preEdit.nbchan - length(remove);
    coh = deriveIcaCoherenceMatrix(preEdit);
    summary = summarizeCoherenceMatrix(coh);
    outputFolder = '/home/data/EEG/processed/Oregon/flanker';
    save(fullfile(outputFolder, preFile), 'coh', 'summary');
  end
  % 112 no response
  % 113 response
  % 114 correct response
  % 115 incorrect response
  % 116 reponse within window
  % 117 response within correct response window
  % 118 response within incorrect response window
  
  %66 cue congruent
  %88 cue incongruent
  %77 nocue congruent
  %99 nocue incongruent
  
  % 1 following xx code with ~2500-3300 msec interval = fixation
  
  % nocue trials:
  % 1 = prestima + c fixation (3200 without)
  % 2 = prestimb (991)
  % 4 = stim (291)
  % 5 = postStim(2000)
  
  % cue trials:
  % 1 = prestima fixation (2476)
  % 2 = prestimb cue (137)
  % 3 = prestimc presentation (854)
  % 4 = stimulation (291)
  % 5 = postStim (2000)
  
  if(~doAllFiles & ~exist('preToPostChanges', 'var'))
    latency = [pre.event.latency]';
    code = [pre.event.type]';
    duration = [diff(latency); pre.pnts - max(latency)];
    preEvents = table(latency, code, duration);
    tab = tabulate(code);
    tab([tab(:,2)] < 2, :) = [];
    
    %arrange data from big matrix into smaller one for plotting
    congruentInd = [66 77];
    incongruentInd = [88 99];
    stimInd = 4;
    
    congruentTimes = latency(ismember(code, congruentInd)) ./ (preEdit.srate * 60);
    incongruentTimes = latency(ismember(code, incongruentInd)) ./ (preEdit.srate * 60);
    codeTimes = latency(ismember(code, [congruentInd incongruentInd])) ./ (preEdit.srate * 60);
    stimTimes = latency(code == stimInd) ./ (preEdit.srate * 60);
    
    preSeconds = 2.4;
    postSeconds = 2;
    preMinutes = preSeconds / 60;
    postMinutes = postSeconds / 60;
    cohSrate = length(coh.timePoints) / (max(coh.timePoints)*60);
    intervalSampleCount = floor((preSeconds + postSeconds) * cohSrate);
    numberOfMeasures = size(coh.matrix,2);
    congruentMat = NaN(intervalSampleCount, length(congruentTimes), numberOfMeasures);
    incongruentMat = NaN(intervalSampleCount, length(incongruentTimes), numberOfMeasures);
    congInd = 1;
    incongInd = 1;
    for i= 1:length(stimTimes)
      preTime = stimTimes(i) - preMinutes;
      postTime = stimTimes(i) + postMinutes;
      inSample = coh.timePoints > preTime & coh.timePoints < postTime;
      data = resample(coh.matrix(inSample, :), intervalSampleCount, sum(inSample));
      isCongruent = false;
      codeTime = codeTimes(max(find(codeTimes < stimTimes(i))));
      if(any(ismember(congruentTimes, codeTime)))
        isCongruent = true;
      end
      if(isCongruent)
        congruentMat(:, congInd, :) = data;
        congInd = congInd + 1;
      else
        incongruentMat(:,incongInd,:) = data;
        incongInd = incongInd + 1;
      end
    end
    if(size(congruentMat, 2) >= congInd)
      congruentMat(:,congInd:end,:) = [];
    end
    if(size(incongruentMat, 2) >= incongInd)
      incongruentMat(:,incongInd:end,:) = [];
    end
    
    %   for i = 1:length(congruentTimes)
    %     preTime = congruentTimes(i) - preMinutes;
    %     postTime = congruentTimes(i) + postMinutes;
    %     inSample = coh.timePoints > preTime & coh.timePoints < postTime;
    %     data = resample(coh.matrix(inSample, :), intervalSampleCount, sum(inSample));
    %     congruentMat(:, i, :) = data;
    %   end
    %   for i = 1:length(incongruentTimes)
    %     preTime = incongruentTimes(i) - preMinutes;
    %     postTime = incongruentTimes(i) + postMinutes;
    %     inSample = coh.timePoints > preTime & coh.timePoints < postTime;
    %     data = resample(coh.matrix(inSample, :), intervalSampleCount, sum(inSample));
    %     incongruentMat(:, i, :) = data;
    %   end
    x = -preSeconds:1/cohSrate:postSeconds;
    x = x(1:intervalSampleCount);
    
    %plot each (?) measure
    i = 1;
    while(i <= numberOfMeasures)'
      fprintf('\nmeasure %d of %d', i, numberOfMeasures);
      %sort by intensity at time = 0
      con = congruentMat(:,:,i);
      incon = incongruentMat(:,:,i);
      if(plotImage)
        zeroPoint = preSeconds * cohSrate;
        [~, conInd] = sort(con(zeroPoint,:));
        con(:, 1:size(con,2)) = con(:, conInd);
        [~, inconInd] = sort(incon(zeroPoint,:));
        incon(:, 1:size(incon,2)) = incon(:, inconInd);
        con = con'; incon = incon';
        buffer = repmat(mean(mean(con)), [10, intervalSampleCount]);
        toPlot = [con; buffer; incon];
        label = coh.labels{i};
        im = imagesc('XData', x, 'CData', toPlot);
        colorbar;
        title(label);
      end
      if(plotMeans)
        preStim = x <= 0;
        postStim = x > 0;
        con = mean(con, 2);
        incon = mean(incon, 2);
        delSum = 0;
        conSum = 0;
        inconSum = 0;
        predelSum = 0;
        preconSum = 0;
        preinconSum = 0;
        postdelSum = 0;
        postconSum = 0;
        postinconSum = 0;
        for j = 1:length(con)
          
          del = con(j) - incon(j);
          co = con(j) - mean(con);
          inc = incon(j) - mean(incon);
          delSum = delSum + del * del;
          conSum = conSum + co * co;
          inconSum = inconSum + inc * inc;
          if(x(j) <= 0)
            predelSum = predelSum + del * del;
            preconSum = preconSum + co * co;
            preinconSum = preinconSum + inc * inc;
          else
            postdelSum = postdelSum + del * del;
            postconSum = postconSum + co * co;
            postinconSum = postinconSum + inc * inc;
          end
        end
        %normalize for length
        postdelSum = postdelSum / sum(postStim);
        postconSum = postconSum / sum(postStim);
        postinconSum = postinconSum / sum(postStim);
        predelSum = predelSum / sum(preStim);
        preconSum = preconSum / sum(preStim);
        preinconSum = preinconSum / sum(preStim);
        
        %normalize for size
        delSums(i) = delSum/mean([conSum inconSum]);
        predelSums(i) = predelSum/mean([conSum inconSum]);
        postdelSums(i) = postdelSum/mean([conSum inconSum]);
        
        %find changes from pre to post
        preToPostChanges(i) = postdelSums(i) - predelSums(i);
      end
      i = i + 1;
    end
    
    [~, delSumInds] = sort(preToPostChanges);
  end
  
  for j = 0:100
    close all;
    i = delSumInds(end-j);
    %plot imagesc showing each trial
    con = congruentMat(:,:,i);
    incon = incongruentMat(:,:,i);
    zeroPoint = round(preSeconds * cohSrate);
    [~, conInd] = sort(con(zeroPoint,:));
    con(:, 1:size(con,2)) = con(:, conInd);
    [~, inconInd] = sort(incon(zeroPoint,:));
    incon(:, 1:size(incon,2)) = incon(:, inconInd);
    con = con'; incon = incon';
    buffer = repmat(mean(mean(con)), [10, intervalSampleCount]);
    toPlot = [con; buffer; incon];
    label = coh.labels{i};
    figure;
    im = imagesc('XData', x, 'CData', toPlot);
    xlabel('time');
    ylabel('congruent    /    incongruent');
    colorbar;
    title(label);
    %plot showing means
    con = congruentMat(:,:,i);
    incon = incongruentMat(:,:,i);
    con = mean(con,2);
    incon = mean(incon,2);
    figure;
    plot(x, [con, incon]);
    legend({'congruent', 'incongruent'});
    title(label);
    xlabel('time');
  end
  if(doAllFiles)
    fileCounter = fileCounter + 1;
  else
    fileCounter = length(preFiles);
  end
end