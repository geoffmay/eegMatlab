performAnalysis = 0;
displayResults = 1;
makeDummy = 1;
dconnFolder = '/media/Seagate/MRI/Analysis/ROBI';
infomapFolder = '/media/Seagate/MRI/subjects';
outputFolder = '/home/data/EEG/processed/Robi/segregation';
dconnFiles = dir(dconnFolder);
dconnPaths = {dconnFiles.name}';
pre = dconnPaths(cellfun(@length, strfind(dconnPaths, '_pre.dconn.nii')) > 0);
post1 = dconnPaths(cellfun(@length, strfind(dconnPaths, '_post101.dconn.nii')) > 0);
post2 = dconnPaths(cellfun(@length, strfind(dconnPaths, '_post102.dconn.nii')) > 0);
post = dconnPaths(cellfun(@length, strfind(dconnPaths, '_post.dconn.nii')) > 0);
allDconn = [pre; post1; post2; post];
filenameFilter = 'ROBI';
remove = cellfun(@length, strfind(allDconn, filenameFilter)) == 0;
allDconn(remove) = [];

if(performAnalysis)
  for fileCounter =1:length(allDconn)
    
    %path = fullfile(dconnFolder, post1{fileCounter});
    dconnFilename = allDconn{fileCounter};
    dconnPath = fullfile(dconnFolder, allDconn{fileCounter});
    subjectId = dconnFilename(1:strfind(dconnFilename, '_')-1);
    sessionId = dconnFilename(1:strfind(dconnFilename, '.')-1);
    outputPath = fullfile(outputFolder, sprintf('%s.mat', sessionId));
    fprintf('\n%s: file %d of %d: %s', char(datetime), fileCounter, length(allDconn), sessionId);
    doThis = true;
    if(makeDummy)
      if(exist(outputPath, 'file'))
        doThis = false;
      else
        placeHolder = sprintf('processing started on %s', char(datetime));
        save(outputPath, 'placeHolder');
      end
    end
    if(doThis)
      infomapPath = fullfile(infomapFolder, subjectId, '/infomap/RSFC/rawassn_minsize400_regularized_recolored.dscalar.nii');
      if(exist(infomapPath, 'file'))
        infomap = ft_read_cifti_mod(infomapPath)
      end
      dconn = ft_read_cifti_mod(dconnPath);
      isStructure = find(dconn.brainstructure > 0); %indexes length 67631, max 73203
      
      cortexLabels = find(cellfun(@length, strfind(dconn.brainstructurelabel, 'CORTEX_')) > 0);
      cortex = ismember(dconn.brainstructure, cortexLabels); %logicals length 73203, 59412 ones, 13791 zeros
      mappedCortex = cortex(isStructure); %logicals length 67631, 59412 ones
      cortexInfomap = infomap.data(mappedCortex);
      cortexDconn = dconn.data(cortex, cortex);
      clear dconn; %it's quite big
      infomapLabels = unique(cortexInfomap);
      corrMatrix = NaN(length(infomapLabels));
      for i = 1:length(infomapLabels)
        %         fprintf('\n%s: row %d of %d', char(datetime), i, length(infomapLabels));
        label1 = infomapLabels(i);
        indexes1 = cortexInfomap == label1;
        for j = i:length(infomapLabels)
          label2 = infomapLabels(j);
          indexes2 = cortexInfomap == label2;
          thisPair = cortexDconn(indexes1, indexes2);
          meanVal = mean(mean(thisPair));
          corrMatrix(i,j) = meanVal;
          corrMatrix(j,i) = meanVal;
        end
      end
      while(exist('outputPath', 'file'))
        outputPath = [outputPath '.mat'];
      end
      save(outputPath, 'corrMatrix', 'infomapLabels', '-v7.3');
      clear cortexDconn;
    end
  end
end

if(displayResults)
  outputFiles = dir(outputFolder);
  
end


