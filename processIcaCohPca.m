function processIcaCohPca( eegDataFilename, icaWeights )
%PROCESSICACOHPCA1 Summary of this function goes here
%   Detailed explanation goes here


if(~exist('icaWeights', 'var'))
  icaCoh = deriveIcaCoherenceMatrix(eegDataFilename);
else
  icaCoh = deriveIcaCoherenceMatrix(eegDataFilename, icaWeights);
end

  
  fprintf('\n(%s) calculating pca', char(datetime));
  a = which('pca');
  if(length(strfind(a, 'eeglab')) > 0)
    rmpath(fileparts(a));
    a = which('pca');
  end
  [cohPca.coeff, cohPca.score, cohPca.latent, cohPca.tsquared, cohPca.explained] = pca(icaCoh.matrix);
  save(fullfile(fileparts(outputFilename), sprintf('%s.mat', subjectName)), 'icaCoh', '-v7.3');
  
  if(true)
    close all;
    x = (1:size(cohPca.score, 1)) ./ 128 ./ 60;
    plot(x, cohPca.score(:,1))
    xlabel('time (minutes)')
    ylabel('pca component 1 coherence')
    title('coherence with eyes open, followed by eyes closed');
    %   figFile = fullfile(fileparts(outputFilename), 'coherenceTimecourse.pdf');
    %   export_fig(figFile);
  end
  
  writetable(icaCoh.icaInfo.spheredWeights, outputFilename);
  
  
  fileId = fopen(pcaFilename, 'w');
  for i = 1:length(icaCoh.labels)
    fprintf(fileId, '%s', icaCoh.labels{i});
    if(i < length(icaCoh.labels))
      fprintf(fileId, '%s', ',');
    end
  end
  fprintf(fileId, '\n');
  for i = 1:size(cohPca.coeff, 1)
    fprintf(fileId, '%f', cohPca.coeff(i));
    if(i < size(cohPca.coeff, 1))
      fprintf(fileId, ',');
    end
  end
  
  
  fprintf('\n(%s) finished processing', char(datetime));
  
  
  
  message = sprintf('attached are ica components and coherence pca 1st component for subject: %s.  Percent variance explained for pca component is %f.', subjectName, cohPca.explained(1));
  message = sprintf('subject: %s.  1st pca comp %%var explained: %f%%.', subjectName, cohPca.explained(1));
  
  try
    % Send an email with the usable parts of our analysis:
    sendFromOutlook('geomay@gmail.com', 'ica components and pca', message, {outputFilename, pcaFilename});
    fprintf('\n(%s) email sent', char(datetime));
  catch ex
  end
  

end

