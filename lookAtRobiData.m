psychFolder = '/home/data/EEG/data/ROBI/psychTests';


if(~exist('rhos', 'var'))
  load(fullfile(psychFolder, 'neuropsych.mat'));
  rhos = NaN(size(masterTable,2));
  ps = NaN(size(masterTable,2));
end
for i = 1:size(masterTable, 2)
  fprintf('%s: %d of %d (%s)\n', char(datetime), i, size(masterTable,2), masterTable.Properties.VariableNames{i});
  rhos(i,i) = 1;
  ps(i,i) = 0;
  x = masterTable{:,i};
  if(isnumeric(x))
    for j = i+1:size(masterTable,2)
      y = masterTable{:,j};
      if(isnumeric(y))
        keep = ~isnan(x) & ~isnan(y);
        if(sum(keep) > 2)
          [rho,p] = corr(x(keep),y(keep));
          rhos(i,j) = rho;
          rhos(j,i) = rho;
          ps(i,j) = p;
          ps(j,i) = p;
        end
      end
    end
  end
end