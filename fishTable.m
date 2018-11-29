function [ tab ] = fishTable( masterTable, toCompare )
%FISHTABLE Summary of this function goes here
%   Detailed explanation goes here


if(~exist('toCompare', 'var'))
  toCompare = size(masterTable, 2);
end

%do pca on changed eeg variables for dimension reduction
cohPowAsyFlags = [1 1 1];  
%1 1 1 yields a 7th component accounting for 3% of eeg signal variance that
%correlates with DMN change p = 0.0186.

changeIndex = cellfun(@length, strfind(masterTable.Properties.VariableNames, 'change')) > 0;
cohIndex = cellfun(@length, strfind(masterTable.Properties.VariableNames, 'coh')) > 0;
powIndex = cellfun(@length, strfind(masterTable.Properties.VariableNames, 'abs')) > 0;
asyIndex = cellfun(@length, strfind(masterTable.Properties.VariableNames, 'asy')) > 0;
cohIndex = cohIndex & cohPowAsyFlags(1);
powIndex = powIndex & cohPowAsyFlags(2);
asyIndex = asyIndex & cohPowAsyFlags(3);
selectedIndex = changeIndex & (cohIndex | asyIndex | powIndex);
changeLabels = masterTable.Properties.VariableNames(selectedIndex);
pcaInput = masterTable{:, selectedIndex};
changeSubjectIndex = ~any(isnan(pcaInput), 2);
changeSubjects = masterTable{changeSubjectIndex, 1};
[changePca.COEFF, changePca.SCORE, changePca.LATENT, changePca.TSQUARED, changePca.EXPLAINED, changePca.MU] = pca(pcaInput(changeSubjectIndex, :));
changeToCompare = masterTable{changeSubjectIndex, toCompare};
clear changeRhos changePs
for i = 1:size(changePca.SCORE, 2)
  x = changePca.SCORE(:, i);
  y = changeToCompare;
  keep = ~(isnan(x) | isnan(y));
  [changeRhos(i), changePs(i)] = corr(x(keep),y(keep));
end
changePca.EXPLAINED


%debug
%end debug


for column1Index = 1:length(toCompare)
  column1 = masterTable{:, toCompare(column1Index)};
  onlyHits = false;
  hitCounter = 1;
  for column2Index = 1:size(masterTable,2)
    column2 = masterTable{:, column2Index};
    if(isnumeric(column2))
      mat = [column1 column2];
      remove = any(isnan(mat),2);
      mat(remove,:) = [];
      if(size(mat,1) > 2)
        [rho(1), p(1)] = corr(mat(:,1), mat(:,2));
        [rho(2), p(2)] = corr(mat(:,1), mat(:,2), 'Type', 'Kendall');
        [rho(3), p(3)] = corr(mat(:,1), mat(:,2), 'Type', 'Spearman');
        if(onlyHits)
          if(any(p < 0.05))
            columns{hitCounter} = masterTable.Properties.VariableNames{column2Index};
            pearsonRhos(hitCounter) = rho(1);
            kendallRhos(hitCounter) = rho(2);
            spearmanRhos(hitCounter) = rho(3);
            pearsonPs(hitCounter) = p(1);
            kendallPs(hitCounter) = p(2);
            spearmanPs(hitCounter) = p(3);
            fprintf('\n%s', columns{hitCounter});
            hitCounter = hitCounter + 1;
          end
        else
          columns{hitCounter} = masterTable.Properties.VariableNames{column2Index};
          pearsonRhos(hitCounter) = rho(1);
          kendallRhos(hitCounter) = rho(2);
          spearmanRhos(hitCounter) = rho(3);
          pearsonPs(hitCounter) = p(1);
          kendallPs(hitCounter) = p(2);
          spearmanPs(hitCounter) = p(3);
          fprintf('\n%s', columns{hitCounter});
          hitCounter = hitCounter + 1;
        end
      end
    end
  end
end

columns = columns';
pearsonRhos = pearsonRhos';
kendallRhos = kendallRhos';
spearmanRhos = spearmanRhos';
pearsonPs = pearsonPs';
kendallPs = kendallPs';
spearmanPs = spearmanPs';

tab = table(columns, pearsonRhos, kendallRhos, spearmanRhos, pearsonPs, kendallPs, spearmanPs);
sortTab = sortrows(tab, 'pearsonPs');

end

