folder = '/home/data/EEG/data/ROBI/psychTests/';

inputFilename = fullfile(folder, 'neuropsych.mat');
outputFilename = fullfile(folder, 'correlations.mat');
if(~exist(inputFilename, 'file'))
  masterTable = openRobiData();
else
  load(inputFilename);
end

tab = fishTable(masterTable);

pearsonTab = sortrows(tab, 'pearsonPs');

save(outputFilename, 'pearsonTab')


%debug
xLabel = 'exit_coh_F8_Oz_1Hz_4Hz';
yLabel = 'mri_dmn_delta';
y = masterTable{:, yLabel};
x = masterTable{:, xLabel};
xLabel = strrep(xLabel, '_', '-');
yLabel = strrep(yLabel, '_', '-');
mat = [x y];
remove = any(isnan(mat),2);
mat(remove, :) = [];
figure;
scatter(mat(:,1), mat(:,2));
xlabel(xLabel);
ylabel(yLabel);
%end debug

isChange = cellfun(@length, strfind(tab{:,1}, 'change')) > 0;
isLag = cellfun(@length, strfind(tab{:,1}, 'lag')) > 0;

lagTitles = tab{(isChange & isLag), 1};
pearsonLagChanges = tab{(isChange & isLag), 'pearsonRhos'};
spearmanLagChanges = tab{(isChange & isLag), 'spearmanRhos'};
kendallLagChanges = tab{(isChange & isLag), 'kendallRhos'};


