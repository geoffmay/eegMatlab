

%infomap
densityThreshold = 0.2;
%                 Run_Infomap_nopar(toShow, ones(size(toShow)), 0, densityThreshold, false, '/home/data/EEG/processed/infomap');%, numreps, numpools,structure_indices)
%                 infoClusters = load('/home/data/EEG/processed/infomap/rawassn.txt');
sessionName = filename(1:strfind(filename,'_63')-1);
sessionName = strrep(sessionName, ' ', '_');
folderName = fullfile('/home/data/EEG/processed/Robi/coherenceInfomap', sessionName);
if(~exist(folderName, 'dir'))
    mkdir(folderName);
end
Run_Infomap_nopar(toShow, ones(size(toShow)), 0, densityThreshold, false, folderName);

infoClusters = load(fullfile(folderName, 'rawassn.txt'));

targetElectrode = 'POz';
%             targetElectrode = 'P3';
targetIndices = cellfun(@length, strfind(pairLabels, targetElectrode));
targetClusters = infoClusters(find(targetIndices));
modeCluster = mode(targetClusters);

[sortCluster, sortIndex] = sort(infoClusters);
cluster2 = pairLabels(sortIndex(find(sortCluster==2)));

sums = sum(toShow,1);
[sortCluster2, sortIndex2] = sort(sums);


f = figure;
colormap('jet');
imagesc(toShow(sortIndex, sortIndex));
title(strrep(sessionName, '_', ' '));
colorbar;
%end infomap

sortedMap = toShow(sortIndex, sortIndex);
sortedLabels = pairLabels(sortIndex);
outputFilename2 = sprintf('/home/data/EEG/processed/Robi/coherenceInfomap/%s.mat', sessionName);

save(outputFilename2, 'sortedMap', 'sortedLabels', 'sortCluster', 'modeCluster', 'filename', 'badReferenceCutoff');
