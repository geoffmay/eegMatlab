rootFolder = 'C:\Vision\Raw Files\eegtest\mri\DICOM';
outputFile = 'C:\Vision\Raw Files\eegtest\mri\info.mat';

contents = dir(rootFolder);
contents(1:2) = [];
folders = contents([contents.isdir]);
files = contents(~[contents.isdir]);
for i = 1:length(files)
    inputPaths{i} = fullfile(rootFolder, files(i).name);
end
pathCounter = length(inputPaths)+1;

for i = 1:length(folders)
    files = dir(fullfile(rootFolder, folders(i).name));
    files([files.isdir]) = [];
    for j = 1:length(files)
        path = fullfile(rootFolder, folders(i).name, files(j).name);
        inputPaths{pathCounter} = path;
        pathCounter = pathCounter + 1;
    end
end
clear files


for i = 1:length(inputPaths)
    fprintf('\n%d of %d', i, length(inputPaths));
    path = inputPaths{i};
    %dicom = dicomread(path);
    info = dicominfo(path);
    infos{i} = info;
    protocols{i} = info.ProtocolName;
    
    if(true)
        if(isfield(info, 'SeriesTime'))
            seriesTimes(i) = str2doubleq(info.SeriesTime);
        else
            seriesTimes(i) = NaN;
        end
        
        if(isfield(info, 'AcquisitionTime'))
            acquisitionTimes(i) = str2doubleq(info.AcquisitionTime);
        else
            acquisitionTimes(i) = NaN;
        end
    end
end

if(false)
    
    isEchoPlanar = find(strcmp(protocols, 'FE_EPI_32chSHC'));
    for i = 1:length(isEchoPlanar)
        info = infos{isEchoPlanar(i)};
        sliceLocation = NaN;
        if(isfield(info, 'SliceLocation'))
            sliceLocation = info.SliceLocation;
        end
        sliceLocations(i) = sliceLocation;
    end
    epTimes = acquisitionTimes(isEchoPlanar);
    normTimes = epTimes - nanmean(epTimes);
    normTimes = normTimes ./ nanstd(normTimes);
    normLocs = sliceLocations - nanmean(sliceLocations);
    normLocs = normLocs ./ nanstd(normLocs);
    figure;
    plot([normTimes', normLocs']);
    isBottom = find(sliceLocations == max(sliceLocations));
    diffBottom = diff(isBottom);
    isBottom(find(diffBottom==1) + 1) = [];
    
    testPath = inputPaths{isEchoPlanar(isBottom(2)-1)};
    dicom = dicomread(testPath);
    figure;
    imagesc(dicom);
    
    tab = tabulate(diffBottom);
    tab(tab(:,2)==0,:)=[];
    
    legend({'times','locs'});
    [sortAcq, sortInd] = sort(epTimes);
    figure;
    plot([normTimes(sortInd)', normLocs(sortInd)']);
    legend({'times','locs'});
    
    tab = table;
    for i = 1:100;
        ind = isEchoPlanar(i);
        path = fullfile(rootFolder, files(ind).name);
        info = dicominfo(path);
        
        names = fieldnames(info);
        isTime = find(cellfun(@length, strfind(lower(names), 'time')) > 0);
        names(isTime)
        
        for j = 1:length(isTime)
            name = names{isTime(j)};
            tab{j, 1} = {name};
            tab{j, 1+i} = {getfield(info, name)};
        end
        img = dicomread(path);
        imageData(:,:,i) = img;
    end
    
    i = 0;
    
    i = i+1;
    path = fullfile(rootFolder, files(isEchoPlanar(sortInd(i))).name);
    dicom = dicomread(path);
    info = dicominfo(path);
    imagesc(dicom);
    title(sprintf('%d', i));
    
end

save(outputFile, '-v7.3');

