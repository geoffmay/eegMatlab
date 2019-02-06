clear
close all

folder = 'C:\Vision\Raw Files\eegtest\export\corrLag';

plotAllIncludingArtifact = 0;
plotIndividual = 1;

files = dir(folder);
files([files.isdir])=[];
startTimes = 63;
endTimes = 187;

for i = 1:length(files)
    try
        data = load(fullfile(files(i).folder, files(i).name));
        if(isfield(data, 'summary'))
            summary = data.summary;
            lagMatrix(:,:,i) = data.summary.lagMatrix;
            pMatrix(:,:,i) = data.summary.pMatrix;
            rhoMatrix(:,:,i) = data.summary.rhoMatrix;
            startTimes(i) = summary.parameters.windowStartFrame;
            endTimes(i) = summary.parameters.windowEndFrame;
        elseif(isfield(data, 'lagMatrix'))
            lagMatrix(:,:,i) = data.lagMatrix;
            pMatrix(:,:,i) = data.pMatrix;
            rhoMatrix(:,:,i) = data.rhoMatrix;
        end
    catch ex
        %file read error; skip the file
    end
end

timePoint = 577;
[permRho, perm] = clusterMatrix(rhoMatrix(:,:,timePoint));
permLab = {summary.parameters.channelLocations.labels};
permLab = permLab(perm);


longCounter = 1;
for i = 1:size(lagMatrix, 3)
    mat = lagMatrix(:,:,i);
    mat = clusterMatrix(mat, perm);
    longLag(longCounter:(longCounter + size(mat,1) - 1), :) = mat;
    squeezedLag(i,:) = reshape(mat, [1, numel(mat)]);
    mat = rhoMatrix(:,:,i);
    mat = clusterMatrix(mat, perm);
    longRho(longCounter:(longCounter + size(mat,1) - 1), :) = mat;
    squeezedRho(i,:) = reshape(mat, [1, numel(mat)]);
    longCounter = longCounter + size(mat,1);
end

if(plotAllIncludingArtifact)
    %linkage
    Z = linkage(longLag');
    figure;
    [H,T,OUTPERM] = dendrogram(Z, 0);
    permLab = {summary.parameters.channelLocations.labels};
    permLab = permLab(OUTPERM);
    
    
    %timecourse
    figure;
    imagesc(squeezedLag);
    a = get(gca, 'ytick');
    set(gca, 'yticklabel', startTimes(a)./250);
    ylabel('time (seconds)');
    xlabel('flattened lag matrix');
    
    figure;
    imagesc(squeezedRho);
    a = get(gca, 'ytick');
    set(gca, 'yticklabel', startTimes(a)./250);
    ylabel('time (seconds)');
    xlabel('flattened rho matrix');
    
    %one slice
    timePoint = 577;
    [permRho, perm] = clusterMatrix(rhoMatrix(:,:,timePoint));
    % permLab = {summary.parameters.channelLocations.labels};
    % permLab = permLab(perm);
    
    %plot rho
    figure;
    imagesc(permRho);
    set(gca, 'xtick', 1:size(permRho,1));
    set(gca, 'xticklabel', permLab);
    set(gca, 'ytick', 1:size(permRho,1));
    set(gca, 'yticklabel', permLab);
    
    
end

plotSecondLimits = [144, 568];
plotFrameStart = min(find((startTimes./250)>plotSecondLimits(1)));
plotFrameEnd = max(find((startTimes./250)<plotSecondLimits(2)));

firstRunLag = squeezedLag(plotFrameStart:plotFrameEnd, :);
firstRunRho = squeezedRho(plotFrameStart:plotFrameEnd, :);


%plot lag
if(plotIndividual)
    figure;
    iRange = plotFrameStart:plotFrameEnd;
    for i = iRange
        permLag = clusterMatrix(lagMatrix(:,:,i), perm);
        imagesc(permLag ./ 250, [-.25,.25]);
        set(gca, 'xtick', 1:size(permLag,1));
        set(gca, 'xticklabel', permLab);
        set(gca, 'ytick', 1:size(permLag,1));
        set(gca, 'yticklabel', permLab);
        colorbar
        title(sprintf('clustered time lag %.02f-%.02f s', startTimes(i)/250, endTimes(i)/250));
    end
end

%timecourse
figure;
imagesc(firstRunLag);
a = get(gca, 'ytick');
set(gca, 'yticklabel', startTimes(a)./250);
ylabel('time (seconds)');
xlabel('flattened lag matrix');

figure;
imagesc(firstRunRho);
a = get(gca, 'ytick');
set(gca, 'yticklabel', startTimes(a)./250);
ylabel('time (seconds)');
xlabel('flattened rho matrix');

scatterLag = tsne(firstRunLag);
scatterRho = tsne(firstRunRho);
figure;
scatter(scatterLag(:,1), scatterLag(:,2));
smallCluster = scatterLag(:,2) > 10;
title('lag tsne');

figure;
scatter(scatterLag(:,1), scatterLag(:,2));
title('rho tsne');





