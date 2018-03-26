function [ topographicMaps, segmentCenters, segmentBoundaries, activationPlot, EEG ] = demoGFP(bdfFileName, numberOfClusters, smoothMaps, plotResults, outputFileName )
%DEMOGFP computes EEG topographic microstates and their timecourses 
%   Segment EEG topographic epochs based on global map dissimilarity;
%   perform hierarchical clustering of segments and output maps and a plot
%   of when each cluster is active.  
%   Based on http://www.hindawi.com/journals/cin/2011/813870/
%   Also thanks http://nemo.nic.uoregon.edu/wiki/images/6/65/JBrain_Topogr_2008_Murray_Topographic_ERP_analyses_a_step_by_step.pdf

argCounter = 1;
if nargin < argCounter
    bdfFileName = '/Users/Geoff/Box Sync/For Geoff Mays/VET MIND EEG files/VM103.1.Tones.bdf';
end
argCounter = argCounter + 1;
if nargin < argCounter
    numberOfClusters=400;
end
argCounter = argCounter + 1;
if nargin < argCounter
    smoothMaps = true;
end
argCounter = argCounter + 1;
if nargin < argCounter
    plotResults = true;
end
argCounter = argCounter + 1;
if nargin < argCounter
    outputFileName = '';
end
argCounter = argCounter + 1;
%check to make sure octave isnt' linked.
scrap = filtfilt([.5, .5],1,[1,2,3,4,5]);
clear scrap;

%read the file
refChan = int32(32);
EEG = pop_readbdf(bdfFileName, {}, 43, refChan, false);
% EEG.data(refChan,:)=[];
% EEG.chanlocs(refChan,:)=[];
% EEG.nbchan = EEG.nbchan-1;

%remove all but EEG data
EEG = removeNonEeg(EEG);
%passband filter the data
EEG.data = eegfilt(EEG.data, EEG.srate , 1, 125, size(EEG.data,2), 3*EEG.srate, 0, 'fir1', 0);

%points of "stable" topography are local minima of global map
%dissimilarity. Local maxima indicate points of transistion.
[gmd, EEG] = globalMapDissimilarity(EEG);
filterSize = int32(EEG.srate/150);
filterVector = (1/double(filterSize))*ones(1,filterSize);
filtGmd = filtfilt(filterVector,1,gmd);
[localMinima, segmentCenters] = findpeaks(1-filtGmd);
[localMaxima, segmentBoundaries] = findpeaks(filtGmd);
display(strvcat([num2str(length(localMinima)), ' segments identified.  Smoothing...']));

%use smoothed EEG to get the average potentials of the GMD local minimum.
smoothEEG = EEG;
if(smoothMaps)
    for i = 1:size(EEG.data,1)
        smoothEEG.data(i,:) = filtfilt(filterVector,1,EEG.data(i,:));
    end
end
mapSegments = ones(size(EEG.data,1),0);
for i = 1:size(segmentCenters, 2)
    map = smoothEEG.data(:, segmentCenters(i));
    map = map - mean(map);
    mapSegments(:, size(mapSegments,2)+1) = map;
end

%cluster the data using hierarchical clustering.
display('clustering...');
try
Z = linkage(mapSegments', 'average', 'euclidean');
fullMatrixUsed = true;
catch
Z = linkage(mapSegments', 'median', 'euclidean', 'savememory', 'on');    
fullMatrixUsed = false;
end

clusters = cluster(Z, numberOfClusters);
%discard clusters that appear rarely.
sortedClusters = sortrows(tabulate(clusters), 2);
percentages = sortedClusters(:,3);
threshold = 0.5;
%display(strvcat(['discarding clusters with fewer than ', num2str(threshold), '% of segments']));
exclude = percentages < threshold;
sortedClusters(exclude,:) = [];
display(strvcat([num2str(size(sortedClusters,1)), ' clusters remain']));
activationPlot = zeros(1, length(localMinima));
for i = 1:length(clusters)
    for j = 1:size(sortedClusters,1)
        if(sortedClusters(j,1) == clusters(i))
            activationPlot(i) = size(sortedClusters,1)-j+1;
        end
    end
end

%compute the mean topography for each cluster
topographicMaps = ones(EEG.nbchan,0);
sdMaps = ones(EEG.nbchan,0);
for i = 1:size(sortedClusters,1);
    clusterMaps = mapSegments(:, activationPlot==(size(sortedClusters,1)-i+1));
    topographicMaps(:, size(topographicMaps,2)+1) = mean(clusterMaps, 2);
    sdMaps(:, size(sdMaps,2)+1) = std(clusterMaps, 0, 2);
end

if(plotResults)
    close all;
    %display results, first showing the time course of clusters
    figure;
    plot(activationPlot);
    zoom xon;
    pan xon;
    % %plot
    % while length(gmd) < length(EEG.times)
    %     gmd(length(gmd)+1) = 0;
    %     filtGmd(length(filtGmd)+1) = 0;
    % end
    % colors = zeros(1, length(gmd));
    % colorStartIndex = 0;
    % maxIndexStart = 1;
    % if(localMinimaLocations(1) > localMaximaLocations(1))
    %     colorStartIndex = localMaximaLocations(1);
    %     maxIndexStart = 2;
    % end
    % for i = maxIndexStart:length(clusters)
    %     nextIndex = localMinimaLocations(i);
    %     colors(colorStartIndex:nextIndex)=clusters(i);
    %     colorStartIndex = nextIndex;
    % end
    % plotTime = EEG.times ./ 1000;
    % figure;
    % hold on;
    %  plot(plotTime, gmd, 'y');
    % plot(plotTime, filtGmd, 'r');
    % plot(plotTime(localMinimaLocations), 1-localMinima, 'k');
    % plot(plotTime(localMaximaLocations), localMaxima, 'b');
    %patch(plotTime(localMaximaLocations), localMaxima, clusters, clusters, 'edgecolor', 'flat');
    % plot(plotTime(localMaximaLocations), clusters);
    % zoom xon;
    % pan xon;
    
    %plot the topographic maps
    for i = 1:size(topographicMaps,2)
        figure;
        topoplot(topographicMaps(:,i),EEG.chanlocs);
    end
    tilefigs;
end
if length(outputFileName) > 0
    save(strcat(outputFileName, 'clusters.mat'));
end
end

