function [fig, coeff, allLabels] = plotCoherencePca(cohPca, allLabels, coefficientIndex, titlePrefix, componentLocations)
%PLOTCOHERENCEPCA Computes a timeseries of coherence for a given frequency.
%   Detailed explanation goes here


%Todo: use component locations to plot ica components

fixedMinMax = 1;
plotRelative = 0;
plotSeparate = 1;
splitColormap = 1;

if(plotRelative)
    powType = 'rel';
else
    powType = 'abs';
end

if(fixedMinMax)
    minWeight = 0;
    maxWeight = 1;
end
if(~exist('coefficientIndex', 'var'))
    coefficientIndex = 1;
end
if(~exist('titlePrefix', 'var'))
    titlePrefix = '';
end
freqLabels1 = {'delta','theta','alpha','beta','hibeta'};
freqLabels2 = {'1Hz-4Hz', '5Hz-8Hz', '9Hz-12Hz', '13Hz-24Hz', '25Hz-30Hz'};

if(isfield(cohPca, 'COEFF'))
    magnitude = std(cohPca.SCORE(:, coefficientIndex));
    coeff = cohPca.COEFF(:,coefficientIndex);
elseif(isfield(cohPca,'coherencePlot'))
    chanCount = length(cohPca.powerPlot);
    chanPairCount = length(cohPca.coherencePlot);
    freqCount = size(cohPca.coherencePlot(1).coherence, 2);
    coeffCount = freqCount * (chanCount + chanPairCount);
    
    %coeff = NaN(1, size(cohPca.coherencePlot,2)*(size(cohPca.coherencePlot,3)+size(cohPca.powerPlot,3)));
    counter = 1;
    for i = 1:chanPairCount
        for j = 1:freqCount
            coeff(counter) = mean(cohPca.coherencePlot(i).coherence(:,j),1);
            a(counter,:) = [i j];
            allLabels{counter} = sprintf('%s %s', cohPca.coherencePlot(i).label, freqLabels1{j});
            counter = counter + 1;
        end
    end
    for i = 1:chanCount
        if(size(cohPca.powerPlot(i).absolutePower, 2) == freqCount)
            for j = 1:freqCount
%                 fprintf('\n%d %d', i, j);
                coeff(counter) = mean(cohPca.powerPlot(i).absolutePower(:,j),1);
                a(counter,:) = [i j];
                allLabels{counter} = sprintf('%s %s %s', cohPca.powerPlot(i).label, freqLabels1{j}, powType);
                counter = counter + 1;
            end
        end
    end
    magnitude = 1;
    
elseif (isfield(cohPca, 'channelPairLabels'))
%    coeff = NaN(size(cohPca.coherence,1), size(cohPca.coherence,2)*(size(cohPca.coherence,3)+size(cohPca.absPower,3)));
    coeff = NaN(1, size(cohPca.coherence,2)*(size(cohPca.coherence,3)+size(cohPca.absPower,3)));
    counter = 1;
    for i = 1:size(cohPca.coherence, 3)
        for j = 1:size(cohPca.coherence, 2)
            coeff(counter) = mean(cohPca.coherence(:,j,i),1);
            counter = counter + 1;
        end
    end
    for i = 1:size(cohPca.absPower, 3)
        for j = 1:size(cohPca.coherence, 2)
%             fprintf('\n%d %d', i, j);
            if(plotRelative)
                coeff(counter) = mean(cohPca.relPower(:,j,i),1);
            else
                coeff(counter) = mean(cohPca.absPower(:,j,i),1);
            end
            counter = counter + 1;
        end
    end
    magnitude = 1;
else
    magnitude = 1;
    coeff = cohPca;
end

%magnitude = std(cohPca.SCORE(:, coefficientIndex));
addpath('/home/data/EEG/scripts/eeglab13_4_4b/functions/sigprocfunc/');
addpath('/home/data/EEG/scripts/eeglab13_4_4b/functions/guifunc/');



if(~exist('allLabels','var'))
    if(length(coeff) == 2325)
        load('/home/data/EEG/processed/Robi/antChanlocs.mat');
        chanlocs([13 19 33 34]) = [];
    elseif(length(coeff) == 2480)
        load('/home/data/EEG/processed/Robi/antChanlocs.mat');
        chanlocs([13 19 34]) = [];
    elseif(length(coeff) == 2805)
        load('/home/data/EEG/processed/Robi/antChanlocs.mat');
        chanlocs([34]) = [];
    elseif(length(coeff) == 950)
        load('/home/data/EEG/processed/Robi/aniChanlocs.mat');
    else
        error('unhandled channel count');
    end
    freqLabels1 = {'delta','theta','alpha','beta','hibeta'};
    allLabels = cell(1,length(cohPca));
    counter = 1;
    for i = 1:length(chanlocs)
        for j = i+1:length(chanlocs)
            for k = 1:length(freqLabels1)
                allLabels{counter} = sprintf('%s-%s %s', chanlocs(i).labels, chanlocs(j).labels, freqLabels1{k});
                counter = counter + 1;
            end
        end
    end
    for i = 1:length(chanlocs)
        for k = 1:length(freqLabels1)
                allLabels{counter} = sprintf('%s %s abs', chanlocs(i).labels, freqLabels1{k});
            counter = counter + 1;
        end
    end
else
    if(length(allLabels) == 950)
        load('/home/data/EEG/processed/Robi/aniChanlocs.mat');
    else
        load('/home/data/EEG/processed/Robi/antChanlocs.mat');
        chanlocs([13 19 34]) = [];
    end
end
if(plotRelative)
    filterSequence = {'rel'};
else
    filterSequence = {'abs'};
end
weightCounter = 1;
% close all;
clear weights allWeights;
if(~plotSeparate)
    fig = figure;
    title(titlePrefix);
    %     set(fig, 'title', titlePrefix);
    scrdim = get( 0, 'ScreenSize' ); % Get screen size.
    scrwid = scrdim(3); % Screen width.
    scrhgt = scrdim(4); % Screen height.
    fig.Position = [0 0 scrwid scrhgt];
else
  fig = cell(0);
end

%plot power
for filterCounter = 1:length(filterSequence)
    filter2 = filterSequence{filterCounter};
    for freqCounter = 1:length(freqLabels1)
        filter1 = sprintf(' %s',freqLabels1{freqCounter});
        filter2 = powType;
        hit1 = cellfun(@length, strfind(allLabels, filter1));
        hit2 = cellfun(@length, strfind(allLabels, filter2));
        ind = find(hit1&hit2);
        %         weights(:,freqCounter) = cohPca.COEFF(ind, coefficientIndex) * magnitude;
        if(length(coeff) > 0)
          weights(:,freqCounter) = coeff(ind) * magnitude;
        else
          weights(:,freqCounter) = NaN(0,1);
        end
        chanLabels{freqCounter} = allLabels(ind);
    end
    if(~fixedMinMax)
        minWeight = min(min(weights));
        maxWeight = max(max(weights));
    end
    for freqCounter = 1:length(freqLabels1)
        filter1 = sprintf(' %s',freqLabels1{freqCounter});
        %filter2 = 'abs';
        if(plotSeparate)
            thisFig = figure;
            topoplot(weights(:,freqCounter), chanlocs, 'maplimits', [minWeight, maxWeight]);
            title(sprintf('%s %s %s',titlePrefix, filter1,filter2));
            colorbar();
            fig{end+1} = thisFig;
        else
            %             subplot(2, length(freqLabels), freqCounter);
            a = 10;
            subplot(a,a,a*a);
            %             subplot(1,1,1);
            xWid = .9/length(freqLabels1);
            pos = [(freqCounter-1) * xWid, .5, xWid, .5];
            set(gca, 'Position', pos);
            
            tempLabels = chanLabels{freqCounter};
            tempWeights = weights(:,freqCounter);
            for tempLabelCounter = 1:length(tempLabels)
                tl = tempLabels{tempLabelCounter};
                tl1 = tl(1:strfind(tl, freqLabels1{freqCounter})-2);
                tl1 = strrep(tl1, 'T3', 'T7');
                tl1 = strrep(tl1, 'T4', 'T8');
                tl1 = strrep(tl1, 'T5', 'P7');
                tl1 = strrep(tl1, 'T6', 'T8');
                %                 fprintf('\n%s', tl1);
                if(strcmp(tl1,'M1'))
                    dummy = 1;
                end
                tempIndex = find(strcmp({chanlocs.labels}, tl1));
                if(length(tempIndex) > 0)
                    tempLocs(tempLabelCounter) = chanlocs(tempIndex);
                else
                    tempWeights(tempLabelCounter) = [];
                end
            end
            if(~exist('tempLocs','var'))
                tempLocs = cell(0);
            end
            if(length(weights) > 0)
                top = topoplot(weights(:,freqCounter), tempLocs, 'maplimits', [minWeight, maxWeight]);
                if(freqCounter == 1)
                    cbar = colorbar();
                    set(cbar,'Position',[.93, .6, .015, .3]);
                end
                
                clear tempLocs;
                title(sprintf('%s %s',filter1,filter2));
            end
        end
    end
    
    allWeights{weightCounter} = weights;
    weightCounter = weightCounter + 1;
end
clear weights;

%plot coherence
for freqCounter = 1:length(freqLabels1)
    filter1 = sprintf(' %s',freqLabels1{freqCounter});
    filter2 = '-';
    hit1 = cellfun(@length, strfind(allLabels, filter1));
    hit2 = cellfun(@length, strfind(allLabels, filter2));
    ind = find(hit1&hit2);
    %     weights(:,freqCounter) = cohPca.COEFF(ind, coefficientIndex) * magnitude;
    weights{freqCounter} = coeff(ind) * magnitude;
    if(length(weights{freqCounter}) > 0)
        minWeights(freqCounter) = min(weights{freqCounter});
        maxWeights(freqCounter) = max(weights{freqCounter});
    else
        minWeights(freqCounter) = 0;
        maxWeights(freqCounter) = 1;
    end
    
    minWeight;
end
if(~fixedMinMax)
    %     minWeight = min(min(weights));
    %     maxWeight = max(max(weights));
    minWeight = min(minWeights);
    maxWeight = max(maxWeights);
end
for freqCounter = 1:length(freqLabels1)
    filter1 = sprintf(' %s',freqLabels1{freqCounter});
    filter2 = '-';
    hit1 = cellfun(@length, strfind(allLabels, filter1));
    hit2 = cellfun(@length, strfind(allLabels, filter2));
    ind = find(hit1&hit2);
    
    if(plotSeparate)
%         plotChannelPairs_PF(allLabels(ind), weights{freqCounter}, [minWeight maxWeight],true);
%function plotChannelPairs_PF(labelPairs, pairWeights, figname, figtype, missingval, colorLimits, drawColorbar, fig)
        figs = plotChannelPairs_PF(allLabels(ind), weights{freqCounter}, 'reliability', [0 0], realmin, [minWeight maxWeight], false);
        figs.Position(3) = 300;
        figs.Position(4) = figs.Position(3);
        figFolder = '/home/data/EEG/processed/Robi/Figures';
        figTitle = sprintf('%s coherence', filter1);
        title(figTitle);
        export_fig(fullfile(figFolder, sprintf('%s %s', titlePrefix, figTitle)));
    else
        a = 10;
        subplot(a,a,a*a);
        xWid = .9/length(freqLabels1);
        margin = 0.005;
        pos = [(freqCounter-1) * xWid + margin, 0, xWid - 2*margin, .5];
        set(gca, 'Position', pos);
        set(gca, 'color', 'white');
        plotChannelPairs_PF(allLabels(ind), weights{freqCounter}, [minWeight maxWeight],false,fig);
        title(sprintf('%s coherence', filter1));
        if(freqCounter == 1)
            cmap = colormap('jet');
            if(splitColormap)
                jetmap = colormap;
                jetmap_top = jetmap(1:floor(size(jetmap,1)/2),:);
                jetmap_bottom = jetmap(ceil(size(jetmap,1)/2):end,:);
                jetmap_middle = jetmap(ceil(size(jetmap,1)/2),:);
                spaces = size(jetmap,1) - ( size(jetmap_top, 1) + size(jetmap_bottom,1));
                new_jetmap = [jetmap_top ; repmat(jetmap_middle,size(jetmap_top,1)*20,1); jetmap_bottom];
                colormap(new_jetmap)
            end
            cbar = colorbar();
            set(cbar,'Position',[.93, .1, .015, .3]);
        end        
    end
end
allWeights{weightCounter} = weights;
weightCounter = weightCounter + 1;
if(plotSeparate)
    if(plotRelative)
        tilefigs([5,3]);
    else
        tilefigs([5,2]);
    end
end
end