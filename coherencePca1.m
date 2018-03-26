function [ chanCoefficients ] = coherencePca1( cohPca, doPlot, doTopo, minCoefficient, maxCoefficient )
%COHERENCEPCA1 Summary of this function goes here
%   Detailed explanation goes here

%display image of pair coefficients
if(false)
    imagesc(cohPca.COEFF(:, 1:10));
    colorbar;
end

%check top values
if(~exist('doPlot','var'))
    doPlot = 0;
end
if(~exist('doTopo','var'))
    doTopo = 0;
end
maxPair = length(cohPca.labelChannelOnly);
if(~exist('maxCoefficient','var'))
    maxCoefficient = 20;
end
if(~exist('minCoefficient','var'))
    minCoefficient = 20;
end

if(doTopo)
    [chanLabs, chanlocs] = antChannelLocs;
else
    [chanLabs] = antChannelLocs;
end
chanCoefficients = zeros(length(chanLabs), maxCoefficient);

M1 = strcmp(chanLabs, 'M1');
M2 = strcmp(chanLabs, 'M2');
CPz = strcmp(chanLabs, 'CPz');
remove = M1 | M2 | CPz;

for i = minCoefficient:maxCoefficient
    a = cohPca.COEFF(:, i);
    [a1, ind] = sort(a, 1, 'descend');
    if(doPlot)
        if(~doTopo)
            figure;
            plot(a1);
        end
    end
    sortLabels = cohPca.labelChannelOnly(ind);
    maxWeight = max(cohPca.COEFF(:, i));
    alpha = cohPca.COEFF(:, i) ./ maxWeight;
    i
    if(true)
        %plot sums from top channels
        for chanPairCounter = 1:maxPair
            chans = strsplit(sortLabels{chanPairCounter}, '-');
            ind1 = find(strcmp(chanLabs, chans{1}));
            ind2 = find(strcmp(chanLabs, chans{2}));
            
            chanCoefficients(ind1, i) = chanCoefficients(ind1, i) + cohPca.COEFF(chanPairCounter, i);
            chanCoefficients(ind2, i) = chanCoefficients(ind2, i) + cohPca.COEFF(chanPairCounter, i);
        end
        if(doPlot)
            figure;
            topoplot(chanCoefficients(~remove,i), chanlocs(~remove), 'maplimits', [-3 3]);
            colorbar;
        end
    end
    sortLabels(1:20)
    if(doPlot)
        if(doTopo)
            plotChannelPairs([{sortLabels(1:100)}, {sortLabels(end-100:end)}]);
            title(sprintf('component %d %s', i, cohPca.filename));
        end
    end
end

end

