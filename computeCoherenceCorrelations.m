function [tab, rho, p, pairs] = computeCoherenceCorrelations(coherenceTimeCourseFilename)

load(coherenceTimeCourseFilename);
labels = antChannelLocs;
labels = labels(1:32);

freqLabels = [{'delta'},{'theta'},{'alpha'},{'beta'},{'hibeta'},{'gamma'}];

%labels = [{'Fp1'},{'AF3'},{'F7'},{'F3'},{'FC1'},{'FC5'},{'T7'},{'C3'},{'CP1'},{'CP5'},{'P7'},{'P3'},{'Pz'},{'PO3'},{'O1'},{'Oz'},{'O2'},{'PO4'},{'P4'},{'P8'},{'CP6'},{'CP2'},{'C4'},{'T8'},{'FC6'},{'FC2'},{'F4'},{'F8'},{'AF4'},{'Fp2'},{'Fz'},{'Cz'},{'EXG1'},{'EXG2'},{'EXG3'},{'EXG4'},{'EXG5'},{'EXG6'},{'EXG7'},{'EXG8'},{'Ana1'},{'Ana3'},{'Status'}];
counter = 1;
for i = 1:length(labels)
    for j = i+1:length(labels)
        pairs{counter} = sprintf('%s-%s', labels{i},labels{j});
        counter = counter + 1;
    end
end


if(false)
    aLabel = 'F8-P8';
    bLabel = 'FC6-Cz';
    ai = find(strcmp(pairs, aLabel));
    bi = find(strcmp(pairs, bLabel));
    freq = 2;
    a = channelPairs(ai).coherence(:,freq);
    b = channelPairs(bi).coherence(:,freq);
    figure;
    scatter(a,b,'.');
    xlabel(aLabel);
    ylabel(bLabel);
end
coh = mean(channelPairs(31).coherence, 2);
windowSize = 1024;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
coh1 = filtfilt(b, a, coh);
threshold1 = 0.99;
supraThreshold = coh1 > threshold1;
firstFilteredSupraThreshold = find(supraThreshold);
if(length(firstFilteredSupraThreshold) == 0)
    badReferenceCutoff = length(coh1);
else
    firstFilteredSupraThreshold = firstFilteredSupraThreshold(1);
    threshold2 = 0.9;
    subThreshold = coh < threshold2;
    subThreshold(firstFilteredSupraThreshold:end) = [];
    lastUnfilteredSubThreshold = find(subThreshold);
    lastUnfilteredSubThreshold = lastUnfilteredSubThreshold(end);
    badReferenceCutoff = lastUnfilteredSubThreshold;
end

%         close all;
%         figure;
%         hold on;
%         x = ((1:length(coh1)) ./ 128)';
%         plot(x, coh1);
%         plot(x, coh, 'r');
%         pan xon;
%         zoom xon;

pairPairCount = length(channelPairs);
freqCount = size(channelPairs(1).coherence,2);
rho = zeros(pairPairCount,pairPairCount,freqCount);
p = ones(size(rho));
counter = 1;
for ai = 1:pairPairCount
    for bi = ai+1:pairPairCount
        for freq=1:freqCount
            if(mod(counter, 1000)==0)
                fprintf('.');
                if(mod(counter, 100000)==0)
                    fprintf('\n%d of %d',counter,numel(p)/2);
                end
            end
            a = channelPairs(ai).coherence(1:badReferenceCutoff,freq);
            b = channelPairs(bi).coherence(1:badReferenceCutoff,freq);
            [thisRho, thisP] = corr(a,b);
            rho(ai,bi,freq) = thisRho;
            rho(bi,ai,freq) = thisRho;
            %                     p(ai,bi,freq) = thisP;
            %                     p(bi,ai,freq) = thisP;
            counter = counter + 1;
        end
    end
end


pairPairs = cell(numel(rho),1);
flatRhos = NaN(numel(rho),1);
flatPs = NaN(numel(rho),1);

counter = 1;
for i = 1:size(rho,1)
    for j = i+1:size(rho,2)
        for k = 1:size(rho,3)
            pairPairLabel = sprintf('%s vs. %s %s', pairs{i}, pairs{j}, freqLabels{k});
            pairPairs{counter} = pairPairLabel;
            flatRhos(counter) = rho(i,j,k);
            flatPs(counter) = p(i,j,k);
            counter = counter + 1;
        end
    end
end
tab = table(pairPairs, flatRhos, flatPs);
remove = isnan(tab{:,2});
tab(remove,:) = [];
remove = zeros(size(tab,1),1);
M1rows = cellfun(@length, strfind(tab{:,1}, 'M1'));
M2rows = cellfun(@length, strfind(tab{:,1}, 'M2'));
remove = M1rows | M2rows;
tab(remove,:) = [];

end