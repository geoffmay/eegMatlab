function [ dataOut, labelsOut ] = pruneCoherenceMeasures( data, labels )
%PRUNECOHERENCEMEASURES Summary of this function goes here
%   Detailed explanation goes here

keepLabels = [{'Fp1'}, {'Fp2'}, ...
    {'F3'}, {'F4'}, {'C3'}, {'C4'}, ...
    {'P3'}, {'P4'}, {'O1'}, {'O2'}, ...
    {'F7'}, {'F8'}, {'T3'}, {'T4'}, {'T5'}, {'T6'}, ...
    {'Fz'}, {'Cz'}, {'Pz'}...
    {'P7'}, {'P8'}, {'T7'}, {'T8'},... 
    ];

keepInd = false(size(labels));
for i = 1:length(labels)
    label = labels{i};
    if(~strcmp(label(1:3), 'coh'))
        keepInd(i) = true;
    else
        spaces = strfind(label, ' ');
        chanPair = label(spaces(1)+1:spaces(2)-1);
        chans = strsplit(chanPair, '-');
        firstMatch = false;
        secondMatch = false;
        for j = 1:length(keepLabels)
            if(strcmp(chans{1}, keepLabels{j}))
                firstMatch = true;
            end
            if(strcmp(chans{2}, keepLabels{j}))
                secondMatch = true;
            end
        end
        if(firstMatch && secondMatch)            
            keepInd(i) = true;            
        end
    end
end

dataOut = data(keepInd);
labelsOut = labels(keepInd);

end

