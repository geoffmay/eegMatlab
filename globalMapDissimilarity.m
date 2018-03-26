function [ gmd, EEG ] = globalMapDissimilarity( EEG )
%GLOBALFIELDPOWER Calculates a timeseries of global field power
%   Detailed explanation goes here

globalFieldPower = squeeze(std(EEG.data,1));
meanVoltage = mean(EEG.data,1);

for i = 1:(size(EEG.data,2)-1)
    map1 = (EEG.data(:,i)-meanVoltage(i))/globalFieldPower(i);
    map2 = (EEG.data(:,i+1)-meanVoltage(i+1))/globalFieldPower(i+1);
    for j = 1:length(map1)
        sqDif = (map1-map2).*(map1-map2);
    end
    gmd(i) = sqrt(sum(sqDif) / length(sqDif));
end

end

