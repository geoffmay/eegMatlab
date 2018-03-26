function testBatch1()

[tones, ant, other] = getBdfFiles();

    fileID = fopen('event.txt', 'w');
fclose(fileID)

for i = 1:length(tones)
    try
        EEG = loadBdf(tones{i}, '44');
    catch exception
        EEG = loadBdf(tones{i}, '209');
    end
    fclose('all');
%     channelNumber = 1;
%     [latency, height, area] = maxErpSample(EEG, channelNumber);
    fileID = fopen('event.txt', 'a');
    formatSpec = '%d,%d,%d,%s\r\n';
    for j = 1:length(EEG.event)
        fprintf(fileID, formatSpec, EEG.event(j).latency, EEG.event(j).type, EEG.event(j).urevent, tones{i});
    end
    fclose(fileID);
end

end
