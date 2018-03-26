function [ output_args ] = testBatch( )
%TESTBATCH Summary of this function goes here
%   Detailed explanation goes here
folder = dir('/Users/Geoff/Box Sync/For Geoff Mays/VET MIND EEG files');
outputDir = '/Users/Geoff/Box Sync/For Geoff Mays/ICA/';
tones = cell(0);
ant = cell(0);
other = cell(0);
for i = 1:length(folder)
    fileName = folder(i).name;
    if(strfind(fileName, 'Tones'))
        tones{length(tones)+1} = fileName;
    elseif(strfind(fileName, 'ANT'))
        ant{length(ant)+1}=fileName;
    elseif(strfind(fileName, '_T1'))
        tones{length(tones)+1} = fileName;
    elseif(strfind(fileName, '_T2'))
        tones{length(tones)+1} = fileName;
    elseif(strfind(fileName, 'tones'))
        tones{length(tones)+1} = fileName;
    else
        other{length(other)+1}=fileName;            
    end
end

eventDefinition = '44';
i = 1;

icaEEG = bdf2ica(tones{i}, eventDefinition);
icaEEG = setStandardLocations(icaEEG);
% loadedEEG = load(strcat(outputDir, 'VM101.1.ANT.bdfPca.mat'));
% icaEEG = loadedEEG.icaEEG;
save(strcat(outputDir, 'pipeline.mat'));
[signal, state] = siftpipeline(icaEEG);
icaEEG.CAT = signal.CAT;
icaEEG.filepath = outputDir;
icaEEG.filename = tones{i};
icaEEG.filename = strcat(icaEEG.filename, 'conn1.set');
pop_saveset(icaEEG, 'filename', icaEEG.filename, 'filepath', icaEEG.filepath);
save(strcat(outputDir, 'pipeline.mat'));

end

