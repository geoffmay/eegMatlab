function [ output_args ] = icaMat2Set( input_args )
%ICAMAT2SET Summary of this function goes here
%   Detailed explanation goes here

filepath = '/Users/Geoff/Box Sync/For Geoff Mays/ICA/'; 
filename = 'VM101.1.ANT.bdfICA.mat';

loadedEEG = load(strcat(filepath, filename));
filename = strcat(filename, '.set');
icaEEG = loadedEEG.icaEEG;
%pop_saveset(icaEEG1, 'filename', filename, 'filepath', filepath);
[signal, state]=siftpipeline(icaEEG);
clear('icaEEG');
saveFile = strcat(filepath, filename, 'sift.mat');
save(saveFile);

end