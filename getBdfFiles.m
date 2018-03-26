function [ tones, ant, other ] = getBdfFiles( )
%GETBDFFILES Summary of this function goes here
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


end

