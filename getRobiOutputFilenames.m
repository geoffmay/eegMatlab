function [ outputFilenames ] = getRobiOutputFilenames( inputFilenames, outputFolder )
%GETOUTPUTFILENAMES Summary of this function goes here
%   Detailed explanation goes here

% outputFolder = '/home/data/EEG/processed/Robi/phaseSlope/';

for i = 1:length(inputFilenames)
  filename= inputFilenames{i};
  outputStart = strfind(filename, 'ROBI_');
  if(length(outputStart) < 1)
    error('invalid filename');
  end
  outFile = filename(outputStart(1):end);
  outFile = strrep(outFile,'/','_');
  outFile = strrep(outFile,'.eegData','.mat');
  outputPath = fullfile(outputFolder, outFile);
  outputFilenames{i} = outputPath;
end
end

