function [ testEEG ] = loadBdf( bdfFileName, epochStartValue )
%LOADBDF Summary of this function goes here
%   Detailed explanation goes here

% inputFolder = '/Users/Geoff/Box Sync/For Geoff Mays/VET MIND EEG files/';
% outputFolder = '/Users/Geoff/Box Sync/For Geoff Mays/Sets/';
inputFolder = '';
outputFolder = '';


refChan = int32(32);
if(~exist('pop_readbdf', 'file'))
  eeglab;
  close all;
end
testEEG = pop_readbdf(strcat(inputFolder,bdfFileName), {}, 43, refChan, false);

if(exist('epochStartValue', 'var'))
  testEEG.data(refChan,:)=[];
  testEEG.chanlocs(refChan,:)=[];
  testEEG.nbchan = testEEG.nbchan-1;
  disp('...reference channel removed.');
  
  
  testEEG = pop_eegfiltnew(testEEG, 1, 15, 3390, 0, [], 0, 0);
  testCell{1}=epochStartValue;
  testEEG = pop_epoch(testEEG, testCell, [-1,2]);
  testEEG = pop_rmbase(testEEG, [-1000 0]);
end

testEEG = setStandardLocations(testEEG);
% pop_saveset(testEEG, strcat(outputFolder, bdfFileName));
end

