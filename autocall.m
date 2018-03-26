function [ output_args ] = autocall( )
%Calls a batch of sift fuctions
%   Detailed explanation goes here
% clear;
% eeglab;
ALLEEG = [];
% OLDSET = [];
% eeglab;
[EEG LASTCOM] = pop_loadset('eb79RespCorr.set','/Users/Geoff/Documents/MATLAB/SIFT sample/');
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG);
siftpipeline(EEG);

% todo: get the selection of multiple datasets to work
% [EEG LASTCOM] = pop_loadset('eb79RespWrong.set','/Users/Geoff/Documents/MATLAB/SIFT sample/');
% [ALLEEG, EEG] = eeg_store(ALLEEG, EEG);
% [EEG LASTCOM] = pop_loadset('eb79RespWrong.set','/Users/Geoff/Documents/MATLAB/SIFT sample/');
% [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, OLDSET, 2);
% setSelectionArgs = cell(1,4);
% setSelectionArgs{1} = 'retrieve';
% setSelectionArgs{2} = 1:2;
% setSelectionArgs{3} = 'study';
% setSelectionArgs{4} = 0;
% popnewset(ALLEEG, EEG, OLDSET, setSelectionArgs);



eeglab('redraw');


output_args = 5;

end

