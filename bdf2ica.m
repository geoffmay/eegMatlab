function icaEEG = bdf2ica( bdfFileName, epochStartValue)%, pcaComponentCount)%, inputFolder, outputFolder )
%BDF2ICA Is specialized to process Wahbeh's files
%   Detailed explanation goes here


% if(~isnumeric(pcaComponentCount))
%     error('pcaComponentCount must be numeric'); 
% end
inputFolder = '/Users/Geoff/Box Sync/For Geoff Mays/VET MIND EEG files/';
outputFolder = '/Users/Geoff/Box Sync/For Geoff Mays/ICA/';
refChan = int32(32);
testEEG = pop_readbdf(strcat(inputFolder,bdfFileName), {}, 43, refChan, false);
testEEG.data(refChan,:)=[];
testEEG.chanlocs(refChan,:)=[];
testEEG.nbchan = testEEG.nbchan-1;
disp('...and now its as though that channel never existed.');
testCell{1}=epochStartValue;
testEEG = pop_epoch(testEEG, testCell, [-1,2]);
lastOption = {'extended' 1};
%icaEEG = pop_runica(testEEG, 'icatype', 'runica', 'dataset', 1, 'pca', pcaComponentCount, 'options', lastOption);
icaEEG = pop_runica(testEEG, 'icatype', 'runica', 'pca', 10, 'dataset', 1, 'options', lastOption);
save(strcat(outputFolder,bdfFileName,'Pca.mat'), 'icaEEG');
end
