function output = loadSplitConvolution(edfFilename)

edfFilename = 'GeoffTestEEG2-edf.edf';

rootFolder = 'C:\Vision\Raw Files\Geoff EEG test\export\fourierHrf\';
folder = fullfile(rootFolder, edfFilename);


files = dir(folder);
files([files.isdir]) = [];


for i = 1:length(files)
    fprintf('consolidating file %d of %d\n', i, length(files));
    filePath = fullfile(folder, files(i).name);
    data = load(filePath);
    if( i == 1)
        output.matrix = NaN(length(data.eegFourier.convolvedValues), length(files));
        output.timeSeconds = [data.eegFourier.convolvedValues.timeSeconds];
        output.labels = cell(1, length(files));
    end
    sig = [data.eegFourier.convolvedValues.signal];
    output.matrix(:, i) = sig';
    label = sprintf('%s %d-%d Hz', data.eegFourier.label, data.eegFourier.freqInfo.lowFreq, data.eegFourier.freqInfo.highFreq);
    output.labels{i} = label;
end
output.freqInfo = data.eegFourier.freqInfo;