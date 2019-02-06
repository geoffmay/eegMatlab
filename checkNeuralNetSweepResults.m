folder = 'C:\Users\Neuro\source\repos\geoffmay\NeuralNet\Testing\bin\Debug\netcoreapp2.0\Hidden Sweep';

files = dir(folder);
counter = 1;
for i = 3:length(files)
    fprintf('\nfile %d of %d', i, length(files));
    fileId = fopen(fullfile(folder, files(i).name));
    text = fscanf(fileId, '%c');
    fclose(fileId);
    json = jsondecode(text);
    hidden(counter) = json.hiddenNeuronCount;
    trainErrors(counter) = json.trainError;
    testErrors(counter) = json.testError;
    trainingSeconds(counter) = json.trainingSeconds;
    counter = counter + 1;
end

eegFile = 'C:\Users\Neuro\Downloads\GhermanPhiliastides\sub-GhermanPhiliastides01\EEG\EEG_data_sub-01_run-01.mat';
eeg = load(eegFile);
close all;
srate = eeg.EEGdata.
x = 
plot(eeg.EEGdata.Y(1,:));



