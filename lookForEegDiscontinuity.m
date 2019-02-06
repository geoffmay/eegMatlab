filename = 'C:\Users\Neuro\Downloads\GhermanPhiliastides\sub-GhermanPhiliastides05\EEG\EEG_data_sub-05_run-02.mat';

rawEeg = load(filename);
eegDiff = diff(rawEeg.EEGdata.Y, 1, 2);
avgDiff = mean(abs(eegDiff));