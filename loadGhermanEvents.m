folder = 'C:\Users\Neuro\Documents\MATLAB\data\GhermanPhilastides\EEG_events_volumemarkers';
file = fullfile(folder, 'efMRIcap64.elp');
chanlocs = readlocs(file);

folder = 'C:\Users\Neuro\Documents\MATLAB\data\GhermanPhilastides\EEG_events_volumemarkers\new_EEG_events_volumemarkers';
file = fullfile(folder, 'EEG_events_sub-01_run-01.mat');
data = load(file);


