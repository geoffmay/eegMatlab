load('/home/data/EEG/data/Oregon/wahbehVariables.mat');

filename = '/home/data/EEG/data/Oregon/Auditory1/VM101.1.Tones.bdf';

eeg = loadBdf(filename);
keepChannels = 1:31;
longEeg = eeg;

startIndex = 174000;
endIndex = 178000;

eeg.data = eeg.data(keepChannels, startIndex:endIndex);
phaseSlopeTimecourse(eeg)

latency = ([eeg.event.latency] ./ eeg.srate)';
type = [eeg.event.type]';
duration = diff([latency; (eeg.pnts / eeg.srate)]);
eventTab = table(latency, type, duration);

eventTab(eventTab(:,2) ==0,:) = [];

