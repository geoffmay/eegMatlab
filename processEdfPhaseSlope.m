% filename = 'C:\Vision\Raw Files\eegtest\export\eegtest_0001-edf.edf';
% markerName = 'C:\Vision\Raw Files\eegtest\export\eegtest_0001_Pulse Artifact Correction.Markers';
% outputFolder = 'C:\Vision\Raw Files\eegtest\export\phaseSlope\eegTest';

filename = 'C:\Vision\Raw Files\Geoff EEG test\export\GeoffTestEEG2-edf.edf';
markerName = 'C:\Vision\Raw Files\Geoff EEG test\export\GeoffTestEEG2_Pulse Artifact Correction.Markers';
outputFolder = 'C:\Users\Neuro\Documents\MATLAB\processed\MrEegLink\';

windowDuration = 1;
epochDuration = 10;


eeg = loadBrainvisionEdf(filename);
positions = loadVolumeMarkers(markerName);

startIndex = 1;
endIndex = startIndex + (windowDuration + epochDuration) * eeg.srate - 1;
slice = eeg;
sliceCounter = 1;
while(endIndex <= eeg.pnts)
    fprintf('\n%d of %d   ', endIndex, eeg.pnts);
    outputFilename = sprintf('%sphase-%08d.mat', outputFolder, sliceCounter);
    if(~exist(outputFilename, 'file'))
        placeholder = sprintf('started on %s', char(datetime));
        save(outputFilename, 'placeholder');
        slice.data = eeg.data(:, startIndex:endIndex);
        slice.pnts = size(slice.data, 2);
        slice.times = (startIndex:endIndex) ./ eeg.srate;
        topography = phaseSlopeTimecourse1(slice);
        topography.times = slice.times;
        save(outputFilename, 'topography', '-v7.3');
    end
    startIndex = startIndex + epochDuration * eeg.srate;
    endIndex = endIndex + epochDuration * eeg.srate - 1;
    sliceCounter = sliceCounter + 1;
end

% tab=tabulate(diff(positions));
% tab(tab(:,2)==0,:) = [];
% tab

