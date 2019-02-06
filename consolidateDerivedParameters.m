function consolidateDerivedParameters(edfFilename)
%compute phase slope info
% edfFilename = 'GeoffTestEEG2-edf.edf';
inputPath = fullfile('C:\Vision\Raw Files\Geoff EEG test\history\', edfFilename);
eeg = loadBrainvisionEdf(inputPath);
phaseIntermediateFolder = fullfile('C:\Vision\Raw Files\Geoff EEG test\export\phaseSlope\', edfFilename);
phaseSlopeTopographies = phaseSlopeTimecourse2( eeg, phaseIntermediateFolder );

for i = 1:size(phaseSlopeTopographies.estimatedTimeLag, 2)
    fprintf('convolving channel %d of %d (%s)\n', i, size(phaseSlopeTopographies.estimatedTimeLag, 1), char(datetime));
    convPhase(:, i) = convolveHrf(phaseSlopeTopographies.estimatedTimeLag(:, i, 1));
end

% consolidatedPhase = 'C:\Vision\Raw Files\Geoff EEG test\export\phaseSlope\GeoffTestEEG2-phaseSlope.mat';
consolidatedPhase = fullfile('C:\Vision\Raw Files\Geoff EEG test\history', [edfFilename, 'phase.mat']);
save(consolidatedPhase, 'phaseSlopeTopographies', 'convPhase', '-v7.3');

[cohInfo, summary] = deriveRobiCoherenceMatrix(eeg);
consolidatedCohFilename  = fullfile('C:\Vision\Raw Files\Geoff EEG test\history', [edfFilename, 'coh.mat']);
save(consolidatedCohFilename, 'cohInfo', 'summary', '-v7.3');


toPlot = squeeze(phaseSlopeTopographies.signalToNoiseRatio);
plot(toPlot);
title('signal to noise ratio');
legend({'1 (1:125)', '2 (20:60)'});

%test
tic;
x = phaseSlopeTopographies.estimatedTimeLag(:, 1, 1);
tr = 1/eeg.srate;
spmHrf = spm_hrf(tr)';

x = rand(1, 100000)';
tic;
convolveResults = convolve(x, spmHrf);
convolveTime = toc;
n1 = 1:length(x);
n2 = 1:length(spmHrf);
tic;
convolutionResults = convolution(x', n1, spmHrf, n2);
convolutionTime = toc;
convolveTime
convolutionTime

% 
% %consolidate phase slope info
% phaseFiles = dir(outputPhase);
% phaseFiles([phaseFiles.isdir]) = [];
% tic;
% fprintf('---');
% for i = 1:length(phaseFiles)
%     fprintf('\b\b\b%02d%%', floor(i/length(phaseFiles)*100));
% 
%     x1 = strrep(phaseFiles(i).name, 'phaseSlope-', '');
%     x1 = strrep(x1, '.mat', '');
%     x(i) = str2num(x1);
%     a = load(fullfile(outputPhase, phaseFiles(i).name));
%     if(isfield(a, 'phaseSlopeTopographies'))
%         if(i == 1)
%             phaseSlopeTopographies.estimatedTimeLag = NaN(length(phaseFiles), size(a.phaseSlopeTopographies.estimatedTimeLag, 2), size(a.phaseSlopeTopographies.estimatedTimeLag, 3));
%             phaseSlopeTopographies.signalToNoiseRatio = NaN(length(phaseFiles), size(a.phaseSlopeTopographies.signalToNoiseRatio,1), size(a.phaseSlopeTopographies.signalToNoiseRatio,2));
%         end
%         phaseSlopeTopographies.estimatedTimeLag(i, :, :) = a.phaseSlopeTopographies.estimatedTimeLag;
%         phaseSlopeTopographies.signalToNoiseRatio(i, :, :) = a.phaseSlopeTopographies.signalToNoiseRatio;
%     end
% end
% time = toc;
% filesPerSecond = length(phaseFiles) / time;
doPlot = true;
if(doPlot)
    count = size(phaseSlopeTopographies.estimatedTimeLag, 1);
%     x = (1:count) ./ 250;
    figure;
    plot(x, phaseSlopeTopographies.estimatedTimeLag(:, 1, 1));
end


folder = 'C:\Vision\Raw Files\Geoff EEG test\export\corrLag\GeoffTestEEG2-edf.edf';
outputFile = 'C:\Vision\Raw Files\Geoff EEG test\export\corrLag\GeoffTestEEG2-edf.mat';

fprintf('\n(%s) loading correlation matrices from folder: %s', char(datetime), folder);
files = dir(folder);

counter = 1;
for i = 3:length(files)
    fprintf('\nfile %d of %d', i, length(files));
    data = load(fullfile(folder, files(i).name));
    if(isfield(data, 'summary'))
    out.startFrame = data.summary.parameters.windowStartFrame;
    out.endFrame = data.summary.parameters.windowStartFrame;
    out.lagMatrix = data.summary.lagMatrix;
    if(counter == 1)
        output.frames = repmat(out, length(files)-2);
    end
    output.frames(counter) = out;
    output.sampleRate = data.summary.parameters.sampleRate;
    counter = counter + 1;
    end
end

save(outputFile, 'output', '-v7.3');

inputPath = 'C:\Vision\Raw Files\Geoff EEG test\history\GeoffTestEEG2-edf.edf';
eeg = loadBrainvisionEdf(inputPath);
phaseSlopeTopographies = phaseSlopeTopography( eeg );
outputFile = 'C:\Vision\Raw Files\Geoff EEG test\export\GeoffTestEEG2-phaseSlope.mat';
save(outputFile, 'phaseSlopeTopographies', '-v7.3');

inputPath = 'C:\Vision\Raw Files\Geoff EEG test\history\GeoffTestEEG2-edf.edf';
eeg = loadBrainvisionEdf(inputPath);
phaseIntermediateFolder = 'C:\Vision\Raw Files\Geoff EEG test\export\phaseSlope';
phaseSlopeTimecourse2( eeg, phaseIntermediateFolder )



