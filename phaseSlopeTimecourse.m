function [ timeLagSeconds, meanPhaseAngle, polyLine, maxTimeShift ] = phaseSlopeTimecourse( EEG, frequencyLimits )
%PHASESLOPETIMECOURSE Summary of this function goes here
%   Detailed explanation goes here

if(~exist('frequencyLimits', 'var'))
  frequencyLimits = [1, EEG.srate / 2 - 1; ...
    1, 4;
    5, 8;
    8, 12;
    12, 15; 
    15, 20;
    20, 100];
  
end

debug = true;
% plot1 = false;
% sizeToIncrementRatio = 16;


sampleRate = EEG.srate;

windowSize = sampleRate;  %1 second window

% maxFrameShift = floor(EEG.srate/(2*maxFrequency));
% maxTimeShift = maxFrameShift / EEG.srate;
% windowIncrement = windowSize / sizeToIncrementRatio;
windowIncrement = 1;
% frequenciesOfInterest = minFrequency:maxFrequency;
% %fftIndices = frequenciesOfInterest + 1;
fftIndices = 1:windowSize/2;
% maxFft = max(fftIndices);
totalWindows = floor(size(EEG.data,2) / windowIncrement);
chanCount = size(EEG.data, 1);
phaseAngles = NaN(length(fftIndices), 1);
% phasePlot = NaN(totalWindows, length(fftIndices), chanCount);
% meanPhaseDifference = NaN(totalWindows,chanCount,chanCount);


phaseSlopeTopographies.estimatedTimeLag = NaN(totalWindows, chanCount);
phaseSlopeTopographies.estimateResiduals = NaN(totalWindows, chanCount, chanCount);
phaseSlopeTopographies.signalToNoiseRatio = NaN(totalWindows, 1);


% phaseSlopes = NaN(totalWindows,chanCount,chanCount);
phaseSlopes = NaN(chanCount,chanCount);
% phaseMeanSlopeCorr.R = NaN(chanCount,chanCount);
% phaseMeanSlopeCorr.P = NaN(chanCount,chanCount);
% if(debug)
%     allPhase1 = NaN(totalWindows, maxFft);
%     allPhase2 = NaN(totalWindows, maxFft);
% end
sampleCounter = 1;
windowCounter = 1;
while(sampleCounter + windowSize < size(EEG.data,2))
  endIndex = sampleCounter + windowSize - 1;
  fprintf('\n(%s) %d of %d', char(datetime), sampleCounter + windowSize, size(EEG.data,2));
  
  for chanIndex1 = 1:chanCount
    sample1 = EEG.data(chanIndex1, sampleCounter:endIndex);
    fft1 = fft(sample1);
    for chanIndex2 = chanIndex1:chanCount
      sample2 = EEG.data(chanIndex2, sampleCounter:endIndex);
      fft2 = fft(sample2);
      for freqCounter = fftIndices
        phase1 = atan2(imag(fft1(freqCounter)), real(fft1(freqCounter)));
        phase2 = atan2(imag(fft2(freqCounter)), real(fft2(freqCounter)));
        phaseAngle = phase2-phase1;
        phaseAngles(freqCounter) = phaseAngle;
        %positive phase angle means phase2 came first (see debug block for proof)
        %         if(debug)
        %             close all;phaseSlopeTopo
        %             allPhase1(windowCounter,freqCounter) = phase1;
        %             allPhase2(windowCounter,freqCounter) = phase2;
        %             if(plot1)
        %                 if(freqCounter == 3)
        %                     keepI = [freqCounter+1, length(fft1)-freqCounter+1];
        %                     keep = zeros(1,length(fft1));
        %                     keep(keepI) = 1;
        %                     filtered1 = fft1;
        %                     filtered2 = fft2;
        %                     filtered1(find(~keep)) = 0;
        %                     filtered2(find(~keep)) = 0;
        %                     fData1 = ifft(filtered1);
        %                     fData2 = ifft(filtered2);
        %                     figure;
        %                     hold on;
        %                     plot(fData1,'b');
        %                     plot(fData2,'r');
        %                 end %break here to compare plot to variables.
        %             end
        %         end
        
        %todo: took this out
        %         if(phaseAngle < -pi)
        %             phaseAngle = phaseAngle + 2*pi;
        %         elseif(phaseAngle > pi)
        %             phaseAngle = phaseAngle - 2*pi;
        %         endphaseSlopeTopo
        %         phaseAngles(windowCounter, freqCounter, chanIndex1, chanIndex2) = phaseAngle;
        %         phaseAngles(windowCounter, freqCounter, chanIndex2, chanIndex1) = -phaseAngle;
      end
      
      %compute phase slopes
      %       meanPhase = mean(phaseAngles);
      %       meanPhaseDifference(windowCounter,chanIndex1,chanIndex2) = meanPhase;
      %       meanPhaseDifference(windowCounter,chanIndex2,chanIndex1) = -meanPhase;
      for freqBinCounter = 1:size(frequencyLimits, 1)
        startFreq = frequencyLimits(freqBinCounter, 1);
        endFreq = frequencyLimits(freqBinCounter, 2);
        x = (startFreq:endFreq)';
        poly = polyfit(x, phaseAngles(startFreq:endFreq), 1);
        phaseSlopes(chanIndex1,chanIndex2, freqBinCounter) = poly(1);
        phaseSlopes(chanIndex2,chanIndex1, freqBinCounter) = -poly(1);
      end
      
      
    end%chanIndex2 loop
    
  end %chanIndex1 loop
  for freqBinCounter = 1:size(frequencyLimits, 1)
    phaseSlopeTopo = solveLinearCombination(squeeze(phaseSlopes(:,:, freqBinCounter)));
    phaseSlopeTopographies.estimatedTimeLag(windowCounter, :, freqBinCounter) = phaseSlopeTopo.estimatedTimeLag;
    phaseSlopeTopographies.estimateResiduals(windowCounter, :, :, freqBinCounter) = phaseSlopeTopo.estimateResiduals(:,:);
    phaseSlopeTopographies.signalToNoiseRatio(windowCounter, freqBinCounter) = phaseSlopeTopo.signalToNoiseRatio;
  end
  sampleCounter = sampleCounter + windowIncrement;
  windowCounter = windowCounter + 1;
end
%     [phaseMeanSlopeCorr.R(chanIndex1,chanIndex2), phaseMeanSlopeCorr.P(windowCounter,chanIndex1,chanIndex2)] = ...
%       corr(meanPhaseDifference(:,chanIndex1,chanIndex2), phaseSlopes(:,chanIndex1,chanIndex2));

%todo: solve for lags at each time step

save('/home/data/EEG/processed/Oregon/phaseSlopeWork.mat', '-v7.3');

phaseSlopeTopo.chanlocs = EEG.chanlocs;
phaseSlopeTopo.filename = EEG.filename;
details.windowSize = windowSize;
details.sampleRate = EEG.srate;
phaseSlopeTopo.details = details;

plot(phaseSlopeTopographies.estimatedTimeLag);
legend({EEG.chanlocs.labels});


if(windowCounter <= totalWindows)
  phaseAngles(windowCounter:end,:,:,:) = [];
  %     if(debug)
  %         allPhase1(windowCounter:end,:)=[];
  %         allPhase2(windowCounter:end,:)=[];
  %     end
end

if(debug)
  %   startIndex = 175800;
  %   endIndex = 176000;
  %   toPlot = phaseAngles(startIndex:endIndex,:)';
  for chan1 = 1:chanCount
    for chan2 = chan1+1:chanCount
      phaseAngles = phaseAngles(:,:,chan1,chan2);
      phaseAngles = phaseAngles';
      meanValues = mean(phaseAngles,1);
      for i = 1:size(phaseAngles,2)
        slice = phaseAngles(:,i)';
        poly = polyfit(1:length(slice), slice, 1);
        slopes(i) = poly(1);
      end
      [r, p] = corr(slopes', meanValues');
      close all;
      imagesc(phaseAngles);
      title(sprintf('%s-%s', EEG.chanlocs(chan1).labels, EEG.chanlocs(chan2).labels));
      colorbar;
    end
  end
end

%This step converts phase slope to an estimated time lead (only if window size is 1 second).
leadTimes = phaseSlopes ./ (2*pi);
% meanPhaseAngle = meanPhaseAngle ./ (2*pi);
% maxSlopeIndex = min(max(fftIndices),maxFrameShift);



[polyLine, myError] = polyfit(fftIndices, meanPhaseAngle(fftIndices), 1);

%positive slope means that phase 2's slope is greater than phase 1's slope.
%In other words, timecourse 1 lags behind timecourse 2.
timeLagSeconds = polyLine(1);

if(debug)
  meanPhase1 = mean(allPhase1,1);
  meanPhase2 = mean(allPhase2,1);
  close all;
  x = fftIndices;
  [angleFit, myError] = polyfit(x, meanPhaseAngle(x), 1);
  [phase1Fit, myError] = polyfit(x, meanPhase1(x), 1);
  [phase2Fit, myError] = polyfit(x, meanPhase2(x), 1);
  angleY = x .* angleFit(1) + angleFit(2);
  phase1Y = x .* phase1Fit(1) + angleFit(2);
  phase2Y = x .* phase2Fit(1) + angleFit(2);
  figure;
  hold all;
  plot(x, meanPhase1(x), 'b');
  plot(x, meanPhase2(x), 'r');
  plot(x, meanPhaseAngle(x), 'g');
  plot(x, angleY, 'b', 'linewidth', 4);
  plot(x, phase1Y, 'r', 'linewidth', 4);
  plot(x, phase2Y, 'g', 'linewidth', 4);
  legend('phase1', 'phase2', 'phaseAngle', 'p1Fit', 'p2Fit', 'angleFit');
end


end

