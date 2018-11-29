function [ phaseSlopeTopographies ] = phaseSlopeTimecourse1( EEG, frequencyLimits )
%PHASESLOPETIMECOURSE Summary of this function goes here
%   Detailed explanation goes here

if(~exist('frequencyLimits', 'var'))
  frequencyLimits = [1, EEG.srate / 2 - 1; ...
%     1, 4; ...
%     5, 8; ...
%     8, 12; ...
%     12, 15; ...
%     15, 20; ...
    20, 60];
  
end
windowDuration = 1/20;

saveResiduals = false;
debug = false;

sampleRate = EEG.srate;
windowSize = floor(sampleRate * windowDuration);  
frequencyLimits = floor(frequencyLimits .* windowDuration);
frequencyLimits(frequencyLimits == 0) = 1;

windowIncrement = 1;
fftIndices = 1:windowSize/2;
totalWindows = floor(size(EEG.data,2) / windowIncrement);
chanCount = size(EEG.data, 1);
phaseAngles = NaN(length(fftIndices), 1);
phaseSlopeTopographies.estimatedTimeLag = NaN(totalWindows, chanCount);
if(saveResiduals)
  phaseSlopeTopographies.estimateResiduals = NaN(totalWindows, chanCount, chanCount);
end
phaseSlopeTopographies.signalToNoiseRatio = NaN(totalWindows, 1);
phaseSlopes = NaN(chanCount,chanCount);
sampleCounter = 1;
windowCounter = 1;
fprintf('  (------')
while(sampleCounter + windowSize < size(EEG.data,2))
  fprintf('\b\b\b\b\b\b%05d)', windowCounter)
  endIndex = sampleCounter + windowSize - 1;
  if(debug)
    fprintf('\n(%s) %d of %d', char(datetime), sampleCounter + windowSize, size(EEG.data,2));
  end
  
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
    if(saveResiduals)
      phaseSlopeTopographies.estimateResiduals(windowCounter, :, :, freqBinCounter) = phaseSlopeTopo.estimateResiduals(:,:);
    end
    phaseSlopeTopographies.signalToNoiseRatio(windowCounter, freqBinCounter) = phaseSlopeTopo.signalToNoiseRatio;
  end
  sampleCounter = sampleCounter + windowIncrement;
  windowCounter = windowCounter + 1;
end

%save('/home/data/EEG/processed/Oregon/phaseSlopeWork.mat', '-v7.3');

phaseSlopeTopographies.chanlocs = EEG.chanlocs;
phaseSlopeTopographies.filename = EEG.filename;
details.windowSize = windowSize;
details.sampleRate = EEG.srate;
phaseSlopeTopographies.details = details;
phaseSlopeTopographies.frequencyLimits = frequencyLimits;

doPlot = false;
if(doPlot)
  if(false)
    %plot multiple frequencies
    channelIndex = 17;
    toPlot = squeeze(phaseSlopeTopographies.estimatedTimeLag(:,channelIndex,:));
    myLegend = cell(0);
    for i = 1:size(frequencyLimits,1)
      myLegend{i} = sprintf('%d-%d', frequencyLimits(i, 1), frequencyLimits(i, 2));
    end
  elseif(false)
    %plot multiple channels
    frequencyIndex = 3;
    channelIndices = [1 17];
    toPlot = squeeze(phaseSlopeTopographies.estimatedTimeLag(:,channelIndices,frequencyIndex));
    myLegend = cell(0);
    for i = 1:length(channelIndices)
      myLegend{i} = sprintf('%s 5-8Hz', phaseSlopeTopographies.chanlocs(channelIndices(i)).labels);
    end
  elseif(false)
    %cluster
    frequencyIndex = 3;
    toPlot = squeeze(phaseSlopeTopographies.estimatedTimeLag(:,:,frequencyIndex));
    Z = linkage(toPlot);
    figure;
    dendrogram(Z,0)
  end
  
  figure;
  x = ((1:size(toPlot, 1)) ./ phaseSlopeTopographies.details.sampleRate)';
  plot(x, toPlot);
  title(sprintf('phase slope height %s', phaseSlopeTopographies.chanlocs(channelIndex).labels));
  legend(myLegend);
  pan xon;
  zoom xon;
end


end

