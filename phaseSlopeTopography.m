function output = phaseSlopeTopography(eeg)


if(~exist('verbose','var'))
  verbose = false;
end

sampleRate = eeg.srate;
windowIncrement = eeg.srate / 16;
subSampleStartFreq = 20;
subSampleEndFreq = 100;

%compute phase slopes
for chanIndex1 = 1:eeg.nbchan
  for chanIndex2 = (chanIndex1+1):eeg.nbchan
    if(verbose)
      fprintf('\n%s computing phase angle between %s and %s', char(datetime), eeg.chanlocs(chanIndex1).labels, eeg.chanlocs(chanIndex2).labels);
    end
    sampleCounter = 1;
    totalWindows = floor(size(eeg.data,2) / sampleRate);
    windowCounter = 1;
    phaseAngles = NaN(totalWindows, subSampleEndFreq);
    
    while(sampleCounter + sampleRate < size(eeg.data,2))
      endIndex = sampleCounter + sampleRate - 1;
      sample1 = eeg.data(chanIndex1, sampleCounter:endIndex);
      sample2 = eeg.data(chanIndex2, sampleCounter:endIndex);
      fft1 = fft(sample1);
      fft2 = fft(sample2);
      for i = subSampleStartFreq:subSampleEndFreq
        phase1 = atan2(imag(fft1(i)), real(fft1(i)));
        phase2 = atan2(imag(fft2(i)), real(fft2(i)));
        if(phase1 < pi)
          phase1 = phase1 + 2*pi;
        elseif(phase1 > pi)
          phase1 = phase1 - 2*pi;
        end
        if(phase2 < pi)
          phase2 = phase2 + 2*pi;
        elseif(phase2 > pi)
          phase2 = phase2 - 2*pi;
        end
        phaseAngles(windowCounter, i) = phase2 - phase1;
      end
      sampleCounter = sampleCounter + windowIncrement;
      windowCounter = windowCounter + 1;
    end
    meanPhaseAngle = mean(phaseAngles,1);
    [polyLine, error] = polyfit(subSampleStartFreq:subSampleEndFreq, meanPhaseAngle(subSampleStartFreq:subSampleEndFreq), 1);
    slope = polyLine(1);
    output.slopes(chanIndex1, chanIndex2) = slope;
    output.slopes(chanIndex2, chanIndex1) = -slope;
    output.phasePlots{chanIndex1, chanIndex2} = meanPhaseAngle;
  end
end

%estimate topography
if(~exist('epsilon','var'))
  epsilon = 0;
end
nIteration = 1000;

%the number of guesses per iteration
guessesPerIteration = 100;
initialGuessRangeRatio = 10;

errorSize = zeros(eeg.nbchan);
%this contains the estimated time lag for each channel
channelLagEstimates = NaN(1,eeg.nbchan);
%set the first channel a zero to fully constrain
channelLagEstimates(1) = 0;
%set other channels as the mean of the remaining error
for chan=2:eeg.nbchan
  estimateLag = channelLagEstimates(1:chan-1);
  %actual = (chan,1:chan-1, fileNumber);
  actualSlope = output.slopes(chan,1:chan-1);
  estimateError = actualSlope - estimateLag;
  channelLagEstimates(chan) = mean(estimateError);
  for i = 1:(chan-1)
    errorSize(chan,i) = estimateError(i) - channelLagEstimates(i);
  end
end
%normalize initial estimates
channelLagEstimates = channelLagEstimates - mean(channelLagEstimates);
%each channel estimate changes proportional to its error
iterationSize = abs(channelLagEstimates) .* initialGuessRangeRatio;
%track the change in error size from one estimate to the next
delError = realmax;
%track the size of the total error
lastError = realmax;
%if best guess is near either end, we need to expand the search space
outOfBounds = false;
%track number of iterations
iteration = 1;
%track number of consecutive estimations that are within the acceptable
%error
zeroCounter = 0;
%loop until we've hit max iterations or within error tolerance
while(iteration < nIteration && zeroCounter < 3)
  iteration = iteration + 1;
  %store the previous estimate before modifications
  lastEstimate = channelLagEstimates;
  if(iteration > 3)
    delError = abs(totalError - lastError);
  end
  if(iteration > 2)
    lastError = totalError;
  end
  if(delError <= epsilon)
    zeroCounter = zeroCounter + 1;
  else
    zeroCounter = 0;
  end
  totalError = 0;
  totalSignal = 0;
  %loop through the channels
  for chan = 1:eeg.nbchan
    %generate a list of guesses in a range surrounding the last
    %estimate
    lastGuess = lastEstimate(chan);
    minGuess = lastGuess-iterationSize(chan)*guessesPerIteration/2;
    maxGuess = lastGuess+iterationSize(chan)*guessesPerIteration/2;
    guessIncrement = iterationSize(chan);
    guesses = minGuess:guessIncrement:maxGuess;
    %     guesses = (lastGuess-iterationSize(chan)*guessesPerIteration/2):iterationSize(chan):(lastGuess+iterationSize(chan)*guessesPerIteration/2);
    %track the error for each guess
    errors = NaN(1,length(guesses));
    minError = realmax;
    minSignal = realmax;
    minIndex = -1;
    for guessNumber = 1:length(guesses)
      
      %       oldErrorCalculation = true;
      %       if(~oldErrorCalculation)
      %         guess = guesses(guessNumber);
      %         virtual = lastEstimate + lastGuess - guess;
      %         virtual(chan) = 0;
      %         realSlopes = output.slopes(chan, :);
      %         guessError = sum(abs(realSlopes- virtual));
      %         guessSignal = sum(abs(virtual));
      %
      %       else
      guessError = 0;
      guessSignal = 0;
      guessChannelCounter = 1;
      actualSlope = NaN(1, eeg.nbchan - 1);
      estimateLag = NaN(1, eeg.nbchan - 1);
      for guessChannel = 1:eeg.nbchan
        doFlip = false;
        if(doFlip)
          if(guessChannel < chan)
            estimateLag(guessChannelCounter) = -guesses(guessNumber);
          elseif(guessChannel > chan)
            estimateLag(guessChannelCounter) = guesses(guessNumber);
          end
        else
          estimateLag(guessChannelCounter) = guesses(guessNumber);
        end
        if(guessChannel~=chan)
          actualSlope(guessChannelCounter) = output.slopes(guessChannel,chan);
          channelError = abs(actualSlope(guessChannelCounter)-estimateLag(guessChannelCounter));
          channelError = channelError * channelError;
          guessError = guessError + channelError;
          guessSignal = guessSignal + abs(guesses(guessNumber));
        end
        guessChannelCounter = guessChannelCounter + 1;
      end
      %       end
      if(guessError < minError)
        minError = guessError;
        minIndex = guessNumber;
        minSignal = guessSignal;
      end
      errors(guessNumber) = guessError;
    end
    %     %debug
    %     if(chan == 3)
    %       fprintf('\nchannel: %f, minIndex: %f, minError: %f', chan, minIndex, minError);
    %       close all;
    %       plot(guesses, errors);
    %     end
    %     %end debug
    if(false)
      figure;
      hold on;
      plot(errors,'b');
    end
    totalError = totalError + minError;
    totalSignal = totalSignal + minSignal;
    if(minIndex > 0)
      channelLagEstimates(chan) = guesses(minIndex);
    else
      dummy = 1;
    end
    %guesses on the edge mean the territory should be expanded.
    if(minIndex < 3 || minIndex > (length(guesses)-2))
      iterationSize(chan) = iterationSize(chan) * 16;
      outOfBounds = true;
    end
  end
  %if all guesses were in bounds, we can narrow our search next time.
  %if(~outOfBounds)
  for chan=1:eeg.nbchan
    iterationSize(chan) = iterationSize(chan) / 2;
  end
  %end
  outOfBounds = false;
  if(verbose)
    fprintf('\niteration #%d, minguess = %f, error = %f, delError = %d', iteration, minGuess, totalError, delError);
    if(minIndex < 3 || minIndex > (length(guesses)-2))
      fprintf('(OOB)');
    end
  end
end
if(verbose)
  fprintf('\n');
end
signalToNoiseRatio = sqrt(totalSignal / totalError);
if(false)
  figure;
  plot(channelPhaseEstimates');
  legends = cell(1,size(channelPhaseEstimates,1));
  for i = 1:length(legends)
    legends{i} = num2str(i);
  end
  legend(legends);
end
%if(nargout > 2)
if(verbose)
  figure;
end
residual = output.slopes;
originalPower = 0;
residualPower = 0;
for chan1 = 1:eeg.nbchan
  for chan2 = (chan1+1):eeg.nbchan
    %             value = slopes(chan2,chan1,fileNumber);
    value = output.slopes(chan2,chan1);
    originalPower = originalPower + value * value;
    value = value - channelLagEstimates(chan1) + channelLagEstimates(chan2);
    %             residual(chan2,chan1,fileNumber) = value;
    residualPower = residualPower + value * value;
    %             residual(chan1,chan2,fileNumber) = value;
    residual(chan1,chan2) = value;
  end
end
originalPower = sqrt(originalPower);
residualPower = sqrt(residualPower);
powerRatio = originalPower / residualPower;
if(verbose)
  %         imagesc(residual(:,:,fileNumber));
  imagesc(residual(:,:));
  colorbar;
  fprintf('original power: %f; residual power: %f; ratio: %f', originalPower, residualPower, powerRatio);
end
for chan1 = 1:eeg.nbchan
  for chan2 = (chan1+1):eeg.nbchan
    value = output.slopes(chan2,chan1);
    value = value - channelLagEstimates(chan1) + channelLagEstimates(chan2);
    residual(chan2,chan1) = value;
  end
end
output.estimatedTimeLag = channelLagEstimates;
output.estimateResiduals = residual;
output.signalToNoiseRatio = signalToNoiseRatio;
output.chanlocs = eeg.chanlocs;

%end
%end
%
% %debug
% close all;
% output
% figure;
% imagesc(output.slopes);
% colorbar;
% title('actual phase slopes');
%
% figure;
% imagesc(output.estimateResiduals);
% colorbar;
% title('residuals');


