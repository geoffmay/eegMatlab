function output = solveLinearCombination(matrix, verbose)

% %for debug purposes
% matrix = squeeze(phaseSlopes(100,:,:));

%if this is a half-matrix, fill the rest in
for i = 1:size(matrix, 1)
  for j = 1:size(matrix,2)
    if(isnan(matrix(i,j)))
      if(i == j)
        matrix(i,j) = 0;
      else
        matrix(i,j) = -matrix(j,i);
      end
    end
  end
end


if(~exist('verbose','var'))
  verbose = false;
end


%estimate topography
if(~exist('epsilon','var'))
  epsilon = 0;
end
nIteration = 1000;

%the number of guesses per iteration
% guessesPerIteration = 100;
guessesPerIteration = 9;
initialGuessRangeRatio = 10;

errorSize = zeros(size(matrix, 1));
%this contains the estimated time lag for each channel
channelLagEstimates = NaN(1,size(matrix, 1));
%set the first channel a zero to fully constrain
channelLagEstimates(1) = 0;
%set other channels as the mean of the remaining error
for chan=2:size(matrix,1)
  estimateLag = channelLagEstimates(1:chan-1);
  %actual = (chan,1:chan-1, fileNumber);
  actualSlope = matrix(chan,1:chan-1);
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
  for chan = 1:size(matrix, 1)
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
      %         realSlopes = matrix(chan, :);
      %         guessError = sum(abs(realSlopes- virtual));
      %         guessSignal = sum(abs(virtual));
      %
      %       else
      guessError = 0;
      guessSignal = 0;
      guessChannelCounter = 1;
      actualSlope = NaN(1, size(matrix, 1) - 1);
      estimateLag = NaN(1, size(matrix, 1) - 1);
      for guessChannel = 1:size(matrix, 1)
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
          actualSlope(guessChannelCounter) = matrix(guessChannel,chan);
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
  for chan=1:size(matrix, 1)
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
residual = matrix;
originalPower = 0;
residualPower = 0;
for chan1 = 1:size(matrix, 1)
  for chan2 = (chan1+1):size(matrix, 1)
    %             value = slopes(chan2,chan1,fileNumber);
    value = matrix(chan2,chan1);
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
  myTitle = sprintf('original power: %f; residual power: %f; ratio: %f', originalPower, residualPower, powerRatio);
  title(myTitle);
end
for chan1 = 1:size(matrix, 1)
  for chan2 = (chan1+1):size(matrix, 1)
    value = matrix(chan2,chan1);
    value = value - channelLagEstimates(chan1) + channelLagEstimates(chan2);
    residual(chan2,chan1) = value;
  end
end
output.estimatedTimeLag = channelLagEstimates;
output.estimateResiduals = residual;
output.signalToNoiseRatio = signalToNoiseRatio;
% output.chanlocs = eeg.chanlocs;

%end
%end
%
% %debug
% close all;
% output
% figure;
% imagesc(matrix);
% colorbar;
% title('actual phase slopes');
%
% figure;
% imagesc(output.estimateResiduals);
% colorbar;
% title('residuals');


