function [channelPhaseEstimates, signalToNoiseRatio, residual] = phaseSlopeTopography(slopes, epsilon, verbose)
%PHASESLOPETOPOGRAPHY calculates the phase slope between two columns of a matrix
%of timeseries data.
%   Phase slope is an approximation of the way that the phase angle between

if(~exist('epsilon','var'))
    epsilon = 0;
end
if(~exist('verbose','var'))
    verbose = false;
end
nIteration = 1000;

searchStepCount = 100;
initialStepSize = 10;


errorSize = zeros(32,32);
channelPhaseEstimates = NaN(1,32);
channelPhaseEstimates(1) = 0;
fileNumber = 1;
for chan=2:32
    estimated = channelPhaseEstimates(1:chan-1);
    actual = slopes(chan,1:chan-1, fileNumber);
    miss = actual - estimated;
    channelPhaseEstimates(chan) = mean(miss);
    for i = 1:(chan-1)
        errorSize(chan,i) = miss(i) - channelPhaseEstimates(i);
    end
end
channelPhaseEstimates = channelPhaseEstimates - mean(channelPhaseEstimates);
iterationSize = abs(channelPhaseEstimates) .* initialStepSize;
delError = realmax;
lastError = realmax;
outOfBounds = false;
iteration = 1;
zeroCounter = 0;
while(iteration < nIteration && zeroCounter < 3)
    iteration = iteration + 1;
    lastEstimate = channelPhaseEstimates;
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
    for chan = 1:32
        lastGuess = lastEstimate(chan);
        guesses = (lastGuess-iterationSize(chan)*searchStepCount/2):iterationSize(chan):(lastGuess+iterationSize(chan)*searchStepCount/2);
        errors = NaN(1,length(guesses));
        minError = realmax;
        minSignal = realmax;
        minIndex = -1;
        for guessNumber = 1:length(guesses)
            guessError = 0;
            guessSignal = 0;
            counter = 1;
            actual = NaN(1, 32 - 1);
            estimated = NaN(1, 32 - 1);
            for guessChannel = 1:32
%                 if(guessChannel < chan)
%                     estimated(counter) = -guesses(guessNumber);
%                 elseif(guessChannel > chan)
%                     estimated(counter) = guesses(guessNumber);
%                 end
                if(guessChannel~=chan)
                    estimated(counter) = guesses(guessNumber);
                    actual(counter) = slopes(guessChannel,chan,fileNumber);
                    channelError = abs(actual(counter)-estimated(counter));
                    channelError = channelError * channelError;
                    guessError = guessError + channelError;
                    guessSignal = guessSignal + abs(guesses(guessNumber));
                end
                counter = counter + 1;
            end
            if(guessError < minError)
                minError = guessError;
                minIndex = guessNumber;
                minSignal = guessSignal;
            end
            if(false)
                if(iteration > 30)
                    figure;
                    hold on;
                    plot(actual,'b');
                    plot(estimated,'r');
                    title(sprintf('channel %d', chan));
                end
            end
            errors(guessNumber) = guessError;
        end
        if(false)
            figure;
            hold on;
            plot(errors,'b');
        end
        totalError = totalError + minError;
        totalSignal = totalSignal + minSignal;
        if(minIndex > 0)
            channelPhaseEstimates(chan) = guesses(minIndex);
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
    if(~outOfBounds)
        for chan=1:32
            iterationSize(chan) = iterationSize(chan) / 2;
        end
    end
    if(verbose)
        fprintf('\niteration #%d, error = %f, delError = %d', iteration, totalError, delError);
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
if(nargout > 2)
    if(verbose)
        figure;
    end
    residual = slopes;
    originalPower = 0;
    residualPower = 0;
    for chan1 = 1:32
        for chan2 = (chan1+1):32
            value = slopes(chan2,chan1,fileNumber);
            originalPower = originalPower + value * value;
            value = value - channelPhaseEstimates(chan1) + channelPhaseEstimates(chan2);
            %             residual(chan2,chan1,fileNumber) = value;
            residualPower = residualPower + value * value;
            residual(chan1,chan2,fileNumber) = value;
        end
    end
    originalPower = sqrt(originalPower);
    residualPower = sqrt(residualPower);
    powerRatio = originalPower / residualPower;
    if(verbose)
        imagesc(residual(:,:,fileNumber));
        colorbar;
        fprintf('original power: %f; residual power: %f; ratio: %f', originalPower, residualPower, powerRatio);
    end
    for chan1 = 1:32
        for chan2 = (chan1+1):32
            value = slopes(chan2,chan1,fileNumber);
            value = value - channelPhaseEstimates(chan1) + channelPhaseEstimates(chan2);
            residual(chan2,chan1,fileNumber) = value;
        end
    end

end
end
