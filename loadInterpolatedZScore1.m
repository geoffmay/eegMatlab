%function [ interpolated ] = loadInterpolatedZScore( filename )
%LOADINTERPOLATEDZSCORE Loads overall z scores interpolated at 128 Hz.
%   filename is the name of the events file produced by the neurofeedback
%   program.   

if(~exist('text', 'var'))
  filename = '/home/data/EEG/data/RobiPilot/RA/tx 1/630564112186073937Events.txt';
  fileId = fopen(filename);
  text = fscanf(fileId, '%c');
  fclose(fileId);
end

newlines = strfind(text, sprintf('\r\n'));
lineStart = 1;
pcaCount = 0;
compCount = 0;
feedbackCount = 0;
extraLineCount = 0;
for i = 1:length(newlines)
  fprintf('\n%f', i/length(newlines));
  line = text(lineStart:newlines(i)-1);
  if(length(strfind(line, 'comp')) > 0)
    compCount = compCount + 1;
  elseif(length(strfind(line, 'pca elem')) > 0)
    pcaCount = pcaCount + 1;

    items = strsplit(line, ' ');
    %     pcaTimeDouble = str2num(items{1});
    %     pcaTimeInt = str2num(items{2});
    
  elseif(length(strfind(line, 'avgValue')) > 0)
    feedbackCount = feedbackCount + 1;
    items = strsplit(line, ' ');
    timeDouble(feedbackCount) = str2num(items{1});
    timeInt(feedbackCount) = str2num(items{2});
    avgHeader = items{3};
    avgValue(feedbackCount) = str2num(items{4});
    momentHeader = items{5};
    momentValue(feedbackCount) = str2num(items{6});
    sizeHeader = items{7};
    if(strcmp(items{8}, '-Infinity'))
      sizeValue(feedbackCount) = -inf;
    elseif(strcmp(items{8}, 'Infinity'))
      sizeValue(feedbackCount) = inf;
    else
      sizeValue(feedbackCount) = str2num(items{8});
    end
    posHeader = items{9};
    posValue(feedbackCount) = strcmp(items{10}, 'True');
    
  else
    extraLineCount = extraLineCount + 1;
    extraLines{extraLineCount} = line;    
  end
    
  lineStart = newlines(i) + 2;
end

numColumns = 13;
numRows = length(text) / numColumns;
text1 = text;
if(numRows ~= floor(numRows))
    numRows = floor(numRows);
    text1 = text1(1:(numRows * numColumns));
end
text2 = reshape(text1, [numColumns, numRows])';

%parse the text
numeric = zeros(size(text2, 2));
numbers = NaN(size(text2));
% fprintf('parsing......');
for pairCounter = 1:size(text2, 2)
%     fprintf('%d', pairCounter);
    try
        a = str2double(text2(1, pairCounter));
        numeric(pairCounter) = 1;
    catch
    end
    if(numeric(pairCounter))
        for j = 1:size(text2, 1)
            debug = false;
            if(debug)
                quickText = text2{j,pairCounter};
                quickNum1 = parseDecimal(quickText);
                quickNum2 = str2num(quickText);
                if(quickNum1 ~= quickNum2)
                    error(sprintf('inconsistent parsing: \n%s\n%.100f\n%.100f', quickText, quickNum1, quickNum2));
                end
            end
            numbers(j,pairCounter) = parseDecimal(text2{j,pairCounter});
        end
    end
end

rawTimes = numbers(:,1);
smoothZScores = numbers(:, 3);
momentaryZScores = numbers(:, 5);
squareSizes = numbers(:, 7);
positiveFeedback = numbers(:, 9);
powerZScore = numbers(:, 11);
coherenceZScore = numbers(:, 13);


timeCourse.rawTimes = numbers(:,1);
timeCourse.smoothZScores = numbers(:, 3);
timeCourse.momentaryZScores = numbers(:, 5);
timeCourse.squareSizes = numbers(:, 7);
timeCourse.positiveFeedback = numbers(:, 9);
timeCourse.powerZScore = numbers(:, 11);
timeCourse.coherenceZScore = numbers(:, 13);


%unwrap times that have overflowed
% fprintf('processing...');
processorFrequency = 2530957;
processorReciprocal = 1 / processorFrequency;
wrapValue = power(2,32) / processorFrequency;
a = diff(numbers(:,1));
wrapPoints = find(abs(a) > 800);
for pairCounter = 1:length(wrapPoints)
    rawTimes(wrapPoints(pairCounter)+1:end) = rawTimes(wrapPoints(pairCounter)+1:end) + wrapValue;
end



%remove nans
remove = isnan(rawTimes);
a = diff(rawTimes);
decreasing = find(a < 0);
if(length(decreasing) > 0)
    remove(decreasing(1):end) = 1;
end
rawTimes(remove) = [];
smoothZScores(remove) = [];
momentaryZScores(remove) = [];
squareSizes(remove) = [];
positiveFeedback(remove) = [];
powerZScore(remove) = [];
coherenceZScore(remove) = [];

maxTime = max(rawTimes);
sampleCount = floor(maxTime * 128);
cohTimes = .5:(1/128):(.5 + (sampleCount-1) / 128);


%interpolate values
interpTotalZ = interp1(rawTimes, momentaryZScores, cohTimes);
interpSmoothZ = interp1(rawTimes, smoothZScores, cohTimes);
keep = ~isnan(interpTotalZ);

%determine whether references have dropped
% interpolated.badRefIndex = checkForBadReference(cohData.channelPairs(31).coherence);
% keep2 = keep;
% if(interpolated.badRefIndex < length(keep2))
%     keep2(interpolated.badRefIndex:end) = 0;
% end

interpolated.InstantZ = interpTotalZ(keep);
interpolated.SmoothZ = interpSmoothZ(keep);
% interpolated.goodMeanZ = interpTotalZ(keep2);
% interpolated.goodStdZ = interpTotalZ(keep2);

