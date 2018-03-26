function [ hasNoise, ratios ] = checkForMainsNoise( filenames, folder )
%CHECKFORMAINSNOISE Summary of this function goes here
%   Detailed explanation goes here

newVersion = true;
alsoCutShortFiles = true;

if(exist('folder','var'))
    for i = 1:length(filenames)
        filenames{i} = fullfile(folder,filenames{i});
    end
end

if(newVersion)
    sampleRate = 2048;
    minDuration = 60*5*sampleRate;
    channelCount = 34;
    
    windowDuration = 1;
    
    windowCount = 20;
    %debug
    windowCount = 200;
    plot60 = [];
    plotCounter = 1;
    %end debug
    powerThreshold = .5;
    
    for fileCounter = 1:length(filenames)
%         fprintf('\n%d of %d', fileCounter, length(filenames));
        filename = filenames{fileCounter};
        fileInfo = dir(filename);
        fileLength = fileInfo.bytes / 8;
        sampleCount = fileLength / channelCount;
        if(sampleCount < minDuration)
            hasNoise(fileCounter) = 1;
            ratios(fileCounter) = NaN;
        else
            %     if(sampleCount ~= floor(sampleCount))
            %         sampleCount = floor(sampleCount);
            %         fileLength = sampleCount * channelCount;
            %     end
            fileId = fopen(filename, 'r');
            %     contents = fread(fileId, fileLength, 'double');
            %     data = reshape(contents, sampleCount, channelCount);
            windowInterval = floor((sampleCount - windowDuration * sampleRate) / windowCount);
            windowStarts = 1 + (0:windowCount-1).*windowInterval;
            for i = 1:windowCount
%                 fprintf('.');
                windowWidth = windowDuration * sampleRate * channelCount;
                windowStart = ((windowStarts(i)-1) * channelCount) * 8;
                fseek(fileId, windowStart, -1);
                contents = fread(fileId, windowWidth, 'double');
                data = reshape(contents, channelCount, windowDuration*sampleRate)';
                
                fourier = fft(data(:,1));
                power = abs(fourier);
                power = power(2:120);
                %                 %debug
                %                 plot60(plotCounter) = power(60);
                %                 plotCounter = plotCounter + 1;
                %                 %end debug
                %                 close all;
                %                 figure;
                if(i == 1)
                    sums = power;
                else
                    sums = sums + power;
                end
            end
            fclose(fileId);
            sums = sums ./ windowCount;
            if(false)
                close all;
                figure;
                plot(sums);
                myTitle = strrep(filenames{fileCounter}, '_', ' ');
                title(myTitle);
            end
            ratio = sums(60) / sums(1);
            if(ratio > powerThreshold)
                hasNoise(fileCounter) = 1;
            else
                hasNoise(fileCounter) = 0;
            end
            ratios(fileCounter) = ratio;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%

if(~newVersion)
sampleRate = 2048;
channelCount = 34;

windowDuration = 1;
windowCount = 10;
powerThreshold = 10000 * windowCount;

for fileCounter = 1:length(filenames)
%     fprintf('\n%d of %d', fileCounter, length(filenames));
    filename = filenames{fileCounter};
    fileInfo = dir(filename);
    fileLength = fileInfo.bytes / 8;
    sampleCount = fileLength / channelCount;
    if(sampleCount ~= floor(sampleCount))
        sampleCount = floor(sampleCount);
        fileLength = sampleCount * channelCount;
    end
    fileId = fopen(filename, 'r');
    contents = fread(fileId, fileLength, 'double');
    data = reshape(contents, sampleCount, channelCount);
    windowInterval = floor((sampleCount - windowDuration * sampleRate) / windowCount);
    windowStarts = 1 + (0:windowCount-1).*windowInterval;
    for i = 1:windowCount
%         fprintf('.');
        window = windowStarts(i):(windowStarts(i)+windowDuration * sampleRate);        
        fourier = fft(data(window,1));
        power = abs(fourier);
        power = power(2:120);
        close all;
        figure;
        if(i == 1)
            sums = power;
        else
            sums = sums + power;
        end
    end
    if(false)
        close all;
        figure;
        plot(sums);
        myTitle = strrep(filenames{fileCounter}, '_', ' ');
        title(myTitle);
    end
    
    ratio = sums(60) / sums(1);
    if(ratio > powerThreshold)
        hasNoise(fileCounter) = 1;
    else
        hasNoise(fileCounter) = 0;
    end
    ratios(fileCounter) = ratio;
end
end
end

