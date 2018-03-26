filenames = {'/home/data/EEG/data/ROBI/ROBI_003/tx 8/630176882864679227.eegData'};

epochInterval = 2048;
epochLength = 2048*4;

channelCount = 34;
if(~exist('chanlocs','var'))
    [labels chanlocs] = antChannelLocs;
end
for fileCounter = 1:length(filenames)
    filename = filenames{fileCounter};
    if(~exist('data','var'))
        fprintf('%s', filename);
        file = dir(filename);
        fileLength = file.bytes / 8;
        sampleCount = fileLength / channelCount;
        if(sampleCount ~= floor(sampleCount))
            sampleCount = floor(sampleCount);
            fileLength = sampleCount * channelCount;
        end
        summary.filename = filename;
        
        fileId = fopen(filename, 'r');
        contents = fread(fileId, fileLength, 'double');
        fclose(fileId);
        data = reshape(contents, channelCount, fileLength / channelCount)';
    end
    
    epochCount = floor(sampleCount / epochInterval);
    phaseHeight = NaN(epochCount,32);
    snr = NaN(epochCount,1);
    residual = NaN(32,32,epochCount);
    epochCounter = 1;
    while(epochCounter <= epochCount)
        fprintf('.');
        if(mod(epochCounter,100)==0)
            fprintf('\n%d/%d',epochCounter,epochCount);
        end
        epochStart = (epochCounter-1) * epochInterval + 1;
        epochEnd = epochStart + epochLength - 1;
        epoch = data(epochStart:epochEnd,:);
        phaseSlopeMatrix = zeros(32,32);
        for chan1 = 1:32
            for chan2 = chan1+1:32
                [slope, phaseAnglePlot, polyFit] = phaseSlopeIndex(epoch, chan1, chan2);
                phaseSlopeMatrix(chan1,chan2) = slope;
                phaseSlopeMatrix(chan2,chan1) = slope;
            end
        end
        [phaseHeight(epochCounter,:), snr(epochCounter), residual(:,:,epochCounter)]...
            = phaseSlopeTopography(phaseSlopeMatrix, 0.0000001);
        
        epochCounter = epochCounter + 1;
    end
    
%     channelCounter = 1;
%     
%     for chan1 = 1:32
%         for chan2 = chan1+1:32
%             %             fprintf('.');
%             %             if(mod(channelCounter,100) == 0)
%             %                 fprintf('\n%d',size(channelPair,1));
%             %             end
%             channelPair.channel1 = labels{chan1};
%             channelPair.channel2 = labels{chan2};
%             [channelPair.slope, channelPair.phaseAngle, channelPair.poly] = phaseSlopeIndex(data, chan1, chan2);
%             fprintf('\n(%d)%s-%s: %f', channelCounter, channelPair.channel1, channelPair.channel2, channelPair.slope);
%             if(~exist('channelPairs', 'var'))
%                 channelPairs = channelPair;
%             else
%                 channelPairs(end+1) = channelPair;
%             end
%             %debug
%             close all;
%             figure;
%             hold on;
%             plot(channelPair.phaseAngle);
%             polyPredict = (20:100) .* channelPair.poly(1) + channelPair.poly(2);
%             plot(20:100, polyPredict, 'r');
%             title(sprintf('\n(%d)%s-%s: %f', channelCounter, channelPair.channel1, channelPair.channel2, channelPair.slope));
%             %end debug
%             channelCounter = channelCounter + 1;
%         end
%     end
%     
%     
%     summary.channelPairs = channelPairs;
%     clear channelPairs;
%     if(~exist('summaries', 'var'))
%         summaries = summary;
%     else
%         summaries(end+1) = summary;
%     end
%     clear summary;
end

% end
% save(outputFile, 'summaries');