filename = 'C:\Users\Neuro\Downloads\TMS eeg\brainvision\patient001_Edit Channels_vector.dat';

eeg = loadBrainvisionBinary(filename);
% 
% %parse header
% if(strcmp(lower(filename(end-3:end)),'.dat'))
%     headerName = [filename(1:end-3) 'vhdr'];
%     if(~exist(headerName, 'file'))
%         error(sprintf('header file does not exist: %s', headerName));
%     end
%     headerFileId = fopen(headerName);
%     header = fscanf(headerFileId, '%c');
%     fclose(headerFileId);
%     chanCount = parseLineValue(header, 'NumberOfChannels=');
%     dataPoints = parseLineValue(header, 'DataPoints=');
%     sampleIntervalNanoseconds = parseLineValue(header, 'SamplingInterval=');
%     sampleRate = 1000000 / sampleIntervalNanoseconds;
%     dataOrientation = parseLineValue(header, 'DataOrientation=');
%     dataType = parseLineValue(header, 'DataType=');
%     dataFormat = parseLineValue(header, 'DataFormat=');
%     binaryFormat = parseLineValue(header, 'BinaryFormat=');
%     channelLabels = cell(1,chanCount);
%     for i = 1:chanCount
%         chan = sprintf('Ch%d=', i);
%         chanLine = parseLineValue(header, chan);
%         chanItems = strsplit(chanLine, ',');
%         channelLabels{i} = chanItems{1};
%     end
% end
% chanlocs = getChanlocs(channelLabels);
% 
% %read data
% file = dir(filename);
% fileLength = file.bytes / 4;
% fileId = fopen(filename);
% contents = fread(fileId, fileLength, 'single');
% fclose(fileId);
% sampleCount = length(contents) / chanCount;
% if(strcmp(dataOrientation, 'VECTORIZED'))
%     allData = reshape(contents, [sampleCount, chanCount]);
% else
%     error('unsupported data orientation, need to write code for multiplexed here.');
% end
% x = (1:size(allData,1)) ./ sampleRate;

% minSecond = 800;
% maxSecond = 1000;
% minSecond = 350;
% maxSecond = 370;
% keepX = (minSecond*sampleRate):(maxSecond*sampleRate);
% keepData = data(keepX, :);
% close all;
% chanIndexes = [1 5];
% plot(x(keepX), data(keepX, chanIndexes));
% legend(channelLabels(chanIndexes));




%tms function
tmsAnalysis = comparePreVsPostTms(eeg);
% 
% tms = findTmsBookends(allData, 2500, sampleRate);
% epochDuration = 1;
% bufferDuration = 0.0;
% epochLength = epochDuration * sampleRate;
% bufferLength = bufferDuration * sampleRate;
% 
% %epoch x time x channel
% preTms = NaN(length(tms.tmsStart), epochLength, chanCount);
% postTms = NaN(length(tms.tmsStart), epochLength, chanCount);
% preRange = NaN(length(tms.tmsStart), chanCount);
% postRange = NaN(length(tms.tmsStart), chanCount);
% rangeChange = NaN(length(tms.tmsStart), chanCount);
% for i = 1:length(tms.tmsStart)
%     startI = (tms.tmsStart(i) - epochLength - bufferLength + 1):(tms.tmsStart(i) - bufferLength);
%     endI = (tms.tmsFinish(i)+bufferLength):(tms.tmsFinish(i) + epochLength + bufferLength - 1);
%     preTms(i, :, :) = allData(startI, :);
%     postTms(i, :, :) = allData(endI, :);
%     for j = 1:chanCount
%         preRange(i,j) = max(preTms(i,:,j)) - min(preTms(i,:,j));
%         postRange(i,j) = max(postTms(i,:,j)) - min(postTms(i,:,j));
%         rangeChange(i,j) = postRange(i,j) - preRange(i,j);
%     end
% end
% 
% %wavelet
% for channelIndex = 1:chanCount
%     preWav = NaN(size(preTms,1), 94, epochLength);
%     postWav = NaN(size(preTms,1), 94, epochLength);
%     for epochIndex = 1:size(preTms,1)
%         fprintf('\nchannel %d, epoch %d', channelIndex, epochIndex);
%         preData = squeeze(preTms(epochIndex, :, channelIndex));
%         postData = squeeze(postTms(epochIndex, :, channelIndex));
%         [waveletTransform, freqInfo, coneOfInfluence] = cwt(preData, sampleRate);
%         preWav(epochIndex, :, :) = cwt(preData);
%         postWav(epochIndex, :, :) = cwt(postData);
%     end
%     %wavelet t test
%     for freqCounter = 1:size(preWav, 2)
%         fprintf('\nchannel %d freq %d of %d', channelIndex, freqCounter, size(preWav, 2));
%         for timeCounter = 1:size(preWav, 3)
%             x = preWav(:, freqCounter, timeCounter);
%             y = postWav(:, freqCounter, timeCounter);
%             [tTest.H,tTest.P,tTest.CI,tTest.STATS] = ttest(x,y);
%             t.H(channelIndex, freqCounter, timeCounter) = tTest.H;
%             t.P(channelIndex, freqCounter, timeCounter) = tTest.P;
%             t.CiMin(channelIndex, freqCounter, timeCounter) = tTest.CI(1);
%             t.CiMax(channelIndex, freqCounter, timeCounter) = tTest.CI(2);
%             t.STATS(channelIndex, freqCounter, timeCounter) = tTest.STATS;
%         end
%     end
% end
% 
% %image y axis tick labels
% freqIndex = length(freqInfo):-5:1;
% plotFreq = freqInfo(freqIndex);
% for i = 1:length(plotFreq)
%     freqLabel{i} = sprintf('%.2f', plotFreq(i));
% end
% plotFreqLabel = fliplr(freqLabel);
% plotFreqIndex= fliplr(freqIndex);
% 
% meanPostWav = squeeze(mean(postWav, 1));
% diffMeanWav = squeeze(mean(postWav, 1)) - squeeze(mean(preWav, 1));
% meanAbsPostWav = squeeze(mean(abs(postWav), 1));









% 
% figure;
% imagesc(abs(diffMeanWav));
% set(gca, 'ytick', plotFreqIndex);
% set(gca, 'yticklabel', plotFreqLabel);
% title('mean(pre - post)');
% colorbar;
% 
% plotP = tPs;
% plotP(plotP > .2) = .2;
% plotP = 1 - plotP;
% figure;
% imagesc(plotP);
% set(gca, 'ytick', plotFreqIndex);
% set(gca, 'yticklabel', plotFreqLabel);
% title('t test p values');
% colorbar;
% 
% figure;
% imagesc(squeeze(postWav(1,:,:)));
% set(gca, 'ytick', plotFreqIndex);
% set(gca, 'yticklabel', plotFreqLabel);
% title('single post trial');
% colorbar;
% 
% 
% figure;
% imagesc(abs(meanPostWav));
% set(gca, 'ytick', plotFreqIndex);
% set(gca, 'yticklabel', plotFreqLabel);
% title('mean post');
% colorbar;
% 
% figure;
% imagesc(abs(meanAbsPostWav));
% set(gca, 'ytick', plotFreqIndex);
% set(gca, 'yticklabel', plotFreqLabel);
% title('mean abs post');
% colorbar;
% 
% 
% figure;
% plotMin = abs(tCiMins);
% imagesc(plotMin);
% set(gca, 'ytick', plotFreqIndex);
% set(gca, 'yticklabel', plotFreqLabel);
% title('95% CI min');
% colorbar;
% 
% figure;
% plotMax = abs(tCiMaxs);
% imagesc(plotMax);
% set(gca, 'ytick', plotFreqIndex);
% set(gca, 'yticklabel', plotFreqLabel);
% title('95% CI max');
% colorbar;
% % 
% % imagesc(abs(tHs));
% % title('pre-tms');
% % xlabel('time (frames)');
% % ylabel('frequency (Hz)');
% % figure;
% % imagesc(abs(postWav));
% % title('post-tms');
% % xlabel('time (frames)');
% % ylabel('frequency (Hz)');
% % 
% 
% if(false)
%     preData = preTms(:, :, channelIndex)';
%     postData = postTms(:, :, channelIndex)';
%     % [ersp,itc,powbase,times,freqs,erspboot,itcboot,itcphase] = ...
%     %                  timef(data,frames,tlimits,srate,cycles,...
%     %                                   'key1',value1,'key2',value2, ... );
%     toProcess = reshape(preData, [1, numel(preData)]);
%     figure;
%     [preWav.ersp,preWav.itc,preWav.powbase,preWav.times,preWav.freqs,preWav.erspboot,preWav.itcboot,preWav.itcphase] = timef(toProcess, size(preData,1), [0 epochDuration*1000], sampleRate, 3);
%     toProcess = reshape(postData, [1, numel(postData)]);
%     figure;
%     [postWav.ersp,postWav.itc,postWav.powbase,postWav.times,postWav.freqs,postWav.erspboot,postWav.itcboot,postWav.itcphase] = timef(toProcess, size(preData,1), [0 epochDuration*1000], sampleRate, 3);
%     %t-test ranges
%     for i = 1:size(rangeChange,2)
%         x = preTms(:,i);
%         y = postTms(:,j);
%         [tTest.H,tTest.P,tTest.CI,tTest.STATS] = ttest(x,y);
%         ttests(i) = tTest;
%     end
% end
% 
% %
% % close all;
% % imagesc(rangeChange);
% % rangeChange(end+1,:) = max(max(rangeChange));
% % figure;
% % imagesc(rangeChange);
% % xlabel('channel');
% % ylabel('epoch');
% %
% % i = 0;
% %
% % i = i+1;
% % clf;
% % plot(data(tms.tmsStart(i):tms.tmsFinish(i),:));
% % plot(data(tms.tmsFinish(i):tms.tmsStart(i+1),:));
% 
% 
% 
% 
% 
% %do ica
% doIca = false;
% if(doIca)
%     if(~exist('runica', 'file'))
%         eeglab
%         close all;
%     end
%     icaSample = allData();
%     icaSample = icaSample;
%     [ica.weights,ica.sphere,ica.compvars,ica.bias,ica.signs,ica.lrates,ica.activations] = runica(keepData');
%     close all;
%     plot(x(keepX), ica.activations([1 5],:)');
% end
% 
% 
% %compute ffts
% doFft = false;
% if(doFft)
%     for i = 1:size(ica.activations, 1)
%         fprintf('\n%d', i);
%         ffts(i,:) = fft(ica.activations(i,:));
%     end
%     ffts(:, (size(ffts,2)/2):end) = [];
%     ffts(:, 1) = [];
%     fftXMin = 1/(maxSecond-minSecond);
%     fftXMax = sampleRate / 2;
%     fftX = fftXMin:fftXMin:fftXMax;
%     fftX(end) = [];
%     toPlot = find(fftX <= 130);
%     
%     %plot ffts
%     plotChan = plotChan + 1
%     close all;
%     plot(fftX(toPlot), abs(ffts(plotChan,(toPlot))));
%     pan xon;
%     zoom xon;
%     
%     
%     chan = 4;
%     weights = ica.weights(chan,:);
%     topoplot(weights, chanlocs);
%     
%     
%     close all;
%     plotChannels = [1, 2, 4];
%     plot(allData(:, plotChannels));
%     legend(channelLabels(plotChannels));
%     pan xon;
%     zoom xon;
% end
% 


save('C:\Users\Neuro\Documents\MATLAB\ttestWorkspaceWithBuffer.mat', '-v7.3');
