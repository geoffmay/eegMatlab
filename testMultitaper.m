function testMultitaper()
folder = '/home/data/EEG/data/Oregon/Flanker1';
files = dir(folder);
files([files.isdir]) = [];
path = fullfile(folder, files(1).name);
EEG = pop_readbdf(path, {}, 43, int32(32), false);


waveletCycles = 0;
eventValue = 114;
channel = 31;
preMills = 500;
postMills = 2500;
if(preMills < 0) preMills = -preMills;end

indices = find([EEG.event.type]== eventValue);
latencies = [EEG.event.latency];
eventCount = length(EEG.event);
eventLatencies = latencies(indices);
preFrames = (preMills * 0.001 * EEG.srate);
postFrames = (postMills * 0.001 * EEG.srate);
frameDuration = preFrames + postFrames + 1;
data = zeros(1, frameDuration * eventCount);
frameNumbers = int64(zeros(1, length(data)));
writeIndex = int64(1);
for i = 1:length(eventLatencies)
  readIndex = eventLatencies(i)-preFrames;
  if(readIndex > 0 && readIndex+frameDuration <= size(EEG.data, 2))
    nextData = EEG.data(channel, readIndex:(readIndex+frameDuration));
    data(writeIndex:(writeIndex+frameDuration)) = nextData;
    frameNumbers(writeIndex:writeIndex+frameDuration) = readIndex:readIndex+frameDuration;
  else
    data(end-frameDuration+1:end) = [];
    frameNumbers(end-frameDuration+1:end) = [];
  end
  writeIndex = writeIndex + frameDuration;
end

srate = EEG.srate;
tlimits = [-preMills, postMills];


[ersp,itc,allErsp] = ...
  timef(data,frameDuration,tlimits,srate,waveletCycles);


end