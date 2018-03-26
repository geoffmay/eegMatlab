filename = '/home/data/EEG/data/ROBI/ROBI_003/baseline eyes open/630158995692243270.eegData';

if(~exist('data','var'))
data = loadRobiDataFile(filename);
end

chan1 = 1;
chan2 = 2;
corrPlot = coherence(data(:,chan1),data(:,chan2));

for framePower = 1:12
  frameSkip = power(2, framePower);
frequency = 4;
sampleTotal = size(corrPlot,1)-frameSkip;

frame1 = NaN(1,sampleTotal);
frame2 = NaN(1, sampleTotal);
for frameNumber = 1:sampleTotal
  frame1(frameNumber) = corrPlot(frameNumber, frequency);
  frame2(frameNumber) = corrPlot(frameNumber + frameSkip, frequency);
end

tic;
[rho, p] = corr(frame1, frame2);
toc;
