function [boldResponse, boldTimes] = s8convolveHrf(brainActivity, sampleRate)

doPlot = false;
if(~exist('sampleRate', 'var'))
    sampleRate = 4;
end

if(size(brainActivity, 1) == 1)
    %we're good
elseif(size(brainActivity, 2) == 1)
    brainActivity = brainActivity';
else
    error('brainActivity must be a row vector');
end

tr = 1/sampleRate;
spmHrf = spm_hrf(tr)';

if(doPlot)
    close all;
    plot(spmHrfT, spmHrf);
    xlabel('time (seconds)');
    ylabel('BOLD response (arbitrary units)');
end

boldResponse = s9convolution(spmHrf, 1:length(spmHrf), brainActivity, 1:length(brainActivity));
boldTimes = tr:tr:(length(boldResponse)*tr);

if(doPlot)
    close all;
    plot(ht, hx);
    xlabel('time (seconds)');
    ylabel('BOLD response (arbitrary units)');
end
