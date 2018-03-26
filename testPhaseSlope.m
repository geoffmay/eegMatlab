sampleRate = 2048;
sampleLength = 100*sampleRate;
phaseShift = 2048 / 100;

coefficients = 1./(1:2048);

sim1 = zeros(1,sampleLength);
sim2 = zeros(1,sampleLength);

rng(1);
jitter = rand(1,sampleRate);
warp = rand(2,sampleLength) .* .001;

%for i = 1:length(coefficients)
for i = 1:100
    addend = sin((1:sampleLength).*(2*pi*i/sampleRate) + jitter(i)*2*pi);
    shiftAddend = sin(((1+phaseShift):(sampleLength+phaseShift)).*(2*pi*i/sampleRate) + jitter(i)*2*pi);
    addend = addend + warp(1,:);
    shiftAddend = shiftAddend + warp(2,:);
    sim1 = sim1 + addend.*coefficients(i);
    sim2 = sim2 + shiftAddend.*coefficients(i);
    if(false)
        %debug
        close all;
        figure;
        hold on;
        plot(addend);
        plot(shiftAddend, 'r');
        %end debug
    end

end
%debug
close all;
figure;
hold on;
plot(sim1);
plot(sim2, 'r');
fullSim = [sim1;sim2]';

[slope, meanPhaseAngle, polyLine] = phaseSlopeIndex(fullSim,1,2);
close all;
figure;
plot(meanPhaseAngle);
% plot(shiftAddend, 'r');
%end debug
