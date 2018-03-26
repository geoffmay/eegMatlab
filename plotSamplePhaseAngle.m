sampleRate = 2048;
startIndex = 1;
endIndex = startIndex + sampleRate - 1;
subStart = 2;
subEnd = sampleRate / 2 + 1;

simHz = 2;
phaseShift1 = 0*pi/2;
phaseShift2 = 1*pi/2;
sim1Data = sin((simHz*2*pi / sampleRate) .* (startIndex:endIndex) + phaseShift1);
sim2Data = sin((simHz*(2*pi) / sampleRate) .* (startIndex:endIndex) + phaseShift2);

x = (1:1024) ./ 1024;

close all;
figure;
hold on;

plot(x, sim1Data(subStart:subEnd),'b', 'Linewidth', 3);
plot(x, sim2Data(subStart:subEnd),'r', 'Linewidth', 3);

fft1 = fft(sim1Data);
fft2 = fft(sim2Data);
pa = phaseAngle(fft1, fft2, 2, 5);
paDegrees = pa / pi * 180

title(sprintf('%d degree phase angle between two sine waves', paDegrees));
legend('wave 1', 'wave 2');
xlabel('period (arbitrary units)');
ylabel('amplitude (arbitrary units)');

