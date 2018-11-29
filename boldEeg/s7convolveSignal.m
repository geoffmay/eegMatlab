
function convSig = s7convolveSignal(sig, sampleRate, downsampleRate, markerTimes, x)

dSignal = downsample(sig, downsampleRate);
dTimes = downsample(x, downsampleRate);

dRate = sampleRate / downsampleRate;


dHrf = s8convolveHrf(dSignal, dRate);
convSig = NaN(1, length(markerTimes));
for i = 1:length(markerTimes)
  %   targetFrame = dMarkers(i);
  %   targetTime = targetFrame / dRate;
  pieceFrame = min(find(x >= markerTimes(i)));
  dFrame = pieceFrame / downsampleRate;
  convSig(i) = interp1(1:length(dHrf), dHrf, dFrame);
end

