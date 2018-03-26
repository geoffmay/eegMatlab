

for i = 1:length(summary.surfResample)
  x(i,1) = summary.surfResample(i).frameCount / 128 / 60;
  del = summary.surfResample(i).averageDifference;
  meanY(i) = mean(del);
  y(i,:) = del;
end

if(true)
  for i = 1:size(y, 2)
    toPlot = y(:,i);
    keep = x > 2;
    poly = polyfit(x(keep), toPlot(keep), 1);
    slopes(i) = poly(1);
    intercepts(i) = poly(2);
    means(i) = mean(toPlot(keep));
  end
end

percentileThreshold = 0.95;
thresholdStdDev = norminv(percentileThreshold);
for i = 1:length(summary.surfResample)
  resample = summary.surfResample(i);
  conf(i,:) = resample.averageDifference + thresholdStdDev .* resample.stddevDifference;  
end

cohInd = cellfun(@length, strfind(summary.surfLabels, 'coh')) > 0;

[sortMeans, sortInd] = sort(means(cohInd));
sortLabels = summary.surfLabels(sortInd);

[sortConf, sortInd] = sort(conf(end, cohInd));
sortLabels = summary.surfLabels(sortInd);

supraThreshold = conf(end,cohInd) > 0.1;

figure;
plot(x, [conf(:,sortInd(1)) y(:, sortInd(1))])
figure;
plot(x, [conf(:,sortInd(end)) y(:, sortInd(end))])

plotCoherencePca(conf(find(x == 5), :), summary.surfLabels);
plotCoherencePca(conf(find(x == 10), :), summary.surfLabels);


figure;
plot(x,y);
xlabel('time (minutes)');