%to recreate the original data, multiply the mixing matrix times the sig.
%to decompose the data into component timecourses, multiply the separation 
%matrix by the original data.


for i = 1:size(cohIca.icaSig,1)
    maxIca = max(cohIca.icaSig(i,:));
    timecourseIndices(i) = find(cohIca.icaSig(i,:) == maxIca);
end
timeCourseMaxes = tabulate(timecourseIndices);
timeCourseMaxes(timeCourseMaxes(:,2)==0,:) = [];
timeCourseMaxes = sortrows(timeCourseMaxes, 2);
timeCourseMaxes


for i = 1:size(cohIca.icaSig,2)
    maxIca = max(cohIca.icaSig(:,i));
    componentIndices(i) = find(cohIca.icaSig(:,i) == maxIca);
end
componentMaxes = tabulate(componentIndices);
componentMaxes(componentMaxes(:,2)==0,:) = [];
componentMaxes = sortrows(componentMaxes, 2);
componentMaxes


addpath(genpath('/home/data/EEG/scripts/FastICA_25'));

doLoop = false;
if(doLoop)
    
    rates = [ 1/7, 1/31, 1/3, 1/2];
    sampleCount = 1000;
    coeff = reshape(1:(length(rates) * length(rates)), [length(rates), length(rates)]);
    
    for i = 1:length(rates)
        x = 1:sampleCount;
        data(i,:) = sin(x .* rates(i)) .* x;
    end
    
    mash = zeros(length(rates), sampleCount);
    for i = 1:size(coeff, 1)
        for j = 1:size(coeff, 2)
            mash(i,:) = mash(i,:) + data(j,:) .* coeff(i,j);
        end
    end
    
else
    rate1 = 2;
    rate2 = 5;
    rate3 = 3;
    
    data1 = sin((1:sampleCount) .* (rate1*2*pi / sampleCount));
    data2 = sin((1:sampleCount) .* (rate2*2*pi / sampleCount));
    data3 = sin((1:sampleCount) .* (rate3*2*pi / sampleCount));

    data = [data1; data2];
    
    mash = [data1.*coeff(1,1) + data2.*coeff(1,2); data1.*coeff(2,1) + data2.*coeff(2,2)];
end


[sig, mix, sep] = fastica(mash);
remix = mix * sig;
resep = sep * mash;

close all;
figure;
plot(mash');
title('original mix');

figure;
plot(remix');
title('remix');

figure;
plot(data');
title('original data');

figure;
plot(resep');
title('resep');


