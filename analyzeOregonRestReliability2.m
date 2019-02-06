folder = '/Users/Geoff/Documents/reliability5/surface digest';
files = dir(folder);
files([files.isdir]) = [];

mat = NaN(70, length(files));

raw = load('~/Documents/reliability5/PM101.mat')
isCoh = cellfun(@length, strfind(raw.summary.surfLabels, 'coh')) > 0;

cohFitOptions = fitoptions('exp2');
cohFitOptions.StartPoint = [0.0931   -0.0266    0.0948   -0.0009];
cohFitOptions.Lower = [0   -1    0   -1];
cohFitOptions.MaxIter = 2000;
%cohFitOptions.Upper = [1   1    1   1];
cohFitOptions.Upper = [3   0    3   0];


powFitOptions = cohFitOptions;
%cohFitOptions.Upper = [5   1    5   1];
powFitOptions.Upper = [10   0    10   0];

for i = 1:length(files)
    a = load(fullfile(folder, files(i).name));
    ind = find(strcmp({a.summaries.label}, 'coh Fp1-Fp2 5Hz-8Hz'));
    ind = 1;
    row = mat(:, i);
    b = a.plot95Raw{ind};
    x = (1:length(b)) .* 15;
    row(1:length(b)) = b;
    mat(:,i) = row;
    
    clear fits
    for ind = 1:length(a.plot95Raw)        
        fprintf('\nfile %d of %d; measure %d of %d', i, length(files), ind, length(a.plot95Raw));
        b = a.plot95Raw{ind};
        x = (1:length(b)) .* 15;
        row(1:length(b)) = b;
        mat(:,i) = row;
        if(isCoh(ind))
            [f.fun, f.good] = fit(x', b', 'exp2', cohFitOptions);
        else
            [f.fun, f.good] = fit(x', b', 'exp2', powFitOptions);
        end
        if(ind == 1)
            fits = repmat(f, [1, length(a.plot95Raw)]);
            goods = NaN(1, length(a.plot95Raw));
        end
        fits(ind) = f;
        goods(ind) = f.good.adjrsquare;
    end
    outfile = fullfile('/Users/Geoff/Documents/reliability5/surface digest 2', files(i).name);
    save(outfile, 'fits', 'goods');
end

intermediateFolder = '/Users/Geoff/Documents/reliability5/surface digest 2/old';
for i = 1:length(files)
    data = load(fullfile(intermediateFolder, files(i).name));
    for j = 1:length(data.fits)
        fprintf('\nfile %d of %d, measure %d of %d', i, length(files), j, length(data.fits));
        coeffs = coeffvalues(data.fits(j).fun);
        fasterFirst = abs(coeffs(2)) > abs(coeffs(4));
        if(~fasterFirst)
            temp = coeffs(1:2);
            coeffs(1:2) = coeffs(3:4);
            coeffs(3:4) = temp;
        end
        coeffMat(i, j, :) = coeffs;
    end
end
save('exponentialFit.mat', 'coeffMat');
% equation is a*e^(b*x) + c*e^(d*x)
%   'coh icaComp1-icaComp19 9Hz-12Hz' is frickin enormous for a and c
%   (9e12 and -9e12)
maxVal = max(max(max(coeffMat)));
minVal = min(min(min(coeffMat)));
maxInd = find(coeffMat == maxVal);
minInd = find(coeffMat == minVal);
minVal = min(min(min(coeffMat)));

extreme = abs(coeffMat) > .5;
coeffMat(extreme)  = 0;
cohMatA = squeeze(coeffMat(:,isCoh,1));
powMatA = squeeze(coeffMat(:,~isCoh,1));
cohMatB = squeeze(coeffMat(:,isCoh,2));
powMatB = squeeze(coeffMat(:,~isCoh,2));
cohMatC = squeeze(coeffMat(:,isCoh,3));
powMatC = squeeze(coeffMat(:,~isCoh,3));
cohMatD = squeeze(coeffMat(:,isCoh,4));
powMatD = squeeze(coeffMat(:,~isCoh,4));

imagesc(cohMatB);
colorbar;

meanD = mean(cohMatD, 2);
[~, indD] = sort(meanD);
sortCohMatD = cohMatD(indD,:);
imagesc(sortCohMatD);

if(false)
    %find extremes
    %here is a very low d (low non-stationarity?)
    minD = cohMatD == min(min(cohMatD));
    minD1 = find(any(minD,2));
    minD2 = find(any(minD,1));
    i = minD1;
    ind = minD2;
    a = load(fullfile(folder, files(i).name));
    b = a.plot95Raw{ind};
    
    %here is a very low b (low measurement noise?)
    minB = cohMatB == min(min(cohMatB));
    
    minB1 = find(any(minB,2));
    minB2 = find(any(minB,1));
    i = minB1;
    ind = minB2;
    a = load(fullfile(folder, files(i).name));
    b = a.plot95Raw{ind};


end

%check for neuropsych correlations
neuropsychDataFilename = '~/Documents/MATLAB/EEG/PTSD MIND for Geoff.mat';
neuropsychData = load(neuropsychDataFilename);
neuropsychData = neuropsychData.neuropsychData;

for i = 1:length(files)
    [folder, file, ext] = fileparts(files(i).name);
    file = str2num(file(4:5));
    rows(i) = find(neuropsychData{:,1} == file);
end
for i = 1:size(neuropsychData, 2)
    x = meanD;
    y = neuropsychData{rows, i};
    [rhos(i), ps(i)] = corr(x,y);
end
pInd = ps < 0.05


hist(cohMat);
legend({'a', 'b', 'c', 'd'});


