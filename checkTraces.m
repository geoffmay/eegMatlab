folder = 'C:\Users\Neuro\Documents\temp\eegFourier\EEG_data_sub-11_run-01.mat';
files=dir(folder);
a = {files.name};
isCoh = cellfun(@length, strfind(a, 'coherence')) > 0;
isPow = cellfun(@length, strfind(a, 'absolutePower')) > 0;
isAlpha = cellfun(@length, strfind(a, '_3.mat')) > 0;
interesting = a(isPow & isAlpha);

for i = 1:length(interesting)
    parts = strsplit(interesting{i}, '_');
    num = parts{2};
    num = num(3:end);
    ind = str2double(num);
    dat = load(fullfile(folder,interesting{i}));
    allPower(ind, :) = dat.piece.signal;
%     plot(dat.piece.times, dat.piece.signal);
%     myTitle = strrep(dat.piece.label, '_', '-');
%     title(myTitle);
end
rhos = ones(size(allPower, 1));
for i = 1:length(rhos)
    for j = (i+1):length(rhos)
        fprintf('%d, %d\n', i, j);
        [rho, p] = corr(allPower(i,:)', allPower(j,:)');
        rhos(i,j) = rho;
        rhos(j,i) = rho;
    end
end
figure;
imagesc(rhos);
[perRhos, perm] = clusterMatrix(rhos);
