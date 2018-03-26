function benchmarkPca()

maxSize = 10;
base = 2;

pcas = struct();
times = NaN(maxSize, 1);
times2 = NaN(maxSize, 1);
times3 = NaN(maxSize, 1);

if(true)
    for i = 1:maxSize
        dataSize = power(base, i);
        data = rand(dataSize*10, dataSize);
        
        toProcess = data;
        tic;
        [dataPca.COEFF, dataPca.SCORE, dataPca.LATENT, dataPca.TSQUARED] = pca(toProcess, 'Algorithm', 'als');
        times(i) = toc;
        fprintf('\n time to process %d by %d is %f seconds', size(toProcess,1), size(toProcess,2), times(i));
    end
end

if(false)
for i = 1:maxSize
    disp(i);
    dataSize = power(base, i);
    data = rand(dataSize);
    data2 = rand(dataSize*base, dataSize);
    data3 = rand(dataSize, dataSize*base);

    toProcess = data;
    tic;
    [dataPca.COEFF, dataPca.SCORE, dataPca.LATENT, dataPca.TSQUARED] = pca(toProcess, 'Algorithm', 'als');
    times(i) = toc;
    fprintf('\n time to process %d by %d is %f seconds', size(toProcess,1), size(toProcess,2), times(i));

    toProcess = data2;
    tic;
    [dataPca.COEFF, dataPca.SCORE, dataPca.LATENT, dataPca.TSQUARED] = pca(toProcess, 'Algorithm', 'als');
    times2(i) = toc;
    fprintf('\n time to process %d by %d is %f seconds', size(toProcess,1), size(toProcess,2), times2(i));

    toProcess = data3;
    tic;
    [dataPca.COEFF, dataPca.SCORE, dataPca.LATENT, dataPca.TSQUARED] = pca(toProcess, 'Algorithm', 'als');
    times3(i) = toc;
    fprintf('\n time to process %d by %d is %f seconds', size(toProcess,1), size(toProcess,2), times3(i));
end
end


end