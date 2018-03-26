function clusterAndPlotCoherenceCorrelations(corrMatrix, filename)

samePerm = false;
doTest = true;
showPlot = true;
doInfomap = false;

freqLabels = [{'delta'},{'theta'},{'alpha'},{'beta'},{'hibeta'},{'gamma'}];
freqNumber = 3;

if(showPlot)
    
    %remove mastoids
    pairs = corrMatrix.pairs;
    m1 = cellfun(@length, strfind(pairs, 'M1'));
    m2 = cellfun(@length, strfind(pairs, 'M2'));
    mastoid = find(m1|m2);
    corrMatrix.rho(mastoid,:,:) = [];
    corrMatrix.rho(:,mastoid,:) = [];
    corrMatrix.p(mastoid,:,:) = [];
    corrMatrix.p(:,mastoid,:) = [];
    corrMatrix.pairs(mastoid) = [];
    
    if(~doInfomap)
        toShow = corrMatrix.rho(:,:,freqNumber);
        toShowName = filename;
        %hierarchical cluster
        if(~exist('clustPerm', 'var'))
            [clustMat, clustPerm] = clusterMatrix(toShow);
        else
            if(samePerm)
                [clustMat, clustPerm] = clusterMatrix(toShow, clustPerm);
            else
                [clustMat, clustPerm] = clusterMatrix(toShow);
            end
        end
        
        if(false)
            %high pass filter
            filtMat = NaN(size(clustMat, 1)-1);
            for i = 1:length(filtMat)
                for j = 1:length(filtMat)
                    if(abs(i-j) < 2)
                        filtMat(i,j) = 0;
                    else
                        xfilt = clustMat(i,j) - clustMat(i+1,j);
                        yfilt = clustMat(i,j) - clustMat(i,j+1);
                        %filtMat(i,j) = mean([xfilt,yfilt]);
                        filtMat(i,j) = xfilt;
                    end
                end
            end
            %sum
            figure;
            sumMat = mean(clustMat,1);
            plot(sumMat);
        end
        %         figure;
        %         im = imagesc(clustMat(1:10,1:10), [-.2 1]);
        %         set(im, 'ButtonDownFcn',@handleImgClick);
        %         title(sprintf('%s %s',toShowName, freqLabels{freqNumber}));
        %         colorbar;
        %         drawnow;
    end
    
    Z = linkage(corrMatrix.rho(:,:,freqNumber));
    
    
    if(doInfomap)
        performCoherenceMatrixInfomapClustering();
    end
    %plotCoherenceCluster(clustMat, pairs(clustPerm));
    plotCoherenceCluster(clustMat, corrMatrix.pairs(clustPerm), filename);
end


end

