clear;

if(~exist('versions', 'var'))
    inputFolder = '/home/data/EEG/processed/Robi/coherencePca';
    files = dir(inputFolder);
    files([files.isdir]) = [];
    
    figure;
    hold on;
    
    for i = 1:length(files)
        versions{i} = load(fullfile(inputFolder, files(i).name));
        primary{i} = versions{i}.cohPca.SCORE(:,1);
        if(i == 1)
            remove = isnan(primary{i});
        else
            remove = remove | isnan(primary{i});
        end
    end
end

componentCount = size(versions{1}.cohPca.SCORE,2);
componentCount = 100;
referenceCount = length(versions);
rho = NaN(componentCount,referenceCount,componentCount,referenceCount);

for component1 = 1:componentCount
    fprintf('\n%s %d of %d', char(datetime), component1, componentCount);
    for reference1 = 1:referenceCount
        score1 = versions{reference1}.cohPca.SCORE(~remove,component1);
        maxRho = realmin;
        maxRhoIndex = -1;
        for component2 = 1:componentCount
            for reference2 = 1:referenceCount
                if(reference1 ~= reference2)
                    score2 = versions{reference2}.cohPca.SCORE(~remove,component2);
                    [thisRho, thisP] = corr(score1, score2);
                    rho(component1, reference1, component2, reference2) = thisRho;
                    p(component1, reference1, component2, reference2) = thisP;
                    if(abs(thisRho) > maxRho)
                        maxRho = abs(thisRho);
                        maxRhoIndex = [component2, reference2];
                    end
                end
            end
        end
        matchRho(component1, reference1) = maxRho;
        matchRhoIndex{component1, reference1} = maxRhoIndex;
    end
end
close all;
figure;
imagesc(matchRho');
set(gca, 'ytick', [1 2 3]);
set(gca, 'xtick', []);
set(gca, 'yticklabel', {'Avg','Cpz','Mast'});
cbar = colorbar;
title(cbar, 'R^2');
title('PCA component timecourse correlations across EEG references');
ylabel('reference')
xlabel('maximum correlation for each component');


for i = 1:length(primary)
    a = primary{i};
    a(remove) = [];
    primary{i} = a;
end

[rho1, p1] = corr(primary{1}, primary{2});
[rho2, p2] = corr(primary{1}, primary{3});
[rho3, p3] = corr(primary{2}, primary{3});

if(false)
    componentNumber = 1;
    for i = 1:length(versions)
        myTitle = strrep(files(i).name, '_', ' ');
        plotCoherencePca(versions{i}.cohPca, versions{i}.cohPca.labelChannelFreq, componentNumber, myTitle);
    end
end