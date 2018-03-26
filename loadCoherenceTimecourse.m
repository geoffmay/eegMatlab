function output = loadCoherenceTimecourse(inputPath)

filterMainsNoise = 1;
filterBadReference = 1;
doRelative = 0;
maxDuration = 128 * 60 * 30;

output = [];

coherenceData = load(inputPath);

proceed = true;
if(isfield(coherenceData, 'filename'))
    if(filterMainsNoise)
        [baseFolder,baseFile,baseExt]=fileparts(coherenceData.filename);
        [hasNoise, noiseRatio] = checkForMainsNoise({[baseFile,baseExt]}, baseFolder);
        output.noiseRatio = noiseRatio;
        if(hasNoise)
            proceed = false;
        end
    end
    if(proceed)
        if(isfield(coherenceData, 'timeCourse'))
            %determine output filename
            [folder, file, ext] = fileparts(coherenceData.timeCourse.basefile);
            if(filterBadReference)
                badRefCutoff = checkForBadReference(coherenceData.channelPairs(29).coherence);
            end
            %         file = strrep(folder, '/home/data/EEG/data/ROBI/', '');
            %         file = strrep(file, '/', ' ');
            freqLabels = {'delta','theta','alpha','beta','hibeta'};
            singleLabels = antChannelLocs;
            singleLabels = singleLabels(1:33);
            
            %remove mastoids
            m1 = cellfun(@length, strfind(singleLabels, 'M1'));
            m2 = cellfun(@length, strfind(singleLabels, 'M2'));
            cpz = cellfun(@length, strfind(singleLabels, 'CPz'));
            removeSingle = find(m1 | m2 | cpz);
            removeDouble = zeros(1,528);
            doubleLabels = cell(size(removeDouble));
            counter = 1;
            for i = 1:33
                badI = length(find(removeSingle==i)) > 0;
                for j = i+1:33
                    badJ = length(find(removeSingle==j)) > 0;
                    if(badI | badJ)
                        removeDouble(counter) = 1;
                    end
                    doubleLabels{counter} = sprintf('%s-%s', singleLabels{i}, singleLabels{j});
                    counter = counter + 1;
                end
            end
            coherenceData.timeCourse.coherencePlot(:,:,find(removeDouble)) = [];
            coherenceData.timeCourse.powerPlot(:,:,removeSingle) = [];
            doubleLabels(find(removeDouble)) = [];
            
            
            %reshape data with proper labels
            cohSize = size(coherenceData.timeCourse.coherencePlot);
            powSize = size(coherenceData.timeCourse.powerPlot);
            allSize = [cohSize(1), cohSize(2) * cohSize(3) + powSize(2) * powSize(3)];
            allDerivations = NaN(allSize);
            allLabels = cell(1, allSize(2));
            freqLabels = {'delta','theta','alpha','beta','hibeta'};
            
            start = 1;
            finish = start + 4;
            for i = 1:cohSize(3)
                allDerivations(:, start:finish) = coherenceData.timeCourse.coherencePlot(:,:,i);
                for(j = 1:length(freqLabels))
                    allLabels{start + j - 1} = sprintf('%s %s', doubleLabels{i}, freqLabels{j});
                end
                start = start + 5;
                finish = start + 4;
            end
            for i = 1:powSize(3)
                allDerivations(:, start:finish) = coherenceData.timeCourse.coherencePlot(:,:,i);
                for j = 1:length(freqLabels)
                    allLabels{start + j - 1} = sprintf('%s %s abs', singleLabels{i}, freqLabels{j});
                end
                start = start + 5;
                finish = start + 4;
            end
            
            %         if(~concatenated)
            %             if(doFastIca)
            %                 %[icasig, A, W]
            %                 [cohIca.icaSig, cohIca.mixingMatrix, cohIca.separatingMatrix] = fastica(allDerivations');
            %             else
            %                 [cohIca.icaSig,cohIca.sphere,cohIca.meanvar,cohIca.bias,cohIca.signs,cohIca.lrates,cohIca.data,cohIca.y] = runica(allDerivations');
            %             end
            %             close all;
            %
            %             %plot first coefficient
            %             if(plotResult)
            %                 plotCoherencePca(cohIca, allLabels, 1);
            %             end
            %
            %
            %             %get channel lables only (no frequency)
            %             onlyLabels = allLabels;
            %             for i = 1:length(onlyLabels)
            %                 label = onlyLabels{i};
            %                 space = strfind(label, ' ');
            %                 if(length(space) > 0)
            %                     label(space(1):end) = [];
            %                 end
            %                 onlyLabels{i} = label;
            %             end
            %             cohIca.labelChannelOnly = onlyLabels;
            %             cohIca.labelChannelFreq = allLabels;
            %             cohIca.filename = coherenceData.timeCourse.basefile;
            %
            %             %display image of pair coefficients
            %             if(false)
            %                 imagesc(cohPca.COEFF(:, 1:10));
            %                 colorbar;
            %             end
            %
            %             %check top values
            %             doPlot = 0;
            %             doTopo = 0;
            %             maxPair = length(cohIca.labelChannelOnly);
            %             maxCoefficient = 20;
            %
            %             if(doTopo)
            %                 [chanLabs, chanlocs] = antChannelLocs;
            %             else
            %                 [chanLabs] = antChannelLocs;
            %             end
            %
            %             M1 = strcmp(chanLabs, 'M1');
            %             M2 = strcmp(chanLabs, 'M2');
            %             CPz = strcmp(chanLabs, 'CPz');
            %             remove = M1 | M2 | CPz;
            %
            %             cohIca.reference = coherenceData.timeCourse.reference;
            %
            %             if(reduceDiskSize)
            %                 %reduce size to save disk space
            %                 cohIca.COEFF(:, maxCoefficient+1:end) = [];
            %                 cohIca.SCORE(:, maxCoefficient+1:end) = [];
            %                 cohIca.EXPLAINED(:, maxCoefficient+1:end) = [];
            %             end
            %
            %             if(calculateVarianceExplained)
            %                 fprintf('\ncalculating explained variance...');
            %                 actual = allDerivations';
            %                 actualMeans = repmat(mean(actual, 2), [1, size(actual,2)]);
            %                 actualError = actual - actualMeans;
            %                 actualVariancePlot = sqrt(sum(actualError .* actualError, 1));
            %                 originalVariance = sum(actualVariancePlot .* actualVariancePlot) / (size(actual,2)-1);
            %                 explainedVariance = NaN(1, size(cohIca.icaSig,1));
            %                 tic;
            %                 for componentIndex = 1:size(cohIca.icaSig, 1)
            %                     fprintf('\ncalculating explained variance %d of %d (%s)', componentIndex, size(cohIca.icaSig,1), char(datetime));
            %                     %icaSig = 2175 x 84722
            %                     %mixing = 2325 x 2175
            %                     %separating = 2175 x 2325
            %                     %allDerivations' = 2325 x 84722
            %                     explained = cohIca.mixingMatrix(:, componentIndex) * cohIca.icaSig(componentIndex, :);
            %                     residual = actualError - explained;
            %                     residualMeans = repmat(mean(residual, 2), [1, size(residual,2)]);
            %                     residualError = residual - actualMeans;
            %                     residualVariancePlot = sqrt(sum(residualError .* residualError, 1));
            %                     newVariance = sum(residualVariancePlot .* residualVariancePlot) / (size(actual,2)-1);
            %                     explainedVariance(componentIndex) = (originalVariance - newVariance) / originalVariance;
            %                     if(false)
            %                         %debug
            %                         close all;
            %                         figure;
            %                         hold on;
            %                         plot(residualVariancePlot, 'b')
            %                         plot(actualVariancePlot, 'r')
            %                         drawnow;
            %                         %end debug
            %                     end
            %                     fprintf(': %.1f%%', explainedVariance(componentIndex) * 100);
            %                 end
            %                 toc;
            %                 cohIca.explainedVariance = explainedVariance;
            %             end
            %
            %             if(saveResult)
            %                 save(outputFile, 'cohIca', 'allLabels', '-v7.3');
            %             end
            %         else
            %             if(~exist('masterData','var'))
            %                 masterData = allDerivations;
            %             else
            %                 masterData(end+1:end+size(allDerivations,1),:) = allDerivations;
            %             end
            %         end
            
            
            %CHANNELPAIRS
            
        elseif(isfield(coherenceData,'channelPairs'))
            if(filterBadReference)
                badReferenceFrame = checkForBadReference(coherenceData.channelPairs(31).coherence(:,1));
                if(badReferenceFrame < size(coherenceData.channelPairs(31).coherence, 1))
                    for i = 1:length(coherenceData.channelPairs)
                        coherenceData.channelPairs(i).coherence(badReferenceFrame:end,:) = [];
                    end
                end
            end
            
            %         %determine output filename
            %         [folder, file, ext] = fileparts(coherenceData.filename);
            %         file = strrep(folder, '/home/data/EEG/data/ROBI/', '');
            %         file = strrep(file, '/', ' ');
            cohLength = size(coherenceData.channelPairs(1).coherence, 1);
            powLength = size(coherenceData.channels(1).absolutePower, 1);
            if(cohLength > 0 & powLength == cohLength)
                
                allSize = size(coherenceData.channelPairs(1).coherence);
                allSize(2) = allSize(2) * length(coherenceData.channelPairs);
                
                allDerivations = NaN(allSize);
                freqLabels = {'delta','theta','alpha','beta','hibeta'};
                
                start = 1;
                finish = start + 4;
                for i = 1:length(coherenceData.channelPairs)
                    allDerivations(:, start:finish) = coherenceData.channelPairs(i).coherence;
                    for(j = 1:length(freqLabels))
                        allLabels{start + j - 1} = sprintf('%s %s', coherenceData.channelPairs(i).label, freqLabels{j});
                    end
                    
                    start = start + 5;
                    finish = start + 4;
                    
                end
                singleLabels = antChannelLocs;
                singleLabels = singleLabels(1:33);
                for i = 1:length(coherenceData.channels)
                    allDerivations(:, start:finish) = coherenceData.channels(i).absolutePower(:,1:length(freqLabels));
                    for j = 1:length(freqLabels)
                        allLabels{start + j - 1} = sprintf('%s %s abs', singleLabels{i}, freqLabels{j});
                    end
                    start = start + 5;
                    finish = start + 4;
                    if(doRelative)
                        allDerivations(:, start:finish) = coherenceData.channels(i).relativePower(:,1:length(freqLabels));
                        for j = 1:length(freqLabels)
                            allLabels{start + j - 1} = sprintf('%s %s rel', singleLabels{i}, freqLabels{j});
                        end
                        start = start + 5;
                        finish = start + 4;
                    end
                end
                
                %remove mastoids
                m1 = cellfun(@length, strfind(allLabels, 'M1'));
                m2 = cellfun(@length, strfind(allLabels, 'M2'));
                cpz = cellfun(@length, strfind(allLabels, 'CPz'));
                remove = find(m1 | m2 | cpz);
                allDerivations(:,remove) = [];
                allLabels(:,remove) = [];
            end
        end
    end
    if(exist('allDerivations', 'var'))
        if(size(allDerivations,1) > maxDuration)
            allDerivations(maxDuration+1,:) = [];
        end
        output.allDerivations = allDerivations;
        output.allLabels = allLabels;
    end
end

end