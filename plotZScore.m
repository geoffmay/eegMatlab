if(false)
    plotIndividual = 1;
    filter = 'open';
    
    inputFolder = '/media/eegDrive/zScoreTimecourse';
    files = dir(inputFolder);
    files([files.isdir]) = [];
    files(find(cellfun(@length, strfind({files.name}, filter)) == 0)) = [];
    
    inputFile = fullfile(inputFolder, 'ROBI_003_baseline eyes open.mat');
    exitFile = fullfile(inputFolder, 'ROBI_003_outcome eyes open.mat');
    
    allValues = zeros(1,length(files));
    for fileCounter = 1:length(files)
        fprintf('\n%d of %d', fileCounter, length(files));
        outputFile = sprintf('/media/eegDrive/figs/%s', files(fileCounter).name);
        outputFile = strrep(outputFile, '.mat', '.png');
        if(~exist(outputFile,'file'))
            data = load(fullfile(inputFolder, files(fileCounter).name));
            
            freqLabels = {'delta', 'theta', 'alpha', 'beta', 'hibeta'};
            
            counter = 1;
            for i = 1:171
                for j = 1:5
                    values(counter) = mean(data.zScoreTimecourse.coherence(:,j,i));
                    labels{counter} = sprintf('%s %s', data.zScoreTimecourse.channelPairLabels{i}, freqLabels{j});
                    counter = counter + 1;
                end
            end
            for i = 1:19
                for j = 1:5
                    values(counter) = mean(data.zScoreTimecourse.absPower(:,j,i));
                    labels{counter} = sprintf('%s %s abs', data.zScoreTimecourse.channelLabels{i}, freqLabels{j});
                    counter = counter + 1;
                end
            end
            allValues(fileCounter) = mean(abs(values));
            
            if(plotIndividual)
                close all;
                try
                    fprintf('\n%s\n',files(fileCounter).name);
                    fig = plotCoherencePca(values, labels);
                    print(outputFile, '-dpng', '-r0');
                catch ex
                    fprintf('\nerror: %s (%s)', ex.message, files(fileCounter).name);
                end
            end
        end
    end
    save('/home/data/EEG/processed/Robi/zScoreProgressionRobi3.mat','files','allValues');
    
    toPlotNames = {...
        'ROBI_003_baseline eyes open.mat',...
        'ROBI_003_tx 1.mat',...
        'ROBI_003_tx 2.mat',...
        'ROBI_003_tx 3.mat',...
        'ROBI_003_tx 4.mat',...
        'ROBI_003_tx 5.5.mat',...
        'ROBI_003_tx 5.mat',...
        'ROBI_003_tx 6.mat',...
        'ROBI_003_tx 7.mat',...
        'ROBI_003_tx 8.mat',...
        'ROBI_003_tx 9.mat',...
        'ROBI_003_tx 10.mat',...
        'ROBI_003_tx 11.mat',...
        'ROBI_003_tx 12.mat',...
        'ROBI_003_tx 13.mat',...
        'ROBI_003_tx 14.mat',...
        'ROBI_003_tx 15.mat',...
        'ROBI_003_tx 17.mat',...
        'ROBI_003_tx 18.mat',...
        'ROBI_003_tx 19.mat',...
        'ROBI_003_tx 20.mat',...
        'ROBI_003_outcome eyes open.mat',...
        };
    
    threshold = 2;
    clear plotValues;
    counter = 1;
    for i = 1:length(toPlotNames)
        index = find(strcmp({files.name}, toPlotNames{i}));
        val = allValues(index);
        if(~isnan(val) & val < threshold)
            plotValues(counter) = val;
            counter = counter + 1;
        end
    end
    figure;
    plot(plotValues);
    
end

z = false;
if(z)
    folder = '/media/eegDrive/fast/fastZ';
else
    folder = '/media/eegDrive/fast/';
end

subjectId = 3;
clear del;
for subjectId = 1:5
    
    preFile = fullfile(folder, sprintf('ROBI_%03d_baseline eyes open.mat', subjectId));
    postFile = fullfile(folder, sprintf('ROBI_%03d_outcome eyes open.mat', subjectId));
    if(exist(preFile,'file') & exist(postFile, 'file'))
        pre = load(preFile);
        post = load(postFile);
        if(z)
            pre = pre.zScores;
            post = post.zScores;
        else
            pre = pre.timeCourse;
            post = post.timeCourse;
        end
        if(z)
            del.relPower = mean(post.relPower) - mean(pre.relPower);
            del.absPower = mean(post.absPower) - mean(pre.absPower);
            del.coherence = mean(post.coherence) - mean(pre.coherence);
            del.channelPairLabels = post.channelPairLabels;
        else
            preIndex = 1:(length(pre.times) - abs(pre.badReferenceFramesDropped));
            postIndex = 1:(length(post.times) - abs(post.badReferenceFramesDropped));
            for i = 1:length(pre.coherencePlot)
                pre.coherencePlot(i).coherence = mean(pre.coherencePlot(i).coherence(preIndex,:));
                post.coherencePlot(i).coherence = mean(post.coherencePlot(i).coherence(postIndex,:));
                coh.coherence = post.coherencePlot(i).coherence - pre.coherencePlot(i).coherence;
                %                 coh.coherence = mean(post.coherencePlot(i).coherence(postIndex,:)) - mean(pre.coherencePlot(i).coherence(preIndex,:));
                coh.label = post.coherencePlot(i).label;
                del.coherencePlot(i) = coh;
            end
            for i = 1:length(pre.powerPlot)
                pre.powerPlot(i).absolutePower = mean(pre.powerPlot(i).absolutePower(preIndex,:));
                post.powerPlot(i).absolutePower = mean(post.powerPlot(i).absolutePower(postIndex,:));
                coh.absolutePower = post.powerPlot(i).absolutePower - pre.powerPlot(i).absolutePower;
                %                 coh.absolutePower = mean(post.powerPlot(i).absolutePower(postIndex,:)) - mean(pre.powerPlot(i).absolutePower(preIndex,:));
                coh.label = post.powerPlot(i).label;
                del.powerPlot(i) = coh;
                %             del.powerPlot(i) = post.powerPlot(i) - pre.powerPlot(i);
            end
        end
        close all;
        plotCoherencePca(pre);
        plotCoherencePca(post);
        plotCoherencePca(del);
        fprintf('\n%s', preFile);
        dummy = 1;
    end
end
