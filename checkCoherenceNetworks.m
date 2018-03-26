
plotAll = true;

rootFolder = '/media/eegDrive';
allFiles = dir(rootFolder);
allFiles([allFiles.isdir]) = [];

inputFilenames= [{'/media/eegDrive/ROBI_003_outcome eyes open_630230539149591228coherenceStats.mat'},...
    {'/media/eegDrive/ROBI_003_baseline eyes open_630158995692243270coherenceStats.mat'}];

% inputFilename = '/media/eegDrive/ROBI_003_tx 3_630166729229512303coherenceStats.mat';
% load(inputFilename);

for i = 1:length(allFiles)
    if(~any(strcmp(inputFilenames, allFiles(i).name)))
        inputFilenames{end+1} = fullfile(rootFolder, allFiles(i).name);
    end
end


for fileCounter = 1:length(inputFilenames)
    inputFilename = inputFilenames{fileCounter};
    [folder, filename, ext] = fileparts(inputFilename);
    outputFolder = '/home/data/EEG/processed/Robi/coherenceNetwork';
    outputFilename = fullfile(outputFolder, sprintf('%s.mat',filename));
    if(~exist(outputFilename, 'file'))        
        placeholder = char(datetime);
        save(outputFilename, 'placeholder');       
        [tab, rho, p, pairs] = computeCoherenceCorrelations(inputFilename)        
        save(outputFilename, 'tab', 'rho', 'p', 'pairs');
    end
    if(plotAll)
        correlationMatrix = load(outputFilename);
        if(isfield(correlationMatrix, 'rho'))
            clusterAndPlotCoherenceCorrelations(correlationMatrix, filename);
        end
    end
end
tilefigs;
