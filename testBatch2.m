function [ output_args ] = testBatch2( input_args )
%TESTBATCH2 Summary of this function goes here
%   Detailed explanation goes here

%inputFolder = '/Users/Geoff/Box Sync/For Geoff Mays/VET MIND EEG files/';
inputFolder = '/Users/Geoff/Box Sync/For Geoff Mays/PTSD MIND EEG Files';
outputFolder = '/Users/Geoff/Box Sync/For Geoff Mays/Sets/';

% [tones, ant, other] = getBdfFiles();
% allFiles = [ant];
allFiles = cell(0);
files = dir(inputFolder);
for i = 1:length(files)
    file = files(i).name;
    if(strfind(file, '.bdf'))
        allFiles{length(allFiles)+1} = file;
    end
end

for i = 1:length(allFiles)
    bdfFileName = allFiles{i};
    refChan = int32(32);
    testEEG = pop_readbdf(strcat(inputFolder,bdfFileName), {}, 43, refChan, true);
    for(j=1:length(testEEG.chanlocs))
        labels{j} = testEEG.chanlocs(j).labels;

    end

    
    % parameters for PSI-calculation
    segleng=100;epleng=200;
    
    % calculation of PSI. The last argument is empty - meaning that
    % PSI is calculated over all frequencies
    [psi, stdpsi, psisum, stdpsisum]=data2psi(testEEG.data.',segleng,epleng,[]);
    % note, psi, as calculated by data2psi corresponds to \hat{\PSI}
    % in the paper, i.e., it is not normalized. The final is
    % the normalized version given by:
    psi./(stdpsi+eps);
    
    outputFileName = strcat(bdfFileName, 'PhaseShiftIndexAllFrequencies.csv');
    path = strcat(outputFolder,outputFileName);
%     table = array2table(psi, 'VariableNames', labels);
%     writetable(table, path);    
    csvwrite(path, psi);
end

end

