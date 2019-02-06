function convertSvgToEmf(folder)

exe = 'C:\Program Files\Inkscape\inkscape.exe';
if(~exist('folder', 'var'))
    folder = 'C:\Users\Neuro\Documents\reliabilityFigures\';
end
files = dir(folder);
files([files.isdir]) = [];

for i = 1:length(files)
    fprintf('converting image %d of %d (%s)\n', i, length(files), char(datetime));
    inputPath = fullfile(folder, files(i).name);
    outputFile = strrep(inputPath, '.svg', '.emf');
    command = sprintf('"%s" --file "%s" --export-emf "%s"', exe, inputPath, outputFile);
    
    [a,b] = system(command);
end
