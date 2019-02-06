folder = 'C:\Users\Neuro\Documents\reliabilityFigures';
exe = 'C:\Program Files\Inkscape\inkscape.exe';
files = dir(folder);

for i = 1:length(files)
    infile = fullfile(folder, files(i).name);
    if(strcmp(lower(infile(end-3:end)), '.svg'))
        outfile = strrep(infile, '.svg', '.emf');
        cmd = sprintf('"%s" --file "%s" --export-emf "%s"', exe, infile, outfile);
        system(cmd);
    end
end

