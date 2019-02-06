rootFolder = '/home/data/subjects';
outputFolder = '/home/data/EEG/processed/GhermanPhilastides/networkTimecourses';
allSubjects = dir(rootFolder);
GhermanSubjects = allSubjects(cellfun(@length, strfind({allSubjects.name}, 'GhermanPhiliastides')) > 0);

for subjectNumber = 1:length(GhermanSubjects)
    subjectFolder = fullfile(rootFolder, GhermanSubjects(subjectNumber).name);
    fprintf('%s\n', subjectFolder);
    [~, subFile] = fileparts(subjectFolder);
    outputFilename = fullfile(outputFolder, [subFile, '.mat']);
    
    
    infomapFolder = fullfile(subjectFolder, 'infomap', 'main');
    timecourseFolder= fullfile(subjectFolder, 'cifti/cifti_timeseries_normalwall');
    timecourseFilename = dir(timecourseFolder);
    timecourseFilename = timecourseFilename(cellfun(@length, strfind({timecourseFilename.name}, '.nii')) > 0);
    timecourseFilename = timecourseFilename.name;
    timecourse = ft_read_cifti_mod(fullfile(timecourseFolder, timecourseFilename));
    
    tmask = textread(fullfile(subjectFolder, 'fc_processed', 'main_all_tmask.txt'));
    
    %  regularized_ciftifile = fullfile(infomapFolder, 'rawassn_minsize400_regularized.dscalar.nii');
    
    groupnetworksfile = '/home/data/atlases/120_LR_minsize400_recolored_manualconsensus4.dtseries.nii';%groupnetworksfile = '/home/data/atlases/Networks_template.dscalar.nii';
    
    groupfile = ft_read_cifti_mod(groupnetworksfile);
    groupdata = groupfile.data;
    ncortverts = nnz(groupfile.brainstructure==1) + nnz(groupfile.brainstructure==2);
    groupdata = groupdata(1:ncortverts,1);
    
    
    
    %networkNumbers = unique(groupdata);
    networks.numbers = (0:17)';
    
    %potential_colors = [1 2 10 9 3 5 11 16 15 7 8 12 14 13];%[1 2 10 9 3 5 6 11 16 15 7 8 17 12 4 14 13];
    
    networks.colors = (...
        [.5, .5, .5; ... %0: unassigned = gray
        1, 0, 0; ... %1: default = red
        0, 0, .5; ... %2: secondary visual = blue
        1, 1, 0; ... %3: frontoparietal = yellow
        227/255, 172/255, 117/255; ... %4: primary visual = peach
        0, 1, 0; ... %5: dorsal attention = green
        230/255, 156/255, 231/255; ... %6: sensory(?) = lavender
        0, .5, .5; ... %7: ventral attention = teal
        0, 0, 0; ... %8: salience = black
        .5, 0, .5; ... %9: ciguloopercular = violet
        0, 1, 1; ... %10: motor hand = light blue
        227/255, 133/255, 45/255; ... %11: motor mouth = orange
        150/255, 79/255, 219/255; ... %12: opercular = lighter purple
        41/255, 82/255, 123/255; ... %13: ? = pale blue
        77/255, 236/255, 78/255; ... %14: ? = lime green
        42/255, 42/255, 252/255; ... %15: sky blue = parietal memory (?)
        1, 1, 1; ... %16: ? = white
        41/255, 129/255, 42/255 % 1,1,1 ... %17: motor foot = dark green
        ]);
    
    networks.labels = {'unassigned', 'default', 'secondary visual', 'frontoparietal',...
        'primary visual', 'dorsal attention', 'sensory', 'ventral attention', ...
        'salience', 'cingulooperular', 'motor hand', 'motor mouth', 'opercular', ...
        'pale blue', 'lime green', 'parietal memory', 'white', 'motor foot'}';
    
    output = NaN(length(tmask), length(networks.numbers));
    times = (((1:length(tmask)) - 1) .* 2)';
    for i = 1:length(networks.numbers)
        networks.vertexCounts(i, 1) = sum(groupdata == networks.numbers(i));
        networkInd = find(groupdata == networks.numbers(i));
        timeCounter = 0;
        for j = 1:length(tmask)
            if(tmask(j))
                timeCounter = timeCounter+ 1;
                output(j,i) = mean(timecourse.data(networkInd, timeCounter));
            end
        end
    end
    
    networks.times = timecourse.time';
    networks.timecourse = output;
    save(outputFilename, 'networks', '-v7.3');
    
end
