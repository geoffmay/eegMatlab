%% Subjects

%subjects = '/home/data/subjects/MAV_final_list.txt';

subjects = {'sub-GhermanPhiliastides01', 'sub-GhermanPhiliastides02', 'sub-GhermanPhiliastides03', 'sub-GhermanPhiliastides04', 'sub-GhermanPhiliastides05', 'sub-GhermanPhiliastides06', 'sub-GhermanPhiliastides07', 'sub-GhermanPhiliastides08', 'sub-GhermanPhiliastides09', 'sub-GhermanPhiliastides10', 'sub-GhermanPhiliastides11', 'sub-GhermanPhiliastides12', 'sub-GhermanPhiliastides13', 'sub-GhermanPhiliastides14', 'sub-GhermanPhiliastides15', 'sub-GhermanPhiliastides16', 'sub-GhermanPhiliastides17', 'sub-GhermanPhiliastides18', 'sub-GhermanPhiliastides19', 'sub-GhermanPhiliastides20', 'sub-GhermanPhiliastides21', 'sub-GhermanPhiliastides22', 'sub-GhermanPhiliastides23', 'sub-GhermanPhiliastides24'};
%'sub-FillmoreTest1'


%% Steps to run

choose_T1 =             0;
choose_T2 =             0;
T1_preprocess =         0;
T2_preprocess =         0;
fieldmap_process =      0;
B0map_process =         1;
DTI_preprocess =        1;
T1_segment =            0;
T1_hippocampal_subfields = 0;
surface_registration =  0;
BOLD_preprocess =       1;
BOLD_fcprocess =        0;
cifti_creation =        0;
cifti_correlation =     0;
cifti_distances =       0;


%% Parameters

%parallel processing params
parallelize = 0;
maxworkers = 8;

%MNI template
template = '/usr/local/fsl/data/standard/MNI152lin_T1_1mm_brain.nii.gz';

%subjects directory
subsdirectory = '/home/data/subjects/'; %individual subject folders are found here

%data signifiers in BIDS files
T1string = 'T1';
T2string = 'T2';
DTIstring = 'dwi1';
fieldmapstrings = {'FieldMap_BOLD'};
B0mapstrings = {'B0Map'};


%desired resolution of outputs
functional_normalized_voxdim = 2;
DTI_normalized_voxdim = 2;

%functional data params
functional_sequences = {'RSFC1'}; %task names in BIDS files
functional_sequence_TRs = {3.01}; %TRs for each fMRI sequence
fcprocess_sequences = {true}; %whether or not to do fc processing for each sequence
sliceorderinfo = {''}; %leave blank if standard
skipframes = 0; %frames to skip at beginning of each run

%params for functional data field map correction
apply_fieldmap_tofunc = true;
echospacing_func = .66522;%.6652;%.5943;%Acheiva seq %.6056;%Ingenia non-MB seq %ms

%params for DTI data B0 map correction
apply_B0map_todti = true;
readout_time = .048063;% ./ 2.5; %ms
DTI_PE_dir = 'P';
DTI_mb_factor = 3;

%segmentation options
freesurfer_use_T2 = false;

%fc processing params
FDthresh = .15;
mintime = 0; % minutes
lowpassmotionfilt = .1;

%cifti creation params
geodesic_smooth = 2.55;
make_smallwall_ciftis = false; %run for gradient parcellation




%%
if length(functional_sequence_TRs)==1
    functional_sequence_TRs = repmat(functional_sequence_TRs,1,length(functional_sequences));
end
if length(fcprocess_sequences)==1
    fcprocess_sequences = repmat(fcprocess_sequences,1,length(functional_sequences));
end
if length(sliceorderinfo)==1
    sliceorderinfo = repmat(sliceorderinfo,1,length(functional_sequences));
end



%% Get subjects


if iscell(subjects)
    
elseif exist(subjects,'file')
    
    subjects = textread(subjects,'%s');
    
elseif ischar(subjects)
    
    subjects = {subjects};
    
end



subcount = length(subjects);

seq_counter = cell(subcount,1);


%% Choose best T1 image

if choose_T1
    for subnum = 1:subcount
        subject = subjects{subnum};
        
        basefolder = [subsdirectory '/' subject '/'];
        anatomicalfolder = [basefolder '/anatomical/']; mkdir(anatomicalfolder)
        
        sessions = dir([basefolder '/ses-*']);
        
        all_T1files = cell(0,1);
        firstT1 = true;
        
        for sess = 1:length(sessions)
            
            T1files = dir([basefolder '/' sessions(sess).name '/anat/*' T1string '.nii.gz']);
            for file = 1:length(T1files)
                filename = [basefolder '/' sessions(sess).name '/anat/' T1files(file).name];
                [~,sizestr] = system(['fslsize ' filename ' -s']);
                tokens = tokenize(sizestr,' ');
                datasize = [str2num(tokens{3}) str2num(tokens{5}) str2num(tokens{7}) str2num(tokens{9})];
                if firstT1
                    T1size = datasize;
                    all_T1files{end+1} = filename;
                    firstT1 = false;
                else
                    if all(T1size==datasize)
                        all_T1files{end+1} = filename;
                    end
                end
            end
        end
        
        fslviewstring = 'fslview';
        for filenum = 1:length(all_T1files)
            fslviewstring = [fslviewstring ' ' all_T1files{filenum}];
            if length(all_T1files) > 1
                disp(['File ' num2str(filenum) ': ' all_T1files{filenum}])
            end
        end
        [~,~] = system([fslviewstring ' &']);
        
        if length(all_T1files) > 1
            bestnum = input(['Subject ' subject ': input index of best image (counting from bottom of fslview list): ']);
        else
            bestnum = 1;
        end
        copyfile(all_T1files{bestnum},[anatomicalfolder '/T1.nii.gz']);
        dlmwrite([anatomicalfolder '/T1source.txt'],all_T1files{bestnum},'delimiter','')
    end
end


%% Choose best T2 image

if choose_T2
    for subnum = 1:subcount
        subject = subjects{subnum};
        
        basefolder = [subsdirectory '/' subject '/'];
        anatomicalfolder = [basefolder '/anatomical/']; mkdir(anatomicalfolder)
        
        sessions = dir([basefolder '/ses-*']);
        
        all_T2files = cell(0,1);
        firstT2 = true;
        
        for sess = 1:length(sessions)
            
            T2files = dir([basefolder '/' sessions(sess).name '/anat/*' T2string '.nii.gz']);
            for file = 1:length(T2files)
                filename = [basefolder '/' sessions(sess).name '/anat/' T2files(file).name];
                [~,sizestr] = system(['fslsize ' filename ' -s']);
                tokens = tokenize(sizestr,' ');
                datasize = [str2num(tokens{3}) str2num(tokens{5}) str2num(tokens{7}) str2num(tokens{9})];
                if firstT2
                    T2size = datasize;
                    all_T2files{end+1} = filename;
                    firstT2 = false;
                else
                    if all(T2size==datasize)
                        all_T2files{end+1} = filename;
                    end
                end
            end
        end
        
        fslviewstring = 'fslview';
        for filenum = 1:length(all_T2files)
            fslviewstring = [fslviewstring ' ' all_T2files{filenum}];
            if length(all_T2files) > 1
                disp(['File ' num2str(filenum) ': ' all_T2files{filenum}])
            end
        end
        [~,~] = system([fslviewstring ' &']);
        
        if length(all_T2files) > 1
            bestnum = input(['Subject ' subject ': input index of best image (counting from bottom of fslview list): ']);
        else
            bestnum = 1;
        end
        copyfile(all_T2files{bestnum},[anatomicalfolder '/T2.nii.gz']);
        dlmwrite([anatomicalfolder '/T2source.txt'],all_T2files{bestnum},'delimiter','')
    end
end



warning off

%% Structural processing

if any([T1_preprocess T1_segment T2_preprocess surface_registration DTI_preprocess fieldmap_process])
    
    disp('STRUCTURAL PREPROCESSING')
    
    %% Set up parallel processing
    
    if logical(parallelize) && (subcount > 1)
        nworkers = min([subcount maxworkers]);
        if isempty(gcp('nocreate'))
            processingpool = parpool(nworkers);
        end
    else
        nworkers = 0;
    end
    
    %parfor (subnum = [1:subcount],nworkers)
    for subnum = 1:subcount
        
        systemresult = cell(0,2);
        
        subject = subjects{subnum};
        
        warning off
        
        try
            
            %% Set up folders
            
            
            basefolder = [subsdirectory '/' subject '/'];
            DTIfolder = [basefolder '/DTI/'];
            anatomicalfolder = [basefolder '/anatomical/'];
            fieldmapfolder = [basefolder '/fieldmap/'];
            functionalfolder = [basefolder '/functional/'];
            freesurferinitfolder = [basefolder '/freesurfer/'];
            freesurferfolder = [basefolder '/freesurfer/' subject '/'];
            fsLRfolder = [basefolder '/fs_LR/'];
            fc_processfolder = [basefolder '/fc_processed/'];
            ciftifolder = [basefolder '/cifti/'];
            logfolder = [basefolder '/processing_logs/']; mkdir(logfolder)
            
            %Create error tracking file to record issues with processing
            logfile = [logfolder datestr(now) '.txt']; tokens = tokenize(logfile,' '); logfile = [tokens{1} '_' tokens{2}];
            fid = fopen(logfile,'at');
            fclose(fid);
            
            
            
            %% T1 Preprocessing
            
            
            T1file = 'T1';
            
            if T1_preprocess
                
                cd(anatomicalfolder)
                
                dlmwrite(logfile,'preprocessing T1','-append','delimiter','')
                
                disp(['Subject ' subject ': skull-stripping T1, registering to MNI, downsampling to BOLD resolution, and fast segmenting']);
                
                % BIAS CORRECTION: NOT DONE
                %             %reorient, crop, and bias-correct T1 image
                %             [systemresult{end+1,1},systemresult{end+1,2}] = system(['fsl_anat -i ' T1file ' -o struct --clobber --noreg --nononlinreg --noseg --nosubcortseg']);
                %             %copyfile([anatomicalfolder '/struct.anat/' T1file '_orig2roi.mat'],[anatomicalfolder '/' T1file '_orig2biascorr.mat'])
                %             %copyfile([anatomicalfolder '/struct.anat/' T1file '_roi2orig.mat'],[anatomicalfolder '/' T1file '_biascorr2orig.mat'])
                %             copyfile([anatomicalfolder '/struct.anat/T1.nii.gz'],[anatomicalfolder '/' T1file '_orig_biascorrspace.nii.gz']);
                % %
                %
                %              T1file = [T1file '_biascorr'];
                %
                %             %remove folder used for processing
                %             copyfile([anatomicalfolder '/struct.anat/T1_biascorr.nii.gz'],[anatomicalfolder '/' T1file '.nii.gz']);
                %             copyfile([anatomicalfolder '/struct.anat/T1_biascorr_brain.nii.gz'],[anatomicalfolder '/' T1file '_bet.nii.gz']);
                %              [systemresult{end+1,1},systemresult{end+1,2}] = system(['chmod 777 -R ' anatomicalfolder '/struct.anat']);
                % %             try delete([anatomicalfolder '/struct.anat/*']);
                % %             rmdir([anatomicalfolder '/struct.anat'],'s'); catch; end
                %             delete([anatomicalfolder '/struct.anat/*']);
                %             rmdir([anatomicalfolder '/struct.anat']);
                %
                
                
                %brain-extract T1 image
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['bet ' T1file ' ' T1file '_bet -R -n -g -.1 -f .2 -m']);
                
                %register T1 to MNI
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T1file '_bet -interp spline -dof 12 -ref ' template ' -omat T1_2MNI.mat -out ' T1file '_bet_MNI']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T1file ' -interp spline -ref ' template ' -applyxfm -init T1_2MNI.mat -out ' T1file '_MNI']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system('convert_xfm -omat MNI_2T1.mat -inverse T1_2MNI.mat');
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T1file ' -interp spline -ref ' template ' -applyisoxfm ' num2str(functional_normalized_voxdim) ' -init T1_2MNI.mat -out ' T1file '_MNI_' num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim)]);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T1file '_bet_mask -interp nearestneighbour -ref ' template ' -applyisoxfm ' num2str(functional_normalized_voxdim) ' -init T1_2MNI.mat -out ' T1file '_bet_mask_MNI_' num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim)]);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' T1file '_bet_mask_MNI_' num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) ' -thr .5 -bin ' T1file '_bet_mask_MNI_' num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim)]);
                
                if DTI_normalized_voxdim ~= functional_normalized_voxdim
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T1file ' -interp spline -ref ' template ' -applyisoxfm ' num2str(DTI_normalized_voxdim) ' -init T1_2MNI.mat -out ' T1file '_MNI_' num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim)]);
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T1file '_bet_mask -interp nearestneighbour -ref ' template ' -applyisoxfm ' num2str(DTI_normalized_voxdim) ' -init T1_2MNI.mat -out ' T1file '_bet_mask_MNI_' num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim)]);
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' T1file '_bet_mask_MNI_' num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim) ' -thr .5 -bin ' T1file '_bet_mask_MNI_' num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim)]);
                end
                
                %FAST segment
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['fast -t 1 ' T1file '_bet']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' T1file '_bet_seg -thr 3 -bin ' T1file '_bet_wm']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' T1file '_bet_seg -uthr 1.5 -bin ' T1file '_bet_csf']);
                
                if any(cell2mat(systemresult(:,1)))
                    error(['Subject ' subject ': Problem detected with T1 preprocessing'])
                end
                
            else
                
                %track name
                %T1file = [T1file '_biascorr'];
                
            end
            
            
            %% T2 Preprocessing
            
            T2file = 'T2';
            
            if T2_preprocess
                
                %% register T2 to T1
                cd(anatomicalfolder)
                
                dlmwrite(logfile,'preprocessing T2','-append','delimiter','')
                
                disp(['Subject ' subject ': skull stripping T2, registering to T1 and through to MNI']);
                
                % BIAS CORRECTION: NOT DONE
                %                 %reorient, crop, and bias-correct T2 image
                %                 [systemresult{end+1,1},systemresult{end+1,2}] = system(['fsl_anat -i ' T2file ' -o struct -t T2 --clobber --noreg --nononlinreg --noseg --nosubcortseg']);
                %                 copyfile([anatomicalfolder '/struct.anat/T2_biascorr.nii.gz'],[anatomicalfolder '/' T2file '_biascorr.nii.gz'])
                %                 copyfile([anatomicalfolder '/struct.anat/T2_orig2roi.mat'],[anatomicalfolder '/' T2file '_orig2biascorr.mat'])
                %                 T2file = [T2file '_biascorr'];
                %
                %                 delete([anatomicalfolder '/struct.anat/*']);
                %                 rmdir([anatomicalfolder '/struct.anat']);
                
                
                %brain-extract T2 image
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['bet ' T2file ' ' T2file '_bet -R -g -.1']);
                
                %register T2 to T1, and concatenate with T1 to MNI
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T2file ' -ref ' T1file ' -dof 6 -omat T2_2T1_init.mat']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T2file ' -ref ' T1file ' -interp spline -dof 6 -cost bbr -wmseg ' T1file '_bet_wm -init T2_2T1_init.mat -omat T2_2T1.mat -out ' T2file '_T1space -schedule ${FSLDIR}/etc/flirtsch/bbr.sch']);
                
                [systemresult{end+1,1},systemresult{end+1,2}] = system('convert_xfm -omat T2_2MNI.mat -concat T1_2MNI.mat T2_2T1.mat');
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T2file ' -ref ' template ' -applyxfm -init T2_2MNI.mat -out ' T2file '_MNI -interp spline']);
                
                
                
                
                %             else
                %
                %                 %track T2 name
                %                 %T2file = [T2file '_biascorr'];
                
            end
            
            %% Fieldmap preparation
            if fieldmap_process
                mkdir(fieldmapfolder)
                cd(fieldmapfolder)
                dlmwrite(logfile,'calculating fieldmaps','-append','delimiter','')
                disp(['Subject ' subject ': copying field map data, calculating field maps, and averaging within scanning sessions']);
                
                FMcounter = 0;
                
                intendedfor_trackerfile = [fieldmapfolder '/fieldmap_intendedfor.txt'];
                delete(intendedfor_trackerfile);
                fid = fopen(intendedfor_trackerfile,'at'); %open the output file for writing
                fclose(fid);
                
                for fieldmapstringnum = 1:length(fieldmapstrings)
                    fieldmapstring = fieldmapstrings{fieldmapstringnum};
                    
                    %prepare commands to concatenate all field maps
                    allfieldmap_mergestring = ['fslmerge -t ' fieldmapstring '_phasediff_all '];
                    allfieldmapmag_mergestring = ['fslmerge -t ' fieldmapstring '_mag_all '];
                    
                    
                    sessions = dir([basefolder '/ses-*']);
                    
                    for sess = 1:length(sessions)
                        
                        %find field maps in this session
                        FMfiles_thissess = dir([basefolder '/' sessions(sess).name '/fmap/*' fieldmapstring '*magnitude1.nii.gz']);
                        
                        %If there are field map files
                        if ~isempty(FMfiles_thissess)
                            
                            %loop through field maps
                            for i = 1:length(FMfiles_thissess)
                                FMcounter = FMcounter+1;
                                
                                %Define the types of data
                                mag = FMfiles_thissess(i).name(1:end-7);
                                phase1 = [mag(1:end-10) 'phase1'];
                                phase2 = [mag(1:end-10) 'phase2'];
                                phasediff = [mag(1:end-10) 'phasediff'];
                                
                                copyfile([basefolder '/' sessions(sess).name '/fmap/' mag '.nii.gz'],[fieldmapfolder '/' mag '.nii.gz']);
                                copyfile([basefolder '/' sessions(sess).name '/fmap/' phase1 '.nii.gz'],[fieldmapfolder '/' phase1 '.nii.gz']);
                                copyfile([basefolder '/' sessions(sess).name '/fmap/' phase2 '.nii.gz'],[fieldmapfolder '/' phase2 '.nii.gz']);
                                
                                img = load_untouch_nii([fieldmapfolder '/' phase1 '.nii.gz']);
                                img.hdr.dime.scl_slope = 1;
                                img.hdr.dime.scl_inter = 0;
                                save_untouch_nii(img,[fieldmapfolder '/' phase1 '.nii.gz']);
                                
                                img = load_untouch_nii([fieldmapfolder '/' phase2 '.nii.gz']);
                                img.hdr.dime.scl_slope = 1;
                                img.hdr.dime.scl_inter = 0;
                                save_untouch_nii(img,[fieldmapfolder '/' phase2 '.nii.gz']);
                                
                                TE1 = [];
                                TE2 = [];
                                jsontext = textread([basefolder '/' sessions(sess).name '/fmap/' phase1 '.json'],'%s','delimiter','\n');
                                intendedfor_files = cell(0,1);
                                for l = 1:length(jsontext)
                                    if ~isempty(strfind(jsontext{l},'"EchoTime"'))
                                        TE1 = str2num(jsontext{l}((length('"EchoTime"')+2) : end)) .* 1000;
                                        break
                                    end
                                end
                                
                                jsontext = textread([basefolder '/' sessions(sess).name '/fmap/' phase2 '.json'],'%s','delimiter','\n');
                                intendedfor_files = cell(0,1);
                                for l = 1:length(jsontext)
                                    if ~isempty(strfind(jsontext{l},'"EchoTime"'))
                                        TE2 = str2num(jsontext{l}((length('"EchoTime"')+2) : end)) .* 1000;
                                        break
                                    end
                                end
                                
                                fieldmap_deltaTE = TE2 - TE1;
                                
                                [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' fieldmapfolder '/' phase1 ' -mul 3.14159 -div 2048 ' fieldmapfolder '/' phase1 '_rad -odt float']);
                                [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' fieldmapfolder '/' phase2 ' -mul 3.14159 -div 2048 ' fieldmapfolder '/' phase2 '_rad -odt float']);
                                
                                
                                %Unwrap phase images
                                [systemresult{end+1,1},systemresult{end+1,2}] = system(['prelude -a ' fieldmapfolder '/' mag ' -p ' fieldmapfolder '/' phase1 '_rad -o ' fieldmapfolder '/' phase1 '_unwrapped']);
                                [systemresult{end+1,1},systemresult{end+1,2}] = system(['prelude -a ' fieldmapfolder '/' mag ' -p ' fieldmapfolder '/' phase2 '_rad -o ' fieldmapfolder '/' phase2 '_unwrapped']);
                                
                                %Calculate rad/sec fieldmap image
                                [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' fieldmapfolder '/' phase1 '_unwrapped -sub ' fieldmapfolder '/' phase2 '_unwrapped -mul 1000 -div ' num2str(fieldmap_deltaTE) ' ' fieldmapfolder '/' phasediff ' -odt float']);
                                
                                
                                jsontext = textread([basefolder '/' sessions(sess).name '/fmap/' mag '.json'],'%s','delimiter','\n');
                                intendedfor_files = cell(0,1);
                                for l = 1:length(jsontext)
                                    if ~isempty(strfind(jsontext{l},'"IntendedFor"'))
                                        
                                        for l2 = l:(length(jsontext)-1)
                                            quoteinds = strfind(jsontext{l2},'"');
                                            this_filename = [basefolder jsontext{l2}((quoteinds(end-1)+1) : (quoteinds(end)-1))];
                                            [~,file,extension] = fileparts(this_filename);
                                            intendedfor_files{end+1} = [file extension];
                                        end
                                    end
                                end
                                
                                for l = 1:length(intendedfor_files)
                                    dlmwrite(intendedfor_trackerfile,[intendedfor_files{l} ' ' fieldmapfolder '/' mag ' ' num2str(fieldmap_deltaTE)],'-append','delimiter','');%write the data to the output file
                                end
                                
                                
                                
                                %Prepare to concatenate all images
                                allfieldmap_mergestring = [allfieldmap_mergestring fieldmapfolder '/' phasediff ' '];
                                allfieldmapmag_mergestring = [allfieldmapmag_mergestring fieldmapfolder '/' mag ' '];
                                
                                
                            end
                            
                            
                        end
                        
                    end
                    
                    if length(FMcounter)>1
                        
                        %Concatenate magnitude images
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(allfieldmapmag_mergestring);
                        %Register magnitude images to each other
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['mcflirt -in ' fieldmapstring '_mag_all -refvol 0 -mats']);
                        %Average magnitude images
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' fieldmapstring '_mag_all_mcf -Tmean ' fieldmapstring '_mag_mean']);
                        
                        %Concatenate fieldmap images
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(allfieldmap_mergestring);
                        %Apply previously calculated registration
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['applyxfm4D ' fieldmapstring '_phasediff_all ' fieldmapstring '_phasediff_all ' fieldmapstring '_phasediff_all_mcf ' fieldmapstring '_mag_all_mcf.mat -userprefix MAT_']);
                        %Average fieldmap images
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' fieldmapstring '_phasediff_all_mcf -Tmean ' fieldmapstring '_phasediff_mean']);
                        
                        
                        %Remove unneeded files
                        delete([fieldmapfolder '/' fieldmapstring '_mag_all.nii.gz']);
                        delete([fieldmapfolder '/' fieldmapstring '_mag_all_mcf.nii.gz']);
                        delete([fieldmapfolder '/' fieldmapstring '_phasediff_all.nii.gz']);
                        delete([fieldmapfolder '/' fieldmapstring '_phasediff_all_mcf.nii.gz']);
                        delete(['' fieldmapstring '_mag_all_mcf.mat/*'])
                        try rmdir(['' fieldmapstring '_mag_all_mcf.mat'],'s'); catch; end
                        
                        
                    else
                        
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['cp ' fieldmapfolder '/' mag '.nii.gz ' fieldmapstring '_mag_mean.nii.gz']);
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['cp ' fieldmapfolder '/' phasediff '.nii.gz ' fieldmapstring '_phasediff_mean.nii.gz']);
                        
                    end
                    
                end
                %Check for problems
                if any(cell2mat(systemresult(:,1)))
                    error(['Subject ' subject ': Problem detected with fieldmap processing'])
                end
                
            end
            
            %% B0map processing
            if B0map_process
                mkdir(fieldmapfolder)
                cd(fieldmapfolder)
                dlmwrite(logfile,'calculating B0map coefficients','-append','delimiter','')
                disp(['Subject ' subject ': copying B0 map data, calculating B0 map coefficients, and averaging within scanning sessions']);
                
                B0counter = 0;
                
                intendedfor_trackerfile = [fieldmapfolder '/B0map_intendedfor.txt'];
                delete(intendedfor_trackerfile);
                fid = fopen(intendedfor_trackerfile,'at'); %open the output file for writing
                fclose(fid);
                
                for B0mapstringnum = 1:length(B0mapstrings)
                    B0mapstring = B0mapstrings{B0mapstringnum};
                    
                    %prepare commands to concatenate all field maps
                    B0map_mergestring = ['fslmerge -t ' B0mapstring '_fieldcoef_all '];
                    
                    
                    sessions = dir([basefolder '/ses-*']);
                    
                    for sess = 1:length(sessions)
                        
                        %find B0 maps in this session
                        B0files_thissess = dir([basefolder '/' sessions(sess).name '/fmap/*' B0mapstring '*P1.nii.gz']);
                        
                        %If there are B0 map files
                        if ~isempty(B0files_thissess)
                            
                            %loop through B0 maps
                            for i = 1:length(B0files_thissess)
                                B0counter = B0counter+1;
                                
                                this_B0file = [fieldmapfolder '/' B0files_thissess(i).name(1:end-10)];
                                
                                paramsfile = [this_B0file '.params'];
                                
                                delete(paramsfile)
                                fid = fopen(paramsfile,'at'); %open the output file for writing
                                fclose(fid);
                                
                                Pfiles = dir([basefolder '/' sessions(sess).name '/fmap/' B0files_thissess(i).name(1:end-10) '_P*.nii.gz']);
                                Afiles = dir([basefolder '/' sessions(sess).name '/fmap/' B0files_thissess(i).name(1:end-10) '_A*.nii.gz']);
                                
                                fslstr = ['fslmerge -t ' this_B0file];
                                for f = 1:length(Pfiles)
                                    fslstr = [fslstr ' ' basefolder '/' sessions(sess).name '/fmap/' Pfiles(f).name];
                                    dlmwrite(paramsfile,[0 1 0 readout_time],'delimiter',' ','-append');
                                end
                                for f = 1:length(Afiles)
                                    fslstr = [fslstr ' ' basefolder '/' sessions(sess).name '/fmap/' Afiles(f).name];
                                    dlmwrite(paramsfile,[0 -1 0 readout_time],'delimiter',' ','-append');
                                end
                                
                                [systemresult{end+1,1},systemresult{end+1,2}] = system(fslstr);
                                [systemresult{end+1,1},systemresult{end+1,2}] = system(['topup --imain=' this_B0file ' --datain=' paramsfile ' --config=b02b0.cnf --out=' this_B0file '_topup']);
                                
                                
                                jsontext = textread([basefolder '/' sessions(sess).name '/fmap/' B0files_thissess(i).name(1:end-7) '.json'],'%s','delimiter','\n');
                                intendedfor_files = cell(0,1);
                                for l = 1:length(jsontext)
                                    if ~isempty(strfind(jsontext{l},'"IntendedFor"'))
                                        
                                        for l2 = l:(length(jsontext)-1)
                                            quoteinds = strfind(jsontext{l2},'"');
                                            this_filename = [basefolder jsontext{l2}((quoteinds(end-1)+1) : (quoteinds(end)-1))];
                                            [~,file,extension] = fileparts(this_filename);
                                            intendedfor_files{end+1} = [file extension];
                                        end
                                    end
                                end
                                
                                for l = 1:length(intendedfor_files)
                                    dlmwrite(intendedfor_trackerfile,[intendedfor_files{l} ' ' this_B0file '_topup.nii.gz'],'-append','delimiter','');%write the data to the output file
                                end
                                
                            end
                        end
                    end
                end
            end
                                
                                
                                
            %% DTI preprocessing
            if DTI_preprocess
                mkdir(DTIfolder)
                cd(DTIfolder)
                
                dlmwrite(logfile,'preprocessing DTI data','-append','delimiter','')
                disp(['Subject ' subject ': copying diffusion data, fitting tensors, and registering to MNI']);
                
                DTIfiles = cell(0,1);
                DTIcounter = 0;
                TEs = [];
                metrics_touse = {'FA','MD','L1','L2','L3'};
                
                sess_tracker = [DTIfolder '/DTIruns_sessions.txt'];
                delete(sess_tracker);
                fid = fopen(sess_tracker,'at'); %open the output file for writing
                fclose(fid);
                
                sessions = dir([basefolder '/ses-*']);
                
                for sess = 1:length(sessions)
                    DTIfiles_thissess = dir([basefolder '/' sessions(sess).name '/dwi/*' DTIstring '.nii.gz']);
                    for i = 1:length(DTIfiles_thissess)
                        DTIcounter = DTIcounter+1;
                        DTIfiles{DTIcounter} = [DTIfolder '/' DTIfiles_thissess(i).name];
                        copyfile([basefolder '/' sessions(sess).name '/dwi/' DTIfiles_thissess(i).name],DTIfiles{DTIcounter});
                        copyfile([basefolder '/' sessions(sess).name '/dwi/' DTIfiles_thissess(i).name(1:end-7) '.bval'],[DTIfiles{DTIcounter}(1:end-7) '.bval']);
                        copyfile([basefolder '/' sessions(sess).name '/dwi/' DTIfiles_thissess(i).name(1:end-7) '.bvec'],[DTIfiles{DTIcounter}(1:end-7) '.bvec']);
                        
                        DTIfullfile = DTIfiles{DTIcounter}(1:end-7);
                        [DTIpath, DTIfile, DTIext] = fileparts(DTIfullfile);
                        
                        %track session number
                        dlmwrite(sess_tracker,[DTIfile '_ec_FA_MNI.nii.gz ' sessions(sess).name(end-1:end)],'-append','delimiter','');%write the data to the output file
                        
                    end
                end
                
                if logical(DTIcounter)
                    for filenum = 1:DTIcounter
                        
                        DTIfullfile = DTIfiles{filenum}(1:end-7);
                        [DTIpath, DTIfile, DTIext] = fileparts(DTIfullfile);
                        
                        %Figure out which image is the (first) B0
                        bvals = load([DTIfullfile '.bval']);
                        B0inds = find(bvals==0) - 1;
                        refimage_ind = num2str(B0inds(1));
                        
                        if apply_B0map_todti
                            
                            %determine correct B0map file to use
                            [intendedforfiles, B0mapfiles] = textread([fieldmapfolder '/B0map_intendedfor.txt'],'%s%s');
                            index_inlist = strcmp(intendedforfiles,[DTIfile '.nii.gz']);
                            this_B0mapfile = B0mapfiles{index_inlist};
                            this_paramsfile = [this_B0mapfile(1:end-7) '.params'];
                            
                            %Make index file
                            if strcmp(DTI_PE_dir,'P')
                                B0index = 1;
                            elseif strcmp(DTI_PE_dir,'A')
                                params = load(this_paramsfile);
                                inds_withAdir = find(params(:,2)==-1);
                                B0index = inds_withAdir(1);
                            end
                            dlmwrite([DTIfullfile '_index.txt'],repmat(B0index,length(bvals),1),'delimiter',' ');
                            
                            % quick T1 reg and mask creation
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslselectvols -i ' DTIfile ' -o ' DTIfile '_b0 --vols=' refimage_ind]);
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' DTIfile '_b0 -ref ' anatomicalfolder '/' T1file '_bet -dof 6 -omat ' DTIfile '_b0_2T1_init.mat']);
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['convert_xfm -omat T1_2' DTIfile '_b0.mat -inverse ' DTIfile '_b0_2T1_init.mat']);
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' anatomicalfolder '/' T1file '_bet_mask -interp nearestneighbour -ref ' DTIfile '_b0 -applyxfm -init T1_2' DTIfile '_b0.mat -out ' T1file '_' DTIfile '_bet_mask']);
                            delete([DTIfile '_b0.nii.gz']);
                            delete([DTIfile '_b0_2T1_init.mat']);
                            
                            %register B0 coeffs to DTI file
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslselectvols -i ' this_B0mapfile(1:end-13) ' -o ' this_B0mapfile(1:end-13) '_matchDTI --vols=' B0index]);
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' this_B0mapfile(1:end-13) '_matchDTI -ref ' DTIfile '_b0 -dof 6 -omat ' this_B0mapfile(1:end-13) '_matchDTI_2' DTIfile '.mat']);
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' this_B0mapfile(1:end-7) ' -ref ' DTIfile '_b0 -applyxfm -init ' this_B0mapfile(1:end-13) '_matchDTI_2' DTIfile '.mat -interp spline -out ' this_B0mapfile(1:end-7) '_' DTIfile]);
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' this_B0mapfile(1:end-7) '_fieldcoef -ref ' DTIfile '_b0 -applyxfm -init ' this_B0mapfile(1:end-13) '_matchDTI_2' DTIfile '.mat -interp spline -out ' this_B0mapfile(1:end-7) '_' DTIfile '_fieldcoef']);
                            
                            %Improved eddy correction
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['eddy_cuda6.5 --imain=' DTIfullfile ' --mask=' T1file '_' DTIfile '_bet_mask --acqp=' this_paramsfile ' --index=' DTIfullfile '_index.txt --bvecs=' DTIfullfile '.bvec --bvals=' DTIfullfile '.bval --topup=' this_B0mapfile(1:end-7) '_' DTIfile ' --repol --ol_type==mb --mb=' num2str(DTI_mb_factor) ' --out=' DTIfullfile '_ec']);
                            
                        else
                            
                            %Eddy correction
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['eddy_correct ' DTIfile ' ' DTIfile '_ec ' refimage_ind ' spline']);
                            
                        end
                        
                        %Pull out B0 image
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslselectvols -i ' DTIfile '_ec -o ' DTIfile '_ec_b0 --vols=' refimage_ind]);
                        
                        %calculate DTI->T1 registration
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' DTIfile '_ec_b0 -ref ' anatomicalfolder '/' T1file '_bet -dof 6 -omat ' DTIfile '_ec_b0_2T1_init.mat']);
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' DTIfile '_ec_b0 -ref '  anatomicalfolder '/' T1file ' -dof 6 -cost bbr -wmseg ' anatomicalfolder '/' T1file '_bet_wm -init ' DTIfile '_ec_b0_2T1_init.mat -omat ' DTIfile '_ec_b0_2T1.mat -schedule ${FSLDIR}/etc/flirtsch/bbr.sch']);
                        delete([DTIfile '_ec_b0_2T1_init.mat']);
                        
                        %concatenate DTI->T1 registration and T1->MNI transform and apply to DTI data
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['convert_xfm -omat ' DTIfile '_ec_b0_2MNI.mat -concat ' anatomicalfolder '/T1_2MNI.mat ' DTIfile '_ec_b0_2T1.mat']);
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' DTIfile '_ec -interp spline -ref ' anatomicalfolder T1file '_MNI_' num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim) ' -applyxfm -init ' DTIfile '_ec_b0_2MNI.mat -out ' DTIfile '_ec_MNI']);
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' DTIfile '_ec_b0 -interp spline -ref ' anatomicalfolder T1file '_MNI_' num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim) ' -applyxfm -init ' DTIfile '_ec_b0_2MNI.mat -out ' DTIfile '_ec_b0_MNI']);
                        
                        %get T1 mask into DTI space
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['convert_xfm -omat T1_2' DTIfile '_ec_b0.mat -inverse ' DTIfile '_ec_b0_2T1.mat']);
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' anatomicalfolder '/' T1file '_bet_mask -interp nearestneighbour -ref ' DTIfile '_ec_b0 -applyxfm -init T1_2' DTIfile '_ec_b0.mat -out ' T1file '_' DTIfile '_bet_mask']);
                        
                        %Calc DTI metrics
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['dtifit -k ' DTIfile '_ec -m ' T1file '_' DTIfile '_bet_mask -r ' DTIfile '.bvec -b ' DTIfile '.bval -o ' DTIfile '_ec']);
                        
                        for metricnum = 1:length(metrics_touse)
                            this_metric = metrics_touse{metricnum};
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' DTIfile '_ec_' this_metric ' -interp spline -ref ' anatomicalfolder T1file '_MNI_' num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim) ' -applyxfm -init ' DTIfile '_ec_b0_2MNI.mat -out ' DTIfile '_ec_' this_metric '_MNI']);
                            
                            metric_img = load_untouch_nii([DTIfile '_ec_FA_MNI.nii.gz']);
                            if filenum == 1
                                metric_avg{metricnum} = metric_img;
                                metric_avg{metricnum}.img = metric_avg{metricnum}.img / DTIcounter;
                            else
                                metric_avg{metricnum}.img = metric_avg{metricnum}.img + (metric_img.img / DTIcounter);
                            end
                        end
                        
                    end
                    
                    for metricnum = 1:length(metrics_touse)
                        this_metric = metrics_touse{metricnum};
                        save_untouch_nii(metric_avg{metricnum},['DTI_avg_ec_' this_metric '_MNI.nii.gz']);
                    end
                    
                    
                    if any(cell2mat(systemresult(:,1)))
                        error(['Subject ' subject ': Problem detected with DTI processing'])
                    end
                    
                end
            end
            
            %% T1 Segmentation
            
            if T1_segment
                
                disp(['Subject ' subject ': segmenting T1 and making compartment masks']);
                dlmwrite(logfile,'segmenting T1','-append','delimiter','')
                
                try [~,~,~] = rmdir(freesurferinitfolder,'s'); catch; end
                
                mkdir(freesurferinitfolder)
                
                %copy T1 to freesurfer folder
                copyfile([anatomicalfolder '/' T1file '.nii.gz'],[freesurferinitfolder '/' T1file '.nii.gz'])
                cd(freesurferinitfolder)
                
                %transform T1 to .mgz format
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['setenv FREESURFER_HOME /usr/local/freesurfer6/freesurfer;'...
                    'source $FREESURFER_HOME/SetUpFreeSurfer.csh;'...
                    'recon-all -subjid ' subject ' -sd ' freesurferinitfolder ' -i ' freesurferinitfolder '/' T1file '.nii.gz']);
                
                %segment T1
                if freesurfer_use_T2
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['mri_convert ' anatomicalfolder '/' T2file '_T1space.nii.gz ' freesurferinitfolder '/' T2file '_T1space.mgz']);
                    
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['setenv FREESURFER_HOME /usr/local/freesurfer6/freesurfer;'...
                        'source $FREESURFER_HOME/SetUpFreeSurfer.csh;'...
                        'recon-all -subjid ' subject ' -sd ' freesurferinitfolder ' -T2 ' freesurferinitfolder '/' T2file '_T1space.mgz -T2pial -all']);
                    
                else
                    
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['setenv FREESURFER_HOME /usr/local/freesurfer6/freesurfer;'...
                        'source $FREESURFER_HOME/SetUpFreeSurfer.csh;'...
                        'recon-all -subjid ' subject ' -sd ' freesurferinitfolder ' -all']);
                    
                end
                
                %get segmentation masks into nifti format
                systemresult = make_fs_masks_mutualdist(freesurferfolder,functional_normalized_voxdim,[anatomicalfolder '/T1_2MNI.mat'],template,systemresult);
                
                if any(cell2mat(systemresult(:,1)))
                    error(['Subject ' subject ': Problem detected with Freesurfer'])
                end
                
            end
            
            %% Surface-based registration
            
            if surface_registration
                
                disp(['Subject ' subject ': surface-based registration']);
                dlmwrite(logfile,'conducting surface registration','-append','delimiter','')
                
                mkdir(fsLRfolder)
                
                prevsize = size(systemresult,1);
                systemresult = PostFreeSurferPipeline_fsavg2fslr_long_singlesub(subject,T1file,freesurferfolder,fsLRfolder,anatomicalfolder,systemresult);
                if any(cell2mat(systemresult(:,1)))
                    error(['Subject ' subject ': Problem detected with surface registration'])
                end
                
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['csh /home/data/scripts/Cifti_creation/create_ribbon_singlesub.csh ' subject ' ' fsLRfolder '/MNI/ ' T1file '_MNI ' anatomicalfolder ' ' num2str(functional_normalized_voxdim)]);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['csh /home/data/scripts/Cifti_creation/create_ribbon_singlesub.csh ' subject ' ' fsLRfolder '/NativeVol/ ' T1file ' ignore ignore']);
                if any(cell2mat(systemresult(:,1)))
                    error(['Subject ' subject ': Problem detected with ribbon creation'])
                end
                
            end
            
            %%
            cprintf('*blue',['Subject ' subject ' COMPLETED structural preprocessing without error!\n'])
            dlmwrite(logfile,'all processing complete.','-append','delimiter','')
            
        catch errorinfo
            
            dlmwrite(logfile,'ERROR:','-append','delimiter','')
            
            errormessage = [errorinfo.stack(end).file ', line ' num2str(errorinfo.stack(end).line) ': ' errorinfo.message];
            
            probleminds = find(cell2mat(systemresult(:,1)));
            for i = 1:length(probleminds)
                dlmwrite(logfile,systemresult{probleminds(i),2},'-append','delimiter','');
            end
            dlmwrite(logfile,errormessage,'-append','delimiter','');
            
            cprintf('err',['Subject ' subject ' FAILED to complete structural preprocessing! Check error log ' logfile '\n'])
            
        end
        
        
    end
    
    clear systemresult
    
    
    try delete(gcp('nocreate')); catch; end
    
    
end
clear logfile T1file

if BOLD_preprocess
    
    disp('FUNCTIONAL PREPROCESSING')
    
    for subnum = 1:subcount
        
        subject = subjects{subnum};
        
        warning off
        
        try
            
            %% Set up folders
            
            
            basefolder = [subsdirectory '/' subject '/'];
            DTIfolder = [basefolder '/DTI/'];
            anatomicalfolder = [basefolder '/anatomical/'];
            fieldmapfolder = [basefolder '/fieldmap/'];
            functionalfolder = [basefolder '/functional/']; mkdir(functionalfolder)
            freesurferinitfolder = [basefolder '/freesurfer/'];
            freesurferfolder = [basefolder '/freesurfer/' subject '/'];
            fsLRfolder = [basefolder '/fs_LR/'];
            fc_processfolder = [basefolder '/fc_processed/'];
            ciftifolder = [basefolder '/cifti/'];
            logfolder = [basefolder '/processing_logs/'];
            
            %Create error tracking file to record issues with processing
            logfile = [logfolder datestr(now) '.txt']; tokens = tokenize(logfile,' '); logfile = [tokens{1} '_' tokens{2}];
            fid = fopen(logfile,'at');
            fclose(fid);
            
            T1file = 'T1';
            
            
            %% Functional preprocessing
            
            mkdir(functionalfolder)
            cd(functionalfolder)
            
            dlmwrite(logfile,'preprocessing BOLD data','-append','delimiter','')
            disp(['Subject ' subject ': slice time correction, motion correction, T1 registration, atlas registration, and intensity normalization']);
            
            
            
            %% Copy data
            
            brainmaskfile = [anatomicalfolder '/' T1file '_bet_mask_MNI_' num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) '.nii.gz'];
            brainmask = load_untouch_nii(brainmaskfile);
            
            for seq = 1:length(functional_sequences)
                seqname = functional_sequences{seq};
                TR = functional_sequence_TRs{seq};
                
                %make a file that tracks session number
                sess_tracker = [functionalfolder '/BOLDruns_sessions_' functional_sequences{seq} '.txt'];
                delete(sess_tracker);
                fid = fopen(sess_tracker,'at'); %open the output file for writing
                fclose(fid);
                
                seqname = functional_sequences{seq};
                
                allfiles = dir([functionalfolder '/*' seqname '*']);
                for i = 1:length(allfiles)
                    delete([functionalfolder '/' allfiles(i).name]);
                end
                
                BOLDrunnames = cell(0,1);
                sesscounter = [];
                TEs = [];
                sessions = dir([basefolder '/ses-*']);
                for sess = 1:length(sessions)
                    tokens = tokenize(sessions(sess).name,'-');
                    sessnum = str2num(tokens{end});
                    thisseq_files = dir([basefolder '/' sessions(sess).name '/func/*' seqname '*.nii.gz']);
                    
                    for f = 1:length(thisseq_files)
                        
                        [~,sizestr] = system(['fslsize ' basefolder '/' sessions(sess).name '/func/' thisseq_files(f).name ' -s']);
                        
                        tokens = tokenize(sizestr,' ');
                        datasize = [str2num(tokens{3}) str2num(tokens{5}) str2num(tokens{7}) str2num(tokens{9})];
                        
                        if datasize(4) > 20
                            
                            BOLDrunnames{end+1} = [functionalfolder '/' thisseq_files(f).name(1:end-7)];
                            sesscounter(end+1) = sessnum;
                            if skipframes > 0
                                %load data, trim first few frames
                                data = load_untouch_nii([basefolder '/' sessions(sess).name '/func/' thisseq_files(f).name]);
                                data.img = data.img(:,:,:,skipframes+1:end);
                                data.hdr.dime.dim(5) = size(data.img,4);
                                
                                %save to preprocessing folder
                                save_untouch_nii(data,[functionalfolder '/' thisseq_files(f).name])
                                
                            else
                                [~,~] = system(['cp ' basefolder '/' sessions(sess).name '/func/' thisseq_files(f).name ' ' thisseq_files(f).name]);
                            end
                            
                            jsontext = textread([basefolder '/' sessions(sess).name '/func/' thisseq_files(f).name(1:end-7) '.json'],'%s','delimiter','\n');
                            for l = 1:length(jsontext)
                                if ~isempty(strfind(jsontext{l},'"EchoTime":'))
                                    TEs(length(BOLDrunnames)) = str2num(jsontext{l}((length('"EchoTime":')+1) : end)) .* 1000;
                                end
                            end
                        end
                    end
                end
                
                for runnum = 1:length(BOLDrunnames)
                    systemresult{runnum} = cell(0,2);
                end
                
                
                %% Preprocessing data
                
                if logical(parallelize) && (length(BOLDrunnames) > 1)
                    nworkers = min([length(BOLDrunnames) maxworkers]);
                    if isempty(gcp('nocreate'))
                        processingpool = parpool(nworkers);
                    end
                else
                    nworkers = 0;
                end
                
                finalBOLDfilenames = cell(length(BOLDrunnames),1);
                finalBOLDavgfilenames = cell(length(BOLDrunnames),1);
                
                %for runnum = 1:length(BOLDrunnames)
                parfor (runnum = [1:length(BOLDrunnames)],nworkers)
                    BOLDfile = BOLDrunnames{runnum};
                    
                    %get short filename
                    [~,BOLDfilename,~] = fileparts(BOLDfile);
                    
                    
                    %motion-correct raw BOLD images to first image to get "true" motion plots
                    [systemresult{runnum}{end+1,1},systemresult{runnum}{end+1,2}] = system(['mcflirt -in ' BOLDfile ' -refvol 0 -plots']);
                    
                    %delete resulting motion corrected timeseries
                    delete([BOLDfile '_mcf.nii.gz'])
                    
                    %slice-time correct BOLD images
                    [systemresult{runnum}{end+1,1},systemresult{runnum}{end+1,2}] = system(['slicetimer -i ' BOLDfile ' -o ' BOLDfile '_st -r ' num2str(TR) ' ' sliceorderinfo{seq}]);
                    
                    %motion-correct BOLD images to first image to get motion correction registrations
                    [systemresult{runnum}{end+1,1},systemresult{runnum}{end+1,2}] = system(['mcflirt -in ' BOLDfile '_st -refvol 0']);
                    
                    %get avg BOLD for this run
                    [systemresult{runnum}{end+1,1},systemresult{runnum}{end+1,2}] = system(['fslmaths ' BOLDfile '_st_mcf -Tmean ' BOLDfile '_st_mcf_avg']);
                    
                    
                    if apply_fieldmap_tofunc
                        %Figure out which field map file should be used for this run
                        [intendedforfiles, fmapfiles, fmap_deltaTEs] = textread([fieldmapfolder '/fieldmap_intendedfor.txt'],'%s%s%s');
                        index_inlist = strcmp(intendedforfiles,[BOLDfilename '.nii.gz']);
                        this_fieldmapmag = fmapfiles{index_inlist};
                        this_fieldmapphasediff = [this_fieldmapmag(1:end-10) 'phasediff'];
                        deltaTE = fmap_deltaTEs{index_inlist};
                        
                        %calculate BOLD unwarp
                        warp_outfolder = [BOLDfile '_st_mcf_avg_warpdir'];
                        [systemresult{runnum}{end+1,1},systemresult{runnum}{end+1,2}] = system(['csh /home/data/scripts/Processing_Pipeline/epi_unwarp_WashU.csh -o ' warp_outfolder ' -m ' this_fieldmapmag ' -p ' this_fieldmapphasediff ' -e ' BOLDfile '_st_mcf_avg -dwell ' num2str(echospacing_func) ' -te ' num2str(TEs(runnum)) ' -delta ' num2str(deltaTE) ' -nounwrap -meanmap']);
                        copyfile([warp_outfolder '/epi_unwarped.nii.gz'],[BOLDfile '_st_mcf_avg_unwarp.nii.gz'])
                        copyfile([warp_outfolder '/epi_0_unwarp_shift_warp.nii.gz'],[BOLDfile '_st_mcf_avg_warp.nii.gz'])
                        try rmdir(warp_outfolder,'s'); catch; end
                        
                        %calculate unwarped BOLD->T1 registration
                        [systemresult{runnum}{end+1,1},systemresult{runnum}{end+1,2}] = system(['flirt -in ' BOLDfile '_st_mcf_avg_unwarp -ref ' anatomicalfolder T1file '_bet -dof 6 -omat ' BOLDfile '_st_mcf_avg_unwarp_2T1_init.mat']);
                        [systemresult{runnum}{end+1,1},systemresult{runnum}{end+1,2}] = system(['flirt -in ' BOLDfile '_st_mcf_avg_unwarp -ref ' anatomicalfolder T1file '_bet -dof 6 -cost bbr -wmseg ' anatomicalfolder T1file '_bet_wm -init ' BOLDfile '_st_mcf_avg_unwarp_2T1_init.mat -omat ' BOLDfile '_st_mcf_avg_unwarp_2T1.mat -schedule ${FSLDIR}/etc/flirtsch/bbr.sch']);% -out ' BOLDfile '_st_mcf_avg_T1space']);
                        delete([BOLDfile '_st_mcf_avg_unwarp_2T1_init.mat']);
                        
                        %concatenate BOLD->T1 registration and T1->MNI transform
                        [systemresult{runnum}{end+1,1},systemresult{runnum}{end+1,2}] = system(['convert_xfm -omat ' BOLDfile '_st_mcf_avg_unwarp_2MNI.mat -concat ' anatomicalfolder T1file '_2MNI.mat ' BOLDfile '_st_mcf_avg_unwarp_2T1.mat']);
                        
                        %apply BOLD->T1 registration+warp and T1->MNI transform
                        [systemresult{runnum}{end+1,1},systemresult{runnum}{end+1,2}] = system(['applywarp -i ' BOLDfile '_st_mcf_avg -r  ' anatomicalfolder T1file '_MNI_' num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) ' -o ' BOLDfile '_st_mcf_avg_unwarp_MNI -w ' BOLDfile '_st_mcf_avg_warp --postmat=' BOLDfile '_st_mcf_avg_unwarp_2MNI.mat --interp=spline']);
                        [systemresult{runnum}{end+1,1},systemresult{runnum}{end+1,2}] = system(['applywarp -i ' BOLDfile '_st_mcf -r  ' anatomicalfolder T1file '_MNI_' num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) ' -o ' BOLDfile '_st_mcf_unwarp_MNI -w ' BOLDfile '_st_mcf_avg_warp --postmat=' BOLDfile '_st_mcf_avg_unwarp_2MNI.mat --interp=spline']);
                        
                        finalBOLDfilenames{runnum} = [BOLDfile '_st_mcf_unwarp_MNI.nii.gz'];
                        finalBOLDavgfilenames{runnum} = [BOLDfile '_st_mcf_avg_unwarp_MNI.nii.gz'];
                        
                    else
                        
                        %calculate BOLD->T1 registration
                        [systemresult{runnum}{end+1,1},systemresult{runnum}{end+1,2}] = system(['flirt -in ' BOLDfile '_st_mcf_avg -ref ' anatomicalfolder T1file '_bet -dof 6 -omat ' BOLDfile '_st_mcf_avg_2T1_init.mat']);
                        [systemresult{runnum}{end+1,1},systemresult{runnum}{end+1,2}] = system(['flirt -in ' BOLDfile '_st_mcf_avg -ref ' anatomicalfolder T1file '_bet -dof 6 -cost bbr -wmseg ' anatomicalfolder T1file '_bet_wm -init ' BOLDfile '_st_mcf_avg_2T1_init.mat -omat ' BOLDfile '_st_mcf_avg_2T1.mat -schedule ${FSLDIR}/etc/flirtsch/bbr.sch']);% -out ' BOLDfile '_st_mcf_avg_T1space']);
                        delete([BOLDfile '_st_mcf_avg_2T1_init.mat']);
                        
                        %concatenate BOLD->T1 registration and T1->MNI transform
                        [systemresult{runnum}{end+1,1},systemresult{runnum}{end+1,2}] = system(['convert_xfm -omat ' BOLDfile '_st_mcf_avg_2MNI.mat -concat ' anatomicalfolder T1file '_2MNI.mat ' BOLDfile '_st_mcf_avg_2T1.mat']);
                        
                        %apply BOLD->T1 registration + T1->MNI transform
                        [systemresult{runnum}{end+1,1},systemresult{runnum}{end+1,2}] = system(['flirt -in ' BOLDfile '_st_mcf -interp spline -ref ' anatomicalfolder T1file '_MNI_' num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) ' -applyxfm -init '  BOLDfile '_st_mcf_avg_2MNI.mat -out ' BOLDfile '_st_mcf_MNI']);
                        [systemresult{runnum}{end+1,1},systemresult{runnum}{end+1,2}] = system(['flirt -in ' BOLDfile '_st_mcf_avg -interp spline -ref ' anatomicalfolder T1file '_MNI_' num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) ' -applyxfm -init '  BOLDfile '_st_mcf_avg_2MNI.mat -out ' BOLDfile '_st_mcf_avg_MNI']);
                        
                        finalBOLDfilenames{runnum} = [BOLDfile '_st_mcf_MNI.nii.gz'];
                        finalBOLDavgfilenames{runnum} = [BOLDfile '_st_mcf_avg_MNI.nii.gz'];
                        
                    end
                    
                    
                    %Intensity normalization to mode within-brain value of 1000
                    data = load_untouch_nii(finalBOLDfilenames{runnum});
                    data_inmask = data.img .* repmat(brainmask.img,[1,1,1,size(data.img,4)]) .* (data.img > 100);
                    data_inmask = data_inmask(data_inmask>0);
                    [counts,edges] = histcounts(data_inmask,1000);
                    [~,maxind] = max(counts);
                    modeval = mean([edges(maxind) edges(maxind+1)]);
                    data.img = data.img .* (1000 ./ modeval);
                    save_untouch_nii(data,finalBOLDfilenames{runnum});
                    data = []; data_inmask = [];
                    
                    %make pngs for visual check of registration quality
                    if exist([fsLRfolder '/MNI/Native/' subject '.L.pial.native.surf.gii'],'file')
                        Batch_wb_image_capture_volreg(finalBOLDavgfilenames{runnum},[fsLRfolder '/MNI/Native/' subject '.L.pial.native.surf.gii'],[fsLRfolder '/MNI/Native/' subject '.L.white.native.surf.gii'],[fsLRfolder '/MNI/Native/' subject '.R.pial.native.surf.gii'],[fsLRfolder '/MNI/Native/' subject '.R.white.native.surf.gii'],[finalBOLDavgfilenames{runnum}(1:end-7) '_registration'])
                    end
                end
                
                %track session number and show registration pngs
                for runnum = 1:length(BOLDrunnames)
                    dlmwrite(sess_tracker,[finalBOLDfilenames{runnum} ' ' num2str(sesscounter(runnum)) ' ' BOLDrunnames{runnum} '_mcf.par'],'-append','delimiter','');%write the data to the output file
                    figure; imshow([finalBOLDavgfilenames{runnum}(1:end-7) '_registration.png'])
                    titlestring = [finalBOLDavgfilenames{runnum}(1:end-7) ' registration'];
                    titlestring(strfind(titlestring,'_')) = ' ';
                    title(gca,titlestring)
                end
                
                try delete(gcp,'nocreate'); catch; end
            end
            
            brainmask = [];
            
            allsystemresult = cell(0,2);
            for runnum = 1:length(BOLDrunnames)
                allsystemresult = [allsystemresult ; systemresult{runnum}];
            end
            
            if any(cell2mat(allsystemresult(:,1)))
                error(['Subject ' subject ': Problem detected with BOLD preprocessing'])
            end
            
            
            
            cprintf('*blue',['Subject ' subject ' COMPLETED functional preprocessing without error!\n'])
            dlmwrite(logfile,'all processing complete.','-append','delimiter','')
            
        catch errorinfo
            
            dlmwrite(logfile,'ERROR:','-append','delimiter','')
            
            errormessage = [errorinfo.stack(end).file ', line ' num2str(errorinfo.stack(end).line) ': ' errorinfo.message];
            
            allsystemresult = cell(0,2);
            for runnum = 1:length(BOLDrunnames)
                allsystemresult = [allsystemresult ; systemresult{runnum}];
            end
            
            probleminds = find(cell2mat(allsystemresult(:,1)));
            for i = 1:length(probleminds)
                dlmwrite(logfile,allsystemresult{probleminds(i),2},'-append','delimiter','');
            end
            dlmwrite(logfile,errormessage,'-append','delimiter','');
            
            cprintf('err',['Subject ' subject ' FAILED to complete functional preprocessing! Check error log ' logfile '\n'])
            
        end
        
    end
    
end




%% POSTPROCESSING

if any([T1_hippocampal_subfields BOLD_fcprocess cifti_creation cifti_correlation cifti_distances])
    
    enoughtime = true(length(subjects),length(functional_sequences));
    
    disp('POSTPROCESSING')
    
    for subnum = 1:subcount
        
        systemresult = cell(0,2);
        
        subject = subjects{subnum};
        
        warning off
        
        try
            
            %% Set up folders
            
            
            basefolder = [subsdirectory '/' subject '/'];
            DTIfolder = [basefolder '/DTI/'];
            anatomicalfolder = [basefolder '/anatomical/'];
            fieldmapfolder = [basefolder '/fieldmap/'];
            functionalfolder = [basefolder '/functional/'];
            freesurferfolder = [basefolder '/freesurfer/' subject '/'];
            fsLRfolder = [basefolder '/fs_LR/'];
            fc_processfolder = [basefolder '/fc_processed/'];
            ciftifolder = [basefolder '/cifti/'];
            correlationfolder = [basefolder '/connectome/'];
            logfolder = [basefolder '/processing_logs/'];
            
            logfile2 = [logfolder datestr(now) '.txt']; tokens = tokenize(logfile2,' '); logfile2 = [tokens{1} '_' tokens{2}];
            
            %% Hippocampal subfields
            
            if T1_hippocampal_subfields
                
                disp(['Subject ' subject ': segmenting hippocampal subfields']);
                dlmwrite(logfile2,'segmenting hippocampal subfields','-append','delimiter','')
                
                %move files around
                [~,~] = system(['rm ' freesurferfolder '/scripts/IsRunning.lh+rh'])
                %                 mkdir([freesurferfolder '/../temp/']);
                %                 [systemresult{end+1,1},systemresult{end+1,2}] = system(['mv ' freesurferfolder '/* ' freesurferfolder '/../temp/']);
                %                 mkdir([freesurferfolder '/' subject]);
                %                 [systemresult{end+1,1},systemresult{end+1,2}] = system(['mv ' freesurferfolder '/../temp/* ' freesurferfolder '/' subject '/']);
                %                 rmdir([freesurferfolder '/../temp/'])
                
                %segment hippocampal subfields from T1
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['setenv FREESURFER_HOME /usr/local/freesurfer6/freesurfer;'...
                    'source $FREESURFER_HOME/SetUpFreeSurfer.csh;'...
                    'recon-all -subjid ' subject ' -sd ' freesurferfolder ' -hippocampal-subfields-T1']);
                
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['mri_convert -rl ' freesurferfolder '/' subject '/mri/rawavg.mgz -rt nearest ' freesurferfolder '/' subject '/mri/lh.hippoSfLabels-T1.v10.FSvoxelSpace.mgz ' freesurfer6folder '/' subject '/mri/lh.hippoSfLabels-T1.v10.T1voxelSpace.nii.gz']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['mri_convert -rl ' freesurferfolder '/' subject '/mri/rawavg.mgz -rt nearest ' freesurferfolder '/' subject '/mri/rh.hippoSfLabels-T1.v10.FSvoxelSpace.mgz ' freesurfer6folder '/' subject '/mri/rh.hippoSfLabels-T1.v10.T1voxelSpace.nii.gz']);
                
                %move files around
                %                 [systemresult{end+1,1},systemresult{end+1,2}] = system(['mv ' freesurferfolder '/' subject '/* ' freesurferfolder]);
                %                 rmdir([freesurferfolder '/' subject ])
                
                if any(cell2mat(systemresult(:,1)))
                    error(['Subject ' subject ': Problem detected with hippocampal subfield segmentation'])
                end
                
            end
            
            
            %% Fc-processing
            
            
            if BOLD_fcprocess
                
                mkdir(fc_processfolder)
                
                dlmwrite(logfile2,'fc-processing BOLD data','-append','delimiter','')
                
                for seq = 1:length(functional_sequences)
                    if fcprocess_sequences{seq}
                        
                        disp(['Subject ' subject ', ' functional_sequences{seq} ': fc-processing']);
                        FC_Process_filt_Avi_better(subject,[functionalfolder 'BOLDruns_sessions_' functional_sequences{seq} '.txt'],FDthresh,functional_sequence_TRs{seq},lowpassmotionfilt,fc_processfolder,functional_sequences{seq},freesurferfolder,functional_normalized_voxdim)
                        make_grayplots(basefolder,functional_sequences{seq},functional_normalized_voxdim)
                    end
                end
                
            end
            
            
            if any([cifti_creation cifti_correlation])
                for seq = 1:length(functional_sequences)
                    tmask = load([fc_processfolder '/' functional_sequences{seq} '_all_tmask.txt']);
                    if (nnz(tmask) .* functional_sequence_TRs{seq} ./ 60) < mintime
                        enoughtime(subnum,seq) = false;
                        disp(['Subject ' subject ' does not have at least ' num2str(mintime) ' minutes of ' functional_sequences{seq} ' data! Processing will not continue.'])
                    end
                end
            end
            
            
            
            %% Cifti creation
            
            if cifti_creation
                
                mkdir(ciftifolder)
                
                dlmwrite(logfile2,'creating ciftis','-append','delimiter','')
                
                medial_masks = {'/home/data/scripts/Resources/cifti_masks/L.atlasroi.32k_fs_LR.shape.gii','/home/data/scripts/Resources/cifti_masks/R.atlasroi.32k_fs_LR.shape.gii'};
                    
                if make_smallwall_ciftis
                    medial_masks = [medial_masks ; {'/home/data/scripts/Resources/cifti_masks/L.atlasroi_erode3.32k_fs_LR.shape.gii','/home/data/scripts/Resources/cifti_masks/R.atlasroi_erode3.32k_fs_LR.shape.gii'}];
                end
                
                
                for seq = 1:length(functional_sequences)
                    
                    if enoughtime(subnum,seq)
                        
                        disp(['Subject ' subject ', ' functional_sequences{seq} ': mapping functional data to cortical surface']);
                        
                        if fcprocess_sequences{seq}
                            
                            systemresult = post_fc_processing_batch_singlesub(subject,[fc_processfolder '/' functional_sequences{seq} '_fc_processed_tmasked.nii.gz'],ciftifolder,[fsLRfolder '/MNI/'],functional_sequence_TRs{seq},functional_sequences{seq},[fc_processfolder '/' functional_sequences{seq} '_all_tmask.txt'],[functionalfolder '/BOLDruns_sessions_' functional_sequences{seq} '.txt'],[fc_processfolder '/' functional_sequences{seq} '_runs_sessions.txt'],medial_masks,geodesic_smooth,functional_normalized_voxdim,systemresult);
                            
                        else
                            
                            BOLDmergedfile = [functionalfolder '/' functional_sequences{seq} '_merged'];
                            
                            fslmergestr = ['fslmerge -t ' BOLDmergedfile ' '];
                            
                            [BOLDrunnames,~,~] = textread([functionalfolder '/BOLDruns_sessions_' functional_sequences{seq} '.txt'],'%s%s%s');
                            
                            for runnum = 1:length(BOLDrunnames)
                                [~,BOLDfilename,~] = fileparts(BOLDrunnames{runnum});
                                if strcmp(BOLDfilename(1:length(functional_sequences{seq})),functional_sequences{seq});
                                    fslmergestr = [fslmergestr BOLDrunnames{runnum} '_st_mcf_MNI.nii.gz '];
                                end
                            end
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(fslmergestr);
                            
                            systemresult = post_fc_processing_batch_singlesub(subject,BOLDmergedfile,ciftifolder,[fsLRfolder '/MNI/'],functional_sequence_TRs{seq},functional_sequences{seq},[],[],[],medial_masks,geodesic_smooth,systemresult);
                            
                            delete(BOLDmergedfile)
                            
                        end
                        
                        try delete([ciftifolder '/surf_timecourses/*.nii'])
                            rmdir([ciftifolder '/surf_timecourses/'],'s')
                        catch; end
                        
                        if systemresult{end,1}
                            error('Problem detected with Cifti creation')
                        end
                    end
                end
                
                
                
            end
            
            
            
            
            %% Cifti correlation
            
            if cifti_correlation
                
                mkdir(correlationfolder)
                dlmwrite(logfile2,'correlating ciftis','-append','delimiter','')
                
                
                for seq = 1:length(functional_sequences)
                    
                    if enoughtime(subnum,seq)
                        
                        disp(['Subject ' subject ', ' functional_sequences{seq} ': correlating all timecourses']);
                        
                        ciftifile = [ciftifolder '/cifti_timeseries_normalwall/' functional_sequences{seq} '_LR_surf_subcort_' num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) num2str(functional_normalized_voxdim) '_32k_fsLR_smooth' num2str(geodesic_smooth) '.dtseries.nii'];
                        
                        data = ft_read_cifti_mod(ciftifile);
                        
                        corr = paircorr_mod(data.data');
                        corr = FisherTransform(corr);
                        corr(isnan(corr)) = 0;
                        
                        data.dimord = 'pos_pos';
                        data.data = corr;
                        clear corr
                        
                        ft_write_cifti_mod([correlationfolder '/' functional_sequences{seq} '_corr.dconn.nii'],data)
                        
                        data.data = [];
                        
                    end
                    
                end
                
            end
            
            
            
            %% Cifti Distances  NOTE: DOESN'T WORK NOW
            
            if cifti_distances
                
                dlmwrite(logfile2,'conducting surface area and geodesic distance calculation','-append','delimiter','')
                disp(['Subject ' subject ': calculating surface area and getting point-to-point geodesic and euclidean distances'])
                
                surfacetousefolder = [fsLRfolder '/MNI/fsaverage_LR32k/'];
                distfolder = [surfacetousefolder '/distances/'];
                mkdir(distfolder)
                
                surffileL = [surfacetousefolder '/' subject '.L.midthickness.32k_fs_LR.surf.gii'];
                surfaceareafileL = [surfacetousefolder '/' subject '.L.midthickness.32k_fs_LR_surfaceareas.func.gii'];
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['wb_command -surface-vertex-areas ' surffileL ' ' surfaceareafileL]);
                
                surffileR = [surfacetousefolder '/' subject '.R.midthickness.32k_fs_LR.surf.gii'];
                surfaceareafileR = [surfacetousefolder '/' subject '.R.midthickness.32k_fs_LR_surfaceareas.func.gii'];
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['wb_command -surface-vertex-areas ' surffileR ' ' surfaceareafileR]);
                
                
                if logical(parallelize)
                    nworkers = maxworkers;
                    if isempty(gcp('nocreate'))
                        processingpool = parpool(maxworkers);
                    end
                else
                    nworkers = 0;
                end
                
                temp = gifti(surffileL);
                nverts = size(temp.vertices,1);
                
                distances = zeros(nverts,nverts,'uint8');
                distancesLfile = [distfolder '/Surface_distances_L.mat'];
                parfor (v = [1:nverts],nworkers)
                    [~,~] = system(['wb_command -surface-geodesic-distance ' surffileL ' ' num2str(v-1) ' ' distfolder 'temp' num2str(v) '.func.gii -limit 255']);
                    temp = gifti([distfolder 'temp' num2str(v) '.func.gii']);
                    delete([distfolder 'temp' num2str(v) '.func.gii']);
                    temp.cdata(temp.cdata==-1) = 255;
                    distances(:,v) = uint8(temp.cdata);
                end
                save(distancesLfile,'distances','-v7.3')
                
                
                temp = gifti(surffileR);
                nverts = size(temp.vertices,1);
                
                distances = zeros(nverts,nverts,'uint8');
                distancesRfile = [distfolder '/Surface_distances_R.mat'];
                parfor (v = [1:nverts],nworkers)
                    [~,~] = system(['wb_command -surface-geodesic-distance ' surffileR ' ' num2str(v-1) ' ' distfolder 'temp' num2str(v) '.func.gii -limit 255']);
                    temp = gifti([distfolder 'temp' num2str(v) '.func.gii']);
                    delete([distfolder 'temp' num2str(v) '.func.gii']);
                    temp.cdata(temp.cdata==-1) = 255;
                    distances(:,v) = uint8(temp.cdata);
                end
                save(distancesRfile,'distances','-v7.3')
                
                
                
                
                
                ciftifiles = dir([ciftifolder '/cifti_timeseries_normalwall/*.dtseries.nii']);
                ciftifile = [ciftifolder '/cifti_timeseries_normalwall/' ciftifiles(1).name];
                mkdir([ciftifolder '/distances/'])
                outputfile = [ciftifolder '/distances/normalwall_distmat_surf_geodesic_vol_euclidean_xhemlarge_uint8.mat'];
                Make_distmat_32k_fsLR_singlesub(distancesLfile,distancesRfile,surffileL,surffileR,ciftifile,outputfile,true,true)
                
            end
            
            
            
            %%
            
            
            cprintf('*blue',['Subject ' subject ' COMPLETED post-processing without error!\n'])
            dlmwrite(logfile2,'all processing complete.','-append','delimiter','')
            
        catch errorinfo
            
            dlmwrite(logfile2,'ERROR:','-append','delimiter','')
            
            errormessage = [errorinfo.stack(1).file ', line ' num2str(errorinfo.stack(1).line) ': ' errorinfo.message];
            
            probleminds = find(cell2mat(systemresult(:,1)));
            for i = 1:length(probleminds)
                dlmwrite(logfile2,systemresult{probleminds(i),2},'-append','delimiter','');
            end
            dlmwrite(logfile2,errormessage,'-append','delimiter','');
            
            cprintf('err',['Subject ' subject ' FAILED to complete post-processing! Check error log ' logfile2 '\n'])
            
        end
        
        
    end
    
end

