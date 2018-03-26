%% Subjects

%subjects = '/home/data/EEG/processed/Robi/mri/MAV_final_list.txt';
%subjects = '/home/data/EEG/processed/Robi/mri/processing_list_10min.txt';
%subjects = '/home/data/EEG/processed/Robi/mri/processing_list_092916.txt';
%subjects = '/home/data/EEG/processed/Robi/mri/list_forJDP.txt';
%subjects = {'MAV006','MAV024','ROBI006','ROBI004'};
subjects = {'ROBI003'};
%{'MAV029','MAV035'};
%
%{'MAV004','MAV014','MAV033','MAV036','MAV020'};
%'MAV029',

%% Steps to run

choose_T1 = 0;
copy_rawBOLDdata = 1;
T1_preprocess = 0;
T1_segment = 0;
surface_registration = 0;
T2_preprocess = 0;
DTI_preprocess = 0;
apply_fieldmap = 0;
BOLD_preprocess = 1;
BOLD_fcprocess = 1;
cifti_creation = 1;
cifti_correlation = 0;
cifti_distances = 0;


%% Parameters

parallelize = 0;
  maxworkers = 6;

template = '/usr/local/fsl/data/standard/MNI152lin_T1_1mm_brain.nii.gz';

T1string = 'T1';
T2string = 'T2';
DTIstring = 'xDTI';
fieldmapstring = 'PhaseMap';

normalized_voxdim = 3;
DTI_normalized_voxdim = 2;

functional_sequences = {'RSFC'};%{'restingstate'}%{'Movie'};%

functional_sequence_TRs = {3.0};% 2.2

functional_sequence_rawsize = {[80 80 34]};%[64 64 34]

fcprocess_sequences = {true}; %true

sliceorderinfo = {''};%{'--ocustom=/home/data/scripts/Resources/sliceorder_interleavedPhillips.txt'};%'--odd',%

skipframes = 0;

fieldmap_deltaTE = 2.0;
echospacing = .00059;

FDthresh = .04;

mintime = 0; % minutes

geodesic_smooth = 2.55;%4.25;

QCtrackingfile = '/home/data/EEG/processed/Robi/mri/QC/MRI_scan_tracker.xlsx';


if length(functional_sequence_TRs)==1
    functional_sequence_TRs = repmat(functional_sequence_TRs,1,length(functional_sequences));
end
if length(functional_sequence_rawsize)==1
    functional_sequence_rawsize = repmat(functional_sequence_rawsize,1,length(functional_sequences));
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

if(exist(QCtrackingfile, 'file'))
[~,~,QCdata]=xlsread(QCtrackingfile);
else
  QCdata = cell(0);
end

%% Choose best T1 image

if choose_T1
    for subnum = 1:subcount
        subject = subjects{subnum};
        
        rawfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/raw/'];
        preprocfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/preprocessed/']; mkdir(preprocfolder)
        
        sessions = dir([rawfolder '/']); 
        folderindex = struct2cell(sessions); folderindex = logical(cell2mat(folderindex(4,:)));
        sessions = sessions(folderindex);
        sessions = sessions(3:end);
        
        all_T1files = cell(0,1);
        firstT1 = true;
        
        for sess = 1:length(sessions)
                
                T1files = dir([rawfolder '/' sessions(sess).name '/o*' T1string '*.nii.gz']);
                for file = 1:length(T1files)
                    filename = [rawfolder '/' sessions(sess).name '/' T1files(file).name];
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
        end
        [~,~] = system([fslviewstring ' &']);
        
        if length(all_T1files) > 1
            bestnum = input(['Subject ' subject ': input index of best image (counting from bottom of fslview list): ']);
        else
            bestnum = 1;
        end
        copyfile(all_T1files{bestnum},[preprocfolder '/T1.nii.gz']);
    end
end



warning off


if any([copy_rawBOLDdata T1_preprocess T1_segment T2_preprocess surface_registration DTI_preprocess BOLD_preprocess])

disp('PREPROCESSING')

%% Set up parallel processing

if logical(parallelize) && (subcount > 1)
    nworkers = min([subcount maxworkers]);
    processingpool = parpool(nworkers);
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
        
        
        rawfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/raw/'];
        DTIfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/DTI/']; 
        preprocfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/preprocessed/']; mkdir(preprocfolder)
        freesurferfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/freesurfer/']; 
        fsLRfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/fs_LR/']; 
        fc_processfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/fc_processed/']; 
        ciftifolder = ['/home/data/EEG/processed/Robi/mri/' subject '/cifti/']; 
        logfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/processing_logs/']; mkdir(logfolder)
        
        logfile = [logfolder datestr(now) '.txt']; tokens = tokenize(logfile,' '); logfile = [tokens{1} '_' tokens{2}];
        %Create error tracking file to record issues with processing
        fid = fopen(logfile,'at');
        fclose(fid);
       
        
        %% Find nifti data in raw folder and move to preprocessed folder
        
        
        sesscounter = [];
        
        sessions = dir([rawfolder '/']); 
        folderindex = struct2cell(sessions); folderindex = logical(cell2mat(folderindex(4,:)));
        sessions = sessions(folderindex);
        sessions = sessions(3:end);
        sessnums_unsorted = zeros(1,length(sessions));
        for sess = 1:length(sessions)
            tokens = tokenize(sessions(sess).name,'-');
            sessnums_unsorted(sess) = str2num(tokens{end});
        end
        [sessnums, sessi] = sort(sessnums_unsorted,'ascend');
        sessions = sessions(sessi);
        
        
        %sessions = dir([rawfolder]); sessions = sessions(3);
        
        for seq = 1:length(functional_sequences)
            seq_counter{subnum}{seq} = 0;
        end
        BOLDrunnames = cell(0,1);
        
        subQCindices = find(strcmp(QCdata(:,1),subject));
        subQC_sessions = cell2mat(QCdata(subQCindices,2));
        
        if copy_rawBOLDdata
            disp(['Subject ' subject ': finding and copying BOLD data']);
            dlmwrite(logfile,'copying raw BOLD data to preprocessing folder','-append','delimiter','')
        end
            
            
        for sess = 1:length(sessions)
            
            dashloc = strfind(sessions(sess).name,'-');
            sessnum = str2num(sessions(sess).name(dashloc+1:end));
            sessQCindices = subQCindices(subQC_sessions==sessnum);
            sessQC_decision = cell2mat(QCdata(sessQCindices,7));
            sessQC_acquisitionnum = cell2mat(QCdata(sessQCindices,8));
            
            %Find all valid BOLD runs
            for seq = 1:length(functional_sequences)
                
                thisseq_files = dir([rawfolder '/' sessions(sess).name '/*' functional_sequences{seq} '*.nii.gz']);
                for file = 1:length(thisseq_files)
                    filename = [rawfolder '/' sessions(sess).name '/' thisseq_files(file).name];
                    
                    extpoint = strfind(filename,'.nii.gz');
                    acquisitionnum = str2num(filename(extpoint-2:extpoint-1));
                    
                    %check if data is the right size
                    [~,sizestr] = system(['fslsize ' filename ' -s']);
                    tokens = tokenize(sizestr,' ');
                    datasize = [str2num(tokens{3}) str2num(tokens{5}) str2num(tokens{7}) str2num(tokens{9})];
                    correctsize = logical(all(datasize(1:3)==functional_sequence_rawsize{seq}));
                    
                    %check if data has many frames (i.e. wasn't aborted)
                    data_islong = (datasize(4) > 20);
                    
                    %check if data has been flagged by manual QC check (or if it's not in the QC check file)
                    if isempty(acquisitionnum)
                        flaggedasbad = false;
                    else
                        flaggedasbad = (sessQC_decision(sessQC_acquisitionnum==acquisitionnum)==0);
                        if ~any(sessQC_acquisitionnum==acquisitionnum)
                            flaggedasbad = false;
                        end
                    end
                    
                    if correctsize && data_islong && (~flaggedasbad)
                        
                        %keep track of session names, session number for each sequence, and overall session number
                        seq_counter{subnum}{seq} = seq_counter{subnum}{seq} +1;
                        BOLDrunnames(end+1) = {[preprocfolder '/' functional_sequences{seq} '_' num2str(seq_counter{subnum}{seq})]};
                        sesscounter(end+1) = sessnums(sess);
                        
                        if copy_rawBOLDdata
                            
                            %load data, trim first few frames
                            data = load_untouch_nii(filename);
                            data.img = data.img(:,:,:,skipframes+1:end);
                            data.hdr.dime.dim(5) = size(data.img,4);
                            
                            %save to preprocessing folder
                            save_untouch_nii(data,[BOLDrunnames{end} '.nii.gz'])
                        end
                        
                    end
                    
                end
            end
            
        end
            
        
        
        
        
        
        %% T1 Preprocessing
        
        
        T1file = 'T1';
        
        if T1_preprocess
            
            cd(preprocfolder)
            
            dlmwrite(logfile,'preprocessing T1','-append','delimiter','')
            
            disp(['Subject ' subject ': reorienting, cropping, bias-correcting, and skull-stripping T1, and registering to MNI']);
            
            %reorient, crop, and bias-correct T1 image
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fsl_anat -i ' T1file ' -o struct --clobber --noreg --nononlinreg --noseg --nosubcortseg']);
            
             T1file = [T1file '_biascorr'];
            
            %remove folder used for processing
            copyfile([preprocfolder '/struct.anat/T1_biascorr.nii.gz'],[preprocfolder '/' T1file '.nii.gz']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['chmod 777 -R ' preprocfolder '/struct.anat']);
            delete([preprocfolder '/struct.anat/*']);
            rmdir([preprocfolder '/struct.anat']);
            
            %brain-extract T1 image
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['bet ' T1file ' ' T1file '_bet -R -g -.1']);
            
            %register T1 to MNI
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T1file '_bet -interp spline -dof 12 -ref ' template ' -omat T1_2MNI.mat -out ' T1file '_bet_MNI']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T1file ' -interp spline -ref ' template ' -applyxfm -init T1_2MNI.mat -out ' T1file '_MNI']);
            
            if any(cell2mat(systemresult(:,1)))
                error(['Subject ' subject ': Problem detected with T1 preprocessing'])
            end
            
        else
            
            %track name
            T1file = [T1file '_biascorr'];
            
        end
        
        
         %% T1 Segmentation
        
        if T1_segment
            
            disp(['Subject ' subject ': segmenting T1']);
            dlmwrite(logfile,'segmenting T1','-append','delimiter','')
            
            try [~,~,~] = rmdir(freesurferfolder,'s'); catch; end
            
            mkdir(freesurferfolder)
            
            %copy bias-corrected T1 to freesurfer folder
            copyfile([preprocfolder '/' T1file '.nii.gz'],[freesurferfolder '/' T1file '.nii.gz'])
            cd(freesurferfolder)
            
            %transform T1 to .mgz format
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['recon-all -subjid ' subject ' -sd ' freesurferfolder ' -i ' freesurferfolder '/' T1file '.nii.gz']);
            %segment T1
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['recon-all -subjid ' subject ' -sd ' freesurferfolder ' -all']);
            
            %move files around
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['mv ' freesurferfolder '/' subject '/* ' freesurferfolder]);
            rmdir([freesurferfolder '/' subject ])
            cd(preprocfolder)
            
            if any(cell2mat(systemresult(:,1)))
                error(['Subject ' subject ': Problem detected with Freesurfer'])
            end

        end
        
        
        %% Surface-based registration
        
        if surface_registration
            
            disp(['Subject ' subject ': surface-based registration']);
            dlmwrite(logfile,'conducting surface registration','-append','delimiter','')
            
            %get segmentation masks into nifti format
            systemresult = make_fs_masks_mutualdist(freesurferfolder,normalized_voxdim,[preprocfolder '/T1_2MNI.mat'],template,systemresult);
            copyfile([freesurferfolder '/nusmask/aparc+aseg_cerebralwm.nii.gz'],preprocfolder)
            
            mkdir(fsLRfolder)
            
            prevsize = size(systemresult,1);
            systemresult = PostFreeSurferPipeline_fsavg2fslr_long_singlesub(subject,T1file,freesurferfolder,fsLRfolder,preprocfolder,systemresult);
            if any(cell2mat(systemresult(:,1)))
                error(['Subject ' subject ': Problem detected with surface registration'])
            end
            
            
        end
        
        %% T2 Preprocessing
        
        T2file = 'T2';
        
        if T2_preprocess
            
            cd(preprocfolder)
            
            dlmwrite(logfile,'preprocessing T2','-append','delimiter','')
            
            disp(['Subject ' subject ': reorienting, cropping, bias-correcting, and skull stripping T2, and registering to T1 and through to MNI']);
            
            %reorient, crop, and bias-correct T2 image
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fsl_anat -i ' T2file ' -o struct -t T2 --clobber --noreg --nononlinreg --noseg --nosubcortseg']);
            copyfile([preprocfolder '/struct.anat/T2_biascorr.nii.gz'],[preprocfolder '/' T2file '_biascorr.nii.gz'])
            T2file = [T2file '_biascorr'];
            
            delete([preprocfolder '/struct.anat/*']);
            rmdir([preprocfolder '/struct.anat']);
            
            %brain-extract T2 image
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['bet ' T2file ' ' T2file '_bet -R -g -.1']);
            
            %register T2 to T1, and concatenate with T1 to MNI
            %OPTION 1
            %[status(end+1) result{end+1}] = system(['flirt -in ' T2file ' -ref ' T1file ' -interp spline -dof 6 -omat T2_2T1.mat -out T2_T1space']);
            
            
            %OPTION 2
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T2file ' -ref ' T1file ' -dof 6 -omat T2_2T1_init.mat']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T2file ' -ref ' T1file ' -interp spline -dof 6 -cost bbr -wmseg aparc+aseg_cerebralwm -init T2_2T1_init.mat -omat T2_2T1.mat -out T2_T1space -schedule ${FSLDIR}/etc/flirtsch/bbr.sch']);
            
            [systemresult{end+1,1},systemresult{end+1,2}] = system('convert_xfm -omat T1_2T2.mat -inverse T2_2T1.mat');
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in aparc+aseg_cerebralwm -ref ' T2file ' -applyxfm -init T1_2T2.mat -interp nearestneighbour -out aparc+aseg_cerebralwm_T2space']);
            
            [systemresult{end+1,1},systemresult{end+1,2}] = system('convert_xfm -omat T2_2MNI.mat -concat T1_2MNI.mat T2_2T1.mat');
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T2file ' -ref ' template ' -applyxfm -init T2_2MNI.mat -out ' T2file '_MNI']);
            
        else
            
            %track T2 name
            T2file = [T2file '_biascorr'];
            
        end
        
        
        %% DTI preprocessing
        if DTI_preprocess
            mkdir(DTIfolder)
            cd(DTIfolder)
            
            dlmwrite(logfile,'preprocessing DTI data','-append','delimiter','')
            disp(['Subject ' subject ': copying diffusion data, fitting tensors, registering to MNI, and getting FA near the cortex']);
            
%             sessions = dir([rawfolder '/']);
%             folderindex = struct2cell(sessions); folderindex = logical(cell2mat(folderindex(4,:)));
%             sessions = sessions(folderindex);
%             sessions = sessions(3:end);
            
            all_T1files = cell(0,1);
            firstT1 = true;
            
            DTIfiles = cell(0,1);
            DTIcounter = 0;
            
            for sess = 1:length(sessions)
                DTIfiles_thissess = dir([rawfolder '/' sessions(sess).name '/' DTIstring '*.nii.gz']);
                for i = 1:length(DTIfiles_thissess)
                    DTIcounter = DTIcounter+1;
                    DTIfiles{DTIcounter} = [DTIfolder '/DTI_' num2str(DTIcounter) '.nii.gz'];
                    copyfile([rawfolder '/' sessions(sess).name '/' DTIfiles_thissess(i).name],DTIfiles{DTIcounter});
                    copyfile([rawfolder '/' sessions(sess).name '/' DTIfiles_thissess(i).name(1:end-7) '.bval'],[DTIfiles{DTIcounter}(1:end-7) '.bval']);
                    copyfile([rawfolder '/' sessions(sess).name '/' DTIfiles_thissess(i).name(1:end-7) '.bvec'],[DTIfiles{DTIcounter}(1:end-7) '.bvec']);
                end
            end
            
            
            for filenum = 1:DTIcounter
                
                DTIfullfile = DTIfiles{filenum}(1:end-7);
                [DTIpath, DTIfile, DTIext] = fileparts(DTIfullfile);
                
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['eddy_correct ' DTIfile ' ' DTIfile '_ec 0']);
                
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslselectvols -i ' DTIfile '_ec -o ' DTIfile '_ec_b0 --vols=0']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' DTIfile '_ec_b0 -ref ' preprocfolder '/' T1file '_bet -dof 12 -omat ' DTIfile '_2T1_init.mat']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' DTIfile '_ec_b0 -ref '  preprocfolder '/' T1file '_bet -dof 12 -cost bbr -wmseg ' preprocfolder '/aparc+aseg_cerebralwm -init ' DTIfile '_2T1_init.mat -omat ' DTIfile '_2T1.mat -schedule ${FSLDIR}/etc/flirtsch/bbr.sch']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['convert_xfm -omat ' DTIfile '_2MNI.mat -concat ' preprocfolder '/T1_2MNI.mat ' DTIfile '_2T1.mat']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['convert_xfm -omat T1_2' DTIfile '.mat -inverse ' DTIfile '_2T1.mat']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' preprocfolder '/' T1file '_bet -interp spline -ref ' DTIfile '_ec_b0 -applyxfm -init T1_2' DTIfile '.mat -out ' T1file '_bet_' DTIfile]);
                %[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' DTIfile '_ec -interp spline -ref ' template ' -applyisoxfm ' num2str(normalized_voxdim) ' -init ' DTIfile '_2MNI.mat -out ' DTIfile '_ec_MNI']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' T1file '_bet_' DTIfile ' -thr 10 -bin ' T1file '_bet_' DTIfile '_mask']);
                
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['dtifit -k ' DTIfile '_ec -m ' T1file '_bet_' DTIfile '_mask -r ' DTIfile '.bvec -b ' DTIfile '.bval -o ' DTIfile '_ec']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' DTIfile '_ec_FA -interp spline -ref ' template ' -applyisoxfm ' num2str(DTI_normalized_voxdim) ' -init ' DTIfile '_2MNI.mat -out ' DTIfile '_ec_FA_MNI']);
                
                FA = load_untouch_nii([DTIfile '_ec_FA_MNI.nii.gz']);
                if filenum == 1
                    FA_avg = FA;
                    FA_avg.img = FA.img / DTIcounter;
                else
                    FA_avg.img = FA_avg.img + (FA.img / DTIcounter);
                end
                
                
            end
            
            save_untouch_nii(FA_avg,'DTI_avg_ec_FA_MNI.nii.gz')
            
            systemresult = make_fs_masks_mutualdist(freesurferfolder,DTI_normalized_voxdim,[preprocfolder '/T1_2MNI.mat'],template,systemresult);
            copyfile([freesurferfolder '/nusmask/aparc+aseg_cerebralwm_mask_' num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim) num2str(DTI_normalized_voxdim) '.nii.gz'],DTIfolder);
            
            %Map FA near cortex to surface
            distances = 2:8; %mm
            FA_near_surface(subject,distances)
                
            if any(cell2mat(systemresult(:,1)))
                error(['Subject ' subject ': Problem detected with DTI processing'])
            end
            
        end
        
        %% Fieldmap preparation
        if apply_fieldmap
            cd(preprocfolder)
            dlmwrite(logfile,'calculating fieldmaps','-append','delimiter','')
            disp(['Subject ' subject ': copying field map data, calculating field maps, and averaging within scanning sessions']);
            
            %Track field map files
            FMfiles = cell(0,1);
            FMcounter = 0;
            
            %prepare commands to concatenate all field maps
            allfieldmap_mergestring = ['fslmerge -t fieldmap_rads_all '];
            allfieldmapmask_mergestring = ['fslmerge -t fieldmap_mag_bet_mask_all '];
            allfieldmapmag_mergestring = ['fslmerge -t fieldmap_mag_all '];
            
            for sess = 1:length(sessions)
                
                %find field maps in this session
                FMfiles_thissess = dir([rawfolder '/' sessions(sess).name '/' fieldmapstring '*1.nii.gz']);
                
                %If there are field map files
                if ~isempty(FMfiles_thissess)
                    
                    %prepare commands to concatenate field maps in this session
                    fieldmap_mergestring = ['fslmerge -t fieldmap_rads_all_sess' num2str(sessnums(sess)) ' '];
                    fieldmapmask_mergestring = ['fslmerge -t fieldmap_mag_bet_mask_all_sess' num2str(sessnums(sess)) ' '];
                    fieldmapmag_mergestring = ['fslmerge -t fieldmap_mag_all_sess' num2str(sessnums(sess)) ' '];
                    
                    %loop through field maps
                    for i = 1:length(FMfiles_thissess)
                        FMcounter = FMcounter+1;
                        
                        %Define the types of data
                        mag1 = FMfiles_thissess(i).name;
                        phase1 = FMfiles_thissess(i).name; phase1(end-7) = '3';
                        phase2 = FMfiles_thissess(i).name; phase2(end-7) = '4';
                        
                        %Save the magnitude image
                        mag1data = load_untouch_nii([rawfolder '/' sessions(sess).name '/' mag1]);
                        mag1data.hdr.dime.scl_slope = 1;
                        mag1data.hdr.dime.scl_inter = 0;
                        save_untouch_nii(mag1data,[preprocfolder '/fieldmap_' num2str(FMcounter) '_mag.nii.gz']);
                        
                        %Save the phase 1 image
                        phase1data = load_untouch_nii([rawfolder '/' sessions(sess).name '/' phase1]);
                        phase1data.hdr.dime.scl_slope = 1;
                        phase1data.hdr.dime.scl_inter = 0;
                        save_untouch_nii(phase1data,[preprocfolder '/fieldmap_' num2str(FMcounter) '_phase1.nii.gz']);
                        
                        %Save the phase 2 image
                        phase2data = load_untouch_nii([rawfolder '/' sessions(sess).name '/' phase2]);
                        phase2data.hdr.dime.scl_slope = 1;
                        phase2data.hdr.dime.scl_inter = 0;
                        save_untouch_nii(phase2data,[preprocfolder '/fieldmap_' num2str(FMcounter) '_phase2.nii.gz']);
                        
                        %Create the phase difference image
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' preprocfolder '/fieldmap_' num2str(FMcounter) '_phase2.nii.gz -sub ' preprocfolder '/fieldmap_' num2str(FMcounter) '_phase1.nii.gz ' preprocfolder '/fieldmap_' num2str(FMcounter) '_phasediff.nii.gz']);
                        
                        %Brain extract the magnitude image
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['bet fieldmap_' num2str(FMcounter) '_mag fieldmap_' num2str(FMcounter) '_mag_bet -R']);
                        
                        %Delete data on edge of magnitude image in
                        %preparation for erosion
                        betimage = load_untouch_nii(['fieldmap_' num2str(FMcounter) '_mag_bet.nii.gz']);
                        betimage.img(:,:,1) = 0; betimage.img(:,:,end) = 0;
                        save_untouch_nii(betimage,['fieldmap_' num2str(FMcounter) '_mag_bet.nii.gz']);
                        
                        %Erode magnitude image
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths fieldmap_' num2str(FMcounter) '_mag_bet -ero fieldmap_' num2str(FMcounter) '_mag_bet']);
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths fieldmap_' num2str(FMcounter) '_mag_bet -ero fieldmap_' num2str(FMcounter) '_mag_bet']);
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths fieldmap_' num2str(FMcounter) '_mag_bet -ero fieldmap_' num2str(FMcounter) '_mag_bet']);
                        
                        %Make a mask from eroded magnitude image
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths fieldmap_' num2str(FMcounter) '_mag_bet -bin fieldmap_' num2str(FMcounter) '_mag_bet_mask']);
                        
                        %Calculate fieldmap (rad/s) image
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['fsl_prepare_fieldmap SIEMENS fieldmap_' num2str(FMcounter) '_phasediff fieldmap_' num2str(FMcounter) '_mag_bet fieldmap_' num2str(FMcounter) '_rads ' num2str(fieldmap_deltaTE)]);
                        
                        %Prepare to concatenate images within session
                        fieldmap_mergestring = [fieldmap_mergestring 'fieldmap_' num2str(FMcounter) '_rads '];
                        fieldmapmask_mergestring = [fieldmapmask_mergestring 'fieldmap_' num2str(FMcounter) '_mag_bet_mask '];
                        fieldmapmag_mergestring = [fieldmapmag_mergestring 'fieldmap_' num2str(FMcounter) '_mag '];
                        
                        %Prepare to concatenate all images
                        allfieldmap_mergestring = [allfieldmap_mergestring 'fieldmap_' num2str(FMcounter) '_rads '];
                        allfieldmapmask_mergestring = [allfieldmapmask_mergestring 'fieldmap_' num2str(FMcounter) '_mag_bet_mask '];
                        allfieldmapmag_mergestring = [allfieldmapmag_mergestring 'fieldmap_' num2str(FMcounter) '_mag '];
                        
                        %Delete unneeded images
                        delete([preprocfolder '/fieldmap_' num2str(FMcounter) '_phase1.nii.gz']);
                        delete([preprocfolder '/fieldmap_' num2str(FMcounter) '_phase2.nii.gz']);
                        delete([preprocfolder '/fieldmap_' num2str(FMcounter) '_mag_bet.nii.gz']);
                        
                    end
                    
                    %Concatenate magnitude images
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(fieldmapmag_mergestring);
                    %Register magnitude images to each other
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['mcflirt -in fieldmap_mag_all_sess' num2str(sessnums(sess)) ' -refvol 0 -mats']);
                    %Average magnitude images
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths fieldmap_mag_all_sess' num2str(sessnums(sess)) '_mcf -Tmean fieldmap_mag_mean_sess' num2str(sessnums(sess))]);
                    %Brain-extract the average magnitude image
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['bet fieldmap_mag_mean_sess' num2str(sessnums(sess)) ' fieldmap_mag_mean_sess' num2str(sessnums(sess)) '_bet']);
                    
                    %Concatenate mask images
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(fieldmapmask_mergestring);
                    %Apply previously calculated registration
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['applyxfm4D fieldmap_mag_bet_mask_all_sess' num2str(sessnums(sess)) ' fieldmap_mag_bet_mask_all_sess' num2str(sessnums(sess)) ' fieldmap_mag_bet_mask_all_sess' num2str(sessnums(sess)) '_mcf fieldmap_mag_all_sess' num2str(sessnums(sess)) '_mcf.mat -userprefix MAT_']);
                    %Select voxels in final mask
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths fieldmap_mag_bet_mask_all_sess' num2str(sessnums(sess)) '_mcf -Tmean -thr .5 -bin fieldmap_mag_bet_mask_union_sess' num2str(sessnums(sess))]);
                    
                    %Concatenate fieldmap images
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(fieldmap_mergestring);
                    %Apply previously calculated registration
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['applyxfm4D fieldmap_rads_all_sess' num2str(sessnums(sess)) ' fieldmap_rads_all_sess' num2str(sessnums(sess)) ' fieldmap_rads_all_sess' num2str(sessnums(sess)) '_mcf fieldmap_mag_all_sess' num2str(sessnums(sess)) '_mcf.mat -userprefix MAT_']);
                    %Average fieldmap images
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths fieldmap_rads_all_sess' num2str(sessnums(sess)) '_mcf -Tmean fieldmap_rads_mean_sess' num2str(sessnums(sess))]);
                    
                    %Mask averaged fieldmap image
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths fieldmap_rads_mean_sess' num2str(sessnums(sess)) ' -mul fieldmap_mag_bet_mask_union_sess' num2str(sessnums(sess)) ' fieldmap_rads_mean_sess' num2str(sessnums(sess))]);
                    
                    %"Unmask" this image
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['fugue --loadfmap=fieldmap_rads_mean_sess' num2str(sessnums(sess)) ' --mask=fieldmap_mag_bet_mask_union_sess' num2str(sessnums(sess)) ' --unmaskfmap --savefmap=fieldmap_rads_mean_sess' num2str(sessnums(sess))]);
                    
                    %Remove unneeded files
                    delete([preprocfolder '/fieldmap_mag_all_sess' num2str(sessnums(sess)) '.nii.gz']);
                    delete([preprocfolder '/fieldmap_mag_all_sess' num2str(sessnums(sess)) '_mcf.nii.gz']);
                    delete([preprocfolder '/fieldmap_mag_bet_mask_all_sess' num2str(sessnums(sess)) '.nii.gz']);
                    delete([preprocfolder '/fieldmap_mag_bet_mask_all_sess' num2str(sessnums(sess)) '_mcf.nii.gz']);
                    delete([preprocfolder '/fieldmap_rads_all_sess' num2str(sessnums(sess)) '.nii.gz']);
                    delete([preprocfolder '/fieldmap_rads_all_sess' num2str(sessnums(sess)) '_mcf.nii.gz']);
                    delete(['fieldmap_mag_all_sess' num2str(sessnums(sess)) '_mcf.mat/*'])
                    rmdir(['fieldmap_mag_all_sess' num2str(sessnums(sess)) '_mcf.mat'],'s')
                end
                
            end
            
            %Concatenate magnitude images
            [systemresult{end+1,1},systemresult{end+1,2}] = system(allfieldmapmag_mergestring);
            %Register magnitude images to each other
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['mcflirt -in fieldmap_mag_all -refvol 0 -mats']);
            %Average magnitude images
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths fieldmap_mag_all_mcf -Tmean fieldmap_mag_mean']);
            %Brain-extract the average magnitude image
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['bet fieldmap_mag_mean fieldmap_mag_mean_bet']);
            
            %Concatenate mask images
            [systemresult{end+1,1},systemresult{end+1,2}] = system(allfieldmapmask_mergestring);
            %Apply previously calculated registration
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['applyxfm4D fieldmap_mag_bet_mask_all fieldmap_mag_bet_mask_all fieldmap_mag_bet_mask_all_mcf fieldmap_mag_all_mcf.mat -userprefix MAT_']);
            %Select voxels in final mask
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths fieldmap_mag_bet_mask_all_mcf -Tmean -thr .5 -bin fieldmap_mag_bet_mask_union']);
            
            %Concatenate fieldmap images
            [systemresult{end+1,1},systemresult{end+1,2}] = system(allfieldmap_mergestring);
            %Apply previously calculated registration
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['applyxfm4D fieldmap_rads_all fieldmap_rads_all fieldmap_rads_all_mcf fieldmap_mag_all_mcf.mat -userprefix MAT_']);
            %Average fieldmap images
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths fieldmap_rads_all_mcf -Tmean fieldmap_rads_mean']);
            
             %Mask averaged fieldmap image
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths fieldmap_rads_mean -mul fieldmap_mag_bet_mask_union fieldmap_rads_mean']);
            
            %"Unmask" this image
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fugue --loadfmap=fieldmap_rads_mean --mask=fieldmap_mag_bet_mask_union --unmaskfmap --savefmap=fieldmap_rads_mean']);
            
            %Remove unneeded files
            delete([preprocfolder '/fieldmap_mag_all.nii.gz']);
            delete([preprocfolder '/fieldmap_mag_all_mcf.nii.gz']);
            delete([preprocfolder '/fieldmap_mag_bet_mask_all.nii.gz']);
            delete([preprocfolder '/fieldmap_mag_bet_mask_all_mcf.nii.gz']);
            delete([preprocfolder '/fieldmap_rads_all.nii.gz']);
            delete([preprocfolder '/fieldmap_rads_all_mcf.nii.gz']);
            delete(['fieldmap_mag_all_mcf.mat/*'])
            rmdir(['fieldmap_mag_all_mcf.mat'],'s')
            
            %Check for problems
            if any(cell2mat(systemresult(:,1)))
                error(['Subject ' subject ': Problem detected with fieldmap processing'])
            end
            
        end
        
        
        %% Functional preprocessing
        
        
        
        if BOLD_preprocess
            
            cd(preprocfolder)
            
            dlmwrite(logfile,'preprocessing BOLD data','-append','delimiter','')
            disp(['Subject ' subject ': slice time correction, motion correction, atlas registration, and intensity normalization']);
            
            brainmaskfile = [freesurferfolder '/nusmask/aparc+aseg_brainmask_mask_' num2str(normalized_voxdim) num2str(normalized_voxdim) num2str(normalized_voxdim) '.nii.gz'];
            brainmask = load_untouch_nii(brainmaskfile);
            
            fmstring = '';
            
            for seq = 1:length(functional_sequences)
                TR = functional_sequence_TRs{seq};
                
                %make a file that tracks session number
                sess_tracker = [preprocfolder '/BOLDruns_sessions_' functional_sequences{seq} '.txt'];
                delete(sess_tracker);
                fid = fopen(sess_tracker,'at'); %open the output file for writing
                fclose(fid);
                
                avg_mergestr_native = ['fslmerge -t allBOLDavgs_' functional_sequences{seq}];
                avg_mergestr_native_bet = ['fslmerge -t allBOLDavgs_' functional_sequences{seq} '_bet '];
                avg_mergestr_MNI = ['fslmerge -t allBOLDavgs_' functional_sequences{seq} '_MNI '];
                
                runnametracker = cell(0,1);
                
                bestrun_name = cell(length(functional_sequences),1);
                                
                for runnum = 1:length(BOLDrunnames)
                    BOLDfile = BOLDrunnames{runnum};
                    
                    %get short filename
                    [~,BOLDfilename,~] = fileparts(BOLDfile);
                    
                    %If this run is the right sequence
                    if (length(BOLDfilename) >= length(functional_sequences{seq})) && strcmp(BOLDfilename(1:length(functional_sequences{seq})),functional_sequences{seq});
                        
                        runnametracker{end+1} = BOLDfilename;
                        
                        %motion-correct raw BOLD images to first image to get "true" motion plots
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['mcflirt -in ' BOLDfile ' -refvol 0 -plots']);
                        
                        %delete resulting motion corrected timeseries--we will include this in a single-step atlas transformation below
                        delete([BOLDfile '_mcf.nii.gz'])
                        
                        %slice-time correct BOLD images
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['slicetimer -i ' BOLDfile ' -o ' BOLDfile '_st -r ' num2str(TR) ' ' sliceorderinfo{seq}]);
                        
                        try delete([BOLDfile '_st_mcf.mat/*'])
                        rmdir([BOLDfile '_st_mcf.mat'],'s'); catch; end
                    
                        %motion-correct BOLD images to first image to get motion correction registrations
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['mcflirt -in ' BOLDfile '_st -refvol 0 -mats']);
                        
                        %split slice-time corrected timeseries into separate volumes
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslsplit ' BOLDfile '_st ' BOLDfile '_st_mcf.mat/' BOLDfilename '_st_mcf -t']);
                        
                        %get avg BOLD for this run
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' BOLDfile '_st_mcf -Tmean ' BOLDfile '_st_mcf_avg']);
                        
                        %delete motion corrected timeseries--we will include this in a single-step atlas transformation below
                        delete([BOLDfile '_st_mcf.nii.gz']) 
                        
                        %brain-extract average BOLD image
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(['bet ' BOLDfile '_st_mcf_avg ' BOLDfile '_st_mcf_avg_bet -R']);
                        
                        
                        if apply_fieldmap
                            
                            fmstring = '_uwrp';
                            
                            %define fieldmap files. if those files don't exist, use the mean fieldmaps. 
                            this_fieldmaprads = ['fieldmap_rads_mean_sess' num2str(sesscounter(runnum))];
                            this_fieldmapmask = ['fieldmap_mag_bet_mask_union_sess' num2str(sesscounter(runnum))];
                            this_fieldmapmag = ['fieldmap_mag_mean_sess' num2str(sesscounter(runnum)) '_bet'];
                            if ~exist([this_fieldmaprads '.nii.gz'],'file')
                                this_fieldmaprads = ['fieldmap_rads_mean'];
                                this_fieldmapmask = ['fieldmap_mag_bet_mask_union'];
                                this_fieldmapmag = ['fieldmap_mag_mean_bet'];
                                if ~exist([this_fieldmaprads '.nii.gz'],'file')
                                    error(['Subject ' subject ': Field map correction requested but no field maps found!'])
                                end
                            end
                            
                            %Register field map images to this average BOLD image
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' this_fieldmapmag ' -ref ' BOLDfile '_st_mcf_avg_bet -dof 6 -omat ' this_fieldmapmag '_2' BOLDfilename '.mat']);
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' this_fieldmaprads ' -ref ' BOLDfile '_st_mcf_avg_bet -applyxfm -init ' this_fieldmapmag '_2' BOLDfilename '.mat -out ' this_fieldmaprads '_' BOLDfilename 'space']);
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' this_fieldmapmask ' -ref ' BOLDfile '_st_mcf_avg_bet -applyxfm -init ' this_fieldmapmag '_2' BOLDfilename '.mat -interp nearestneighbour -out ' this_fieldmapmask '_' BOLDfilename 'space']);
                            
                            %Warp average BOLD image with field map
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fugue -i ' BOLDfile '_st_mcf_avg_bet --dwell=' num2str(echospacing) ' --loadfmap=' this_fieldmaprads '_' BOLDfilename 'space --mask=' this_fieldmapmask '_' BOLDfilename 'space -u ' BOLDfile '_st_mcf_avg_bet_uwrp --unmaskshift --saveshift=' this_fieldmaprads '_' BOLDfilename 'space_shift']);
                            
                            %Make a "warp" image from field map, for later use
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fugue -i ' BOLDfile '_st_mcf_avg --dwell=' num2str(echospacing) ' --loadfmap=' this_fieldmaprads '_' BOLDfilename 'space --mask=' this_fieldmapmask '_' BOLDfilename 'space -u ' BOLDfile '_st_mcf_avg_uwrp']);
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['convertwarp -s ' this_fieldmaprads '_' BOLDfilename 'space_shift -r ' BOLDfile '_st_mcf_avg_uwrp -o ' this_fieldmaprads '_' BOLDfilename 'space_warp']);
                        end
                        %Add to string to concatenate avg BOLD images
                        avg_mergestr_native_bet = [avg_mergestr_native_bet ' ' BOLDfile '_st_mcf_avg_bet' fmstring];
                        avg_mergestr_native = [avg_mergestr_native ' ' BOLDfile '_st_mcf_avg' fmstring];
                    end
                end
                
                try delete(['allBOLDavgs_' functional_sequences{seq} '_bet_mcf.mat/*'])
                    rmdir(['allBOLDavgs_' functional_sequences{seq} '_bet_mcf.mat'],'s')
                catch; end
                
                [systemresult{end+1,1},systemresult{end+1,2}] = system(avg_mergestr_native_bet);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(avg_mergestr_native);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['mcflirt -in allBOLDavgs_' functional_sequences{seq} '_bet -refvol 0 -bins 10000 -mats']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['applyxfm4D allBOLDavgs_' functional_sequences{seq} ' allBOLDavgs_' functional_sequences{seq} '_bet_mcf allBOLDavgs_' functional_sequences{seq} '_mcf allBOLDavgs_' functional_sequences{seq} '_bet_mcf.mat -userprefix MAT_']);
                alldata = load_untouch_nii_2D(['allBOLDavgs_' functional_sequences{seq} '_mcf.nii.gz']);
                corrs = paircorr_mod(double(alldata.img));
                alldata = [];
                [~,bestind] = max(sum(corrs,2));
                bestrun_name{seq} = runnametracker{bestind};
                try delete(['allBOLDavgs_' functional_sequences{seq} '_bet_mcf.mat/*'])
                    rmdir(['allBOLDavgs_' functional_sequences{seq} '_bet_mcf.mat'],'s')
                catch; end
                
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['mcflirt -in allBOLDavgs_' functional_sequences{seq} '_bet -refvol ' num2str(bestind-1) ' -bins 10000 -mats']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['applyxfm4D allBOLDavgs_' functional_sequences{seq} ' allBOLDavgs_' functional_sequences{seq} '_bet_mcf allBOLDavgs_' functional_sequences{seq} '_mcf allBOLDavgs_' functional_sequences{seq} '_bet_mcf.mat -userprefix MAT_']);
%                 alldata = load_untouch_nii_2D(['allBOLDavgs_' functional_sequences{seq} '_mcf.nii.gz']);
%                 corrs = paircorr_mod(double(alldata.img));
%                 corrs = corrs(:,bestind);
%                 acceptable = corrs >= .98;
%                 alldata = [];
%                 
%                 string = [];
%                 unacceptableinds = find(~acceptable);
%                 for i = 1:length(unacceptableinds)
%                     string = [string runnametracker{unacceptableinds(i)} ' '];
%                 end
%                 if ~isempty(string)
%                     disp(['Subject ' subject ': the following runs failed registration and will not be included: ' string]);
%                 end
                
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' bestrun_name{seq} '_st_mcf_avg_bet' fmstring ' -ref '  T1file '_bet -dof 6 -omat ' bestrun_name{seq} '_avg_2T1_init.mat']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' bestrun_name{seq} '_st_mcf_avg_bet' fmstring ' -ref '  T1file '_bet -dof 12 -cost bbr -wmseg aparc+aseg_cerebralwm -init ' bestrun_name{seq} '_avg_2T1_init.mat -omat ' bestrun_name{seq} '_avg_2T1.mat -schedule ${FSLDIR}/etc/flirtsch/bbr.sch']);
                
                %Make MNI registration
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['convert_xfm -omat ' bestrun_name{seq} '_avg_2MNI.mat -concat T1_2MNI.mat ' bestrun_name{seq} '_avg_2T1.mat']);
                
                runnumtracker = 0;
                for runnum = 1:length(BOLDrunnames)
                    BOLDfile = BOLDrunnames{runnum};
                    
                    %get short filename
                    [~,BOLDfilename,~] = fileparts(BOLDfile);
                    
                    %If this run is the right sequence
                    if (length(BOLDfilename) >= length(functional_sequences{seq})) && strcmp(BOLDfilename(1:length(functional_sequences{seq})),functional_sequences{seq});
                        runnumtracker = runnumtracker+1;
                        
                        this_fieldmaprads = ['fieldmap_rads_mean_sess' num2str(sesscounter(runnum))];
                        this_fieldmapmask = ['fieldmap_mag_bet_mask_union_sess' num2str(sesscounter(runnum))];
                        this_fieldmapmag = ['fieldmap_mag_mean_sess' num2str(sesscounter(runnum)) '_bet'];
                        if ~exist([this_fieldmaprads '.nii.gz'],'file')
                            this_fieldmaprads = ['fieldmap_rads_mean'];
                            this_fieldmapmask = ['fieldmap_mag_bet_mask_union'];
                            this_fieldmapmag = ['fieldmap_mag_mean_bet'];
                        end
                        
                        %if acceptable(runnumtracker)
                            
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['convert_xfm -omat ' BOLDfile '_avg_2MNI.mat -concat ' bestrun_name{seq} '_avg_2MNI.mat allBOLDavgs_' functional_sequences{seq} '_bet_mcf.mat/MAT_' sprintf('%04i',runnumtracker-1)]);
                            
                            %Register session warped average BOLD image to MNI
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' BOLDfile '_st_mcf_avg' fmstring ' -interp spline -ref ' template ' -applyisoxfm ' num2str(normalized_voxdim) ' -init ' BOLDfile '_avg_2MNI.mat -out ' BOLDfile '_st_mcf_avg_MNI']);
                            
                            avg_mergestr_MNI = [avg_mergestr_MNI ' ' BOLDfile '_st_mcf_avg_MNI'];
                            
                            %get volume numbers
                            vols = dir([BOLDfile '_st_mcf.mat/MAT*']);
                            fslmergestr = ['fslmerge -t ' BOLDfile '_st_mcf_MNI'];
                            %for each volume
                            for i = 1:length(vols)
                                %concatenate motion correction transform with run avg to MNI transform
                                [systemresult{end+1,1},systemresult{end+1,2}] = system(['convert_xfm -omat ' BOLDfile '_st_mcf.mat/' vols(i).name '_MNI.mat -concat  ' BOLDfile '_avg_2MNI.mat ' BOLDfile '_st_mcf.mat/' vols(i).name]);
                                
                                %apply concatenated transform using the "warp" image
                                if apply_fieldmap
                                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['applywarp -i ' BOLDfile '_st_mcf.mat/' BOLDfilename '_st_mcf' sprintf('%04i',(i-1)) ' --interp=spline -r  ' BOLDfile '_st_mcf_avg_MNI -o ' BOLDfile '_st_mcf.mat/' BOLDfilename '_st_mcf_MNI' sprintf('%04i',(i-1)) ' -w ' this_fieldmaprads '_' BOLDfilename 'space_warp --postmat=' BOLDfile '_st_mcf.mat/' vols(i).name '_MNI.mat']);
                                else
                                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' BOLDfile '_st_mcf.mat/' BOLDfilename '_st_mcf' sprintf('%04i',(i-1)) ' -interp spline -ref ' template ' -applyisoxfm ' num2str(normalized_voxdim) ' -init '  BOLDfile '_st_mcf.mat/' vols(i).name '_MNI.mat -out ' BOLDfile '_st_mcf.mat/' BOLDfilename '_st_mcf_MNI' sprintf('%04i',(i-1))]);
                                end
                                
                                %prepare to merge all volumes back into a 4D image
                                fslmergestr = [fslmergestr ' ' BOLDfile '_st_mcf.mat/' BOLDfilename '_st_mcf_MNI' sprintf('%04i',(i-1))];
                            end
                            %merge volumes
                            [systemresult{end+1,1},systemresult{end+1,2}] = system(fslmergestr);
                            
                            
                            %Intensity normalization to mode within-brain value of 1000
                            data = load_untouch_nii([BOLDfile '_st_mcf_MNI.nii.gz']);
                            data_inmask = data.img .* repmat(brainmask.img,[1,1,1,size(data.img,4)]) .* (data.img > 100);
                            data_inmask = data_inmask(data_inmask>0);
                            [counts,edges] = histcounts(data_inmask,1000);
                            [~,maxind] = max(counts);
                            modeval = mean([edges(maxind) edges(maxind+1)]);
                            data.img = data.img + (1000 - modeval);
                            save_untouch_nii(data,[BOLDfile '_st_mcf_MNI.nii.gz']);
                            data = []; data_inmask = [];
                            
                            %track session number
                            dlmwrite(sess_tracker,[BOLDfile '_st_mcf_MNI.nii.gz ' num2str(sesscounter(runnum)) ' ' BOLDfile '_mcf.par'],'-append','delimiter','');%write the data to the output file
                            
                        %end
                        
                        %remove unneeded files
                        delete([BOLDfile '_st_mcf.mat/*'])
                        rmdir([BOLDfile '_st_mcf.mat'],'s')
                        delete([BOLDfile '_st.nii.gz'])
                        delete([BOLDfile '_avg_2T1_init.mat'])
                        %delete([BOLDfile '_st_mcf_avg_bet.nii.gz'])
                        
                    end
                    
                end
                
                %merge averages
                [systemresult{end+1,1},systemresult{end+1,2}] = system(avg_mergestr_MNI);
                try delete(['allBOLDavgs_' functional_sequences{seq} '_bet_mcf.mat/*'])
                    rmdir(['allBOLDavgs_' functional_sequences{seq} '_bet_mcf.mat'],'s'); 
                catch; end
            end
            
            brainmask = [];
            
            if any(cell2mat(systemresult(:,1)))
                error(['Subject ' subject ': Problem detected with BOLD preprocessing'])
            end
            
        end
        
        cprintf('*blue',['Subject ' subject ' COMPLETED preprocessing without error!\n'])
        dlmwrite(logfile,'all processing complete.','-append','delimiter','')
        
    catch errorinfo
        
        dlmwrite(logfile,'ERROR:','-append','delimiter','')
        
        errormessage = [errorinfo.stack(end).file ', line ' num2str(errorinfo.stack(end).line) ': ' errorinfo.message];
        
        probleminds = find(cell2mat(systemresult(:,1)));
        for i = 1:length(probleminds)
            dlmwrite(logfile,systemresult{probleminds(i),2},'-append','delimiter','');
        end
        dlmwrite(logfile,errormessage,'-append','delimiter','');
        
        cprintf('err',['Subject ' subject ' FAILED to complete processing! Check error log ' logfile '\n'])
        
    end
    
    
end


try delete(processingpool); catch; end


end






if any([BOLD_fcprocess cifti_creation cifti_correlation])

    enoughtime = true(length(subjects),length(functional_sequences));
    
    disp('POSTPROCESSING')
    
for subnum = 1:subcount
    
    systemresult = cell(0,2);
    
    subject = subjects{subnum};
    
    warning off
    
    try
        
        %% Set up folders
        
        
        rawfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/raw/'];
        preprocfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/preprocessed/'];
        freesurferfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/freesurfer/']; 
        fsLRfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/fs_LR/']; 
        fc_processfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/fc_processed/']; 
        ciftifolder = ['/home/data/EEG/processed/Robi/mri/' subject '/cifti/']; 
        correlationfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/connectome/'];
        infomapfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/infomap/'];
        logfolder = ['/home/data/EEG/processed/Robi/mri/' subject '/processing_logs/']; 
        
        logfile2 = [logfolder datestr(now) '.txt']; tokens = tokenize(logfile2,' '); logfile2 = [tokens{1} '_' tokens{2}];
        
        %% Fc-processing
        
        
        if BOLD_fcprocess
            
            mkdir(fc_processfolder)
            
            dlmwrite(logfile2,'fc-processing BOLD data','-append','delimiter','')
            
            for seq = 1:length(functional_sequences)
                if fcprocess_sequences{seq}
                    
                    disp(['Subject ' subject ', ' functional_sequences{seq} ': fc-processing']);
                    %FC_Process(subject,['/home/data/EEG/processed/Robi/mri/' subject '/preprocessed/BOLDruns_sessions_' functional_sequences{seq} '.txt'],FDthresh,functional_sequence_TRs{seq},['/home/data/EEG/processed/Robi/mri/' subject '/fc_processed/'],functional_sequences{seq})
                    FC_Process_filt_Avi_norunsessregression(subject,['/home/data/EEG/processed/Robi/mri/' subject '/preprocessed/BOLDruns_sessions_' functional_sequences{seq} '.txt'],FDthresh,functional_sequence_TRs{seq},['/home/data/EEG/processed/Robi/mri/' subject '/fc_processed/'],functional_sequences{seq})
                    %FC_Process_filt(subject,['/home/data/EEG/processed/Robi/mri/' subject '/preprocessed/BOLDruns_sessions_' functional_sequences{seq} '.txt'],FDthresh,functional_sequence_TRs{seq},['/home/data/EEG/processed/Robi/mri/' subject '/fc_processed/'],functional_sequences{seq})
                    
                end
            end
            
        end
        
        
        for seq = 1:length(functional_sequences)
            tmask = load([fc_processfolder '/' functional_sequences{seq} '_all_tmask.txt']);
            if (nnz(tmask) .* functional_sequence_TRs{seq} ./ 60) < mintime
                enoughtime(subnum,seq) = false;
                disp(['Subject ' subject ' does not have at least ' num2str(mintime) ' minutes of ' functional_sequences{seq} ' data! Processing will not continue.'])
            end
        end
        
        
        
        %% Cifti creation
        
        if cifti_creation
            
            mkdir(ciftifolder)
            
            dlmwrite(logfile2,'creating ciftis','-append','delimiter','')
            
            medial_masks = {'/home/data/scripts/Resources/cifti_masks/L.atlasroi.32k_fs_LR.shape.gii','/home/data/scripts/Resources/cifti_masks/R.atlasroi.32k_fs_LR.shape.gii';...
                '/home/data/scripts/Resources/cifti_masks/L.atlasroi_erode3.32k_fs_LR.shape.gii','/home/data/scripts/Resources/cifti_masks/R.atlasroi_erode3.32k_fs_LR.shape.gii'};
            
            for seq = 1:length(functional_sequences)
                
                if enoughtime(subnum,seq)
                    
                    disp(['Subject ' subject ', ' functional_sequences{seq} ': mapping functional data to cortical surface']);
                    
                    if fcprocess_sequences{seq}
                        
                        systemresult = post_fc_processing_batch_singlesub(subject,[fc_processfolder '/' functional_sequences{seq} '_fc_processed_tmasked.nii.gz'],ciftifolder,[fsLRfolder '/MNI/'],functional_sequence_TRs{seq},functional_sequences{seq},[fc_processfolder '/' functional_sequences{seq} '_all_tmask.txt'],[preprocfolder '/BOLDruns_sessions_' functional_sequences{seq} '.txt'],[fc_processfolder '/' functional_sequences{seq} '_runs_sessions.txt'],medial_masks,geodesic_smooth,systemresult);
                        
                    else
                        
                        BOLDmergedfile = [preprocfolder '/' functional_sequences{seq} '_merged'];
                        
                        fslmergestr = ['fslmerge -t ' BOLDmergedfile ' '];
                        
                        [BOLDrunnames,~,~] = textread([preprocfolder '/BOLDruns_sessions_' functional_sequences{seq} '.txt'],'%s%s%s');
                        
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
                 
                ciftifile = [ciftifolder '/cifti_timeseries_normalwall/' functional_sequences{seq} '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(geodesic_smooth) '.dtseries.nii'];
                 
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
        
        
        
        %% Cifti Distances
        
        if cifti_distances
            
            dlmwrite(logfile2,'conducting geodesic distance calculation','-append','delimiter','')
            disp(['Subject ' subject ': getting point-to-point geodesic and euclidean distances'])
            
            surfacetousefolder = ['/home/data/EEG/processed/Robi/mri/' subject '/fs_LR/MNI/fsaverage_LR32k/'];
            distfolder = [surfacetousefolder '/distances/'];
            mkdir(distfolder)
            
            surfcoordsfileL = [surfacetousefolder '/' subject '.L.midthickness.32k_fs_LR.coord.gii'];
            surftopofileL = [surfacetousefolder '/' subject '.L.32k_fs_LR.topo.gii'];
            
            delete([distfolder '/Surface_distances_L_noHEAD.func.gii'])
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['caret_command -surface-geodesic ' surfcoordsfileL ' ' surftopofileL ' ' distfolder '/Surface_distances_L.func.gii true']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['caret_command -file-convert -format-convert ASCII ' distfolder '/Surface_distances_L.func.gii']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['awk ''NF > 25'' ' distfolder '/Surface_distances_L.func.gii > ' distfolder '/Surface_distances_L_noHEAD.func.gii']);
            delete([distfolder '/Surface_distances_L.func.gii'])
            
            
            surfcoordsfileR = [surfacetousefolder '/' subject '.R.midthickness.32k_fs_LR.coord.gii'];
            surftopofileR = [surfacetousefolder '/' subject '.R.32k_fs_LR.topo.gii'];
            
            delete([distfolder '/Surface_distances_R_noHEAD.func.gii'])
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['caret_command -surface-geodesic ' surfcoordsfileR ' ' surftopofileR ' ' distfolder '/Surface_distances_R.func.gii true']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['caret_command -file-convert -format-convert ASCII ' distfolder '/Surface_distances_R.func.gii']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['awk ''NF > 25'' ' distfolder '/Surface_distances_R.func.gii > ' distfolder '/Surface_distances_R_noHEAD.func.gii']);
            delete([distfolder '/Surface_distances_R.func.gii'])
            
            
            surffileL = [surfacetousefolder '/' subject '.L.midthickness.32k_fs_LR.surf.gii'];
            surfaceareafileL = [surfacetousefolder '/' subject '.L.midthickness.32k_fs_LR_surfaceareas.func.gii'];
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['wb_command -surface-vertex-areas ' surffileL ' ' surfaceareafileL]);
            
            surffileR = [surfacetousefolder '/' subject '.R.midthickness.32k_fs_LR.surf.gii'];
            surfaceareafileR = [surfacetousefolder '/' subject '.R.midthickness.32k_fs_LR_surfaceareas.func.gii'];
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['wb_command -surface-vertex-areas ' surffileR ' ' surfaceareafileR]);
            
            distancesLfile = [distfolder '/Surface_distances_L.dconn.nii'];
            if ~exist(distancesLfile)
                geo_distances = load([distfolder '/Surface_distances_L_noHEAD.func.gii']);
                geo_distances = single(geo_distances); geo_distances(:,1) = [];
                
                out = data;
                out.data = geo_distances; clear geo_distances
                out.brainstructure = ones(max(find(data.brainstructure==1)),1);
                out.pos = data.pos(1:length(out.brainstructure),:);
                out.brainstructurelabel = data.brainstructurelabel(1);
                out = rmfield(out,{'dim','transform','time'});
                ft_write_cifti_mod(distancesLfile,out);
                
                %save(distancesLfile,'geo_distances','-v7.3');%parsave(distancesLfile,geo_distances,'-v7.3')
                %clear geo_distances
                delete([distfolder '/Surface_distances_L_noHEAD.func.gii'])
            end
            
            distancesRfile = [distfolder '/Surface_distances_R.dconn.nii'];
            if ~exist(distancesRfile)
                geo_distances = load([distfolder '/Surface_distances_R_noHEAD.func.gii']);
                geo_distances = single(geo_distances); geo_distances(:,1) = [];
                
                out = data;
                out.data = geo_distances; clear geo_distances
                out.brainstructure = ones(max(find(data.brainstructure==2)) - max(find(data.brainstructure==1)),1);
                out.pos = data.pos((max(find(data.brainstructure==2)) - max(find(data.brainstructure==1))+1):max(find(data.brainstructure==2)),:);
                out.brainstructurelabel = data.brainstructurelabel(2);
                out = rmfield(out,{'dim','transform','time'});
                ft_write_cifti_mod(distancesRfile,out);
                
                
                %save(distancesRfile,'geo_distances','-v7.3');%parsave(distancesRfile,geo_distances,'-v7.3')
                %clear geo_distances
                delete([distfolder '/Surface_distances_R_noHEAD.func.gii'])
            end
            
            ciftifile = [ciftifolder '/cifti_timeseries_normalwall/' functional_sequences{seq} '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(geodesic_smooth) '.dtseries.nii'];
            mkdir([ciftifolder '/distances/'])
            Make_distmat_32k_fsLR_singlesub(subject,ciftifile,[ciftifolder '/distances/normalwall_distmat'])
            
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

