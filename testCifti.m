folder = '/home/data/subjects/ROBI003/cifti/cifti_timeseries_normalwall';
file = 'RSFC_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii';

a = ft_read_cifti(fullfile(folder,file), 'readdata', 'true');


folder2 = '/home/data/subjects/ROBI003/cifti/parcel_timecourses';
file2 = 'RSFC_Parcels_LR.ptseries.nii';

ptSeries = ft_read_cifti(fullfile(folder2,file2), 'readdata', 'true');

c = load('/home/data/subjects/ROBI003/connectome/RSFC_Parcels_LR_corr.mat');

d = load('/home/data/subjects/ROBI003/connectome/RSFC_parcels_edgethresh_0.5_corr.mat');



e = gifti('/home/data/subjects/ROBI003/cifti/surf_timecourses/RSFC_L_dil10_32k_fs_LR_smooth4.25.func.gii');

f = ft_read_cifti('/home/data/subjects/ROBI003/infomap/RSFC/rawassn_minsize400_regularized_recoloredv4_striped_164.dtseries.nii', 'readdata', 'true');
f = ft_read_cifti_mod('/home/data/subjects/ROBI003/infomap/RSFC/rawassn_minsize400_regularized_recoloredv4.dscalar.nii');

%f.dtseries appears to map each voxel to a network.

f1 = ft_read_cifti('/home/data/subjects/ROBI003/infomap/RSFC/rawassn_minsize400_regularized_allcolumns_recoloredv4.dscalar.nii', 'readdata', 'true');

corticalVoxels = find(a.brainstructure == 1 | a.brainstructure == 2);

timeSeries = a.x55(corticalVoxels, :);






rawfolder = '/home/data/EEG/processed/Robi/mri/ROBI003/raw';
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


functional_sequences = {'RSFC'};


counter = 1;
sums = [];
for sess = 1:length(sessions)
  sums(sess) = 0;
  dashloc = strfind(sessions(sess).name,'-');
  sessnum = str2num(sessions(sess).name(dashloc+1:end));
%   sessQCindices = subQCindices(subQC_sessions==sessnum);
%   sessQC_decision = cell2mat(QCdata(sessQCindices,7));
%   sessQC_acquisitionnum = cell2mat(QCdata(sessQCindices,8));
  for seq = 1:length(functional_sequences)
    thisseq_files = dir([rawfolder '/' sessions(sess).name '/*' functional_sequences{seq} '*.nii.gz']);
    for file = 1:length(thisseq_files)
      filename = [rawfolder '/' sessions(sess).name '/' thisseq_files(file).name];
      [~,sizestr] = system(['fslsize ' filename ' -s']);
      tokens = tokenize(sizestr,' ');
      datasize = [str2num(tokens{3}) str2num(tokens{5}) str2num(tokens{7}) str2num(tokens{9})];
      disp(filename);
      datasize
      allSizes(counter) = datasize(4);
      sums(sess) = sums(sess) + datasize(4);
      allSessions{counter} = sessions(sess).name;
      counter = counter + 1;
      
    end
  end
end
% 
% copied = dir('/home/data/EEG/processed/Robi/mri/ROBI003/preprocessed');
% copied([copied.isdir]) = [];
% counter = 1;
% for i = 1:length(copied)
% 
%       filename = ['/home/data/EEG/processed/Robi/mri/ROBI003/preprocessed/' copied(i).name];
%       [~,sizestr] = system(['fslsize ' filename ' -s']);
%       tokens = tokenize(sizestr,' ');
%       datasize = [str2num(tokens{3}) str2num(tokens{5}) str2num(tokens{7}) str2num(tokens{9})];
%       disp(filename);
%       datasize
%       allSizes1(counter) = datasize(4);
%       counter = counter + 1;
% 
% end


atlas = ft_read_cifti('/home/data/atlases/120_LR_minsize400_recolored_manualconsensus4.dtseries.nii', 'readdata', 'true', 'dataformat', 'raw');

atlas1 = ft_read_cifti_mod('/home/data/atlases/120_LR_minsize400_recolored_manualconsensus4.dtseries.nii');

atlas.dtseries




subjects = {'ROBI003'};
outfolder = '/home/data/EEG/processed/Robi/mri/analysis';%'/home/data/Analysis/ROBI/';

for s = 1:length(subjects)
    
    datafile = ['/home/data/subjects/' subjects{s} '/cifti/cifti_timeseries_normalwall/RSFC_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'];
    data = ft_read_cifti_mod(datafile);
    
    datastruct = data; datastruct.data = [];
    datastruct.dimord = 'pos_pos';
    
    runs_sessionsfile = ['/home/data/subjects/' subjects{s} '/fc_processed/RSFC_runs_sessions.txt'];
    
    runs_sessions = load(runs_sessionsfile);
    sessions = runs_sessions(:,2);
    clear run_sessions
    
    tmask = load(['/home/data/subjects/' subjects{s} '/fc_processed/RSFC_all_tmask.txt']);
    sessions = sessions(logical(tmask));
    
    cortex = find(data.brainstructure == 1 | data.brainstructure == 2a);
    uSess = unique(sessions);
    
    
    corr_pre = paircorr_mod(data.data(:,sessions<100)');
    corr_pre(isnan(corr_pre)) = 0;
    corr_pre = FisherTransform(corr_pre);
    
    datastruct.data = corr_pre;
    ft_write_cifti_mod([outfolder '/' subjects{s} '_pre'],datastruct);
    datastruct.data = [];
    
    clear corr_pre
    
    
    
    corr_post = paircorr_mod(data.data(:,sessions>100)');
    corr_post(isnan(corr_post)) = 0;
    corr_post = FisherTransform(corr_post);
    
    clear data
    
    datastruct.data = corr_post;
    ft_write_cifti_mod([outfolder '/' subjects{s} '_post'],datastruct);
    datastruct.data = [];
    
    corr_pre = ft_read_cifti_mod([outfolder '/' subjects{s} '_pre.dconn.nii']);
    
    corr_postminpre = corr_post - corr_pre.data;
    clear corr_pre corr_post
    
    datastruct.data = corr_postminpre;
    ft_write_cifti_mod([outfolder '/' subjects{s} '_postminpre'],datastruct);
    datastruct.data = [];
    
    postminpre_mean = mean(abs(corr_postminpre),2);
    datastruct.data = postminpre_mean;
    datastruct.dimord = 'pos_time';
    ft_write_cifti_mod([outfolder '/' subjects{s} '_postminpre'],datastruct);
    
    clear corr_postminpre
    
end

