function [ EEG ] = setStandardLocations( EEG )
%SETSTANDARDLOCATIONS Sets EEG.chanlocs using Standard-10-5-Cap385.sfp
%   Sets xyz, polar, and spherical coordinates of all channels that match
%   the label exactly (case-sensitive) of the standard EEG.  No
%   transformations are done.

%chans = readlocs('Standard-10-5-Cap385.sfp');
%chans = readlocs('/Applications/MATLAB_R2014a.app/toolbox/eeglab13_4_4b/plugins/dipfit2.3/standard_BESA/standard_BESA.mat');
%chans = readlocs('/Applications/MATLAB_R2014a.app/toolbox/eeglab13_4_4b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
eeglabMissing = false;
if(~exist('finputcheck'))
    eeglabFolder = '/home/data/EEG/scripts/eeglab13_4_4b';
    if(~exist('eeglab', 'file'))
        eeglabMissing = true;
        addpath(eeglabFolder);
    end
    eeglab;
end
chans = readlocs('/home/data/EEG/scripts/bespoke/standard-10-5-cap385.elp');
totalFound =0;
chanString = cell(0);
for j = 1:length(EEG.chanlocs)
    found = false;
    for i = 1:length(chans)
        if(strcmp(chans(i).labels, EEG.chanlocs(j).labels))
                chanString{j} = chans(i).labels;
            
            found = true;
            totalFound = totalFound +1;
            EEG.chanlocs(j).Y = chans(i).Y;
            EEG.chanlocs(j).X = chans(i).X;
            EEG.chanlocs(j).Z = chans(i).Z;
            EEG.chanlocs(j).sph_theta = chans(i).sph_theta;
            EEG.chanlocs(j).sph_phi = chans(i).sph_phi;
            EEG.chanlocs(j).sph_radius = chans(i).sph_radius;
            EEG.chanlocs(j).theta = chans(i).theta;
            EEG.chanlocs(j).radius = chans(i).radius;
            EEG.chanlocs(j).sph_theta_besa = chans(i).sph_theta_besa;
            EEG.chanlocs(j).sph_phi_besa = chans(i).sph_phi_besa;
            %EEG.chanlocs(j) = chans(j);
        end
    end
    if(~found)
        %here you would remove the channel if that's what needs to be done
    end    
end

% if(clear)
% for j = 1:length(EEG.chanlocs)
%             EEG.chanlocs(j).Y = [];
%             EEG.chanlocs(j).X = [];
%             EEG.chanlocs(j).Z = [];
%             EEG.chanlocs(j).sph_theta = [];
%             EEG.chanlocs(j).sph_phi = [];
%             EEG.chanlocs(j).sph_radius = [];
%             EEG.chanlocs(j).theta = [];
%             EEG.chanlocs(j).radius = [];
%             EEG.chanlocs(j).sph_theta_besa = [];
%             EEG.chanlocs(j).sph_phi_besa= [];
% end
% end

if(eeglabMissing)
    rmpath(genpath(eeglabFolder));
end

end

