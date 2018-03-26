function output = lookAtErp()
folder = '/home/data/EEG/processed/Oregon/ERP/summary4';
files = dir(folder);
files([files.isdir]) =[];
end