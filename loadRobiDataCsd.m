function [ X ] = loadRobiDataCsd( filename )
%LOADROBIDATACSD Summary of this function goes here
%   Detailed explanation goes here

if(~exist('filename', 'var'))
    filename = '/home/data/EEG/data/ROBI/ROBI_003/baseline eyes open/630158995692243270.eegData';
end

rereference = true;

channelCount = 34;
file = dir(filename);
fileLength = file.bytes / 8;
sampleCount = fileLength / channelCount;
if(floor(sampleCount) ~= sampleCount)
    sampleCount = floor(sampleCount);
    fileLenght = sampleCount * channelCount;
end
fileId = fopen(filename);
fprintf('\nreading...');
contents = fread(fileId, fileLength, 'double');
fclose(fileId);

%truncate if necessary
sampleCount = fileLength / channelCount;
if(sampleCount ~= floor(sampleCount))
    sampleCount = floor(sampleCount);
    fileLength = sampleCount * channelCount;
    contents = contents(1:fileLength);
end
[labels, chanlocs] = antChannelLocs;

data = reshape(contents, channelCount, fileLength / channelCount)';
clear contents;

if(rereference)
    fprintf('\nrereferencing...');
    cpzIndex = find(strcmp(labels,'CPz'));
    m1Index = find(strcmp(labels,'M1'));
    m2Index = find(strcmp(labels,'M2'));
    counterIndex = find(strcmp(labels,'counter'));
    
    for i=1:size(data,1)
        sample = data(i,:);
        linkedMastoids = (sample(m1Index) + sample(m2Index)) / 2;
        avg = mean(sample(1:33));
        newSample = sample - linkedMastoids;
        avgSample = sample - avg;
        data(i,1:33) = newSample(1:33);
    end
end

fprintf('\ncomputing current source density...');
remove = [m1Index, m2Index, counterIndex];
data(:, remove) = [];
labels(:, remove) = [];
chanlocs(:, remove) = [];

% ------------ Step 1 -----------------------------------------------------
% understand the spherical coordinate system of the CSD toolbox
% ------------ Step 2 -----------------------------------------------------
E = textread('E31.asc','%s');
chanMontage = ExtractMontage('10-5-System_Mastoids_EGI129.csd',labels');
% ------------ Step 3 -----------------------------------------------------
[G,H] = GetGH(chanMontage);
% ------------ Step 4 -----------------------------------------------------
% D = textread('NR_C66_trr.dat');
% D = D';
data = data';
% ------------ Step 5 -----------------------------------------------------
X = CSD (data, G, H);
X = X';
%debug
t = ((1:length(X)) ./ 2048)';
figure;           % create new figure
plot(t, [X(:,1), data(1,:)']);
legend([{'csd'},{'raw'}]);

%end debug
if(false)
    figure;           % create new figure
    plot(X);
    plot(X(:,[14 24]));
    D = D';
    figure;           % create new figure
    plot(D);
    plot(D(:,[14 24]));
    % ------------ Step 6 -----------------------------------------------------
    WriteMatrix2Text(X,'CSD_C66_trr.dat');
end
fprintf('\n');
