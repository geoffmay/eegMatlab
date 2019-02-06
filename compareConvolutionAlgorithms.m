folder = 'C:\Users\Neuro\Documents\MATLAB\processed\MrEegLink';

a1 = load(fullfile(folder, 'GeoffTestEEG2-edf.edfconvolvedFirstPass.mat'));
a2 = load(fullfile(folder, 'GeoffTestEEG2-edf.withGamma.mat'));
a1 = a1.convolvedEeg;
a2 = a2.convolvedFourier;

t1 = a1.timeSeconds;
t2 = a2.times;

crossInd = NaN(size(t1));
counter = 1;
for i = 1:length(t1)
    for j = 1:length(t2)
        if(t1(i) == t2(j))
            crossInd(i) = j;
            counter = counter + 1;
        end
    end
end

% match = NaN(length(t2), size(a1.matrix, 2));

for j = 1:length(a1.labels)
    l1 = a1.labels{j};
    parts = strsplit(l1, ' ');
    search = [parts{2}, ' ', parts{1}, ' ', parts{3}, 'Hz'];
    j2 = find(strcmp(a2.labels, search));
    if(length(j2) > 0)
        col1 = a1.matrix(:, j);
        col2 = a2.signals(j2,:)';
        col1 = col1(~isnan(crossInd));
% for i = 1:length(t2)
%     fprintf('\n%d of %d', i, length(t2));
%     i1 = find(crossInd == i);
%     row1 = a1.matrix(i1, :);
%     row2 = a2.signals(:, i)';
%     
% end        
    end
end

