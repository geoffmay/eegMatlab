function [ components, counts ] = countComponentOccurences( labels )
%COUNTCOMPONENTOCCURENCES Summary of this function goes here
%   Detailed explanation goes here

components = cell(0);
counts = [];
for i = 1:length(labels)
  line = labels{i};
  items = strsplit(line, ' ');
  theseComps = strsplit(items{2}, '-');  
  for j = 1:length(theseComps)
    found = false;
    for k = 1:length(components)
      if(strcmp(components{k}, theseComps{j}))
        counts(k) = counts(k) + 1;
        found = true;
      end
    end
    if(~found)
      components{length(components) + 1} = theseComps{j};
      counts(length(counts) + 1) = 1;
    end
  end
end
components = components';
counts = counts';

tab = table(components, counts);
tab = sortrows(tab, 2);
tab