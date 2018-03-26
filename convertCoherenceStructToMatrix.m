function [ matrix, labels ] = convertCoherenceStructToMatrix( coh, freqInfo, pow )
%CONVERTCOHERENCESTRUCTTOMATRIX Summary of this function goes here
%   Detailed explanation goes here
matrix = NaN(size(coh(1).coherence, 1), size(coh(1).coherence, 2) * length(coh));
labels = cell(1, size(matrix,2));
counter = 1;
for i = 1:length(coh)
  input = coh(i).coherence;
  for j = 1:size(input,2)
    output = input(:,j);
    labels{counter} = sprintf('coh %s %dHz-%dHz', coh(i).label, freqInfo.lowFrequencies(j), freqInfo.highFrequencies(j));
    if(size(matrix(:,counter), 1) == size(output,1))
      matrix(:, counter) = output;
      counter = counter + 1;
    end
  end
end
if(exist('pow', 'var'))
  for i = 1:length(pow)
%     if(i == length(pow))
%       dummy = 1;
%     end
    input = pow(i).absolutePower;
    for j = 1:size(input,2)
      output = input(:,j);
      labels{counter} = sprintf('abs %s %dHz-%dHz', pow(i).label, freqInfo.lowFrequencies(j), freqInfo.highFrequencies(j));
      if(size(matrix(:,1), 1) == size(output,1))
        if(~any(isnan(output)))
          matrix(:, counter) = output;
          counter = counter + 1;
        end
      end
    end
  end
  if(isfield(pow, 'relativePower'))
    for i = 1:length(pow)
      input = pow(i).relativePower;
      for j = 1:size(input,2)
        output = input(:,j);
        labels{counter} = sprintf('rel %s %dHz-%dHz', pow(i).label, freqInfo.lowFrequencies(j), freqInfo.highFrequencies(j));
        if(size(matrix(:,1), 1) == size(output,1))
          if(~any(isnan(output)))
            matrix(:, counter) = output;
            counter = counter + 1;
          end
        end
      end
    end
  end
end
matrix(:, counter:end) = [];
labels(counter:end) = [];

end

