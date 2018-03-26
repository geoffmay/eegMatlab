function [ newMat, perm ] = clusterMatrix( mat, perm, verbosity )
%CLUSTERMATRIX Summary of this function goes here
%   Detailed explanation goes here

%pass in text like 'none' to set verbosity without setting perm.
if(exist('perm','var'))
  if(ischar(perm))
    clear perm;
  end
end
if(~exist('verbosity', 'var'))
  verbosity = 0;
end
meanMat = mean(mat, 3);
if(~exist('perm', 'var'))
  if(verbosity > 0)
    fprintf('\n(%s) calculating linkage hierarchy', char(datetime));
  end
  perm = clusterPermutation(meanMat);
end
newMat = ones(size(meanMat));
if(verbosity > 0)
  fprintf('\n(%s) reordering matrix', char(datetime));
end
for dest1 = 1:size(newMat,1)
  source1 = perm(dest1);
  if(verbosity > 1)
    fprintf('\n(%s) row %d of %d', char(datetime), dest1, size(newMat,1));
  end
  for dest2 = 1:size(newMat,2)
    source2 = perm(dest2);
    for k = 1:size(mat,3)
      newMat(dest1,dest2,k) = mat(source1,source2,k);
    end
  end
end

end

