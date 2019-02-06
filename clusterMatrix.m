function [ newMat, perm ] = clusterMatrix( mat, perm )
%CLUSTERMATRIX Summary of this function goes here
%   Detailed explanation goes here

meanMat = mean(mat, 3);
if(~exist('perm', 'var'))
    perm = clusterPermutation(meanMat);
end
newMat = ones(size(meanMat));
for dest1 = 1:size(newMat,1)
    source1 = perm(dest1);
    for dest2 = 1:size(newMat,2)
        source2 = perm(dest2);
        for k = 1:size(mat,3)
            newMat(dest1,dest2,k) = mat(source1,source2,k);
        end
    end
end

end

