function [ outputMatrix ] = removeEmptyRowsAndColumns( inputMatrix )
%REMOVESOLIDNANROWSANDCOLUMNS Summary of this function goes here
%   Detailed explanation goes here
outputMatrix = inputMatrix(~all(isnan(inputMatrix),2),:); %remove solid NaN rows
outputMatrix = outputMatrix(:,~all(isnan(outputMatrix))); %remove solid NaN columns

outputMatrix(~any(outputMatrix,2),:) = []; %remove solid 0 rows
outputMatrix(:,~any(outputMatrix,1)) = []; %remove solid 0 columns

end
