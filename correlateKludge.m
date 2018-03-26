function [ pvalM ] = correlateKludge(data )
%CORRELATE Summary of this function goes here
%   Detailed explanation goes here


for(dependent = 1:135)
    subtable = data(:,dependent);
    sample = subtable{1,1};
    if(isnumeric(sample))
        [rho, pval] = corr(table2array(subtable),table2array(data(:,135:end)), 'rows', 'pairwise');
        pval(pval==0)=NaN;
        minP = nanmin(pval)*1848*135;
        if(minP < 0.05)
            pvalM(:,dependent) = pval;
        end
    end
end

end

