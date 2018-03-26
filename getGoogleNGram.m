for gramNumber = 2:3


rootUrl = 'http://storage.googleapis.com/books/ngrams/books/datasetsv2.html';

page = webread(rootUrl);
linkStarts = strfind(page, '<a href=''');
linkPseudoEnds = strfind(page, '</a>');
linkPseudoEnds(find(linkPseudoEnds < linkStarts(1))) = [];

links = cell(1,min(length(linkStarts),length(linkPseudoEnds)));
quote = '''';
endCounter = 1;
for i = 1:length(links)
    linkPseudoEnd = linkPseudoEnds(endCounter);
    while(linkPseudoEnd < linkStarts(i))
        endCounter = endCounter + 1;
        linkPseudoEnd = linkPseudoEnds(endCounter);
    end
    link = page(linkStarts(i):linkPseudoEnd);
    quotes = strfind(link, quote);    
    links{i} = link(quotes(1)+1:quotes(2)-1);
    endCounter = endCounter + 1;
end

searchText = sprintf('googlebooks-eng-all-%dgram', gramNumber);
excludeText = '.csv.zip';
goodLinks = find(cellfun(@length, strfind(links, searchText)));
goodLinks(find(cellfun(@length, strfind(links(goodLinks), excludeText)))) = [];

destinationFolder = sprintf('/media/eegDrive/googleNGram/%dgram/',gramNumber);
if(~exist(destinationFolder, 'dir'))
    mkdir(destinationFolder);
end
for i = 1:length(goodLinks)
    fprintf('\ndownloading %d of %d', i, length(goodLinks));
    link = links{goodLinks(i)};
    [folder, file, extension] = fileparts(link);
    destination = sprintf('%s%s%s', destinationFolder, file, extension);
    if(~exist(destination(1:end-3), 'file'))
        websave(destination, link);
        fprintf('...decompressing');
        [a, b] = system(sprintf('gunzip %s', destination));
        fprintf('...done');
    else
        fprintf('...already done');
    end
end
fprintf('\n');

end