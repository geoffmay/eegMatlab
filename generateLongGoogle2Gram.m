minLength = 5;
minFrequency = 10000;

folder = '/media/eegDrive/googleNGram/2gram';
hitCounter = 1;


files = dir(folder);
files([files.isdir]) = [];
for fileCounter = 22:length(files)
    fileId = fopen(fullfile(folder,files(fileCounter).name));
    counter = 1;
    fprintf('|');
    while(~feof(fileId))
        if(mod(counter, 10000) == 0)
            tick = counter / 10000;
            fprintf('\b');
            if(mod(tick, 4) == 0)
                fprintf('/');
            elseif(mod(tick,4)==1)
                fprintf('-');
            elseif(mod(tick,4)==2)
                fprintf('\\');
            elseif(mod(tick,4)==3)
                fprintf('|');
            end
        end
        line = fgets(fileId);
        tabs = strfind(line, sprintf('\t'));
        %ngram year occurences uniqueBooks
        ngram = line(1:tabs(1)-1);
        freqText = line(tabs(2)+1:tabs(3)-1);
        frequency = parseInt(freqText);

        
        if(counter == 1)
            oldNGram = ngram;
            totalFrequency = 0;
        else
            if(~strcmp(ngram, oldNGram))
                if(totalFrequency >= minFrequency && length(oldNGram) >= minLength)
                    words = strsplit(oldNGram, ' ');
                    word1 = words{1};
                    word2 = words{2};
                    skip = false;
                    if(word1(1) == '_' || word2(1) == '_')
                        skip = true;
                    end
                    if(~skip)
                        fprintf('\b\n%s: %d|',oldNGram, totalFrequency);
                        ngrams{hitCounter} = oldNGram;
                        frequencies(hitCounter) = totalFrequency;
                        oldNGram = ngram;
                        hitCounter = hitCounter + 1;
                    end
                end
                totalFrequency = 0;
            end
        end
        totalFrequency = totalFrequency + frequency;
        counter = counter + 1;
    end
    fclose(fileId);
    outputFilename = sprintf('/media/eegDrive/googleNGram/%s.mat', files(fileCounter).name);
    save(outputFilename, 'ngrams', 'frequencies');
    clear ngrams frequencies;

end