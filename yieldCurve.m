filename = 'C:\Users\Neuro\Downloads\DailyTreasuryYieldCurveRateData.txt';
fileId = fopen(filename);
text = fscanf(fileId, '%c');
fclose(fileId);

dayStart = '<content type="application/xml">';
dayEnd = '</content>'
lineStart = '<d:';
lineEnd = '</d';

starts = strfind(text, dayStart);
ends = strfind(text, dayEnd);

for i = 1:length(starts)
    item = text((starts(i) + length(dayStart)):(ends(i) - 1))
    newlines = strfind(item, sprintf('\n'));
    lineStarts = strfind(item, lineStart);
    lineEnds = strfind(item, lineEnd);
end