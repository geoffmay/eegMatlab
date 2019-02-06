%sampleMinutes = 0.25;
sampleMinutes = 24;

recordingMinutes = 60;
framesPerMinute = 60 * 128;

firstTerm = (recordingMinutes - sampleMinutes)*framesPerMinute;
secondTerm = (recordingMinutes - sampleMinutes*2)*framesPerMinute;

totalCombinations = firstTerm*secondTerm;


fprintf('\n%f (%d) combos \n',totalCombinations,totalCombinations)

