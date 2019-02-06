
testFile = '/home/data/EEG/data/ROBI/ROBI_003/baseline eyes open/630158995692243270.eegData';
testIca1 = deriveIcaCoherenceMatrix(testFile);
testIca2 = deriveIcaCoherenceMatrix(testFile, testIca1.icaInfo.spheredWeights);

folder = '/home/data/EEG/processed/Robi/icaCohPca';

comb = load(fullfile(folder, 'ROBI_002combined.mat'))
pre = load(fullfile(folder, 'ROBI_002pre.mat'));
post = load(fullfile(folder, 'ROBI_002post.mat'));



comLoc = comb.icaCoh.icaInfo.spheredLocations;
preLoc = pre.icaCoh.icaInfo.spheredLocations;
posLoc = post.icaCoh.icaInfo.spheredLocations;

xs = [[preLoc.X], [posLoc.X], [comLoc.X]]';
ys = [[preLoc.Y], [posLoc.Y], [comLoc.Y]]';
zs = [[preLoc.Z], [posLoc.Z], [comLoc.Z]]';
gs = [repmat(1, [length(preLoc),1]); repmat(2, [length(preLoc),1]); repmat(3, [length(preLoc),1])];

xs = [[preLoc.X], [posLoc.X]]';
ys = [[preLoc.Y], [posLoc.Y]]';
zs = [[preLoc.Z], [posLoc.Z]]';
gs = [repmat(1, [length(preLoc),1]); repmat(2, [length(preLoc),1])];

figure;
gscatter(xs, ys, gs);
legend({'pre', 'post', 'combined'});
xlabel('x');
ylabel('y');

figure;
gscatter(xs, zs, gs);
legend({'pre', 'post', 'combined'});
xlabel('x');
ylabel('z');

figure;
gscatter(ys, zs, gs);
legend({'pre', 'post', 'combined'});
xlabel('y');
ylabel('z');




