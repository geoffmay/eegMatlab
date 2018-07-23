preserveFiletimes;

oldfolder = pwd;
cd /home/data/EEG/scripts/bespoke;
[out1, out2] = system('git add --all');
[out3, out4] = system('git commit -am automated');
[out5, out6] = system('git push');
cd(oldfolder);
