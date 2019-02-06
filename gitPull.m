%this script pulls from github and updates file modified times.

cd(fileparts(which('gitPull')));
[out1, out2] = system('git pull');

applyFiletimes