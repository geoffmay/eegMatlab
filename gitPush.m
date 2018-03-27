message = sprintf('gitPush.m-%s', getComputerName);

preserveFiletimes;
oldfolder = pwd;
cd(fileparts(which('gitPush')));
[error1, message1] = system('git add --all');
if(error1)
    error(sprintf('error using gitAdd: %s', message1));
end
a = sprintf('git commit -am %s', message);
[out3, out4] = system(a);
[out5, out6] = system('git push');
cd(oldfolder);
