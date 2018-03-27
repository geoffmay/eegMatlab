message = sprintf('gitPush.m-%s', getComputerName);

preserveFiletimes;
oldfolder = pwd;
cd(fileparts(which('gitPush')));
[error1, message1] = system('git add --all');
if(error1)
    error(sprintf('error using "git add": %s', message1));
end
a = sprintf('git commit -am %s', message);
[error2, message2] = system(a);
if(error2)
    error(sprintf('error using "git commit": %s', message2));
end

[out5, out6] = system('git push');
cd(oldfolder);
