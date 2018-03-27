function gitPush(commitMessage)

if(~exist('commitMessage'))
    commitMessage = sprintf('gitPush.m-%s', getComputerName);
end

preserveFiletimes;

oldfolder = pwd;
cd(fileparts(which('gitPush')));
[error1, message1] = system('git add --all');
if(error1)
    error(sprintf('error using "git add": %s', message1));
else
    fprintf('\n%s', message1);
end

commitCommand = sprintf('git commit -am %s', commitMessage);
[error2, message2] = system(commitCommand);
if(error2)
    error(sprintf('error using "git commit": %s', message2));
else
    fprintf('\n%s', message2);
end

[error3, message3] = system('git push');
if(error2)
    error(sprintf('error using "git push": %s', message3));
else
    fprintf('\n%s', message3);
end

cd(oldfolder);
end