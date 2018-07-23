function gitPush(commitMessage)

if(~exist('commitMessage'))
    commitMessage = sprintf('gitPush.m-%s', getComputerName);
end

preserveFiletimes;

oldfolder = pwd;
<<<<<<< HEAD
cd /home/data/EEG/scripts/bespoke;
[out1, out2] = system('git add --all');
[out3, out4] = system('git commit -am automated');
[out5, out6] = system('git push');
=======
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

>>>>>>> 13ca5b7edd1213a84493b6d2641fc2cff415801d
cd(oldfolder);
end
