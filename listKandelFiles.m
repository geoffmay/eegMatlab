function [ files ] = listKandelFiles( folder )
%LISTKANDELFILES Summary of this function goes here
%   Detailed explanation goes here

username = 'admin';
password = 'coeadmin';

[code,fileBlock] = system(['curl -l -u ', username, ':', password, ' ftp://192.168.1.196' folder]);
if(code == 0)
  files = strsplit(fileBlock, sprintf('\n'));
else
  error(fileBlock);
end

end

