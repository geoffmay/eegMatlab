hashinfoFilename = '/Users/Geoff/Documents/MATLAB/filehash/dirs.txt';

ipAddr = '192.168.1.196';
usr = 'admin';
pass = 'coeadmin';
ftpObj = ftp(ipAddr, usr, pass);

hashFolder(ftpObj, '/', 0);
