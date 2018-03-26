tempDir = '/tmp/ftp';
if(~exist(tempDir, 'dir'))
  mkdir(tempDir);
end

ipAddr = '192.168.1.196';
usr = 'admin';
pass = 'coeadmin';

ftpObj = ftp(ipAddr, usr, pass);
a = cd(ftpObj, '/Backup');
contents = dir(ftpObj);


ngdc = ftp('ftp.ngdc.noaa.gov');
cd(ngdc,'pub');
a = dir(ngdc)

ftpObj.mget(


