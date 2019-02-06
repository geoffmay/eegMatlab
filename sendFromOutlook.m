function sendFromOutlook( toAddress, subject, message, attachments )
%SENDFROMOUTLOOK Summary of this function goes here
%   Detailed explanation goes here


fromAddress = 'geoffrey.may@outlook.com'; 
setpref('Internet','SMTP_Server','smtp-mail.outlook.com');
cipher1 = [161    31   239    24   125    26   250    30    88    37    80   177     3   134     3   113];
cipher2 = [178   208   222   195   133   156   177   109   169   107     3   191   218   176    88   205];
load('/home/gmay/Documents/cipher.mat');
if(~exist('inv_ciper', 'file'))
  addpath('/home/data/EEG/scripts/AES');
end
plain1 = inv_cipher(cipher1, w1, inv_s_box, inv_poly_mat);
plain2 = inv_cipher(cipher2, w1, inv_s_box, inv_poly_mat);
plain1(plain1 == 0) = [];
plain2(plain2 == 0) = [];
plain = [plain1 plain2];
password = char(plain);
setpref('Internet','E_mail',fromAddress);
setpref('Internet','SMTP_Username',fromAddress);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
props.setProperty('mail.smtp.starttls.enable','true');

% Send the email.
sendmail1(toAddress, subject, message, attachments);

end

