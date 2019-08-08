 function [err,ErrMessage,Info]=WriteBrikWrap(pathToData,img,info,prefix,view)
InfoWrite=info;
opt.Scale=0;
opt.Prefix=prefix;
opt.View=view;
opt.Verbose=1;
opt.AppendHistory=1;
opt.NoCheck=0;
opt.Overwrite=1;
origPath=pwd;
cd(pathToData);
brikFilename=[prefix '+' view '.BRIK'];
headFilename=[prefix '+' view '.HEAD'];
if exist(brikFilename,'file')
    eval(['delete ' brikFilename]);
    eval(['delete ' headFilename]);
end
[err, ErrMessage, Info] = WriteBrik (img, InfoWrite,opt);
cd(origPath);