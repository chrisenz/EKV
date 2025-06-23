function save_model(mod,name);
fid=fopen(name,'w');
if fid==-1
    error('File doesn''t exist!');
end
if mod.type==1
    type='nmos';
else
    type='pmos';
end
s=['.model ',type(1),' ',type,'\n'];
fprintf(fid,'*EKV model\n*Created with Matlab\n*by Nicola Scolari\n\n');
fprintf(fid,s);
nam=fieldnames(mod);
for i=1:length(nam);
    s=['+',nam{i}, ' = ',num2str(getfield(mod,nam{i})),'\n'];
    fprintf(fid,s);
end
fclose(fid);