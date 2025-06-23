function mod=read_model(fname);
fid=fopen(fname);
if fid==-1
    error('File doesn''t exist!');
end
tl=lower(fgetl(fid));
t=tokenize(tl);

%Search the model
while strcmp(t{1},'.model')==0&feof(fid)==0
    tl=lower(fgetl(fid));
    if isempty(tl)==0
        t=tokenize(tl);
    end
end
if feof(fid)
    fclose(fid);
    mod=0;
    return;
end
mod=default_ekv;
if strcmp(t{3},'nmos')
    mod.type=1;
else
    mod.type=0;
end
if length(t)>3
    t=t{4:length(t)};
end
tl=lower(fgetl(fid));
t=tokenize(tl)
fin=0;
while (t{1}(1)=='+'|t{1}(1)=='*'|feof(fid)==0)&fin==0
    if t{1}(1)=='*'
        tl=lower(fgetl(fid));
        t=tokenize(tl);
        continue;
    end
    if t{1}(1)~='+'
        fin=1;
        break;
    end
    t{1}=t{1}(2:length(t{1}));
    if isempty(t{1})
        t={t{2:length(t)}};
    end
    ll=length(t);
    if (round(ll/2)*2)~=ll
        warning('Strange line. Pass away!');
    else
        for i=1:(ll/2)
            let=isletter(t{i*2-1});
            val=str2num(t{i*2-1}(find(let==0)));
            if sum(1-let)==0|isempty(val)==0
                mod=setfield(mod,t{i*2-1},str2num_s(t{i*2}));
            else
                warning('Strange parameters encountred.');
            end
        end
    end
    tl=lower(fgetl(fid));
    t=tokenize(tl);
end
fclose(fid);

function tok=tokenize(tl);
c=1;
while isempty(tl)==0
    [s,tl]=strtok(tl);
    if isempty(find('='==s))
        if isempty(s)==0
            tok{c}=s;
            c=c+1;
        end
    else
        p=find('='==s);
        if p~=1
            tok{c}=s(1:(p-1));
            c=c+1;
        end
        if length(s)>p
            tok{c}=s((p+1):length(s));
            c=c+1;    
        end
    end
end        

function n=str2num_s(s);
s=s(find((s=='{'|s=='}')==0));
p=find(s=='+'|s=='-');
if isempty(p)==0
    if p(1)~=1
        if s(p(1)-1)~='e'
            s=s(1:(p(1)-1));
        end
    else
        if length(p)>1
            if s(p(2)-1)~='e'
                s=s(1:(p(2)-1));
            end
       end
   end
end
n=str2num(s);
if isempty(n)==0
    n
    return;
end    
l=isletter(s);
n=str2num(s(find(l==0)));
ex=lower(s(find(l==1)));
if strcmp(ex,'m')
    n=n*1e-3;
elseif strcmp(ex,'u')
    n=n*1e-6;
elseif strcmp(ex,'n')
    n=n*1e-9;
elseif strcmp(ex,'p')
    n=n*1e-12;
elseif strcmp(ex,'a')
    n=n*1e-15;
elseif strcmp(ex,'k')
    n=n*1e3;
elseif strcmp(ex,'meg')
    n=n*1e6;
elseif strcmp(ex,'g')
    n=n*1e9;
end
