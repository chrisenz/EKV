function varargout = ekvint(varargin)
% EKVINT Application M-file for ekvint.fig
%    FIG = EKVINT launch ekvint GUI.
%    EKVINT('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 31-May-2002 14:41:46

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    handles.ekvint_ver=0.3;
	guidata(fig, handles);
    
    ekvint_ver=0.3;
    global ver;
    ekvint_init(handles);
	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.






function ekvint_init(h);
h.min=struct('cox',1e-4,'xj',1e-8,'vto',0,'gamma',0,'phi',.3,'kp',1e-5,...
    'e0',1e6,'ucrit',1e6,'dl',-1e-7,'dw',-1e-7,'lambda',0,'leta',0,'weta',0,...
    'q0',0,'lk',.05e-6,'iba',0,'ibb',1.8e8,'ibn',.4,'rsh',0,'hdif',0);
h.max=struct('cox',1e-2,'xj',1e-6,'vto',2,'gamma',2,'phi',2,'kp',1e-3,...
    'e0',1e12,'ucrit',25e6,'dl',1e-7,'dw',1e-7,'lambda',3,'leta',2,'weta',2,...
    'q0',1e-3,'lk',2e-6,'iba',5e8,'ibb',4e8,'ibn',1,'rsh',1e3,'hdif',1e-6);
h.model=struct('cox',3.45e-3,'xj',.15e-6,'vto',.7,'gamma',.7,'phi',.5,'kp',150e-6,...
    'e0',200e6,'ucrit',2.3e6,'dl',-7.5e-8,'dw',-7.5e-8,'lambda',0.8,'leta',0.3,'weta',0.2,...
    'q0',230e-6,'lk',.4e-6,'iba',2.0e8,'ibb',2e8,'ibn',.6,'theta',0,'ekvint',0,...
    'tcv',1e-3,'bex',-1.5,'ucex',.8,'ibbt',9e-4,'avto',0,'akp',0,'agamma',0,...
    'kf',0,'af',0,'nqs',0,'satlim',exp(4),'xqc',0.4,'ekvdyn',0,'type',0,'update',26,...
    'rsh',0,'hdif',0);
h.device=struct('l',1e-6,'w',1e-6,'ns',1,'m',1);
[id,gm,gms,gmd,o]=ekv(h.model,h.device,1,0,1,25);
opt_names{1,1}='IDS';
opt_names{2,1}='gm';
opt_names{3,1}='gms';
opt_names{4,1}='gmd';
names=fieldnames(o);
for k=5:(4+length(names))
    opt_names{k,1}=names{k-4};
end
set(h.graphtype,'String',opt_names);
guidata(h.ekvint,h);
fix_sliders(h);
set(h.version,'String',num2str(h.ekvint_ver));

function fix_sliders(h)
p={'cox','xj','vto','gamma','phi','kp','e0','ucrit','dl','dw','lambda',...
    'leta','weta','q0','lk','iba','ibb','ibn','rsh','hdif'};
set(h.minmax,'String',upper(p));
for i=1:length(p)
    ma=getfield(h.max,p{i});
    mi=getfield(h.min,p{i});
    val=getfield(h.model,p{i});
    if i==7 | i==8
        set(getfield(h,p{i}),'Value',log(val/mi)/log(ma/mi));
    else
        set(getfield(h,p{i}),'Value',(val-mi)/(ma-mi));
    end
end
set(h.w,'String',num2str(h.device.w/1e-6));
set(h.l,'String',num2str(h.device.l/1e-6));
set(h.m,'String',num2str(h.device.m));
set(h.ns,'String',num2str(h.device.ns));
plotcall(h);
minmax_Callback(0,0,h,0);

% --------------------------------------------------------------------
function plotcall(h)
tt=get(h.Vsweep,'Value');
v=linspace(str2num(get(h.xmin,'String')),str2num(get(h.xmax,'String')),100)';
h.device.w=str2num(get(h.w,'String'))*1e-6;
h.device.l=str2num(get(h.l,'String'))*1e-6;
h.device.ns=str2num(get(h.ns,'String'));
h.device.m=str2num(get(h.m,'String'));
rs=h.model.rsh.*h.model.hdif./(h.device.m.*h.device.w);
[v,h.device.l]=matrix_voltages(v,h.device.l);
[v,h.device.w]=matrix_voltages(v,h.device.w);
[v,h.device.m]=matrix_voltages(v,h.device.m);
switch tt
case 1
    vs=str2num(get(h.vs,'String'));
    [v,vs]=matrix_voltages(v,vs);
    vd=str2num(get(h.vd,'String'));
    [v,vd]=matrix_voltages(v,vd);
    [id,gm,gms,gmd,o]=ekv(h.model,h.device,v,vs,vd,25);
case 2
    vs=str2num(get(h.vs,'String'));
    [v,vs]=matrix_voltages(v,vs);
    vg=str2num(get(h.vg,'String'));
    [v,vg]=matrix_voltages(v,vg);
    [id,gm,gms,gmd,o]=ekv(h.model,h.device,vg,vs,v,25);
case 3
    vd=str2num(get(h.vd,'String'));
    [v,vd]=matrix_voltages(v,vd);
    vg=str2num(get(h.vg,'String'));
    [v,vg]=matrix_voltages(v,vg);
    [id,gm,gms,gmd,o]=ekv(h.model,h.device,vg,v,vd,25);
end
axes(h.graph);
tt=get(h.graphtype,'Value');
if tt==1
    y=id;
elseif tt==2
    y=gm;
elseif tt==3
    y=gms;
elseif tt==4
    y=gmd;
else
    par_str=get(h.graphtype,'String');
    y=getfield(o,par_str{tt});
end
compare=get(h.compare,'Value'); %compare values
if compare==1
    fname=char(get(h.filename,'String'));
    if exist(fname)
        lo=load(fname);
    else
        warndlg(['File ''',fname,''' doesn''t exist!']);
        return;
    end
end
if get(h.diff,'Value')==1
    y2=interp1(lo(:,1),lo(:,2),v);
    if get(h.log,'Value')==1
        semilogy(v,abs(y-y2));
    else
        plot(v,(y-y2));
    end
else
    if get(h.log,'Value')==1
        semilogy(v,y);
    else
        plot(v,y);
    end
    if compare==1
        hold on
        plot(lo(:,1),lo(:,2),'r')
        hold off
    end
end
if get(h.grid,'Value')==1
    grid on;
end
guidata(h.ekvint,h)
str=get(h.minmax,'String');
n=get(h.minmax,'Value');
set(h.curr,'String',num2str(getfield(h.model,lower(str{n}))));
if get(h.fix,'Value')==1
    a=zeros(1,4);
    a(1)=str2num(get(h.xmin,'String'));
    a(2)=str2num(get(h.xmax,'String'));
    a(3)=str2num(get(h.ymin,'String'));
    a(4)=str2num(get(h.ymax,'String'));
    axis(a);
else
    a=axis;
    %set(h.xmin,'String',num2str(a(1)));
    %set(h.xmax,'String',num2str(a(2)));
    a(1)=str2num(get(h.xmin,'String'));
    a(2)=str2num(get(h.xmax,'String'));
    set(h.ymin,'String',num2str(a(3)));
    set(h.ymax,'String',num2str(a(4)));
    axis(a);
end
set(h.rs,'String',num2str(rs));



% --------------------------------------------------------------------
function varargout = cox_Callback(h, eventdata, handles, varargin)
handles.model.cox=get(handles.cox,'Value')*(handles.max.cox-handles.min.cox)+handles.min.cox;
plotcall(handles);

function varargout = xj_Callback(h, eventdata, handles, varargin)
handles.model.xj=get(handles.xj,'Value')*(handles.max.xj-handles.min.xj)+handles.min.xj;
plotcall(handles)

function varargout = vto_Callback(h, eventdata, handles, varargin)
handles.model.vto=get(handles.vto,'Value')*(handles.max.vto-handles.min.vto)+handles.min.vto;
plotcall(handles)

function varargout = gamma_Callback(h, eventdata, handles, varargin)
handles.model.gamma=get(handles.gamma,'Value')*(handles.max.gamma-handles.min.gamma)+handles.min.gamma;
plotcall(handles)

function varargout = phi_Callback(h, eventdata, handles, varargin)
handles.model.phi=get(handles.phi,'Value')*(handles.max.phi-handles.min.phi)+handles.min.phi;
plotcall(handles)

function varargout = kp_Callback(h, eventdata, handles, varargin)
handles.model.kp=get(handles.kp,'Value')*(handles.max.kp-handles.min.kp)+handles.min.kp;
plotcall(handles)

function varargout = e0_Callback(h, eventdata, handles, varargin)
handles.model.e0=exp(get(handles.e0,'Value')*(log(handles.max.e0)-log(handles.min.e0))+log(handles.min.e0));
plotcall(handles)

function varargout = ucrit_Callback(h, eventdata, handles, varargin)
handles.model.ucrit=exp(get(handles.ucrit,'Value')*(log(handles.max.ucrit)-log(handles.min.ucrit))+log(handles.min.ucrit));
plotcall(handles)

function varargout = dl_Callback(h, eventdata, handles, varargin)
handles.model.dl=get(handles.dl,'Value')*(handles.max.dl-handles.min.dl)+handles.min.dl;
plotcall(handles)

function varargout = dw_Callback(h, eventdata, handles, varargin)
handles.model.dw=get(handles.dw,'Value')*(handles.max.dw-handles.min.dw)+handles.min.dw;
plotcall(handles)

function varargout = lambda_Callback(h, eventdata, handles, varargin)
handles.model.lambda=get(handles.lambda,'Value')*(handles.max.lambda-handles.min.lambda)+handles.min.lambda;
plotcall(handles)

function varargout = leta_Callback(h, eventdata, handles, varargin)
handles.model.leta=get(handles.leta,'Value')*(handles.max.leta-handles.min.leta)+handles.min.leta;
plotcall(handles)

function varargout = weta_Callback(h, eventdata, handles, varargin)
handles.model.weta=get(handles.weta,'Value')*(handles.max.weta-handles.min.weta)+handles.min.weta;
plotcall(handles)

function varargout = q0_Callback(h, eventdata, handles, varargin)
handles.model.q0=get(handles.q0,'Value')*(handles.max.q0-handles.min.q0)+handles.min.q0;
plotcall(handles)

function varargout = lk_Callback(h, eventdata, handles, varargin)
handles.model.lk=get(handles.lk,'Value')*(handles.max.lk-handles.min.lk)+handles.min.lk;
plotcall(handles)

function varargout = iba_Callback(h, eventdata, handles, varargin)
handles.model.iba=get(handles.iba,'Value')*(handles.max.iba-handles.min.iba)+handles.min.iba;
plotcall(handles)

function varargout = ibb_Callback(h, eventdata, handles, varargin)
handles.model.ibb=get(handles.ibb,'Value')*(handles.max.ibb-handles.min.ibb)+handles.min.ibb;
plotcall(handles)

function varargout = ibn_Callback(h, eventdata, handles, varargin)
handles.model.ibn=get(handles.ibn,'Value')*(handles.max.ibn-handles.min.ibn)+handles.min.ibn;
plotcall(handles)

function varargout = rsh_Callback(h, eventdata, handles, varargin)
handles.model.rsh=get(handles.rsh,'Value')*(handles.max.rsh-handles.min.rsh)+handles.min.rsh;
plotcall(handles)

function varargout = hdif_Callback(h, eventdata, handles, varargin)
handles.model.hdif=get(handles.hdif,'Value')*(handles.max.hdif-handles.min.hdif)+handles.min.hdif;
plotcall(handles)


% --------------------------------------------------------------------
function varargout = minmax_Callback(h, eventdata, handles, varargin)
str=get(handles.minmax,'String');
n=get(handles.minmax,'Value');
set(handles.minb,'String',num2str(getfield(handles.min,lower(str{n}))));
set(handles.maxb,'String',num2str(getfield(handles.max,lower(str{n}))));
set(handles.curr,'String',num2str(getfield(handles.model,lower(str{n}))));



% --------------------------------------------------------------------
function varargout = mostype_Callback(h, eventdata, handles, varargin)
val=get(handles.mostype,'Value');
if val==1
    handles.model.type=1;
else
    handles.model.type=2;
end
plotcall(handles);

% --------------------------------------------------------------------
function varargout = chgminmax_Callback(hcc, eventdata, h, varargin)
str=get(h.minmax,'String');
n=get(h.minmax,'Value');
mi=str2num(get(h.minb,'String'));
ma=str2num(get(h.maxb,'String'));
h.min=setfield(h.min,lower(str{n}),min(mi,getfield(h.model,lower(str{n}))));
h.max=setfield(h.max,lower(str{n}),max(ma,getfield(h.model,lower(str{n}))));
p=lower(str{n});
mi=getfield(h.min,p);
ma=getfield(h.max,p);
val=getfield(h.model,p);
if n==7|n==8
    set(getfield(h,p),'Value',log(val/mi)/log(ma/mi));
else
    set(getfield(h,p),'Value',(val-mi)/(ma-mi));
end
guidata(h.ekvint,h)
minmax_Callback(hcc, eventdata, h, varargin)


% --------------------------------------------------------------------
function varargout = curr_Callback(hcc, eventdata, h, varargin)
str=get(h.minmax,'String');
n=get(h.minmax,'Value');
v=str2num(get(h.curr,'String'));
p=lower(str{n});
h.model=setfield(h.model,p,max(min(v,getfield(h.max,p)),getfield(h.min,p)));
mi=getfield(h.min,p);
ma=getfield(h.max,p);
val=getfield(h.model,p);
if n==7|n==8
    set(getfield(h,p),'Value',log(val/mi)/log(ma/mi));
else
    set(getfield(h,p),'Value',(val-mi)/(ma-mi));
end
plotcall(h);



% --------------------------------------------------------------------
function varargout = graphtype_Callback(h, eventdata, handles, varargin)
n=get(handles.Vsweep,'Value');
switch n
case 1
    set(handles.vg,'Enable','Off');
    set(handles.vd,'Enable','On');
    set(handles.vs,'Enable','On');
case 2
    set(handles.vd,'Enable','Off');
    set(handles.vg,'Enable','On');
    set(handles.vs,'Enable','On');
case 3
    set(handles.vs,'Enable','Off');
    set(handles.vg,'Enable','On');
    set(handles.vd,'Enable','On');   
end
plotcall(handles);

% --------------------------------------------------------------------
function varargout = w_Callback(h, eventdata, handles, varargin)
plotcall(handles);


% --------------------------------------------------------------------
function varargout = vgsd_Callback(h, eventdata, handles, varargin)
plotcall(handles);



% --------------------------------------------------------------------
function varargout = compare_Callback(h, eventdata, handles, varargin)
c=get(handles.compare,'Value');
if c==1
    set(handles.filename,'Enable','On');
    set(handles.diff,'Enable','On');
    set(handles.fit,'Enable','On');
    set(handles.comp_button,'Enable','On');
else
    set(handles.filename,'Enable','Off');
    set(handles.fit,'Enable','Off');
    set(handles.diff,'Enable','Off');
    set(handles.diff,'Value',0);
    set(handles.comp_button,'Enable','Off');
end
plotcall(handles);


% --------------------------------------------------------------------
function varargout = fix_Callback(h, eventdata, handles, varargin)
if get(handles.fix,'Value')==1
    %set(handles.xmin,'Enable','On');
    %set(handles.xmax,'Enable','On');
    set(handles.ymin,'Enable','On');
    set(handles.ymax,'Enable','On');
else
    %set(handles.xmin,'Enable','Off');
    %set(handles.xmax,'Enable','Off');
    set(handles.ymin,'Enable','Off');
    set(handles.ymax,'Enable','Off');
end
plotcall(handles);

function [v,vb]=matrix_voltages(v,vb);
[h,w]=size(v);
[h2,w2]=size(vb);
if w>1
    if w2==w
        if h2==1
            vb=ones(h,1)*vb;
        end
    else
        if w2>1
            errordlg('Voltage vectors must have the same length!');
        end
    end
else
    if w2>1
        if w==1
            v=v*ones(1,w2);
        end
        if h2==1
            vb=ones(h,1)*vb;
        end
    end
end




% --------------------------------------------------------------------
function varargout = load_Callback(h, eventdata, handles, varargin)
[f,p]=uigetfile('*.lib','Select a model file');
if f==0
    return
end
mod=read_model([p,f]);
if isstruct(mod)==0
    warndlg('Error reading file!');
    return;
end
mod
handles.model=mod;
fix_sliders(handles);

% --------------------------------------------------------------------
function varargout = save_Callback(h, eventdata, handles, varargin)
[f,p]=uiputfile('*.lib','Select a file name to save model');
if f==0
    return;
end
if length(f)<4
    f=[f,'.lib'];
else
    if strcmp(f((length(f)-3):length(f)),'.lib')==0
       f=[f,'.lib'];
   end
end
save_model(handles.model,[p,f]);




% --------------------------------------------------------------------
function varargout = comp_button_Callback(h, eventdata, handles, varargin)
[f,p]=uigetfile('*.dat','Select a data file');
if f~=0
    set(handles.filename,'String',[p,f]);
    plotcall(handles);
end




% --------------------------------------------------------------------
function varargout = fit_Callback(hprimo, eventdata, h, varargin)
s=get(h.minmax,'String');
s=lower(s{get(h.minmax,'Value')});
val=getfield(h.model,s);
var=linspace(val*.9,val*1.1,100);
tt=get(h.Vsweep,'Value');
v=linspace(str2num(get(h.xmin,'String')),str2num(get(h.xmax,'String')),100)';
v=v*ones(1,size(var,2));
var=ones(size(v,1),1)*var;
h.model=setfield(h.model,s,var);
switch tt
case 1
    vs=str2num(get(h.vs,'String'));
    vd=str2num(get(h.vd,'String'));
    [id,gm,gms,gmd,o]=ekv(h.model,h.device,v,vs,vd,25);
case 2
    vs=str2num(get(h.vs,'String'));
    vg=str2num(get(h.vg,'String'));
    [id,gm,gms,gmd,o]=ekv(h.model,h.device,vg,vs,v,25);
case 3
    vd=str2num(get(h.vd,'String'));
    vg=str2num(get(h.vg,'String'));
    [id,gm,gms,gmd,o]=ekv(h.model,h.device,vg,v,vd,25);
end
tt=get(h.graphtype,'Value');
if tt==1
    y=id;
elseif tt==2
    y=gm;
elseif tt==3
    y=gms;
elseif tt==4
    y=gmd;
else
    par_str=get(h.graphtype,'String');
    y=getfield(o,par_str{tt});
end
compare=get(h.compare,'Value'); %compare values
if compare==1
    fname=char(get(h.filename,'String'));
    if exist(fname)
        lo=load(fname);
    else
        warndlg(['File ''',fname,''' doesn''t exist!']);
        return;
    end
end
y2=interp1(lo(:,1),lo(:,2),v(:,1),'cubic','extrap');
y2=y2*ones(1,size(var,2));
if get(h.logfit,'Value')==0
    diff=sum((y2-y).^2);
else
    diff=sum((log(y2+1e-16)-log(y+1e-16)).^2);
end
val=var(1,find(diff==min(diff)));
h.model=setfield(h.model,s,val(1));
fix_sliders(h);    
    




% --------------------------------------------------------------------
function varargout = logfit_Callback(h, eventdata, handles, varargin)
return;




