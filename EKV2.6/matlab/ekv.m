function [IDS,gm,gms,gmd,other]=ekv(ekvmodel,device,VG,VS,VD,TC);
%EKV v2.6
%input: model,device,vg,vs,vd,tc
%output: id,gm,gms,gmd,other

T0C=273;
T=T0C+TC;
k=1.381e-23;
q=1.609e-19;


% Device geometry
W=device.w;
L=device.l;
M=device.m;
NS=device.ns;
Rd=0;
Rs=0;


% Model parameters
UPDATE=ekvmodel.update;
NQS=ekvmodel.nqs;
SATLIM=ekvmodel.satlim;
XQC=ekvmodel.xqc;
EKVINT=ekvmodel.ekvint;
EKVDYN=ekvmodel.ekvdyn;
TYPE=ekvmodel.type;
COX=ekvmodel.cox;
XJ=ekvmodel.xj;
DW=ekvmodel.dw;
DL=ekvmodel.dl;
VTO=ekvmodel.vto;
GAMMA=ekvmodel.gamma;
PHI=ekvmodel.phi;
KP=ekvmodel.kp;
E0=ekvmodel.e0;
THETA=ekvmodel.theta;
UCRIT=ekvmodel.ucrit;
LAMBDA=ekvmodel.lambda;
WETA=ekvmodel.weta;
LETA=ekvmodel.leta;
Q0=ekvmodel.q0;
LK=ekvmodel.lk;
IBA=ekvmodel.iba;
IBB=ekvmodel.ibb;
IBN=ekvmodel.ibn;
TCV=ekvmodel.tcv;
BEX=ekvmodel.bex;
UCEX=ekvmodel.ucex;
IBBT=ekvmodel.ibbt;
AVTO=ekvmodel.avto;
AKP=ekvmodel.akp;
AGAMMA=ekvmodel.agamma;
KF=ekvmodel.kf;
AF=ekvmodel.af;
if isfield(ekvmodel,'hdif')&isfield(ekvmodel,'rsh')
    Rs=ekvmodel.hdif.*ekvmodel.rsh./(M.*W);
    Rd=Rs;
end

epsilonSi=104.5e-12;
Cepsilon=4*(22e-3)^2;
CA=0.028;
Vt=k*T/q;
Tref=300;
EgTr=1.16-.000702*Tref^2/(Tref+1108);
EgT=1.16-.000702*T^2/(T+1108);

%Temperature effects
VTO=VTO-TCV*(T-Tref);
KP=KP*(T/Tref)^BEX;
UCRIT=UCRIT*(T/Tref)^UCEX;
PHI=PHI*T/Tref-3*Vt*log(T/Tref)-EgTr*T/Tref+EgT;
IBB=IBB*(1+IBBT*(T-Tref));


if TYPE==1.0
    eta=0.5;
else
    eta=1/3;
end

%Effective channel length & width

Weff=W+DW;
Leff=L+DL;

area=sqrt(Weff.*Leff.*NS.*M);
VTOa=VTO+AVTO./area;
KPa=KP*(1+AKP./area);
GAMMAa=GAMMA+AGAMMA./area;

%define some variables
qB0=GAMMAa.*sqrt(PHI);
COX_eSi=COX/epsilonSi;
LETA_L=LETA./Leff;
WETA_W=3*WETA./Weff;

%reverse short channel

csi=CA*(10*Leff/LK-1);
DeltaVRSCE=2*Q0./(COX.*(1+0.5*(csi+sqrt(csi.^2+Cepsilon))).^2);

%Effective gate voltage including RSCE

VG_p=VG-VTOa-DeltaVRSCE+PHI+qB0;

%pinch off voltage for narrow channel

%if VG_p>0
%    VP0=VG_p-PHI-GAMMAa.*(sqrt(VG_p+(GAMMAa/2)^2)-GAMMAa/2);
%else
%    VP0=-PHI;
%end
ch=VG_p>0;
VP0=ch.*(VG_p-PHI-GAMMAa.*(sqrt(VG_p+(GAMMAa/2).^2)-GAMMAa/2));
VP0=VP0+(1-ch)*(-PHI);

%Effective substrate factor accounting for charge sharing

VS_p=0.5*(VS+PHI+sqrt((VS+PHI).^2+(4*Vt)^2));
VD_p=0.5*(VD+PHI+sqrt((VD+PHI).^2+(4*Vt)^2));

gamma_circ=GAMMAa-1./COX_eSi.*(LETA_L.*(sqrt(VS_p)+sqrt(VD_p))-WETA_W.*sqrt(VP0+PHI+0.1*Vt));

%Effective substrate factor, including charge-sharing for short & narrow channels

gamma_p=0.5*(gamma_circ+sqrt(gamma_circ.^2+0.1*Vt));

%pinch-off voltage, including short & narrow channels effects

%if VG_p>0
%    VP_sqrt=sqrt(VG_p+(gamma_p/2).^2);
%    VP=VG_p-PHI-gamma_p.*(VP_sqrt-gamma_p/2);
%else
%    VP_sqrt=(gamma_p/2);
%    VP=-PHI;
%end
VP_sqrt=ch.*sqrt(VG_p+(gamma_p/2).^2);
VP=ch.*(VG_p-PHI-gamma_p.*(VP_sqrt-gamma_p/2));
VP_sqrt=VP_sqrt+(1-ch).*(gamma_p/2);
VP=VP+(1-ch).*(-PHI);

%forward normalized current

[i_f,D_F_if]=F((VP-VS)/Vt);
[ir,D_F_ir]=F((VP-VD)/Vt);

%Velocity saturation voltage

VC=UCRIT.*NS.*Leff;
VDSS=VC.*(sqrt(1/4+Vt*sqrt(i_f)./VC)-1/2);

%channel length modulation

DeltaV=4*Vt*sqrt(LAMBDA*(sqrt(i_f)-VDSS/Vt)+1/64);
Vds=(VD-VS)/2;
Vip=sqrt(VDSS.^2+DeltaV.^2)-sqrt((Vds-VDSS).^2+DeltaV.^2);
LC=sqrt(XJ./COX_eSi);
DeltaL_log=1+(Vds-Vip)./(LC.*UCRIT);
DeltaL=LAMBDA.*LC.*log(DeltaL_log);

%Equivalent channel length including channel length modulation ans velocity saturation

L_p=NS.*Leff-DeltaL+(Vds+Vip)./UCRIT;
Lmin=NS.*Leff/10;
Leq=(L_p+sqrt(L_p.^2+Lmin.^2))/2;

%Drain-to-source saturation voltage for reverse normalized current

VDSS_p=VC.*(sqrt(0.25+Vt./VC.*(sqrt(i_f)-3/4*log(i_f)))-1/2)+Vt*(log(VC/(2*Vt))-.6);

%reverse normalized current

irp_sqrt1=sqrt(VDSS_p.^2+DeltaV.^2);
irp_sqrt2=sqrt((Vds-VDSS_p).^2+DeltaV.^2);
[ir_p,D_F_ir_p]=F((VP-Vds-VS-irp_sqrt1+irp_sqrt2)/Vt);

%normalized intrinsic node charges

nq=1+GAMMAa./(2*sqrt(VP+PHI+1e-6));
xf=sqrt(1/4+i_f);
xr=sqrt(1/4+ir);
xs=xf+xr;
qI=-nq.*(4/3*(xs-xf.*xr./xs)-1);
%if VG_p>0
%    qB=(-GAMMAa*sqrt(VP+PHI+1e-6))*1/Vt-((nq-1)./nq).*qI;
%else
%    qB=-VG_p/Vt;
%end
qB=ch.*((-GAMMAa.*sqrt(VP+PHI+1e-6))*1/Vt-((nq-1)./nq).*qI);
qB=qB+(1-ch).*(-VG_p/Vt);
qBqI=qB+eta*qI;

%Slope factor

VP_PHI_4Vt=VP+PHI+4*Vt;
n=1+GAMMAa./(2*sqrt(VP_PHI_4Vt));

%transconductance factor and mobility reduction due to vertical field

beta0=KPa.*Weff.*M./Leq;
if THETA==0
    beta0_p=beta0.*(1+COX_eSi.*qB0./E0);
    beta=beta0_p./(1+COX_eSi./E0.*Vt.*abs(qBqI));
else
    VP_p=1/2*(VP+sqrt(VP.^2+2*Vt^2));
    beta=beta0./(1+THETA*VP_p);
end

%Specific current

IS=2*n.*beta*Vt^2;

%Drain-to-source current

IDS=IS.*(i_f-ir_p);

%Impact ionization current

Vib=VD-VS-IBN*2*VDSS;
ch=Vib>0;
other.IDB=ch.*IDS*IBA/IBB.*Vib.*exp(-IBB*LC./Vib);

%evaluate gm: d(IDS)/d(Vg)
%         gms:d(IDS)/d(Vs)
%         gmd:d(IDS)/d(Vd)
for i=1:3
    D_VG = i==1;
    D_VS = i==2;
    D_VD = i==3;
    D_VG_p = D_VG;
    D_VD_p = D_VD/2+1/4*(2*D_VD*VD+2*D_VD*PHI)./sqrt((VD+PHI).^2+16*Vt^2);
    if VG_p>0
        D_VP0 = D_VG_p-GAMMAa*D_VG_p./sqrt(4*VG_p+GAMMAa.^2);
    else
        D_VP0 = 0;
    end
    D_VS_p = 1/2*D_VS+1/4*(2*D_VS.*(VS+PHI))./sqrt((VS+PHI).^2+16*Vt^2);
    D_gamma_circ = -1./COX_eSi.*(LETA_L.*(1/2*D_VS_p./sqrt(VS_p)+1/2*D_VD_p./sqrt(VD_p))-WETA_W/2.*D_VP0./sqrt(VP0+PHI+Vt/10));
    D_gamma_p = D_gamma_circ.*gamma_p./(2*gamma_p-gamma_circ);
    if VG_p>0
        D_VP = D_VG_p-D_gamma_p.*(VP_sqrt-gamma_p/2)-gamma_p.*(1/2*(D_VG_p+D_gamma_p.*gamma_p/2)./VP_sqrt-D_gamma_p/2);
    else
        D_VP=0;
    end
    D_n = D_VP.*(1-n)./(2*VP_PHI_4Vt);
    D_if = D_F_if.*(D_VP-D_VS)/Vt;
    D_VDSS = 1/2*VDSS.*(VDSS+VC).*D_if./(i_f.*(2*VDSS+VC));
    D_VDSS_p = Vt*(1/2*D_if./sqrt(i_f)-3/4*1./i_f.*D_if)./sqrt(1+4*Vt./VC.*(sqrt(i_f)-3/4*log(i_f)));
    D_DeltaV = 16*Vt*LAMBDA*(1/2*D_if./sqrt(i_f)-D_VDSS/Vt)./sqrt(64*LAMBDA*(sqrt(i_f)-VDSS/Vt)+1);
    D_Vds = (D_VD-D_VS)/2;
    D_Vip = (D_VDSS.*VDSS+D_DeltaV.*DeltaV)./sqrt(VDSS.^2+DeltaV.^2)-((D_Vds-D_VDSS).*(Vds-VDSS)+D_DeltaV.*DeltaV)./sqrt((Vds-VDSS).^2+DeltaV.^2);
    D_DeltaL = LAMBDA./DeltaL_log.*(D_Vds-D_Vip)./UCRIT; 
    D_L_p = -D_DeltaL+(D_Vds+D_Vip)./UCRIT;
    D_Leq = 1/2*D_L_p+5*D_L_p.*L_p./sqrt(100*L_p.^2+Leff.^2);
    D_ir_p=D_F_ir_p.*(D_VP-D_Vds-D_VS-(D_VDSS_p.*VDSS_p+D_DeltaV.*DeltaV)./irp_sqrt1+((D_Vds-D_VDSS_p).*(Vds-VDSS_p)+D_DeltaV.*DeltaV)./irp_sqrt2)/Vt;
    D_ir = D_F_ir.*(D_VP-D_VD)/Vt;
    D_beta0 = -beta0.*D_Leq./Leq;
    
    D_nq = -1/4*GAMMAa.*D_VP./(VP+PHI+.1e-5).^(3/2);
    D_xf = D_if./(2*xf);
    D_xr = D_ir./(2*xr);
    D_xs = D_xf+D_xr;
    D_qI = -D_nq.*(4/3*(xf.^2+xf.*xr+xr.^2)./(xf+xr)-1)-nq.*(4/3*(D_xs.*(1+xr.*xf./xs.^2)-(D_xr.*xf+D_xf.*xr)./xs));
    if VG_p>0
        D_qB = -1/2*GAMMAa.*D_VP./(sqrt(VP+PHI+.1e-5)*Vt)-D_nq./nq.*qI+(nq-1).*D_nq./nq.^2.*qI-(nq-1)./nq.*D_qI;
    else
        D_qB = 0;
    end
    D_qBqI=sign(qBqI).*(D_qB+eta*D_qI);
    
    if THETA==0
        D_beta0_p = beta0_p./beta0.*D_beta0;
        D_beta = beta./beta0.*(D_beta0+Vt*(beta-beta0).*D_qBqI./(Vt*abs(qBqI)-qB0));
    else
        D_VP_p=.5+D_VP/sqrt(VP^2+2*Vt^2);
        D_beta=D_beta0./(1+THETA*VP_p)-beta0*THETA.*D_VP_p./(1+THETA*VP_p).^2;
    end
    D_IS = 2*D_n.*beta*Vt^2+2*n.*D_beta*Vt^2;
    D_IDS = D_IS.*(i_f-ir_p)+IS.*(D_if-D_ir_p);
    switch i
    case 1
        gm = D_IDS;
    case 2
        gms = -D_IDS;
    case 3
        gmd = D_IDS;
    end
end
other.IDS0=IDS;
other.gm0=gm;
other.gms0=gms;
other.gmd0=gmd;
other.gmbs=gms-gm-gmd;


%Intrinsic noise model equations
%thermal noise

Cox=COX.*Weff.*Leff;
QI=Cox.*Vt.*qI;
%Sthermal=4*k*T*KPa./(COX*NS.*Leff.^2).*abs(QI);
Sthermal=4*k*T*beta.*abs(qI);

%flicker noise

%Sflicker=KF*gm.^2./(Weff.*M.*NS.*Leff*COX.*f.^AF);
other.SINDS=Sthermal;

%assign other informations

other.n=n;
other.mode=(i_f./ir)>SATLIM;
other.IS=IS;
other.VM=IDS./gmd;
other.if_=i_f;
other.ir=ir;
other.ir_p=ir_p;
other.VDSAT=2*VDSS+4*Vt;
other.Vp=VP;
other.VTH=VTOa+DeltaVRSCE+gamma_p.*sqrt(VS_p)-qB0;
other.beta=beta;
other.qI=qI;
other.qB=qB;

%Correction for gate access resistor (Cf. PhD thesis of Matthias Bucher)
krser=1+Rs.*gms+Rd.*gmd;
IDS=IDS./krser;
gm=gm./krser;
gms=gms./krser;
gmd=gmd./krser;

%rest
%other.cgs=2/3*(1-(xr.^2+xr+xf/2)./xs.^2);
%other.cgd=2/3*(1-(xf.^2+xf+xr/2)./xs.^2);
%other.cgb=(nq-1)./nq.*(1-other.cgs-other.cgd);
%other.csb=(nq-1).*other.cgs;
%other.cdb=(nq-1).*other.cgd;
%Cox=COX.*Weff.*Leff.*M.*NS;
%tau0=Cox./(2.*Vt.*beta);
%other.tau0=tau0;
%other.tau=tau0*4/15.*(xf.^2+3*xf.*xr+xr.^2)./xs.^3;

%Large signal interpolation function

function [i,d]=F(v);
ch1=v>(-15);
ch2=v<(-.35);
z0=ch2.*(1.55*exp(-v));
z0=z0+(1-ch2).*(2./(1.3+v-log(v+1.6)));
D_z0=ch2.*(-exp(-v));
D_z0=D_z0+(1-ch2).*(-2./(1.3+v-log(v+1.6)).^2.*(1-1./(v+1.6)));
z1=(2+z0)./(1+v+log(z0));
D_z1=D_z0./(1+v+log(z0))-(2+z0)./(1+v+log(z0)).^2.*(1+1./z0.*D_z0);
y=ch1.*((1+v+log(z1))./(2+z1));
D_y=ch1.*((1+1./z1.*D_z1)./(2+z1)-(1+v+log(z1))./(2+z1).^2.*D_z1);
y=y+(1-ch1).*(1./(2+exp(-v)));
D_y=D_y+(1-ch1).*(1./(2+exp(-v)).^2.*exp(-v)); 
%if v>(-15)
%    if v<(-.35)
%        z0=1.55*exp(-v);
%        D_z0=-exp(-v);
%    else
%        z0=2./(1.3+v-log(v+1.6));
%        D_z0=-2./(1.3+v-log(v+1.6)).^2.*(1-1./(v+1.6));
%    end
%    z1=(2+z0)./(1+v+log(z0));
%    D_z1=D_z0./(1+v+log(z0))-(2+z0)./(1+v+log(z0)).^2.*(1+1./z0.*D_z0);
%    y=(1+v+log(z1))./(2+z1);
%    D_y=(1+1./z1.*D_z1)./(2+z1)-(1+v+log(z1))./(2+z1).^2.*D_z1;
%else
%    y=1./(2+exp(-v));
%    D_y=1./(2+exp(-v)).^2.*exp(-v); 
%end
i=y.*(1+y);
d=D_y.*(1+y)+y.*D_y;
