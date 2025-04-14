clear all
close all

% Constants
F = 96485.3415;     % coulomb_per_mole (in model_parameters)
R = 8.314472;       % joule_per_mole_kelvin (in model_parameters)
T = 310.0;          % kelvin (in model_parameters) %37°C

% Cell geometry
V_SR = 465.199;        % micrometre_cube (in model_parameters)
Vc   = 7012.0;        % micrometre_cube (in model_parameters)
Cm   = 7.8667e-11;   % farad (in model_parameters)

% Extracellular concentrations
Nao = 154.0; % millimolar (in model_parameters)
%Ko  = 5.4;   % millimolar (in model_parameters)
Ko  = 4;   % millimolar (in model_parameters)
Cao = 2.0;   % millimolar (in model_parameters)

% Intracellular concentrations
Ki = 150.0;   % millimolar (in model_parameters)


% PER PASSARE DA PACED A UNPACED CAMBIARE
% * gK1
% * condizioni iniziali
% * istim

%g_K1 = 0.06; % unpaced
g_K1 = 0.169261259039458; % paced




g_Na        = 9800.86868852711;
myCoefTauM  = 1;
tauINaL     = 200; 
GNaLmax     = 13.3509157535054; 
Vh_hLate    = 87.61;
g_f         = 25.9541335109086;
fNa         = 0.37;
fK          = 1 - fNa;
g_CaL       = 7.44995532741049e-05;
constf2     = 1.0;
g_to        = 59.4564976228642;  
g_Ks        = 2.08556136889656;
L0          = 0.025;  
Q           = 2.3;    
g_Kr        = 16.8722;   
V_half      = 1000.0*(-R*T/(F*Q)*log((1.0+Cao/2.6)^4.0/(L0*(1.0+Cao/0.58)^4.0))-0.019);
%g_K1       = 28.1492;  
KmCa        = 1.38;  
KmNai       = 87.5;  
Ksat        = 0.1;    
gamma       = 0.35;  
alpha       = 2.16659;
kNaCa       = 3450.73331808885; 
Km_K        = 1.0;  
Km_Na       = 40.0; 
PNaK        = 2.27178916268546;
KPCa        = 0.0005;   
g_PCa       = 0.456955096510698;  
g_b_Na      = 1.14;
g_b_Ca      = 0.8727264; 
VmaxUp		= 0.3226;
Kup			= 4.40435e-4;
V_leak		= 4.48209e-4;
g_irel_max	= 75.4190;
RyRa1       = 0.1034;
RyRa2       = 0.050001;
RyRahalf    = 0.02632;
RyRohalf    = 0.00944;
RyRchalf    = 0.00167;
RyRtauadapt = 1; 
Buf_C       = 0.25;   % millimolar (in calcium_dynamics)
Buf_SR      = 10.0;   % millimolar (in calcium_dynamics)
Kbuf_C      = 0.001;   % millimolar (in calcium_dynamics)
Kbuf_SR     = 0.3;   % millimolar (in calcium_dynamics)
KQ10 = 3.0;
coeff_Kur = 3.5;
g_KCa = 0.0753851103617082;
KCa_on = 47.0e6;
KCa_off = 13.0;



% condizioni iniziali paced
U_rest=-0.0908810000000000;
CaSR0=0.0374480000000000;
Cai0=3.65880000000000e-05;
g0=0;
d0=5.33930000e-06;
f10=0.925250000000000;
f20=1;
fCa0=0.997840000000000;
Xr10=4.10050000000000e-07;
Xr20=0.514850000000000;
Xs0=0.0117060000000000;
h0=0.973370000000000;
j0=0.974170000000000;
m0=0.00958940000000000;
Xf0=0.512630000000000;
q0=0.948970000000000;
r0=0.00237220000000000;
Nai0=14.9650000000000;
mL0=0.000108850000000000;
hL0=0.625550000000000;
RyRa0=0.0981820000000000;
RyRo0=3.62930000000000e-11;
RyRc0=0.984730000000000;
ua0=0.297490000000000;
ui0=0.995440000000000;
O0=0.00499310000000000;


% condizioni iniziali unpaced
% U_rest=-0.081651;
% CaSR0=0.039624;
% Cai0=3.4412e-05;
% g0=0;
% d0=2.021e-05;
% f10=0.63335;
% f20=0.94148;
% fCa0=0.99796;
% Xr10=0.017269;
% Xr20=0.4683;
% Xs0=0.02077;
% h0=0.90275;
% j0=0.84352;
% m0=0.021708;
% Xf0=0.14923;
% q0=0.90058;
% r0=0.0038953;
% Nai0=8.6397;
% mL0=0.00062885;
% hL0=0.22058;
% RyRa0=0.091555;
% RyRo0=2.4012e-10;
% RyRc0=0.87108;
% ua0=0.37035;
% ui0=0.99279;
% O0=0.0068298;



% Nernst potential
E_Na =@(Nai) R*T/F*log(Nao/Nai);
E_Ca =@(Cai) 0.5*R*T/F*log(Cao/Cai);
E_K  = R*T/F*log(Ko/Ki);
PkNa = 0.03;   % dimensionless (in electric_potentials)
E_Ks =@(Nai) R*T/F*log((Ko+PkNa*Nao)/(Ki+PkNa*Nai));


% funzioni
m_inf       =@(V) 1 / (1 + exp((V*1000 + 39)/-11.2));
tau_m       =@(V) (0.00001 + 0.00013*exp(-((V*1000 + 48)/15)^2) + 0.000045 / (1 + exp((V*1000 + 42)/-5)));
h_inf       =@(V) 1 / (1 + exp((V*1000 + 66.5)/6.8));
tau_h       =@(V) (0.00007 + 0.034 / (1 + exp((V*1000 + 41)/5.5) + exp(-(V*1000 + 41)/14)) + 0.0002 / (1 + exp(-(V*1000 + 79)/14)));
j_inf       = h_inf;
tau_j       =@(V) 10*(0.0007 + 0.15 / (1 + exp((V*1000 + 41)/5.5) + exp(-(V*1000 + 41)/14)) + 0.002 / (1 + exp(-(V*1000 + 79)/14)));

i_Na        =@(V,m,h,j,Nai)  g_Na*m^3.0*h*j*(V - E_Na(Nai));


mL_inf     =@(V) 1/(1+exp(-(V*1000+42.85)/(5.264)));
alpha_m_L   =@(V) 1/(1+exp((-60-V*1000)/5));
beta_m_L    =@(V) 0.1/(1+exp((V*1000+35)/5))+0.1/(1+exp((V*1000-50)/200));
tau_mL     =@(V) 1/1000 * myCoefTauM*alpha_m_L(V)*beta_m_L(V);
hL_inf     =@(V) 1/(1+exp((V*1000+Vh_hLate)/(7.488)));
tau_hL     = 1/1000 * tauINaL;

i_NaL       =@(V,mL,hL,Nai) GNaLmax* mL^(3)*hL*(V-E_Na(Nai));

Xf_inf      =@(V) 1.0/(1.0 + exp((V*1000 + 69)/8));
tau_Xf      =@(V) (5600 / (1 + exp((V*1000 + 65)/7) + exp(-(V*1000 + 65)/19)))/1000;

i_fK        =@(V,Xf) fK*g_f*Xf*(V - E_K);
i_fNa       =@(V,Xf,Nai) fNa*g_f*Xf*(V - E_Na(Nai));
i_f         =@(V,Xf,Nai) i_fK(V,Xf) + i_fNa(V,Xf,Nai);

d_inf       =@(V) 1.0/(1.0+exp(-(V*1000.0+5.9860)/7.0));
%d_inf       =@(V) 1.0/(1.0+exp(-(V*1000.0+9.1)/7.0));
alpha_d     =@(V) 0.25+1.4/(1.0+exp((-V*1000.0-35.0)/13.0));
beta_d      =@(V) 1.4/(1.0+exp((V*1000.0+5.0)/5.0));
gamma_d     =@(V) 1.0/(1.0+exp((-V*1000.0+50.0)/20.0));
tau_d       =@(V) (alpha_d(V)*beta_d(V)+gamma_d(V))*1.0/1000.0;

f1_inf      =@(V) 1.0/(1.0+exp((V*1000.0+25.2260)/3.0));
constf1     =@(V,Cai,f1) ((1.0+1433.0*(Cai-50.0*1.0e-6))*1.35)*(f1_inf(V)-f1 > 0.0) + 1.0*(f1_inf(V)-f1 < 0.0);
tau_f1      =@(V,Cai,f1) (20.0+1102.5*exp(-((V*1000.0+50.0)/15.0)^2.0)+200.0/(1.0+exp((13.0-V*1000.0)/10.0))+280.0/(1.0+exp((30.0+V*1000.0)/10.0)))*constf1(V,Cai,f1)/1000.0;

f2_inf      =@(V) 0.33+0.67/(1.0+exp((V*1000.0+31.2260)/4.0));
tau_f2      =@(V) ((600.0*exp(-(V*1000.0+50.0)^2.0/400.0)+31.0/(1.0+exp((25.0-V*1000.0)/10.0))+1.0/(1.0+exp((30.0+V*1000.0)/10.0)))*constf2/1000.0)*2.0;

alpha_fCa   =@(Cai) 1.0/(1.0+(Cai/0.0006)^8.0);
beta_fCa    =@(Cai) 0.1/(1.0+exp((Cai-0.0009)/0.0001));
gamma_fCa   =@(Cai) 0.3/(1.0+exp((Cai-0.00075)/0.0008));
fCa_inf     =@(Cai) (alpha_fCa(Cai)+beta_fCa(Cai)+gamma_fCa(Cai))/1.3156;
constfCa    =@(V,Cai,fCa) 0.0*((V > -0.06) && (fCa_inf(Cai) > fCa)) + 1.0*((V < -0.06) || (fCa_inf(Cai) > fCa));
tau_fCa     =@(V,Cai,fCa) 0.002/constfCa(V,Cai,fCa);   % second (in i_CaL_fCa_gate)

i_CaL       =@(V,Cai,d,f1,f2,fCa) g_CaL*4.0*V*F^2.0/(R*T)*(Cai*exp(2.0*V*F/(R*T))-0.341*Cao)/(exp(2.0*V*F/(R*T))-1.0)*d*f1*f2*fCa;

q_inf       =@(V) 1.0/(1.0+exp((V*1000.0+53.0)/13.0));
tau_q       =@(V) (6.06+39.102/(0.57*exp(-0.08*(V*1000.0+44.0))+0.065*exp(0.1*(V*1000.0+45.93))))/1000.0;
r_inf       =@(V) 1.0/(1.0+exp(-(V*1000.0-22.3)/18.75));
tau_r       =@(V) (2.75352+14.40516/(1.037*exp(0.09*(V*1000.0+30.61))+0.369*exp(-0.12*(V*1000.0+23.84))))/1000.0;

i_to        =@(V,q,r) g_to*(V-E_K)*q*r;

Xs_inf      =@(V) 1.0/(1.0+exp((-V*1000.0-20.0)/16.0));
alpha_Xs    =@(V) 1100.0/sqrt(1.0+exp((-10.0-V*1000.0)/6.0));
beta_Xs     =@(V) 1.0/(1.0+exp((-60.0+V*1000.0)/20.0));
tau_Xs      =@(V) 1.0*alpha_Xs(V)*beta_Xs(V)/1000.0;

i_Ks        =@(V,Cai,Xs,Nai) g_Ks*(V-E_Ks(Nai))*Xs^2.0*(1.0+0.6/(1.0+(3.8*0.00001/Cai)^1.4));

Xr1_inf      =@(V) 1.0/(1.0+exp((V_half-V*1000.0)/4.9));
alpha_Xr1    =@(V) 450.0/(1.0+exp((-45.0-V*1000.0)/10.0));
beta_Xr1     =@(V) 6.0/(1.0+exp((30.0+V*1000.0)/11.5));
tau_Xr1      =@(V) 1.0*alpha_Xr1(V)*beta_Xr1(V)/1000.0;

Xr2_inf      =@(V) 1.0/(1.0+exp((V*1000.0+88.0)/50.0));
alpha_Xr2    =@(V) 3.0/(1.0+exp((-60.0-V*1000.0)/20.0));
beta_Xr2     =@(V) 1.12/(1.0+exp((-60.0+V*1000.0)/20.0));
tau_Xr2      =@(V) 1.0*alpha_Xr2(V)*beta_Xr2(V)/1000.0;

i_Kr         =@(V,Xr1,Xr2) g_Kr*(V-E_K)*Xr1*Xr2*sqrt(Ko/5.4);

% alpha_K1    =@(V) 3.91/(1.0+exp(0.5942*(V*1000.0-E_K*1000.0-200.0)));
% beta_K1     =@(V) (-1.509*exp(0.0002*(V*1000.0-E_K*1000.0+100.0))+exp(0.5886*(V*1000.0-E_K*1000.0-10.0)))/(1.0+exp(0.4547*(V*1000.0-E_K*1000.0)));
% XK1_inf     =@(V) alpha_K1(V)/(alpha_K1(V)+beta_K1(V));

% i_K1        =@(V) g_K1*XK1_inf(V)*(V-E_K)*sqrt(Ko/5.4);
i_K1        =@(V) g_K1*Ko^(0.4457)*1000.0*(V-E_K)/(1.0+exp(1.5*(V*1000.0-E_K*1000.0+3.6)*F/(R*1000.0*T)));

i_NaCa      =@(V,Nai,Cai) kNaCa*(exp(gamma*V*F/(R*T))*Nai^3.0*Cao-exp((gamma-1.0)*V*F/(R*T))*Nao^3.0*Cai*alpha)/((KmNai^3.0+Nao^3.0)*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*V*F/(R*T))));

i_NaK       =@(V,Nai) 1.3*(PNaK*Ko/(Ko+Km_K)*Nai/(Nai+Km_Na)/(1.0+0.1245*exp(-0.1*V*F/(R*T))+0.0353*exp(-V*F/(R*T))));

i_PCa       =@(Cai) g_PCa*Cai/(Cai+KPCa);

i_b_Na      =@(V,Nai) 0.1*(g_b_Na*(V-E_Na(Nai)));

i_b_Ca      =@(V,Cai) g_b_Ca*(V-E_Ca(Cai));

i_up        =@(Cai) VmaxUp/(1.0+Kup^2.0/Cai^2.0);

i_leak      =@(CaSR,Cai) (CaSR-Cai)*V_leak;

RyRainfss   =@(Cai) RyRa1-RyRa2/(1 + exp((1000.0*Cai-(RyRahalf))/0.0082));
RyRoinfss   =@(RyRa,Cai) (1 - 1/(1 +  exp((1000.0*Cai-(RyRa+ RyRohalf))/0.003)));
RyRtauact   =@(RyRo,Cai,RyRa) (18.75e-3)*(RyRoinfss(RyRa,Cai)>= RyRo) + (0.1*18.75e-3)*(RyRoinfss(RyRa,Cai)< RyRo);

RyRcinfss   =@(RyRa,Cai) (1/(1 + exp((1000.0*Cai-(RyRa+RyRchalf))/0.001)));
RyRtauinact =@(RyRc,Cai,RyRa) (2*87.5e-3)*(RyRcinfss(RyRa,Cai)>= RyRc) + (87.5e-3)*(RyRcinfss(RyRa,Cai)< RyRc);

RyRSRCass   =@(CaSR) (1 - 1/(1 +  exp((CaSR-0.3)/0.1)));
i_rel       =@(CaSR,Cai,RyRo,RyRc) g_irel_max*RyRSRCass(CaSR)*RyRo*RyRc*(CaSR-Cai);

Cai_bufc    =@(Cai) 1.0/(1.0+Buf_C*Kbuf_C/(Cai+Kbuf_C)^2.0);
Ca_SR_bufSR =@(CaSR) 1.0/(1.0+Buf_SR*Kbuf_SR/(CaSR+Kbuf_SR)^2.0);

g_Kur       =@(V) coeff_Kur * (0.005 + 0.05/(1+exp((V*1000.0-15.0)/-13.0))); % nS/pF
ua_inf     =@(V) 1.0 / (1.0+exp(-(1000.0*V+30.3)/9.6));
ui_inf     =@(V) 1.0 / (1.0+exp((1000.0*V-99.45)/27.48));
alpha_ua    =@(V) 0.65 / (exp(-(1000.0*V+10.0)/8.5)+exp(-(1000.0*V-30.0)/59.0));
beta_ua     =@(V) 0.65 / (exp(-(1000.0*V+82.0)/17.0)+2.5);
alpha_ui    =@(V) 1.0 / (21.0+exp((1000.0*V-185.0)/28.0));
beta_ui     =@(V) exp((V*1000.0-158.0)/16.0);
tau_ua      =@(V) KQ10/(alpha_ua(V)+beta_ua(V));
tau_ui      =@(V) KQ10/(alpha_ui(V)+beta_ui(V));

i_Kur       =@(V,ua,ui)  g_Kur(V) * ua^3.0 * ui * (1000.0*V-1000.0*E_K);

dO          =@(O,Cai) (1.0-O) * KCa_on * Cai^2.0 - O * KCa_off;
i_KCa       =@(V,O) g_KCa * O *(1.0/(1.0 + exp((V*1000.0 - E_K*1000.0 + 120.0)/45.0))) * (V*1000.0 - E_K*1000.0);


Tf=2.5; %500;
nt=5000;
dt=Tf/nt;
T=linspace(0,Tf,nt+1);

V=zeros(nt+1,1);
CaSR=zeros(nt+1,1);
Cai=zeros(nt+1,1);
g=zeros(nt+1,1);
d=zeros(nt+1,1);
f1=zeros(nt+1,1);
f2=zeros(nt+1,1);
fCa=zeros(nt+1,1);
Xr1=zeros(nt+1,1);
Xr2=zeros(nt+1,1);
Xs=zeros(nt+1,1);
h=zeros(nt+1,1);
j=zeros(nt+1,1);
m=zeros(nt+1,1);
Xf=zeros(nt+1,1);
q=zeros(nt+1,1);
r=zeros(nt+1,1);
Nai=zeros(nt+1,1);
mL=zeros(nt+1,1);
hL=zeros(nt+1,1);
RyRa=zeros(nt+1,1);
RyRo=zeros(nt+1,1);
RyRc=zeros(nt+1,1);
ua=zeros(nt+1,1);
ui=zeros(nt+1,1);
O=zeros(nt+1,1);


V(1)=U_rest;
CaSR(1)=CaSR0;
Cai(1)=Cai0;
g(1)=g0;
d(1)=d0;
f1(1)=f10;
f2(1)=f20;
fCa(1)=fCa0;
Xr1(1)=Xr10;
Xr2(1)=Xr20;
Xs(1)=Xs0;
h(1)=h0;
j(1)=j0;
m(1)=m0;
Xf(1)=Xf0;
q(1)=q0;
r(1)=r0;
Nai(1)=Nai0;
mL(1)=mL0;
hL(1)=hL0;
RyRa(1)=RyRa0;
RyRo(1)=RyRo0;
RyRc(1)=RyRc0;
ua(1)=ua0;
ui(1)=ui0;
O(1)=O0;

for i=2:nt+1
    t=T(i);

    % if mod(i,1000)
    %     te=num2str(t);
    %     Vma=num2str(V(i-1));
    %     dma=num2str(d(i-1));
    %     f2ma=num2str(f2(i-1));
    %     fCassma=num2str(fCass(i-1));
    %     pfma=num2str(pf(i-1));
    %     Casrma=num2str(Ca_sr(i-1));
    %     Caima=num2str(Ca_i(i-1));
    %     Cassma=num2str(Ca_ss(i-1));
    %     Relma=num2str(Rel(i-1));
    %     hma=num2str(h(i-1));
    %     jma=num2str(j(i-1));
    %     mma=num2str(m(i-1));
    %     Kima=num2str(K_i(i-1));
    %     xr1ma=num2str(xr1(i-1));
    %     xr2ma=num2str(xr2(i-1));
    %     xsma=num2str(xs(i-1));
    %     Naima=num2str(Na_i(i-1));
    %     prma=num2str(pr(i-1));
    %     sma=num2str(s(i-1));
    % 
    % 
    %     %disp(['t=',te,' Vmax=',Vma,' dma=',dma,' f2ma=',f2ma,' fCassma=',fCassma,' pfma=',pfma,' Casrma=',Casrma,' Caima=',Caima,' Cassma=',Cassma,' Relma=',Relma,' hma=',hma,' jma=',jma,' mma=',mma,' Kima=',Kima,' xr1ma=',xr1ma,' xr2ma=',xr2ma,' xsma=',xsma,' Naima=',Naima,' prma=',prma,' sma=',sma])
    % end

    Vold=V(i-1)
    CaSRold=CaSR(i-1);
    Caiold=Cai(i-1);
    gold=g(i-1);
    dold=d(i-1);
    f1old=f1(i-1);
    f2old=f2(i-1);
    fCaold=fCa(i-1);
    Xr1old=Xr1(i-1);
    Xr2old=Xr2(i-1);
    Xsold=Xs(i-1);
    hold=h(i-1);
    jold=j(i-1);
    mold=m(i-1);
    Xfold=Xf(i-1);
    qold=q(i-1);
    rold=r(i-1);
    Naiold=Nai(i-1);
    mLold=mL(i-1);
    hLold=hL(i-1);
    RyRaold=RyRa(i-1);
    RyRoold=RyRo(i-1);
    RyRcold=RyRc(i-1);
    uaold=ua(i-1);
    uiold=ui(i-1);
    Oold=O(i-1);
    
    m_inf_eval = m_inf(Vold);
    tau_m_eval = tau_m(Vold);
    m(i)=m_inf_eval+(mold-m_inf_eval)*exp(-dt/tau_m_eval);

    h_inf_eval = h_inf(Vold);
    tau_h_eval = tau_h(Vold);
    h(i)=h_inf_eval+(hold-h_inf_eval)*exp(-dt/tau_h_eval);

    j_inf_eval = j_inf(Vold);
    tau_j_eval = tau_j(Vold);
    j(i)=j_inf_eval+(jold-j_inf_eval)*exp(-dt/tau_j_eval);

    mL_inf_eval = mL_inf(Vold);
    tau_mL_eval = tau_mL(Vold);
    mL(i)=mL_inf_eval+(mLold-mL_inf_eval)*exp(-dt/tau_mL_eval);

    hL_inf_eval = hL_inf(Vold);
    tau_hL_eval = tau_hL;
    hL(i)=hL_inf_eval+(hLold-hL_inf_eval)*exp(-dt/tau_hL_eval);

    Xf_inf_eval = Xf_inf(Vold);
    tau_Xf_eval = tau_Xf(Vold);
    Xf(i)=Xf_inf_eval+(Xfold-Xf_inf_eval)*exp(-dt/tau_Xf_eval);
   
    d_inf_eval = d_inf(Vold);
    tau_d_eval = tau_d(Vold);
    d(i)=d_inf_eval+(dold-d_inf_eval)*exp(-dt/tau_d_eval);

    f1_inf_eval = f1_inf(Vold);
    tau_f1_eval = tau_f1(Vold,Caiold,f1old);
    f1(i)=f1_inf_eval+(f1old-f1_inf_eval)*exp(-dt/tau_f1_eval);

    f2_inf_eval = f2_inf(Vold);
    tau_f2_eval = tau_f2(Vold);
    f2(i)=f2_inf_eval+(f2old-f2_inf_eval)*exp(-dt/tau_f2_eval);

    fCa_inf_eval = fCa_inf(Caiold);
    tau_fCa_eval = tau_fCa(Vold,Caiold,fCaold);
    fCa(i)=fCa_inf_eval+(fCaold-fCa_inf_eval)*exp(-dt/tau_fCa_eval);
   
    q_inf_eval = q_inf(Vold);
    tau_q_eval = tau_q(Vold);
    q(i)=q_inf_eval+(qold-q_inf_eval)*exp(-dt/tau_q_eval);

    r_inf_eval = r_inf(Vold);
    tau_r_eval = tau_r(Vold);
    r(i)=r_inf_eval+(rold-r_inf_eval)*exp(-dt/tau_r_eval);

    Xs_inf_eval = Xs_inf(Vold);
    tau_Xs_eval = tau_Xs(Vold);
    Xs(i)=Xs_inf_eval+(Xsold-Xs_inf_eval)*exp(-dt/tau_Xs_eval);

    Xr1_inf_eval = Xr1_inf(Vold);
    tau_Xr1_eval = tau_Xr1(Vold);
    Xr1(i)=Xr1_inf_eval+(Xr1old-Xr1_inf_eval)*exp(-dt/tau_Xr1_eval);

    Xr2_inf_eval = Xr2_inf(Vold);
    tau_Xr2_eval = tau_Xr2(Vold);
    Xr2(i)=Xr2_inf_eval+(Xr2old-Xr2_inf_eval)*exp(-dt/tau_Xr2_eval);

    RyRainfss_eval=RyRainfss(Caiold);
    dRyRa_eval=(RyRainfss_eval- RyRaold)/RyRtauadapt;
    RyRa(i)= RyRaold+dt*dRyRa_eval;

    RyRoinfss_eval=RyRoinfss(RyRaold,Caiold);
    RyRtauact_eval=RyRtauact(RyRoold,Caiold,RyRaold);
    dRyRo_eval=(RyRoinfss_eval- RyRoold)/(RyRtauact_eval);
    RyRo(i)= RyRoold+dt*dRyRo_eval;

    RyRcinfss_eval   =RyRcinfss(RyRaold,Caiold);
    RyRtauinact_eval =RyRtauinact(RyRcold,Caiold,RyRaold);
    dRyRc_eval = (RyRcinfss_eval- RyRcold)/(RyRtauinact_eval);
    RyRc(i)= RyRcold+dt*dRyRc_eval;

    ua_inf_eval = ua_inf(Vold);
    tau_ua_eval = tau_ua(Vold);
    ua(i)=ua_inf_eval+(uaold-ua_inf_eval)*exp(-dt/tau_ua_eval);
    
    ui_inf_eval = ui_inf(Vold);
    tau_ui_eval = tau_ui(Vold);
    ui(i)=ui_inf_eval+(uiold-ui_inf_eval)*exp(-dt/tau_ui_eval);

    dO_eval=dO(Oold,Caiold);
    O(i)= Oold+dt*dO_eval;
    
    
    i_Na_eval=i_Na(Vold,m(i),h(i),j(i),Naiold);
    i_NaL_eval=i_NaL(Vold,mL(i),hL(i),Naiold);
    i_f_eval=i_f(Vold,Xf(i),Naiold);
    i_CaL_eval=i_CaL(Vold,Caiold,d(i),f1(i),f2(i),fCa(i));
    i_to_eval=i_to(Vold,q(i),r(i));
    i_Ks_eval=i_Ks(Vold,Caiold,Xs(i),Naiold);
    i_Kr_eval=i_Kr(Vold,Xr1(i),Xr2(i));
    i_K1_eval=i_K1(Vold);
    i_NaCa_eval=i_NaCa(Vold,Naiold,Caiold);
    i_NaK_eval=i_NaK(Vold,Naiold);
    i_PCa_eval=i_PCa(Caiold);
    i_b_Na_eval=i_b_Na(Vold,Naiold);
    i_b_Ca_eval=i_b_Ca(Vold,Caiold);
    i_up_eval=i_up(Caiold);
    i_leak_eval=i_leak(CaSRold,Caiold);
    i_fNa_eval=i_fNa(Vold,Xf(i),Naiold);
    i_rel_eval=i_rel(CaSRold,Caiold,RyRo(i),RyRc(i));
    i_Kur_eval=i_Kur(Vold,ua(i),ui(i));
    i_KCa_eval=i_KCa(Vold,O(i));

    dNai_eval   = -Cm*(i_Na_eval+i_NaL_eval+i_b_Na_eval+3.0*i_NaK_eval+3.0*i_NaCa_eval+i_fNa_eval)/(F*Vc*1.0e-18);
    Nai(i)      = Naiold+dt*dNai_eval;

    dCai_eval   =@(Cai) Cai_bufc(Cai)*(i_leak_eval-i_up_eval+i_rel_eval-(i_CaL_eval+i_b_Ca_eval+i_PCa_eval-2.0*i_NaCa_eval)*Cm/(2.0*Vc*F*1.0e-18));
    Cai(i)      = Caiold+dt*dCai_eval(Caiold);

    dCaSR_eval   =@(CaSR) Ca_SR_bufSR(CaSR)*Vc/V_SR*(i_up_eval-(i_rel_eval+i_leak_eval));
    CaSR(i)      = CaSRold+dt*dCaSR_eval(CaSRold);
    
    
    if (0<t && t<0.003)
        Iapp=14.1e-10/Cm; %7.5e-10;
    else
        Iapp=0;
    end

    Itot=Iapp-(i_KCa_eval+i_Kur_eval+i_K1_eval+i_to_eval+i_Kr_eval+i_Ks_eval+i_CaL_eval+i_NaK_eval+i_Na_eval+i_NaL_eval+i_NaCa_eval+i_PCa_eval+i_f_eval+i_b_Na_eval+i_b_Ca_eval);

    
    V(i)=Vold+dt*Itot;

end

plot(T,V)
