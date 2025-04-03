clear all
close all

% parametri
g_CaL = 0.0000398d0;
g_bca = 0.000592d0;
Buf_c = 0.2d0;
Buf_sr = 10.0d0;
Buf_ss = 0.4d0;
Ca_o = 2.0d0;
EC = 1.5d0;
K_buf_c = 0.001d0;
K_buf_sr = 0.3d0;
K_buf_ss = 0.00025d0;
K_up = 0.00025d0;
V_leak = 0.00036d0;
V_rel = 0.102d0;
V_sr = 0.001094d0;
V_ss = 0.00005468d0;
V_xfer = 0.0038d0;
Vmax_up = 0.006375d0;
k1_prime = 0.15d0;
k2_prime = 0.045d0;
k3 = 0.06d0;
k4 = 0.005d0;
max_sr = 2.5d0;
min_sr = 1.0d0;
K_pCa = 0.0005d0;
g_pCa = 0.1238d0;
g_K1 = 5.405d0;
Cm = 0.185d0;
F = 96485.3415d0;
R = 8314.472d0;
T = 310.0d0;
V_c = 0.016404d0;

K_o=5.4;


g_pK = 0.0146d0;
g_Kr = 0.153d0;
P_kna = 0.03d0;
g_Ks = 0.392d0;
g_bna = 0.00029d0;
K_NaCa = 1000.0d0;
K_sat = 0.1d0;
Km_Ca = 1.38d0;
Km_Nai = 87.5d0;
alpha = 2.5d0;
gamma = 0.35d0;
Na_o = 140.0d0;
K_mNa = 40.0d0;
K_mk = 1.0d0;
P_NaK = 2.724d0;
g_to = 0.294d0;
g_Na = 14.838d0;



% condizioni iniziali
U_rest=-8.523000000000000e+01;
d0=3.373000000000000e-05;
f20=9.755000000000000e-01;
fCass0=9.953000000000000e-01;
pf0=7.887999999999999e-01;
Ca_sr0=3.640000000000000e+00;
Ca_i0=1.260000000000000e-04;
Ca_ss0=3.600000000000000e-04;
Rel0=9.073000000000000e-01;
h0=7.444000000000000e-01;
j0=7.045000000000000e-01;
m0=1.720000000000000e-03;
K_i0=1.368900000000000e+02;
xr10=6.210000000000000e-03;
xr20=4.712000000000000e-01;
xs0=9.500000000000000e-03;
Na_i0=8.603999999999999e+00;
pr0=2.420000000000000e-08;
s0=9.999980000000001e-01;

E_Ca =@(Ca_i) 0.5d0*R*T/F*log(Ca_o/Ca_i);
E_K =@(K_i) R*T/F*log(K_o/K_i);

i_CaL =@(V,d,pf,f2,fCass,Ca_ss) g_CaL*d*pf*f2*fCass*4.d0*(V-15.d0)*F^2/...
          (R*T)*(0.25d0*Ca_ss*exp(2.d0*(V-15.d0)*F/(R*T))-Ca_o)/...
          (exp(2.d0*(V-15.d0)*F/(R*T))-1.d0);
i_b_Ca =@(V,Ca_i) g_bca*(V-E_Ca(Ca_i));

d_inf =@(V) 1.d0/(1.d0+exp((-8.d0-V)/7.5d0));
alpha_d =@(V) 1.4d0/(1.d0+exp((-35.d0-V)/13.d0))+0.25d0;
beta_d =@(V) 1.4d0/(1.d0+exp((V+5.d0)/5.d0));
gamma_d =@(V) 1.d0/(1.d0+exp((50.d0-V)/20.d0));
tau_d =@(V) 1.d0*alpha_d(V)*beta_d(V)+gamma_d(V);

f2_inf =@(V) 0.67d0/(1.d0+exp((V+35.d0)/7.d0))+0.33d0;
tau_f2 =@(V) 562.d0*exp(-(V+27.d0)^2/240.d0)+31.d0/...
       (1.d0+exp((25.d0-V)/10.d0))+80.d0/...
       (1.d0+exp((V+30.d0)/10.d0)); 

fCass_inf =@(V,Ca_ss) 0.6d0/(1.d0+(Ca_ss/0.05d0)^2)+0.4d0;
tau_fCass =@(V,Ca_ss) 80.d0/(1.d0+(Ca_ss/0.05d0)^2)+2.d0;

f_inf =@(V) 1.d0/(1.d0+exp((V+20.d0)/7.d0));
tau_f =@(V) 1102.5d0*exp(-(V+27.d0)^2/225.d0)+200.d0/...
      (1.d0+exp((13.d0-V)/10.d0))+180.d0/...
      (1.d0+exp((V+30.d0)/10.d0))+20.d0;
   
kcasr =@(Ca_sr) max_sr-(max_sr-min_sr)/(1.d0+(EC/Ca_sr)^2);
k1 =@(Ca_sr) k1_prime/kcasr(Ca_sr);
O =@(Rel, Ca_sr, Ca_ss) k1(Ca_sr)*Ca_ss^2*Rel/(k3+k1(Ca_sr)*Ca_ss^2);
i_rel =@(Rel, Ca_sr, Ca_ss) V_rel*O(Rel, Ca_sr, Ca_ss)*(Ca_sr-Ca_ss);
i_up =@(Ca_i) Vmax_up/(1.d0+(K_up/Ca_i)^2);
i_leak =@(Ca_sr, Ca_i) V_leak*(Ca_sr-Ca_i);
i_xfer =@(Ca_ss, Ca_i) V_xfer*(Ca_ss-Ca_i);

k2 =@(Ca_sr) k2_prime*kcasr(Ca_sr);
dRel =@(Ca_ss, Rel, Ca_sr) (-k2(Ca_sr)*Ca_ss*Rel+k4*(-Rel+1.0));

Ca_i_bufc =@(Ca_i) 1.d0/(1.d0+Buf_c*K_buf_c/(Ca_i+K_buf_c)^2);
Ca_sr_bufsr =@(Ca_sr) 1.d0/(1.d0+Buf_sr*K_buf_sr/(Ca_sr+K_buf_sr)^2);
Ca_ss_bufss =@(Ca_ss) 1.d0/(1.d0+Buf_ss*K_buf_ss/(Ca_ss+K_buf_ss)^2);
i_p_Ca =@(Ca_i) g_pCa*Ca_i/(Ca_i+K_pCa);

i_NaCa =@(V, Na_i, Ca_i) K_NaCa*(exp(gamma*V*F/(R*T))*Na_i^3*Ca_o-...
        exp((gamma-1.d0)*V*F/(R*T))*Na_o^3*Ca_i*alpha)/...
        ((Km_Nai^3+Na_o^3)*(Km_Ca+Ca_o)*...
        (1.d0+K_sat*exp((gamma-1.d0)*V*F/(R*T))));

dCa_i =@(Na_i, Ca_ss, Ca_sr, Ca_i, V) Ca_i_bufc(Ca_i)*((i_leak(Ca_sr, Ca_i)-i_up(Ca_i))*V_sr/V_c+i_xfer(Ca_ss, Ca_i)-...
       1.d0*(i_b_Ca(V,Ca_i)+i_p_Ca(Ca_i)-2.d0*i_NaCa(V, Na_i, Ca_i))*Cm/(2.d0*1.d0*V_c*F));

dCa_sr =@(Rel, Ca_sr, Ca_i, Ca_ss) Ca_sr_bufsr(Ca_sr)*(i_up(Ca_i)-(i_rel(Rel, Ca_sr, Ca_ss)+i_leak(Ca_sr, Ca_i)));

dCa_ss =@(Rel,V,d,pf,f2,fCass,Ca_ss,Ca_sr,Ca_i) Ca_ss_bufss(Ca_ss)*(-1.d0*i_CaL(V,d,pf,f2,fCass,Ca_ss)*Cm/(2.d0*1.d0*V_ss*F)+...
           i_rel(Rel, Ca_sr, Ca_ss)*V_sr/V_ss-i_xfer(Ca_ss, Ca_i)*V_c/V_ss);

E_Na =@(Na_i) R*T/F*log(Na_o/Na_i);
i_Na =@(V, m, h, j, Na_i) g_Na*m^3*h*j*(V-E_Na(Na_i));

h_inf =@(V) 1.d0/(1.d0+exp((V+71.55d0)/7.43d0))^2;
alpha_h_l=@(V) 0.057d0*exp(-(V+80.d0)/6.8d0);
alpha_h =@(V) (alpha_h_l(V)*(V<-40.0)+0.d0*(V>-40.0));
beta_h_l=@(V) (2.7d0*exp(0.079d0*V)+3.1d+5*exp(0.3485d0*V) );
beta_h_g=@(V) (0.77d0/(0.13d0*(1.d0+exp((V+10.66d0)/(-11.1d0)))) );
beta_h =@(V) (beta_h_l(V)*(V<-40.0)+beta_h_g(V)*(V>-40.0));
tau_h =@(V) 1.d0/(alpha_h(V) + beta_h(V));

j_inf =@(V) 1.d0/(1.d0+exp((V+71.55d0)/7.43d0))^2;
alpha_j_l=@(V) (-25428.d0*exp(0.2444d0*V)-6.948d-6*...
               exp(-0.04391d0*V))*(V+37.78d0)/...
               1.d0/(1.d0+exp(0.311d0*(V+79.23d0)));
alpha_j_g=@(V) (0.d0);
alpha_j =@(V) (alpha_j_l(V)*(V<-40.0)+alpha_j_g(V)*(V>-40.0));
beta_j_l=@(V) (0.02424d0*exp(-0.01052d0*V)/...
              (1.d0+exp(-0.1378d0*(V+40.14d0))) );
beta_j_g=@(V) (0.6d0*exp(0.057d0*V)/...
             (1.d0+exp(-0.1d0*(V+32.d0))) );
beta_j =@(V) (beta_j_l(V)*(V<-40.0)+beta_j_g(V)*(V>-40.0));
tau_j =@(V) 1.d0/(alpha_j(V)+beta_j(V));

m_inf =@(V) 1.d0/(1.d0+exp((-56.86d0-V)/9.03d0))^2;
alpha_m =@(V) 1.d0/(1.d0+exp((-60.d0-V)/5.d0));
beta_m =@(V) 0.1d0/(1.d0+exp((V+35.d0)/5.d0))+...
        0.1d0/(1.d0+exp((V-50.d0)/200.d0));
tau_m =@(V) 1.d0*alpha_m(V)*beta_m(V);

xr1_inf =@(V) 1.d0/(1.d0+exp((-26.d0-V)/7.d0));
alpha_xr1 =@(V) 450.d0/(1.d0+exp((-45.d0-V)/10.d0));
beta_xr1 =@(V) 6.d0/(1.d0+exp((V+30.d0)/11.5d0));
tau_xr1 =@(V) 1.d0*alpha_xr1(V)*beta_xr1(V);

xr2_inf =@(V) 1.d0/(1.d0+exp((V+88.d0)/24.d0));
alpha_xr2 =@(V) 3.d0/(1.d0+exp((-60.d0-V)/20.d0));
beta_xr2 =@(V) 1.12d0/(1.d0+exp((V-60.d0)/20.d0));
tau_xr2 =@(V) 1.d0*alpha_xr2(V)*beta_xr2(V);

xs_inf =@(V) 1.d0/(1.d0+exp((-5.d0-V)/14.d0));
alpha_xs =@(V) 1400.d0/sqrt(1.d0+exp((5.d0-V)/6.d0));
beta_xs =@(V) 1.d0/(1.d0+exp((V-35.d0)/15.d0));
tau_xs =@(V) 1.d0*alpha_xs(V)*beta_xs(V)+80.d0;

r_inf =@(V) 1.d0/(1.d0+exp((20.d0-V)/6.d0));
tau_r =@(V) 9.5d0*exp(-(V+40.d0)^2/1800.d0)+0.8d0;

s_inf =@(V) 1.d0/(1.d0+exp((V+20.d0)/5.d0));
tau_s =@(V) 85.d0*exp(-(V+45.d0)^2/320.d0)+...
       5.d0/(1.d0+exp((V-20.d0)/5.d0))+3.d0;

alpha_K1 =@(V,K_i) 0.1d0/(1.d0+exp(0.06d0*(V-E_K(K_i)-200.d0)));
beta_K1 =@(V,K_i) (3.d0*exp(0.0002d0*(V-E_K(K_i)+100.d0))+...
         exp(0.1d0*(V-E_K(K_i)-10.d0)))/...
         (1.d0+exp(-0.5d0*(V-E_K(K_i))));
xK1_inf =@(V,K_i) alpha_K1(V,K_i)/(alpha_K1(V,K_i)+beta_K1(V,K_i));
i_K1 =@(V,K_i) g_K1*xK1_inf(V,K_i)*sqrt(K_o/5.4d0)*(V-E_K(K_i));

i_to =@(V, pr, s, K_i) g_to*pr*s*(V-E_K(K_i));
i_Kr =@(V, xr1, xr2, K_i) g_Kr*sqrt(K_o/5.4d0)*xr1*xr2*(V-E_K(K_i));

E_Ks =@(Na_i,K_i) R*T/F*log((K_o+P_kna*Na_o)/(K_i+P_kna*Na_i));
i_Ks =@(V, Na_i, xs,K_i) g_Ks*xs^2*(V-E_Ks(Na_i,K_i));

i_NaK =@(V,Na_i) P_NaK*K_o/(K_o+K_mk)*Na_i/(Na_i+K_mNa)/...
       (1.d0+0.1245d0*exp(-0.1d0*V*F/(R*T))+...
       0.0353d0*exp(-V*F/(R*T)));
i_b_Na =@(V,Na_i) g_bna*(V-E_Na(Na_i));
i_p_K =@(V,K_i) g_pK*(V-E_K(K_i))/(1.d0+exp((25.d0-V)/5.98d0));





Tf=500;
nt=5000;
dt=Tf/nt;
T=linspace(0,Tf,nt+1);

V=zeros(nt+1,1);
d=zeros(nt+1,1);
f2=zeros(nt+1,1);
fCass=zeros(nt+1,1);
pf=zeros(nt+1,1);
Ca_sr=zeros(nt+1,1);
Ca_i=zeros(nt+1,1);
Ca_ss=zeros(nt+1,1);
Rel=zeros(nt+1,1);
h=zeros(nt+1,1);
j=zeros(nt+1,1);
m=zeros(nt+1,1);
K_i=zeros(nt+1,1);
xr1=zeros(nt+1,1);
xr2=zeros(nt+1,1);
xs=zeros(nt+1,1);
Na_i=zeros(nt+1,1);
pr=zeros(nt+1,1);
s=zeros(nt+1,1);


V(1)=U_rest;
d(1)=d0;
f2(1)=f20;
fCass(1)=fCass0;
pf(1)=pf0;
Ca_sr(1)=Ca_sr0;
Ca_i(1)=Ca_i0;
Ca_ss(1)=Ca_ss0;
Rel(1)=Rel0;
h(1)=h0;
j(1)=j0;
m(1)=m0;
K_i(1)=K_i0;
xr1(1)=xr10;
xr2(1)=xr20;
xs(1)=xs0;
Na_i(1)=Na_i0;
pr(1)=pr0;
s(1)=s0;

for i=2:nt+1
    t=T(i);

    if mod(i,1000)
        te=num2str(t);
        Vma=num2str(V(i-1));
        dma=num2str(d(i-1));
        f2ma=num2str(f2(i-1));
        fCassma=num2str(fCass(i-1));
        pfma=num2str(pf(i-1));
        Casrma=num2str(Ca_sr(i-1));
        Caima=num2str(Ca_i(i-1));
        Cassma=num2str(Ca_ss(i-1));
        Relma=num2str(Rel(i-1));
        hma=num2str(h(i-1));
        jma=num2str(j(i-1));
        mma=num2str(m(i-1));
        Kima=num2str(K_i(i-1));
        xr1ma=num2str(xr1(i-1));
        xr2ma=num2str(xr2(i-1));
        xsma=num2str(xs(i-1));
        Naima=num2str(Na_i(i-1));
        prma=num2str(pr(i-1));
        sma=num2str(s(i-1));
        

        %disp(['t=',te,' Vmax=',Vma,' dma=',dma,' f2ma=',f2ma,' fCassma=',fCassma,' pfma=',pfma,' Casrma=',Casrma,' Caima=',Caima,' Cassma=',Cassma,' Relma=',Relma,' hma=',hma,' jma=',jma,' mma=',mma,' Kima=',Kima,' xr1ma=',xr1ma,' xr2ma=',xr2ma,' xsma=',xsma,' Naima=',Naima,' prma=',prma,' sma=',sma])
    end

    Vold=V(i-1)
    dold=d(i-1);
    f2old=f2(i-1);
    fCassold=fCass(i-1);
    pfold=pf(i-1);
    Ca_srold=Ca_sr(i-1);
    Ca_iold=Ca_i(i-1);
    Ca_ssold=Ca_ss(i-1);
    Relold=Rel(i-1);
    hold=h(i-1);
    jold=j(i-1);
    mold=m(i-1);
    K_iold=K_i(i-1);
    xr1old=xr1(i-1);
    xr2old=xr2(i-1);
    xsold=xs(i-1);
    Na_iold=Na_i(i-1);
    prold=pr(i-1);
    sold=s(i-1);
    
    d_inf_eval = d_inf(Vold);
    tau_d_eval = tau_d(Vold);
    d(i)=d_inf_eval+(dold-d_inf_eval)*exp(-dt/tau_d_eval);

    f2_inf_eval = f2_inf(Vold);
    tau_f2_eval = tau_f2(Vold);
    f2(i)=f2_inf_eval+(f2old-f2_inf_eval)*exp(-dt/tau_f2_eval);

    fCass_inf_eval = fCass_inf(Vold,Ca_ssold);
    tau_fCass_eval = tau_fCass(Vold,Ca_ssold);
    fCass(i)=fCass_inf_eval+(fCassold-fCass_inf_eval)*exp(-dt/tau_fCass_eval);
    
    f_inf_eval = f_inf(Vold);
    tau_f_eval = tau_f(Vold);
    pf(i)=f_inf_eval+(pfold-f_inf_eval)*exp(-dt/tau_f_eval);

    dRel_eval = dRel(Ca_ssold, Relold, Ca_srold);
    Rel(i)=Relold+dt*dRel_eval;

    dCa_i_eval = dCa_i(Na_iold, Ca_ssold, Ca_srold, Ca_iold, Vold);
    Ca_i(i)=Ca_iold+dt*dCa_i_eval;

    dCa_sr_eval=dCa_sr(Relold, Ca_srold, Ca_iold, Ca_ssold);
    Ca_sr(i)=Ca_srold+dt*dCa_sr_eval;

    dCa_ss_eval=dCa_ss(Relold,Vold,dold,pfold,f2old,fCassold,Ca_ssold,Ca_srold,Ca_iold);
    Ca_ss(i)=Ca_ssold+dt*dCa_ss_eval;

    h_inf_eval=h_inf(Vold);
    tau_h_eval=tau_h(Vold);
    h(i)=h_inf_eval+(hold-h_inf_eval)*exp(-dt/tau_h_eval);

    j_inf_eval=j_inf(Vold);
    tau_j_eval=tau_j(Vold);
    j(i)=j_inf_eval+(hold-j_inf_eval)*exp(-dt/tau_j_eval);

    m_inf_eval=m_inf(Vold);
    tau_m_eval=tau_m(Vold);
    m(i)=m_inf_eval+(mold-m_inf_eval)*exp(-dt/tau_m_eval);

    xr1_inf_eval=xr1_inf(Vold);
    tau_xr1_eval=tau_xr1(Vold);
    xr1(i)=xr1_inf_eval+(xr1old-xr1_inf_eval)*exp(-dt/tau_xr1_eval);

    xr2_inf_eval=xr2_inf(Vold);
    tau_xr2_eval=tau_xr2(Vold);
    xr2(i)=xr2_inf_eval+(xr2old-xr2_inf_eval)*exp(-dt/tau_xr2_eval);

    xs_inf_eval=xs_inf(Vold);
    tau_xs_eval=tau_xs(Vold);
    xs(i)=xs_inf_eval+(xsold-xs_inf_eval)*exp(-dt/tau_xs_eval);

    r_inf_eval=r_inf(Vold);
    tau_r_eval=tau_r(Vold);
    pr(i)=r_inf_eval+(prold-r_inf_eval)*exp(-dt/tau_r_eval);

    s_inf_eval=s_inf(Vold);
    tau_s_eval=tau_s(Vold);
    s(i)=s_inf_eval+(sold-s_inf_eval)*exp(-dt/tau_s_eval);

    i_K1_eval=i_K1(Vold,K_iold);
    i_to_eval=i_to(Vold, pr(i), s(i),K_iold);
    i_Kr_eval=i_Kr(Vold, xr1(i), xr2(i),K_iold);
    i_Ks_eval=i_Ks(Vold, Na_iold, xs(i),K_iold);
    i_CaL_eval=i_CaL(Vold,d(i),pf(i),f2(i),fCass(i),Ca_ss(i));
    i_NaK_eval=i_NaK(Vold,Na_iold);
    i_Na_eval=i_Na(Vold,m(i),h(i),j(i),Na_iold);
    i_b_Na_eval=i_b_Na(Vold,Na_iold);
    i_NaCa_eval=i_NaCa(Vold,Na_iold,Ca_i(i));
    i_b_Ca_eval=i_b_Ca(Vold,Ca_i(i));
    i_p_Ca_eval=i_p_Ca(Ca_i(i));
    i_p_K_eval=i_p_K(Vold,K_iold);  
    
    dK_i_eval = -1.d0*(i_K1_eval+i_to_eval+i_Kr_eval+i_Ks_eval+i_p_K_eval-2.d0*i_NaK_eval)/...
               (1.d0*V_c*F)*Cm;
    K_i(i) = K_iold+dt*dK_i_eval;
    
    dNa_i_eval = -1.d0*(i_Na_eval+i_b_Na_eval+3.d0*i_NaK_eval+3.d0*i_NaCa_eval)/...
                (1.d0*V_c*F)*Cm;
    Na_i(i)=Na_iold+dt*dNa_i_eval;
    
    
    if (2<t && t<2.3)
        Iapp=280;
    else
        Iapp=0;
    end

    Itot=Iapp-(i_K1_eval+i_to_eval+i_Kr_eval+i_Ks_eval+i_CaL_eval+i_NaK_eval+i_Na_eval+...
           i_b_Na_eval+i_NaCa_eval+i_b_Ca_eval+i_p_K_eval+i_p_Ca_eval);

    V(i)=Vold+dt*Itot;

end

plot(T,V)
