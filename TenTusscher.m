classdef TenTusscher < handle
    properties
        % material properties
        g_CaL
        g_bca
        Buf_c 
        Buf_sr
        Buf_ss 
        Ca_o 
        EC 
        K_buf_c 
        K_buf_sr 
        K_buf_ss
        K_up
        V_leak
        V_rel 
        V_sr 
        V_ss
        V_xfer 
        Vmax_up
        k1_prime
        k2_prime 
        k3
        k4 
        max_sr 
        min_sr 
        K_pCa 
        g_pCa 
        g_K1 
        Cm 
        F
        R 
        T 
        V_c
        
        K_o      
        
        g_pK 
        g_Kr 
        P_kna
        g_Ks 
        g_bna 
        K_NaCa
        K_sat 
        Km_Ca
        Km_Nai
        alpha 
        gamma
        Na_o 
        K_mNa
        K_mk 
        P_NaK 
        g_to 
        g_Na 

        % initial conditions
        d0
        f20
        fCass0
        pf0
        Ca_sr0
        Ca_i0
        Ca_ss0
        Rel0
        h0
        j0
        m0
        K_i0
        xr10
        xr20
        xs0
        Na_i0
        pr0
        s0

        % TT functions
        E_Ca 
        E_K 
        i_CaL 
        i_b_Ca 
        d_inf 
        alpha_d 
        beta_d 
        gamma_d 
        tau_d 
        
        f2_inf 
        tau_f2 
        
        fCass_inf 
        tau_fCass 
        
        f_inf 
        tau_f 
           
        kcasr 
        k1 
        O 
        i_rel 
        i_up 
        i_leak 
        i_xfer 
        
        k2 
        dRel 
        
        Ca_i_bufc 
        Ca_sr_bufsr 
        Ca_ss_bufss 
        i_p_Ca 
        
        i_NaCa 
        
        dCa_i 
        
        dCa_sr 
        
        dCa_ss 
        
        E_Na 
        i_Na 
        
        h_inf 
        alpha_h_l
        alpha_h 
        beta_h_l
        beta_h_g
        beta_h 
        tau_h 
        
        j_inf 
        alpha_j_l
        alpha_j_g
        alpha_j 
        beta_j_l
        beta_j_g
        beta_j 
        tau_j 
        
        m_inf 
        alpha_m 
        beta_m 
        tau_m 
        
        xr1_inf 
        alpha_xr1 
        beta_xr1 
        tau_xr1 
        
        xr2_inf 
        alpha_xr2 
        beta_xr2 
        tau_xr2 
        
        xs_inf 
        alpha_xs 
        beta_xs 
        tau_xs 
        
        r_inf 
        tau_r 
        
        s_inf 
        tau_s 
        
        alpha_K1 
        beta_K1 
        xK1_inf 
        i_K1 
        
        i_to 
        i_Kr 
        
        E_Ks 
        i_Ks 
        
        i_NaK 
        i_b_Na 
        i_p_K 


        % variables()
        Do
        F2o
        FCASSo
        PFo
        CASRo
        CAIo
        CASSo
        RELo
        Ho
        Jo
        Mo
        KIo
        XR1o
        XR2o
        XSo
        NAIo
        PRo
        So
        Dn
        F2n
        FCASSn
        PFn
        CASRn
        CAIn
        CASSn
        RELn
        Hn
        Jn
        Mn
        KIn
        XR1n
        XR2n
        XSn
        NAIn
        PRn
        Sn
        
        %
        dt

    end
    methods
        function obj=TenTusscher(V,dt)
            obj.g_CaL = 0.0000398d0;
            obj.g_bca = 0.000592d0;
            obj.Buf_c = 0.2d0;
            obj.Buf_sr = 10.0d0;
            obj.Buf_ss = 0.4d0;
            obj.Ca_o = 2.0d0;
            obj.EC = 1.5d0;
            obj.K_buf_c = 0.001d0;
            obj.K_buf_sr = 0.3d0;
            obj.K_buf_ss = 0.00025d0;
            obj.K_up = 0.00025d0;
            obj.V_leak = 0.00036d0;
            obj.V_rel = 0.102d0;
            obj.V_sr = 0.001094d0;
            obj.V_ss = 0.00005468d0;
            obj.V_xfer = 0.0038d0;
            obj.Vmax_up = 0.006375d0;
            obj.k1_prime = 0.15d0;
            obj.k2_prime = 0.045d0;
            obj.k3 = 0.06d0;
            obj.k4 = 0.005d0;
            obj.max_sr = 2.5d0;
            obj.min_sr = 1.0d0;
            obj.K_pCa = 0.0005d0;
            obj.g_pCa = 0.1238d0;
            obj.g_K1 = 5.405d0;
            obj.Cm = 0.185d0;
            obj.F = 96485.3415d0;
            obj.R = 8314.472d0;
            obj.T = 310.0d0;
            obj.V_c = 0.016404d0;
            
            obj.K_o=5.4;
            
            obj.g_pK = 0.0146d0;
            obj.g_Kr = 0.153d0;
            obj.P_kna = 0.03d0;
            obj.g_Ks = 0.392d0;
            obj.g_bna = 0.00029d0;
            obj.K_NaCa = 1000.0d0;
            obj.K_sat = 0.1d0;
            obj.Km_Ca = 1.38d0;
            obj.Km_Nai = 87.5d0;
            obj.alpha = 2.5d0;
            obj.gamma = 0.35d0;
            obj.Na_o = 140.0d0;
            obj.K_mNa = 40.0d0;
            obj.K_mk = 1.0d0;
            obj.P_NaK = 2.724d0;
            obj.g_to = 0.294d0;
            obj.g_Na = 14.838d0;

            obj.d0=3.373000000000000e-05;
            obj.f20=9.755000000000000e-01;
            obj.fCass0=9.953000000000000e-01;
            obj.pf0=7.887999999999999e-01;
            obj.Ca_sr0=3.640000000000000e+00;
            obj.Ca_i0=1.260000000000000e-04;
            obj.Ca_ss0=3.600000000000000e-04;
            obj.Rel0=9.073000000000000e-01;
            obj.h0=7.444000000000000e-01;
            obj.j0=7.045000000000000e-01;
            obj.m0=1.720000000000000e-03;
            obj.K_i0=1.368900000000000e+02;
            obj.xr10=6.210000000000000e-03;
            obj.xr20=4.712000000000000e-01;
            obj.xs0=9.500000000000000e-03;
            obj.Na_i0=8.603999999999999e+00;
            obj.pr0=2.420000000000000e-08;
            obj.s0=9.999980000000001e-01;

            
            obj.E_Ca =@(Ca_i) 0.5d0.*obj.R.*obj.T./obj.F.*log(obj.Ca_o./Ca_i);
            obj.E_K =@(K_i) obj.R.*obj.T./obj.F.*log(obj.K_o./K_i);

            obj.i_CaL =@(V,d,pf,f2,fCass,Ca_ss) obj.g_CaL.*d.*pf.*f2.*fCass.*4.d0.*(V-15.d0).*obj.F.^2./(obj.R.*obj.T).*(0.25d0.*Ca_ss.*exp(2.d0.*(V-15.d0).*obj.F./(obj.R.*obj.T))-obj.Ca_o)./(exp(2.d0.*(V-15.d0).*obj.F./(obj.R.*obj.T))-1.d0);
            obj.i_b_Ca =@(V,Ca_i) obj.g_bca.*(V-obj.E_Ca(Ca_i));
 
            obj.d_inf =@(V) 1.d0./(1.d0+exp((-8.d0-V)./7.5d0));
            obj.alpha_d =@(V) 1.4d0./(1.d0+exp((-35.d0-V)./13.d0))+0.25d0;
            obj.beta_d =@(V) 1.4d0./(1.d0+exp((V+5.d0)./5.d0));
            obj.gamma_d =@(V) 1.d0./(1.d0+exp((50.d0-V)./20.d0));
            obj.tau_d =@(V) 1.d0.*obj.alpha_d(V).*obj.beta_d(V)+obj.gamma_d(V);

            obj.f2_inf =@(V) 0.67d0./(1.d0+exp((V+35.d0)./7.d0))+0.33d0;
            obj.tau_f2 =@(V) 562.d0.*exp(-(V+27.d0).^2./240.d0)+31.d0./(1.d0+exp((25.d0-V)./10.d0))+80.d0./(1.d0+exp((V+30.d0)./10.d0)); 

            obj.fCass_inf =@(V,Ca_ss) 0.6d0./(1.d0+(Ca_ss./0.05d0).^2)+0.4d0;
            obj.tau_fCass =@(V,Ca_ss) 80.d0./(1.d0+(Ca_ss./0.05d0).^2)+2.d0;
            
            obj.f_inf =@(V) 1.d0./(1.d0+exp((V+20.d0)./7.d0));
            obj.tau_f =@(V) 1102.5d0.*exp(-(V+27.d0).^2./225.d0)+200.d0./(1.d0+exp((13.d0-V)./10.d0))+180.d0./(1.d0+exp((V+30.d0)./10.d0))+20.d0;
 
            obj.kcasr =@(Ca_sr) obj.max_sr-(obj.max_sr-obj.min_sr)./(1.d0+(obj.EC./Ca_sr).^2);
            obj.k1 =@(Ca_sr) obj.k1_prime./obj.kcasr(Ca_sr);
            obj.O =@(Rel, Ca_sr, Ca_ss) obj.k1(Ca_sr).*Ca_ss.^2.*Rel./(obj.k3+obj.k1(Ca_sr).*Ca_ss.^2);
            obj.i_rel =@(Rel, Ca_sr, Ca_ss) obj.V_rel.*obj.O(Rel, Ca_sr, Ca_ss).*(Ca_sr-Ca_ss);
            obj.i_up =@(Ca_i) obj.Vmax_up./(1.d0+(obj.K_up./Ca_i).^2);
            obj.i_leak =@(Ca_sr, Ca_i) obj.V_leak.*(Ca_sr-Ca_i);
            obj.i_xfer =@(Ca_ss, Ca_i) obj.V_xfer.*(Ca_ss-Ca_i);
 
            obj.k2 =@(Ca_sr) obj.k2_prime.*obj.kcasr(Ca_sr);
            obj.dRel =@(Ca_ss, Rel, Ca_sr) (-obj.k2(Ca_sr).*Ca_ss.*Rel+obj.k4.*(-Rel+1.0));

            obj.Ca_i_bufc =@(Ca_i) 1.d0./(1.d0+obj.Buf_c.*obj.K_buf_c./(Ca_i+obj.K_buf_c).^2);
            obj.Ca_sr_bufsr =@(Ca_sr) 1.d0./(1.d0+obj.Buf_sr.*obj.K_buf_sr./(Ca_sr+obj.K_buf_sr).^2);
            obj.Ca_ss_bufss =@(Ca_ss) 1.d0./(1.d0+obj.Buf_ss.*obj.K_buf_ss./(Ca_ss+obj.K_buf_ss).^2);
            obj.i_p_Ca =@(Ca_i) obj.g_pCa.*Ca_i./(Ca_i+obj.K_pCa);

            obj.i_NaCa =@(V, Na_i, Ca_i) obj.K_NaCa.*(exp(obj.gamma.*V.*obj.F./(obj.R.*obj.T)).*Na_i.^3.*obj.Ca_o-exp((obj.gamma-1.d0).*V.*obj.F./(obj.R.*obj.T)).*obj.Na_o.^3.*Ca_i.*obj.alpha)./((obj.Km_Nai.^3+obj.Na_o.^3).*(obj.Km_Ca+obj.Ca_o).*(1.d0+obj.K_sat.*exp((obj.gamma-1.d0).*V.*obj.F./(obj.R.*obj.T))));

            obj.dCa_i =@(Na_i, Ca_ss, Ca_sr, Ca_i, V) obj.Ca_i_bufc(Ca_i).*((obj.i_leak(Ca_sr, Ca_i)-obj.i_up(Ca_i)).*obj.V_sr./obj.V_c+obj.i_xfer(Ca_ss, Ca_i)-1.d0.*(obj.i_b_Ca(V,Ca_i)+obj.i_p_Ca(Ca_i)-2.d0.*obj.i_NaCa(V, Na_i, Ca_i)).*obj.Cm./(2.d0.*1.d0.*obj.V_c.*obj.F));
            obj.dCa_sr =@(Rel, Ca_sr, Ca_i, Ca_ss) obj.Ca_sr_bufsr(Ca_sr).*(obj.i_up(Ca_i)-(obj.i_rel(Rel, Ca_sr, Ca_ss)+obj.i_leak(Ca_sr, Ca_i)));
            obj.dCa_ss =@(Rel,V,d,pf,f2,fCass,Ca_ss,Ca_sr,Ca_i) obj.Ca_ss_bufss(Ca_ss).*(-1.d0.*obj.i_CaL(V,d,pf,f2,fCass,Ca_ss).*obj.Cm./(2.d0.*1.d0.*obj.V_ss.*obj.F)+obj.i_rel(Rel, Ca_sr, Ca_ss).*obj.V_sr./obj.V_ss-obj.i_xfer(Ca_ss, Ca_i).*obj.V_c./obj.V_ss);

            obj.E_Na =@(Na_i) obj.R.*obj.T./obj.F.*log(obj.Na_o./Na_i);
            obj.i_Na =@(V, m, h, j, Na_i) obj.g_Na.*m.^3.*h.*j.*(V-obj.E_Na(Na_i));

            obj.h_inf =@(V) 1.d0./(1.d0+exp((V+71.55d0)./7.43d0)).^2;
            obj.alpha_h_l=@(V) 0.057d0.*exp(-(V+80.d0)./6.8d0);
            obj.alpha_h =@(V) (obj.alpha_h_l(V).*(V<-40.0)+0.d0.*(V>-40.0));
            obj.beta_h_l=@(V) (2.7d0.*exp(0.079d0.*V)+3.1d+5.*exp(0.3485d0.*V) );
            obj.beta_h_g=@(V) (0.77d0./(0.13d0.*(1.d0+exp((V+10.66d0)./(-11.1d0)))) );
            obj.beta_h =@(V) (obj.beta_h_l(V).*(V<-40.0)+obj.beta_h_g(V).*(V>-40.0));
            obj.tau_h =@(V) 1.d0./(obj.alpha_h(V) + obj.beta_h(V));

            obj.j_inf =@(V) 1.d0./(1.d0+exp((V+71.55d0)./7.43d0)).^2;
            obj.alpha_j_l=@(V) (-25428.d0.*exp(0.2444d0.*V)-6.948d-6.*exp(-0.04391d0.*V)).*(V+37.78d0)./1.d0./(1.d0+exp(0.311d0.*(V+79.23d0)));
            obj.alpha_j_g=@(V) (0.d0);
            obj.alpha_j =@(V) (obj.alpha_j_l(V).*(V<-40.0)+obj.alpha_j_g(V).*(V>-40.0));
            obj.beta_j_l=@(V) (0.02424d0.*exp(-0.01052d0.*V)./ (1.d0+exp(-0.1378d0.*(V+40.14d0))) );
            obj.beta_j_g=@(V) (0.6d0.*exp(0.057d0.*V)./(1.d0+exp(-0.1d0.*(V+32.d0))) );
            obj.beta_j =@(V) (obj.beta_j_l(V).*(V<-40.0)+obj.beta_j_g(V).*(V>-40.0));
            obj.tau_j =@(V) 1.d0./(obj.alpha_j(V)+obj.beta_j(V));

            obj.m_inf =@(V) 1.d0./(1.d0+exp((-56.86d0-V)./9.03d0)).^2;
            obj.alpha_m =@(V) 1.d0./(1.d0+exp((-60.d0-V)./5.d0));
            obj.beta_m =@(V) 0.1d0./(1.d0+exp((V+35.d0)./5.d0))+0.1d0./(1.d0+exp((V-50.d0)./200.d0));
            obj.tau_m =@(V) 1.d0.*obj.alpha_m(V).*obj.beta_m(V);

            obj.xr1_inf =@(V) 1.d0./(1.d0+exp((-26.d0-V)./7.d0));
            obj.alpha_xr1 =@(V) 450.d0./(1.d0+exp((-45.d0-V)./10.d0));
            obj.beta_xr1 =@(V) 6.d0./(1.d0+exp((V+30.d0)./11.5d0));
            obj.tau_xr1 =@(V) 1.d0.*obj.alpha_xr1(V).*obj.beta_xr1(V);

            obj.xr2_inf =@(V) 1.d0./(1.d0+exp((V+88.d0)./24.d0));
            obj.alpha_xr2 =@(V) 3.d0./(1.d0+exp((-60.d0-V)./20.d0));
            obj.beta_xr2 =@(V) 1.12d0./(1.d0+exp((V-60.d0)./20.d0));
            obj.tau_xr2 =@(V) 1.d0.*obj.alpha_xr2(V).*obj.beta_xr2(V);
  
            obj.xs_inf =@(V) 1.d0./(1.d0+exp((-5.d0-V)./14.d0));
            obj.alpha_xs =@(V) 1400.d0./sqrt(1.d0+exp((5.d0-V)./6.d0));
            obj.beta_xs =@(V) 1.d0./(1.d0+exp((V-35.d0)./15.d0));
            obj.tau_xs =@(V) 1.d0.*obj.alpha_xs(V).*obj.beta_xs(V)+80.d0;

            obj.r_inf =@(V) 1.d0./(1.d0+exp((20.d0-V)./6.d0));
            obj.tau_r =@(V) 9.5d0.*exp(-(V+40.d0).^2./1800.d0)+0.8d0;

            obj.s_inf =@(V) 1.d0./(1.d0+exp((V+20.d0)./5.d0));
            obj.tau_s =@(V) 85.d0.*exp(-(V+45.d0).^2./320.d0)+5.d0./(1.d0+exp((V-20.d0)./5.d0))+3.d0;

            obj.alpha_K1 =@(V,K_i) 0.1d0./(1.d0+exp(0.06d0.*(V-obj.E_K(K_i)-200.d0)));
            obj.beta_K1 =@(V,K_i) (3.d0.*exp(0.0002d0.*(V-obj.E_K(K_i)+100.d0))+exp(0.1d0.*(V-obj.E_K(K_i)-10.d0)))./(1.d0+exp(-0.5d0.*(V-obj.E_K(K_i))));
            obj.xK1_inf =@(V,K_i) obj.alpha_K1(V,K_i)./(obj.alpha_K1(V,K_i)+obj.beta_K1(V,K_i));
            obj.i_K1 =@(V,K_i) obj.g_K1.*obj.xK1_inf(V,K_i).*sqrt(obj.K_o./5.4d0).*(V-obj.E_K(K_i));

            obj.i_to =@(V, pr, s, K_i) obj.g_to.*pr.*s.*(V-obj.E_K(K_i));
            obj.i_Kr =@(V, xr1, xr2, K_i) obj.g_Kr.*sqrt(obj.K_o./5.4d0).*xr1.*xr2.*(V-obj.E_K(K_i));

            obj.E_Ks =@(Na_i,K_i) obj.R.*obj.T./obj.F.*log((obj.K_o+obj.P_kna.*obj.Na_o)./(K_i+obj.P_kna.*Na_i));
            obj.i_Ks =@(V, Na_i, xs,K_i) obj.g_Ks.*xs.^2.*(V-obj.E_Ks(Na_i,K_i));

            obj.i_NaK =@(V,Na_i) obj.P_NaK.*obj.K_o./(obj.K_o+obj.K_mk).*Na_i./(Na_i+obj.K_mNa)./(1.d0+0.1245d0.*exp(-0.1d0.*V.*obj.F./(obj.R.*obj.T))+0.0353d0.*exp(-V.*obj.F./(obj.R.*obj.T)));  
            obj.i_b_Na =@(V,Na_i) obj.g_bna.*(V-obj.E_Na(Na_i));
            obj.i_p_K =@(V,K_i) obj.g_pK.*(V-obj.E_K(K_i))./(1.d0+exp((25.d0-V)./5.98d0));

            obj.Do=zeros(size(V));
            obj.F2o=zeros(size(V));
            obj.FCASSo=zeros(size(V));
            obj.PFo=zeros(size(V));
            obj.CASRo=zeros(size(V));
            obj.CAIo=zeros(size(V));
            obj.CASSo=zeros(size(V));
            obj.RELo=zeros(size(V));
            obj.Ho=zeros(size(V));
            obj.Jo=zeros(size(V));
            obj.Mo=zeros(size(V));
            obj.KIo=zeros(size(V));
            obj.XR1o=zeros(size(V));
            obj.XR2o=zeros(size(V));
            obj.XSo=zeros(size(V));
            obj.NAIo=zeros(size(V));
            obj.PRo=zeros(size(V));
            obj.So=zeros(size(V));

            obj.Dn=zeros(size(V));
            obj.F2n=zeros(size(V));
            obj.FCASSn=zeros(size(V));
            obj.PFn=zeros(size(V));
            obj.CASRn=zeros(size(V));
            obj.CAIn=zeros(size(V));
            obj.CASSn=zeros(size(V));
            obj.RELn=zeros(size(V));
            obj.Hn=zeros(size(V));
            obj.Jn=zeros(size(V));
            obj.Mn=zeros(size(V));
            obj.KIn=zeros(size(V));
            obj.XR1n=zeros(size(V));
            obj.XR2n=zeros(size(V));
            obj.XSn=zeros(size(V));
            obj.NAIn=zeros(size(V));
            obj.PRn=zeros(size(V));
            obj.Sn=zeros(size(V));

            obj.Dn(:,1)=obj.d0;
            obj.F2n(:,1)=obj.f20;
            obj.FCASSn(:,1)=obj.fCass0;
            obj.PFn(:,1)=obj.pf0;
            obj.CASRn(:,1)=obj.Ca_sr0;
            obj.CAIn(:,1)=obj.Ca_i0;
            obj.CASSn(:,1)=obj.Ca_ss0;
            obj.RELn(:,1)=obj.Rel0;
            obj.Hn(:,1)=obj.h0;
            obj.Jn(:,1)=obj.j0;
            obj.Mn(:,1)=obj.m0;
            obj.KIn(:,1)=obj.K_i0;
            obj.XR1n(:,1)=obj.xr10;
            obj.XR2n(:,1)=obj.xr20;
            obj.XSn(:,1)=obj.xs0;
            obj.NAIn(:,1)=obj.Na_i0;
            obj.PRn(:,1)=obj.pr0;
            obj.Sn(:,1)=obj.s0;
            % dt
            obj.dt=dt;

        end


        function Iion=solveTimestep(obj,Vold,i)
            d_inf_eval = obj.d_inf(Vold);
            tau_d_eval = obj.tau_d(Vold);
            obj.Dn=d_inf_eval+(obj.Do-d_inf_eval).*exp(-obj.dt./tau_d_eval);
        
            f2_inf_eval = obj.f2_inf(Vold);
            tau_f2_eval = obj.tau_f2(Vold);
            obj.F2n=f2_inf_eval+(obj.F2o-f2_inf_eval).*exp(-obj.dt./tau_f2_eval);
        
            fCass_inf_eval = obj.fCass_inf(Vold,obj.CASSo);
            tau_fCass_eval = obj.tau_fCass(Vold,obj.CASSo);
            obj.FCASSn=fCass_inf_eval+(obj.FCASSo-fCass_inf_eval).*exp(-obj.dt./tau_fCass_eval);
            
            f_inf_eval = obj.f_inf(Vold);
            tau_f_eval = obj.tau_f(Vold);
            obj.PFn=f_inf_eval+(obj.PFo-f_inf_eval).*exp(-obj.dt./tau_f_eval);
        
            dRel_eval = obj.dRel(obj.CASSo, obj.RELo, obj.CASRo);
            obj.RELn=obj.RELo+obj.dt.*dRel_eval;
        
            dCa_i_eval = obj.dCa_i(obj.NAIo, obj.CASSo, obj.CASRo, obj.CAIo, Vold);
            obj.CAIn=obj.CAIo+obj.dt.*dCa_i_eval;
        
            dCa_sr_eval = obj.dCa_sr(obj.RELo, obj.CASRo, obj.CAIo, obj.CASSo);
            obj.CASRn=obj.CASRo+obj.dt.*dCa_sr_eval;
        
            dCa_ss_eval=obj.dCa_ss(obj.RELo,Vold,obj.Do,obj.PFo,obj.F2o,obj.FCASSo,obj.CASSo,obj.CASRo,obj.CAIo);
            obj.CASSn=obj.CASSo+obj.dt.*dCa_ss_eval;
        
            h_inf_eval=obj.h_inf(Vold);
            tau_h_eval=obj.tau_h(Vold);
            obj.Hn=h_inf_eval+(obj.Ho-h_inf_eval).*exp(-obj.dt./tau_h_eval);
        
            j_inf_eval=obj.j_inf(Vold);
            tau_j_eval=obj.tau_j(Vold);
            obj.Jn=j_inf_eval+(obj.Jo-j_inf_eval).*exp(-obj.dt./tau_j_eval);
        
            m_inf_eval=obj.m_inf(Vold);
            tau_m_eval=obj.tau_m(Vold);
            obj.Mn=m_inf_eval+(obj.Mo-m_inf_eval).*exp(-obj.dt./tau_m_eval);
        
            xr1_inf_eval=obj.xr1_inf(Vold);
            tau_xr1_eval=obj.tau_xr1(Vold);
            obj.XR1n=xr1_inf_eval+(obj.XR1o-xr1_inf_eval).*exp(-obj.dt./tau_xr1_eval);
        
            xr2_inf_eval=obj.xr2_inf(Vold);
            tau_xr2_eval=obj.tau_xr2(Vold);
            obj.XR2n=xr2_inf_eval+(obj.XR2o-xr2_inf_eval).*exp(-obj.dt./tau_xr2_eval);
        
            xs_inf_eval=obj.xs_inf(Vold);
            tau_xs_eval=obj.tau_xs(Vold);
            obj.XSn=xs_inf_eval+(obj.XSo-xs_inf_eval).*exp(-obj.dt./tau_xs_eval);
        
            r_inf_eval=obj.r_inf(Vold);
            tau_r_eval=obj.tau_r(Vold);
            obj.PRn=r_inf_eval+(obj.PRo-r_inf_eval).*exp(-obj.dt./tau_r_eval);
        
            s_inf_eval=obj.s_inf(Vold);
            tau_s_eval=obj.tau_s(Vold);
            obj.Sn=s_inf_eval+(obj.So-s_inf_eval).*exp(-obj.dt./tau_s_eval);


            i_K1_eval=obj.i_K1(Vold,obj.KIo);
            i_to_eval=obj.i_to(Vold, obj.PRn, obj.Sn,obj.KIo);
            i_Kr_eval=obj.i_Kr(Vold,obj.XR1n, obj.XR2n,obj.KIo);
            i_Ks_eval=obj.i_Ks(Vold, obj.NAIo, obj.XSn,obj.KIo);
            i_CaL_eval=obj.i_CaL(Vold,obj.Dn,obj.PFn,obj.F2n,obj.FCASSn,obj.CASSn);
            i_NaK_eval=obj.i_NaK(Vold,obj.NAIo);
            i_Na_eval=obj.i_Na(Vold,obj.Mn,obj.Hn,obj.Jn,obj.NAIo);
            i_b_Na_eval=obj.i_b_Na(Vold,obj.NAIo);
            i_NaCa_eval=obj.i_NaCa(Vold,obj.NAIo,obj.CAIn);
            i_b_Ca_eval=obj.i_b_Ca(Vold,obj.CAIn);
            i_p_Ca_eval=obj.i_p_Ca(obj.CAIn);
            i_p_K_eval=obj.i_p_K(Vold,obj.KIo);   
            
            dK_i_eval = -1.d0.*(i_K1_eval+i_to_eval+i_Kr_eval+i_Ks_eval+i_p_K_eval-2.d0.*i_NaK_eval)./(1.d0.*obj.V_c.*obj.F).*obj.Cm;
            obj.KIn = obj.KIo+obj.dt.*dK_i_eval;
            
            dNa_i_eval = -1.d0*(i_Na_eval+i_b_Na_eval+3.d0.*i_NaK_eval+3.d0.*i_NaCa_eval)./(1.d0*obj.V_c*obj.F)*obj.Cm;
            obj.NAIn=obj.NAIo+obj.dt.*dNa_i_eval;

            Iion=i_K1_eval+i_to_eval+i_Kr_eval+i_Ks_eval+i_CaL_eval+i_NaK_eval+i_Na_eval+i_b_Na_eval+i_NaCa_eval+i_b_Ca_eval+i_p_K_eval+i_p_Ca_eval;
        end
        
        
        function Iion=getCurr(obj,Vold,i)
            obj.Do=obj.Dn;
            obj.F2o=obj.F2n;
            obj.FCASSo=obj.FCASSn;
            obj.PFo=obj.PFn;
            obj.CASRo=obj.CASRn;
            obj.CAIo=obj.CAIn;
            obj.CASSo=obj.CASSn;
            obj.RELo=obj.RELn;
            obj.Ho=obj.Hn;
            obj.Jo=obj.Jn;
            obj.Mo=obj.Mn;
            obj.KIo=obj.KIn;
            obj.XR1o=obj.XR1n;
            obj.XR2o=obj.XR2n;
            obj.XSo=obj.XSn;
            obj.NAIo=obj.NAIn;
            obj.PRo=obj.PRn;
            obj.So=obj.Sn;

            Iion=obj.solveTimestep(Vold,i);
        end
    end
end
