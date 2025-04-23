classdef Paci < handle
    properties
        
        % Constants
        F 
        R 
        T 
        
        V_SR 
        Vc   
        Cm   
        
        Nao 
        Ko  
        Cao 
        
        Ki 
        
        g_Na        
        myCoefTauM  
        tauINaL     
        GNaLmax     
        Vh_hLate    
        g_f         
        fNa         
        fK          
        g_CaL       
        constf2     
        g_to        
        g_Ks        
        L0           
        Q            
        g_Kr         
        V_half       
        g_K1        
        KmCa        
        KmNai       
        Ksat        
        gamma       
        alpha       
        kNaCa       
        Km_K        
        Km_Na       
        PNaK        
        KPCa        
        g_PCa       
        g_b_Na      
        g_b_Ca      
        VmaxUp		
        Kup			
        V_leak		
        g_irel_max	
        RyRa1       
        RyRa2       
        RyRahalf    
        RyRohalf   
        RyRchalf   
        RyRtauadapt
        Buf_C      
        Buf_SR     
        Kbuf_C     
        Kbuf_SR    


        % initial conditions
        CaSR0
        Cai0
        g0
        d0
        f10
        f20
        fCa0
        Xr10
        Xr20
        Xs0
        h0
        j0
        m0
        Xf0
        q0
        r0
        Nai0
        mL0
        hL0
        RyRa0
        RyRo0
        RyRc0

        % Nernst potential
        E_Na 
        E_Ca 
        E_K  
        PkNa 
        E_Ks 


        % Paci functions
        m_inf       
        tau_m       
        h_inf       
        tau_h       
        j_inf       
        tau_j       
        
        i_Na        
        
        
        mL_inf     
        alpha_m_L   
        beta_m_L    
        tau_mL     
        hL_inf     
        tau_hL     
        
        i_NaL       
        
        Xf_inf      
        tau_Xf      
        
        i_fK        
        i_fNa       
        i_f         
        
        d_inf       
        alpha_d     
        beta_d      
        gamma_d     
        tau_d       
        
        f1_inf      
        constf1     
        tau_f1      
        
        f2_inf      
        tau_f2      
        
        alpha_fCa   
        beta_fCa    
        gamma_fCa   
        fCa_inf     
        constfCa    
        tau_fCa     
        
        i_CaL       
        
        q_inf       
        tau_q       
        r_inf       
        tau_r       
        
        i_to        
        
        Xs_inf      
        alpha_Xs    
        beta_Xs     
        tau_Xs      
        
        i_Ks        
        
        Xr1_inf      
        alpha_Xr1    
        beta_Xr1     
        tau_Xr1      
        
        Xr2_inf      
        alpha_Xr2    
        beta_Xr2     
        tau_Xr2      
        
        i_Kr         
        
        alpha_K1    
        beta_K1     
        XK1_inf     
        
        i_K1        
        
        i_NaCa      
        
        i_NaK       
        
        i_PCa       
        
        i_b_Na      
        
        i_b_Ca      
        
        i_up        
        
        i_leak      
        
        RyRainfss   
        RyRoinfss   
        RyRtauact   
        
        RyRcinfss   
        RyRtauinact 
        
        RyRSRCass   
        i_rel       
        
        Cai_bufc    
        Ca_SR_bufSR 



        % variables()
        CaSR_o
        Cai_o
        g_o
        d_o
        f1_o
        f2_o
        fCa_o
        Xr1_o
        Xr2_o
        Xs_o
        h_o
        j_o
        m_o
        Xf_o
        q_o
        r_o
        Nai_o
        mL_o
        hL_o
        RyRa_o
        RyRo_o
        RyRc_o

        CaSR_n
        Cai_n
        g_n
        d_n
        f1_n
        f2_n
        fCa_n
        Xr1_n
        Xr2_n
        Xs_n
        h_n
        j_n
        m_n
        Xf_n
        q_n
        r_n
        Nai_n
        mL_n
        hL_n
        RyRa_n
        RyRo_n
        RyRc_n
        
        %
        dt

    end
    methods
        function obj=Paci(V,dt)
            % Constants
            obj.F = 96485.3415;     % coulomb_per_mole (in model_parameters)
            obj.R = 8.314472;       % joule_per_mole_kelvin (in model_parameters)
            obj.T = 310.0;          % kelvin (in model_parameters) %37°C
            
            obj.V_SR = 583.73;        % micrometre_cube (in model_parameters)
            obj.Vc   = 8800.0;        % micrometre_cube (in model_parameters)
            obj.Cm   = 9.87109e-11;   % farad (in model_parameters)
            
            obj.Nao = 151.0; % millimolar (in model_parameters)
            obj.Ko  = 5.4;   % millimolar (in model_parameters)
            obj.Cao = 1.8;   % millimolar (in model_parameters)
            
            obj.Ki = 150.0;   % millimolar (in model_parameters)
            
            obj.g_Na        = 6447.1896;
            obj.myCoefTauM  = 1;
            obj.tauINaL     = 200; 
            obj.GNaLmax     = 2.3.*7.5; 
            obj.Vh_hLate    = 87.61;
            obj.g_f         = 22.2763088;
            obj.fNa         = 0.37;
            obj.fK          = 1 - obj.fNa;
            obj.g_CaL       = 8.635702e-5;
            obj.constf2     = 1.0;
            obj.g_to        = 29.9038;  
            obj.g_Ks        = 2.041;
            obj.L0           = 0.025;  
            obj.Q            = 2.3;    
            obj.g_Kr         = 29.8667;   
            obj.V_half       = 1000.0.*(-obj.R.*obj.T./(obj.F.*obj.Q).*log((1.0+obj.Cao./2.6).^4.0./(obj.L0.*(1.0+obj.Cao./0.58).^4.0))-0.019);
            obj.g_K1        = 28.1492;  
            obj.KmCa        = 1.38;  
            obj.KmNai       = 87.5;  
            obj.Ksat        = 0.1;    
            obj.gamma       = 0.35;  
            obj.alpha       = 2.16659;
            obj.kNaCa       =  6514.47574; 
            obj.Km_K        = 1.0;  
            obj.Km_Na       = 40.0; 
            obj.PNaK        = 2.74240;
            obj.KPCa        = 0.0005;   
            obj.g_PCa       = 0.4125;  
            obj.g_b_Na      = 1.14;
            obj.g_b_Ca      = 0.8727264; 
            obj.VmaxUp		= 0.82205;
            obj.Kup			=  4.40435e-4;
            obj.V_leak		= 4.48209e-4;
            obj.g_irel_max	= 55.808061;
            obj.RyRa1       = 0.05169;
            obj.RyRa2       = 0.050001;
            obj.RyRahalf    = 0.02632;
            obj.RyRohalf    = 0.00944;
            obj.RyRchalf    = 0.00167;
            obj.RyRtauadapt = 1; 
            obj.Buf_C       = 0.25;   % millimolar (in calcium_dynamics)
            obj.Buf_SR      = 10.0;   % millimolar (in calcium_dynamics)
            obj.Kbuf_C      = 0.001;   % millimolar (in calcium_dynamics)
            obj.Kbuf_SR     = 0.3;   % millimolar (in calcium_dynamics)
                
            % initial conditions
            obj.CaSR0=0.115107956531682;
            obj.Cai0=1.76731736262123e-05;
            obj.g0=0;
            obj.d0=0.000101671059827320;
            obj.f10=0.979415292087528;
            obj.f20=0.999979398897691;
            obj.fCa0=0.998948743326411;
            obj.Xr10=0.00894672801754690;
            obj.Xr20=0.427809818209464;
            obj.Xs0=0.0341200913071620;
            obj.h0=0.736080688718430;
            obj.j0=0.742052141479516;
            obj.m0=0.0441025126145443;
            obj.Xf0=0.192854555974207;
            obj.q0=0.829179626325527;
            obj.r0=0.00601271295586856;
            obj.Nai0=9.00417069274294;
            obj.mL0=0.00297784296969314;
            obj.hL0=0.135864730293044;
            obj.RyRa0=0.0297153296103769;
            obj.RyRo0=0.000710450936345816;
            obj.RyRc0=0.948119074119825;

            % Nernst potential
            obj.E_Na =@(Nai) obj.R.*obj.T./obj.F.*log(obj.Nao./Nai);
            obj.E_Ca =@(Cai) 0.5.*obj.R.*obj.T./obj.F.*log(obj.Cao./Cai);
            obj.E_K  = obj.R.*obj.T./obj.F.*log(obj.Ko./obj.Ki);
            obj.PkNa = 0.03;   % dimensionless (in electric_potentials)
            obj.E_Ks =@(Nai) obj.R.*obj.T./obj.F.*log((obj.Ko+obj.PkNa.*obj.Nao)./(obj.Ki+obj.PkNa.*Nai));
            
            
            % funzioni
            obj.m_inf       =@(V) 1 ./ (1 + exp((V.*1000 + 39)./-11.2));
            obj.tau_m       =@(V) (0.00001 + 0.00013.*exp(-((V.*1000 + 48)./15).^2) + 0.000045 ./ (1 + exp((V.*1000 + 42)./-5)));
            obj.h_inf       =@(V) 1 ./ (1 + exp((V.*1000 + 66.5)./6.8));
            obj.tau_h       =@(V) (0.00007 + 0.034 ./ (1 + exp((V.*1000 + 41)./5.5) + exp(-(V.*1000 + 41)./14)) + 0.0002 ./ (1 + exp(-(V.*1000 + 79)./14)));
            obj.j_inf       = obj.h_inf;
            obj.tau_j       =@(V) 10.*(0.0007 + 0.15 ./ (1 + exp((V.*1000 + 41)./5.5) + exp(-(V.*1000 + 41)./14)) + 0.002 ./ (1 + exp(-(V.*1000 + 79)./14)));
            
            obj.i_Na        =@(V,m,h,j,Nai)  obj.g_Na.*m.^3.0.*h.*j.*(V - obj.E_Na(Nai));
            
            
            obj.mL_inf     =@(V) 1./(1+exp(-(V.*1000+42.85)./(5.264)));
            obj.alpha_m_L   =@(V) 1./(1+exp((-60-V.*1000)./5));
            obj.beta_m_L    =@(V) 0.1./(1+exp((V.*1000+35)./5))+0.1./(1+exp((V.*1000-50)./200));
            obj.tau_mL     =@(V) 1./1000 .* obj.myCoefTauM.*obj.alpha_m_L(V).*obj.beta_m_L(V);
            obj.hL_inf     =@(V) 1./(1+exp((V.*1000+obj.Vh_hLate)./(7.488)));
            obj.tau_hL     = 1./1000 .* obj.tauINaL;
            
            obj.i_NaL       =@(V,mL,hL,Nai) obj.GNaLmax.* mL.^(3).*hL.*(V-obj.E_Na(Nai));
            
            obj.Xf_inf      =@(V) 1.0./(1.0 + exp((V.*1000 + 69)./8));
            obj.tau_Xf      =@(V) (5600 ./ (1 + exp((V.*1000 + 65)./7) + exp(-(V.*1000 + 65)./19)))./1000;
            
            obj.i_fK        =@(V,Xf) obj.fK.*obj.g_f.*Xf.*(V - obj.E_K);
            obj.i_fNa       =@(V,Xf,Nai) obj.fNa.*obj.g_f.*Xf.*(V - obj.E_Na(Nai));
            obj.i_f         =@(V,Xf,Nai) obj.i_fK(V,Xf) + obj.i_fNa(V,Xf,Nai);
            
            obj.d_inf       =@(V) 1.0./(1.0+exp(-(V.*1000.0+9.1)./7.0));
            obj.alpha_d     =@(V) 0.25+1.4./(1.0+exp((-V.*1000.0-35.0)./13.0));
            obj.beta_d      =@(V) 1.4./(1.0+exp((V.*1000.0+5.0)./5.0));
            obj.gamma_d     =@(V) 1.0./(1.0+exp((-V.*1000.0+50.0)./20.0));
            obj.tau_d       =@(V) (obj.alpha_d(V).*obj.beta_d(V)+obj.gamma_d(V)).*1.0./1000.0;
            
            obj.f1_inf      =@(V) 1.0./(1.0+exp((V.*1000.0+26.0)./3.0));
            obj.constf1     =@(V,Cai,f1) (1.0+1433.0.*(Cai-50.0.*1.0e-6)).*(obj.f1_inf(V)-f1 > 0.0) + 1.0.*(obj.f1_inf(V)-f1 < 0.0);
            obj.tau_f1      =@(V,Cai,f1) (20.0+1102.5.*exp(-((V.*1000.0+27.0)./15.0).^2.0)+200.0./(1.0+exp((13.0-V.*1000.0)./10.0))+180.0./(1.0+exp((30.0+V.*1000.0)./10.0))).*obj.constf1(V,Cai,f1)./1000.0;
            
            obj.f2_inf      =@(V) 0.33+0.67./(1.0+exp((V.*1000.0+32.0)./4.0));
            obj.tau_f2      =@(V) (600.0.*exp(-(V.*1000.0+25.0).^2.0./170.0)+31.0./(1.0+exp((25.0-V.*1000.0)./10.0))+16.0./(1.0+exp((30.0+V.*1000.0)./10.0))).*obj.constf2./1000.0;
            
            obj.alpha_fCa   =@(Cai) 1.0./(1.0+(Cai./0.0006).^8.0);
            obj.beta_fCa    =@(Cai) 0.1./(1.0+exp((Cai-0.0009)./0.0001));
            obj.gamma_fCa   =@(Cai) 0.3./(1.0+exp((Cai-0.00075)./0.0008));
            obj.fCa_inf     =@(Cai) (obj.alpha_fCa(Cai)+obj.beta_fCa(Cai)+obj.gamma_fCa(Cai))./1.3156;
            obj.constfCa    =@(V,Cai,fCa) 0.0.*((V > -0.06) & (obj.fCa_inf(Cai) > fCa)) + 1.0.*((V < -0.06) | (obj.fCa_inf(Cai) > fCa));
            obj.tau_fCa     =@(V,Cai,fCa) 0.002./obj.constfCa(V,Cai,fCa);   % second (in i_CaL_fCa_gate)
            
            obj.i_CaL       =@(V,Cai,d,f1,f2,fCa) obj.g_CaL.*4.0.*V.*obj.F.^2.0./(obj.R.*obj.T).*(Cai.*exp(2.0.*V.*obj.F./(obj.R.*obj.T))-0.341.*obj.Cao)./(exp(2.0.*V.*obj.F./(obj.R.*obj.T))-1.0).*d.*f1.*f2.*fCa;
            
            obj.q_inf       =@(V) 1.0./(1.0+exp((V.*1000.0+53.0)./13.0));
            obj.tau_q       =@(V) (6.06+39.102./(0.57.*exp(-0.08.*(V.*1000.0+44.0))+0.065.*exp(0.1.*(V.*1000.0+45.93))))./1000.0;
            obj.r_inf       =@(V) 1.0./(1.0+exp(-(V.*1000.0-22.3)./18.75));
            obj.tau_r       =@(V) (2.75352+14.40516./(1.037.*exp(0.09.*(V.*1000.0+30.61))+0.369.*exp(-0.12.*(V.*1000.0+23.84))))./1000.0;
            
            obj.i_to        =@(V,q,r) obj.g_to.*(V-obj.E_K).*q.*r;
            
            obj.Xs_inf      =@(V) 1.0./(1.0+exp((-V.*1000.0-20.0)./16.0));
            obj.alpha_Xs    =@(V) 1100.0./sqrt(1.0+exp((-10.0-V.*1000.0)./6.0));
            obj.beta_Xs     =@(V) 1.0./(1.0+exp((-60.0+V.*1000.0)./20.0));
            obj.tau_Xs      =@(V) 1.0.*obj.alpha_Xs(V).*obj.beta_Xs(V)./1000.0;
            
            obj.i_Ks        =@(V,Cai,Xs,Nai) obj.g_Ks.*(V-obj.E_Ks(Nai)).*Xs.^2.0.*(1.0+0.6./(1.0+(3.8.*0.00001./Cai).^1.4));
            
            obj.Xr1_inf      =@(V) 1.0./(1.0+exp((obj.V_half-V.*1000.0)./4.9));
            obj.alpha_Xr1    =@(V) 450.0./(1.0+exp((-45.0-V.*1000.0)./10.0));
            obj.beta_Xr1     =@(V) 6.0./(1.0+exp((30.0+V.*1000.0)./11.5));
            obj.tau_Xr1      =@(V) 1.0.*obj.alpha_Xr1(V).*obj.beta_Xr1(V)./1000.0;
            
            obj.Xr2_inf      =@(V) 1.0./(1.0+exp((V.*1000.0+88.0)./50.0));
            obj.alpha_Xr2    =@(V) 3.0./(1.0+exp((-60.0-V.*1000.0)./20.0));
            obj.beta_Xr2     =@(V) 1.12./(1.0+exp((-60.0+V.*1000.0)./20.0));
            obj.tau_Xr2      =@(V) 1.0.*obj.alpha_Xr2(V).*obj.beta_Xr2(V)./1000.0;
            
            obj.i_Kr         =@(V,Xr1,Xr2) obj.g_Kr.*(V-obj.E_K).*Xr1.*Xr2.*sqrt(obj.Ko./5.4);
            
            obj.alpha_K1    =@(V) 3.91./(1.0+exp(0.5942.*(V.*1000.0-obj.E_K.*1000.0-200.0)));
            obj.beta_K1     =@(V) (-1.509.*exp(0.0002.*(V.*1000.0-obj.E_K.*1000.0+100.0))+exp(0.5886.*(V.*1000.0-obj.E_K.*1000.0-10.0)))./(1.0+exp(0.4547.*(V.*1000.0-obj.E_K.*1000.0)));
            obj.XK1_inf     =@(V) obj.alpha_K1(V)./(obj.alpha_K1(V)+obj.beta_K1(V));
            
            obj.i_K1        =@(V) obj.g_K1.*obj.XK1_inf(V).*(V-obj.E_K).*sqrt(obj.Ko./5.4);
            
            obj.i_NaCa      =@(V,Nai,Cai) obj.kNaCa.*(exp(obj.gamma.*V.*obj.F./(obj.R.*obj.T)).*Nai.^3.0.*obj.Cao-exp((obj.gamma-1.0).*V.*obj.F./(obj.R.*obj.T)).*obj.Nao.^3.0.*Cai.*obj.alpha)./((obj.KmNai.^3.0+obj.Nao.^3.0).*(obj.KmCa+obj.Cao).*(1.0+obj.Ksat.*exp((obj.gamma-1.0).*V.*obj.F./(obj.R.*obj.T))));
            
            obj.i_NaK       =@(V,Nai) obj.PNaK.*obj.Ko./(obj.Ko+obj.Km_K).*Nai./(Nai+obj.Km_Na)./(1.0+0.1245.*exp(-0.1.*V.*obj.F./(obj.R.*obj.T))+0.0353.*exp(-V.*obj.F./(obj.R.*obj.T)));
            
            obj.i_PCa       =@(Cai) obj.g_PCa.*Cai./(Cai+obj.KPCa);
            
            obj.i_b_Na      =@(V,Nai) obj.g_b_Na.*(V-obj.E_Na(Nai));
            
            obj.i_b_Ca      =@(V,Cai) obj.g_b_Ca.*(V-obj.E_Ca(Cai));
            
            obj.i_up        =@(Cai) obj.VmaxUp./(1.0+obj.Kup.^2.0./Cai.^2.0);
            
            obj.i_leak      =@(CaSR,Cai) (CaSR-Cai).*obj.V_leak;
            
            obj.RyRainfss   =@(Cai) obj.RyRa1-obj.RyRa2./(1 + exp((1000.0.*Cai-(obj.RyRahalf))./0.0082));
            obj.RyRoinfss   =@(RyRa,Cai) (1 - 1./(1 +  exp((1000.0.*Cai-(RyRa+ obj.RyRohalf))./0.003)));
            obj.RyRtauact   =@(RyRo,Cai,RyRa) (18.75e-3).*(obj.RyRoinfss(RyRa,Cai)>= RyRo) + (0.1.*18.75e-3).*(obj.RyRoinfss(RyRa,Cai)< RyRo);
            
            obj.RyRcinfss   =@(RyRa,Cai) (1./(1 + exp((1000.0.*Cai-(RyRa+obj.RyRchalf))./0.001)));
            obj.RyRtauinact =@(RyRc,Cai,RyRa) (2.*87.5e-3).*(obj.RyRcinfss(RyRa,Cai)>= RyRc) + (87.5e-3).*(obj.RyRcinfss(RyRa,Cai)< RyRc);
            
            obj.RyRSRCass   =@(CaSR) (1 - 1./(1 +  exp((CaSR-0.3)./0.1)));
            obj.i_rel       =@(CaSR,Cai,RyRo,RyRc) obj.g_irel_max.*obj.RyRSRCass(CaSR).*RyRo.*RyRc.*(CaSR-Cai);
            
            obj.Cai_bufc    =@(Cai) 1.0./(1.0+obj.Buf_C.*obj.Kbuf_C./(Cai+obj.Kbuf_C).^2.0);
            obj.Ca_SR_bufSR =@(CaSR) 1.0./(1.0+obj.Buf_SR.*obj.Kbuf_SR./(CaSR+obj.Kbuf_SR).^2.0);
            
            
            
            obj.CaSR_o=zeros(size(V));
            obj.Cai_o=zeros(size(V));
            obj.g_o=zeros(size(V));
            obj.d_o=zeros(size(V));
            obj.f1_o=zeros(size(V));
            obj.f2_o=zeros(size(V));
            obj.fCa_o=zeros(size(V));
            obj.Xr1_o=zeros(size(V));
            obj.Xr2_o=zeros(size(V));
            obj.Xs_o=zeros(size(V));
            obj.h_o=zeros(size(V));
            obj.j_o=zeros(size(V));
            obj.m_o=zeros(size(V));
            obj.Xf_o=zeros(size(V));
            obj.q_o=zeros(size(V));
            obj.r_o=zeros(size(V));
            obj.Nai_o=zeros(size(V));
            obj.mL_o=zeros(size(V));
            obj.hL_o=zeros(size(V));
            obj.RyRa_o=zeros(size(V));
            obj.RyRo_o=zeros(size(V));
            obj.RyRc_o=zeros(size(V));

            obj.CaSR_n=zeros(size(V));
            obj.Cai_n=zeros(size(V));
            obj.g_n=zeros(size(V));
            obj.d_n=zeros(size(V));
            obj.f1_n=zeros(size(V));
            obj.f2_n=zeros(size(V));
            obj.fCa_n=zeros(size(V));
            obj.Xr1_n=zeros(size(V));
            obj.Xr2_n=zeros(size(V));
            obj.Xs_n=zeros(size(V));
            obj.h_n=zeros(size(V));
            obj.j_n=zeros(size(V));
            obj.m_n=zeros(size(V));
            obj.Xf_n=zeros(size(V));
            obj.q_n=zeros(size(V));
            obj.r_n=zeros(size(V));
            obj.Nai_n=zeros(size(V));
            obj.mL_n=zeros(size(V));
            obj.hL_n=zeros(size(V));
            obj.RyRa_n=zeros(size(V));
            obj.RyRo_n=zeros(size(V));
            obj.RyRc_n=zeros(size(V));

            obj.CaSR_n(:,1)=obj.CaSR0;
            obj.Cai_n(:,1)=obj.Cai0;
            obj.g_n(:,1)=obj.g0;
            obj.d_n(:,1)=obj.d0;
            obj.f1_n(:,1)=obj.f10;
            obj.f2_n(:,1)=obj.f20;
            obj.fCa_n(:,1)=obj.fCa0;
            obj.Xr1_n(:,1)=obj.Xr10;
            obj.Xr2_n(:,1)=obj.Xr20;
            obj.Xs_n(:,1)=obj.Xs0;
            obj.h_n(:,1)=obj.h0;
            obj.j_n(:,1)=obj.j0;
            obj.m_n(:,1)=obj.m0;
            obj.Xf_n(:,1)=obj.Xf0;
            obj.q_n(:,1)=obj.q0;
            obj.r_n(:,1)=obj.r0;
            obj.Nai_n(:,1)=obj.Nai0;
            obj.mL_n(:,1)=obj.mL0;
            obj.hL_n(:,1)=obj.hL0;
            obj.RyRa_n(:,1)=obj.RyRa0;
            obj.RyRo_n(:,1)=obj.RyRo0;
            obj.RyRc_n(:,1)=obj.RyRc0;
            
            % dt
            obj.dt=dt;

        end


        function Iion=solveTimestep(obj,Vold,i)
            m_inf_eval = obj.m_inf(Vold);
            tau_m_eval = obj.tau_m(Vold);
            obj.m_n=m_inf_eval+(obj.m_o-m_inf_eval).*exp(-obj.dt./tau_m_eval);
        
            h_inf_eval = obj.h_inf(Vold);
            tau_h_eval = obj.tau_h(Vold);
            obj.h_n=h_inf_eval+(obj.h_o-h_inf_eval).*exp(-obj.dt./tau_h_eval);
        
            j_inf_eval = obj.j_inf(Vold);
            tau_j_eval = obj.tau_j(Vold);
            obj.j_n=j_inf_eval+(obj.j_o-j_inf_eval).*exp(-obj.dt./tau_j_eval);
        
            mL_inf_eval = obj.mL_inf(Vold);
            tau_mL_eval = obj.tau_mL(Vold);
            obj.mL_n=mL_inf_eval+(obj.mL_o-mL_inf_eval).*exp(-obj.dt./tau_mL_eval);
        
            hL_inf_eval = obj.hL_inf(Vold);
            tau_hL_eval = obj.tau_hL;
            obj.hL_n=hL_inf_eval+(obj.hL_o-hL_inf_eval).*exp(-obj.dt./tau_hL_eval);
        
            Xf_inf_eval = obj.Xf_inf(Vold);
            tau_Xf_eval = obj.tau_Xf(Vold);
            obj.Xf_n=Xf_inf_eval+(obj.Xf_o-Xf_inf_eval).*exp(-obj.dt./tau_Xf_eval);
           
            d_inf_eval = obj.d_inf(Vold);
            tau_d_eval = obj.tau_d(Vold);
            obj.d_n=d_inf_eval+(obj.d_o-d_inf_eval).*exp(-obj.dt./tau_d_eval);
        
            f1_inf_eval = obj.f1_inf(Vold);
            tau_f1_eval = obj.tau_f1(Vold,obj.Cai_o,obj.f1_o);
            obj.f1_n=f1_inf_eval+(obj.f1_o-f1_inf_eval).*exp(-obj.dt./tau_f1_eval);
        
            f2_inf_eval = obj.f2_inf(Vold);
            tau_f2_eval = obj.tau_f2(Vold);
            obj.f2_n=f2_inf_eval+(obj.f2_o-f2_inf_eval).*exp(-obj.dt./tau_f2_eval);
        
            fCa_inf_eval = obj.fCa_inf(obj.Cai_o);
            tau_fCa_eval = obj.tau_fCa(Vold,obj.Cai_o,obj.fCa_o);
            obj.fCa_n=fCa_inf_eval+(obj.fCa_o-fCa_inf_eval).*exp(-obj.dt./tau_fCa_eval);
           
            q_inf_eval = obj.q_inf(Vold);
            tau_q_eval = obj.tau_q(Vold);
            obj.q_n=q_inf_eval+(obj.q_o-q_inf_eval).*exp(-obj.dt./tau_q_eval);
        
            r_inf_eval = obj.r_inf(Vold);
            tau_r_eval = obj.tau_r(Vold);
            obj.r_n=r_inf_eval+(obj.r_o-r_inf_eval).*exp(-obj.dt./tau_r_eval);
        
            Xs_inf_eval = obj.Xs_inf(Vold);
            tau_Xs_eval = obj.tau_Xs(Vold);
            obj.Xs_n=Xs_inf_eval+(obj.Xs_o-Xs_inf_eval).*exp(-obj.dt./tau_Xs_eval);
        
            Xr1_inf_eval = obj.Xr1_inf(Vold);
            tau_Xr1_eval = obj.tau_Xr1(Vold);
            obj.Xr1_n=Xr1_inf_eval+(obj.Xr1_o-Xr1_inf_eval).*exp(-obj.dt./tau_Xr1_eval);
        
            Xr2_inf_eval = obj.Xr2_inf(Vold);
            tau_Xr2_eval = obj.tau_Xr2(Vold);
            obj.Xr2_n=Xr2_inf_eval+(obj.Xr2_o-Xr2_inf_eval).*exp(-obj.dt./tau_Xr2_eval);
        
            RyRainfss_eval=obj.RyRainfss(obj.Cai_o);
            dRyRa_eval=(RyRainfss_eval- obj.RyRa_o)./obj.RyRtauadapt;
            obj.RyRa_n= obj.RyRa_o+obj.dt.*dRyRa_eval;
        
            RyRoinfss_eval=obj.RyRoinfss(obj.RyRa_o,obj.Cai_o);
            RyRtauact_eval=obj.RyRtauact(obj.RyRo_o,obj.Cai_o,obj.RyRa_o);
            dRyRo_eval=(RyRoinfss_eval- obj.RyRo_o)./(RyRtauact_eval);
            obj.RyRo_n= obj.RyRo_o+obj.dt.*dRyRo_eval;
        
            RyRcinfss_eval   =obj.RyRcinfss(obj.RyRa_o,obj.Cai_o);
            RyRtauinact_eval =obj.RyRtauinact(obj.RyRc_o,obj.Cai_o,obj.RyRa_o);
            dRyRc_eval = (RyRcinfss_eval- obj.RyRc_o)./(RyRtauinact_eval);
            obj.RyRc_n= obj.RyRc_o+obj.dt.*dRyRc_eval;
            
         
        
        
            i_Na_eval=obj.i_Na(Vold,obj.m_n,obj.h_n,obj.j_n,obj.Nai_o);
            i_NaL_eval=obj.i_NaL(Vold,obj.mL_n,obj.hL_n,obj.Nai_o);
            i_f_eval=obj.i_f(Vold,obj.Xf_n,obj.Nai_o);
            i_CaL_eval=obj.i_CaL(Vold,obj.Cai_o,obj.d_n,obj.f1_n,obj.f2_n,obj.fCa_n);
            i_to_eval=obj.i_to(Vold,obj.q_n,obj.r_n);
            i_Ks_eval=obj.i_Ks(Vold,obj.Cai_o,obj.Xs_n,obj.Nai_o);
            i_Kr_eval=obj.i_Kr(Vold,obj.Xr1_n,obj.Xr2_n);
            i_K1_eval=obj.i_K1(Vold);
            i_NaCa_eval=obj.i_NaCa(Vold,obj.Nai_o,obj.Cai_o);
            i_NaK_eval=obj.i_NaK(Vold,obj.Nai_o);
            i_PCa_eval=obj.i_PCa(obj.Cai_o);
            i_b_Na_eval=obj.i_b_Na(Vold,obj.Nai_o);
            i_b_Ca_eval=obj.i_b_Ca(Vold,obj.Cai_o);
            i_up_eval=obj.i_up(obj.Cai_o);
            i_leak_eval=obj.i_leak(obj.CaSR_o,obj.Cai_o);
            i_fNa_eval=obj.i_fNa(Vold,obj.Xf_n,obj.Nai_o);
            i_rel_eval=obj.i_rel(obj.CaSR_o,obj.Cai_o,obj.RyRo_n,obj.RyRc_n);
        
            dNai_eval   = -obj.Cm.*(i_Na_eval+i_NaL_eval+i_b_Na_eval+3.0.*i_NaK_eval+3.0.*i_NaCa_eval+i_fNa_eval)./(obj.F.*obj.Vc.*1.0e-18);
            obj.Nai_n      = obj.Nai_o+obj.dt.*dNai_eval;
        
            dCai_eval   =@(Cai) obj.Cai_bufc(Cai).*(i_leak_eval-i_up_eval+i_rel_eval-(i_CaL_eval+i_b_Ca_eval+i_PCa_eval-2.0.*i_NaCa_eval).*obj.Cm./(2.0.*obj.Vc.*obj.F.*1.0e-18));
            obj.Cai_n      = obj.Cai_o+obj.dt.*dCai_eval(obj.Cai_o);
        
            dCaSR_eval   =@(CaSR) obj.Ca_SR_bufSR(CaSR).*obj.Vc./obj.V_SR.*(i_up_eval-(i_rel_eval+i_leak_eval));
            obj.CaSR_n      = obj.CaSR_o+obj.dt.*dCaSR_eval(obj.CaSR_o);
            
            
            Iion=i_K1_eval+i_to_eval+i_Kr_eval+i_Ks_eval+i_CaL_eval+i_NaK_eval+i_Na_eval+i_NaL_eval+i_NaCa_eval+i_PCa_eval+i_f_eval+i_b_Na_eval+i_b_Ca_eval;




        end
        
        
        function Iion=getCurr(obj,Vold,i)
            obj.CaSR_o=obj.CaSR_n;
            obj.Cai_o=obj.Cai_n;
            obj.g_o=obj.g_n;
            obj.d_o=obj.d_n;
            obj.f1_o=obj.f1_n;
            obj.f2_o=obj.f2_n;
            obj.fCa_o=obj.fCa_n;
            obj.Xr1_o=obj.Xr1_n;
            obj.Xr2_o=obj.Xr2_n;
            obj.Xs_o=obj.Xs_n;
            obj.h_o=obj.h_n;
            obj.j_o=obj.j_n;
            obj.m_o=obj.m_n;
            obj.Xf_o=obj.Xf_n;
            obj.q_o=obj.q_n;
            obj.r_o=obj.r_n;
            obj.Nai_o=obj.Nai_n;
            obj.mL_o=obj.mL_n;
            obj.hL_o=obj.hL_n;
            obj.RyRa_o=obj.RyRa_n;
            obj.RyRo_o=obj.RyRo_n;
            obj.RyRc_o=obj.RyRc_n;

            Iion=obj.solveTimestep(Vold,i);
            %Iion = Iion*1e3;
        end
    end
end
