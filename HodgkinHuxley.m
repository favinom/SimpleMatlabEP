classdef HodgkinHuxley < handle
    properties
        % material properties
        g_Na
        g_K
        g_L
        U_Na
        U_K
        U_L
        % initial conditions
        m0
        h0
        n0
        % hh functions
        am
        bm
        ah
        bh
        an
        bn
        % rhs
        rhs1
        rhs2
        rhs3
        % currents
        I_Na
        I_K
        I_L
        % variables()
        M
        H
        N
        %
        dt

    end
    methods
        function obj=HodgkinHuxley(V,dt)
            obj.g_Na = 120.0;
            obj.g_K = 36.0;
            obj.g_L = 0.3;
            obj.U_Na = 120;
            obj.U_K = -77.0;
            obj.U_L = -54.387;

            obj.m0 = 0.0000;
            obj.h0 = 0.9998;
            obj.n0 = 0.0040;

            obj.am =@(V)( 0.1 * (25 - V) ./ (exp((25 - V) / 10) - 1) );
            obj.bm =@(V)( 4.0 * exp(-V / 18) );
            obj.ah =@(V)( 0.07 * exp(-V / 20) );
            obj.bh =@(V)( 1 ./ (exp((30 - V) / 10) + 1) );
            obj.an =@(V)( 0.01 * (10 - V) ./ (exp((10 - V) / 10) - 1) );
            obj.bn =@(V)( 0.125 * exp(-V / 80) );

            obj.rhs1=@(m,V)(obj.am(V).*(1-m)-obj.bm(V).*m);
            obj.rhs2=@(h,V)(obj.ah(V).*(1-h)-obj.bh(V).*h);
            obj.rhs3=@(n,V)(obj.an(V).*(1-n)-obj.bn(V).*n);

            obj.I_Na=@(V,m,h)( obj.g_Na * m.^3 .* h .* (V - obj.U_Na));
            obj.I_K =@(V,n)  ( obj.g_K * n.^4 .* (V - obj.U_K)    );
            obj.I_L =@(V)    ( obj.g_L * (V - obj.U_L)          );

            obj.M=zeros(size(V));
            obj.H=zeros(size(V));
            obj.N=zeros(size(V));

            % supponiamo sia 2D
            obj.M(:,1)=obj.m0;
            obj.H(:,1)=obj.h0;
            obj.N(:,1)=obj.n0;
            % dt
            obj.dt=dt;

        end
        function solveTimestep(obj,Vold,i)
            mold=obj.M(:,i-1);
            hold=obj.H(:,i-1);
            nold=obj.N(:,i-1);
            rhs1eval=obj.rhs1(mold,Vold);
            rhs2eval=obj.rhs2(hold,Vold);
            rhs3eval=obj.rhs3(nold,Vold);
            obj.M(:,i)=mold+obj.dt*rhs1eval;
            obj.H(:,i)=hold+obj.dt*rhs2eval;
            obj.N(:,i)=nold+obj.dt*rhs3eval;
        end
        function Iion=evaluateCurr(obj,Vold,i)
                I_Na_eval=obj.I_Na(Vold,obj.M(:,i),obj.H(:,i));
                I_K_eval=obj.I_K(Vold,obj.N(:,i));
                I_L_eval=obj.I_L(Vold);
                Iion=I_Na_eval+I_K_eval+I_L_eval;
        end
        function Iion=getCurr(obj,Vold,i)
            obj.solveTimestep(Vold,i);
            Iion=obj.evaluateCurr(Vold,i);
        end
    end
end
