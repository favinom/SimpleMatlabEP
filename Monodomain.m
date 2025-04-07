classdef Monodomain < handle
    properties
        pg

        M
        L
        Mat
        H
        
        T
        dt
        
        V
        Vn
        Vo
        
        factorize
        ionicModelType

        ionicModel

        diff
        Iapp
        IappStartTime
        IappStopTime
        %Vsave
        %saveStep
        exportStep

        U_rest

    end
    methods
        function obj=Monodomain(pg,M,L,T,ionicModelType,factorize,U_rest)
            obj.pg=pg;

            obj.U_rest=U_rest;

            obj.T=T;
            obj.dt=T(2)-T(1);
            obj.diff=1e-3;

            obj.M=M;
            obj.L=L;
            obj.Mat=M+obj.dt*1e-3*L;
            if factorize
                obj.H=chol(obj.Mat);
            end
            
            obj.V=zeros(size(L,1),length(T));
            obj.V(:,1)=U_rest;

            obj.ionicModelType=ionicModelType;
            if ionicModelType==1
                obj.ionicModel=HodgkinHuxley(obj.V,obj.dt);
            end
            obj.exportStep=1;
        end
        function run(obj)
            nameCounter=0;
            for i=2:length(obj.T)
                t=obj.T(i);
                if mod(i,100)==0
                    disp(['t=', num2str(t),'Vmax=',num2str(max(obj.V(:,i-1)))])
                end
                Vold=obj.V(:,i-1);
                Iion=obj.ionicModel.getCurr(Vold,i);
                if obj.IappStartTime<t && t<obj.IappStopTime
                    Iapp2=280*obj.Iapp;
                else
                    Iapp2=0*obj.Iapp;
                end
                
                Itot=Iapp2-Iion;

                rhs=Vold+obj.dt*Itot;
                rhs=obj.M*rhs;

                if obj.factorize
                    y=obj.H'\rhs;
                    Vn=obj.H\y;
                else
                    [Vn,~,~]=pcg(obj.Mat,rhs,[],[],[],[],Vold);
                end

                % PROTOCOL FOR SPIRAL WAVES
                if (i==3334)
                    [X,Y,Z]=obj.pg.getCoo;
                    which=find(Y<0.5 | X<0.5);
                    Vn(which)=obj.U_rest;
                end

                obj.V(:,i)=Vn;

                if mod(i,obj.exportStep)==0
                    nameCounter=nameCounter+1;
                    exportVTK(Vn,obj.pg,nameCounter);
                end
            end
        end
    end
end
