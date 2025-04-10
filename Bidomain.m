classdef Bidomain < handle
    properties
        pg

        Di
        De
        Li
        Lie
        
        M
        L
        Mat
        Mat2
        H
        
        T
        dt

        V
        u
        Vn
        un
        Vo
        uo
        
        factorize
        ionicModelType

        ionicModel

        %diff
        Iapp
        IappStartTime
        IappStopTime
        %Vsave
        %saveStep
        exportStep

        U_rest

    end
    methods
        function obj=Bidomain(pg,M,L,T,ionicModelType,factorize,U_rest,Di,De)
            obj.pg=pg;
            
            obj.factorize=factorize;

            obj.U_rest=U_rest;

            obj.T=T;
            obj.dt=T(2)-T(1);
            %obj.diff=1e-3;

            obj.Di=Di;
            obj.De=De;

            obj.M=M;
            obj.L=L;
            obj.Li=obj.Di*obj.L;
            obj.Lie=(obj.Di+obj.De)*obj.L;

            obj.Mat=M+obj.dt*obj.Li;
            %if factorize
            obj.H=chol(obj.Mat);
            %end
            obj.Mat2=obj.Lie;
            
            obj.Vn=zeros(obj.pg.nv,1);
            obj.Vo=zeros(obj.pg.nv,1);
            obj.un=zeros(obj.pg.nv,1);
            obj.uo=zeros(obj.pg.nv,1);
            obj.Vn(:,1)=U_rest;

            %obj.V=zeros(size(L,1),length(T));
            %obj.V(:,1)=obj.Vn;

            obj.ionicModelType=ionicModelType;
            if ionicModelType==1
                obj.ionicModel=HodgkinHuxley(obj.Vn,obj.dt);
            end
            if ionicModelType==2
                obj.ionicModel=TenTusscher(obj.V,obj.dt);
            end
            if ionicModelType==3
                obj.ionicModel=Paci(obj.V,obj.dt);
            end
            obj.exportStep=1;
            
            if obj.ionicModelType==3
                exportVTK(obj.Vn*1e3,obj.un,obj.pg,0,1);
            else
                exportVTK(obj.Vn,obj.un,obj.pg,0,1);
            end
        end
        function run(obj)
            nameCounter=0;
            for i=2:length(obj.T)

                t=obj.T(i);
                
                obj.Vo=obj.Vn;
                obj.uo=obj.un;

                if mod(i,100)==0
                    disp(['t=', num2str(t),'Vmax=',num2str(max(obj.Vo))])
                end
                
                Iion=obj.ionicModel.getCurr(obj.Vo,i);
                if obj.IappStartTime<t && t<obj.IappStopTime
                    Iapp2=obj.Iapp;
                else
                    Iapp2=0*obj.Iapp;
                end
                
                Itot=Iapp2-Iion;

                rhs=obj.Vo+obj.dt*Itot;
                rhs=obj.M*rhs;
                rhs=rhs-obj.dt*obj.Li*obj.uo;

                if obj.factorize
                    y=obj.H'\rhs;
                    obj.Vn=obj.H\y;
                else
                    [obj.Vn,~,~]=pcg(obj.Mat,rhs,1e-7,1000,obj.H,obj.H',obj.Vo);
                end

                clear rhs

                rhs=-obj.Li*obj.Vn;
                [obj.un,~,~]=pcg(obj.Mat2,rhs,1e-7,1000,[],[],obj.uo);

                clear rhs


                % PROTOCOL FOR SPIRAL WAVES
                %if (i==3334)
                %    [X,Y,Z]=obj.pg.getCoo;
                %    which=find(Y<0.5 | X<0.5);
                %    Vn(which)=obj.U_rest;
                %end

                %obj.V(:,i)=Vn;

                if mod(i,obj.exportStep)==0
                    nameCounter=nameCounter+1;
                    if obj.ionicModelType==3
                        exportVTK(obj.Vn*1e3,obj.un,obj.pg,nameCounter,1);
                    else
                        exportVTK(obj.Vn,obj.un,obj.pg,nameCounter,1);
                    end
                    
                end
            end
        end
    end
end
