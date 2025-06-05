classdef MEA  < handle
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
        H2
        
        T
        dt

        V
        u
        U
        I
        Vn
        un
        Un
        In
        Vo
        uo
        Uo
        Io

        tau
        Ri
        Cel
        thick
        
        ionicModelType

        ionicModel

        %diff
        Iapp
        IappStartTime
        IappStopTime
        Vsave
        usave
        Usave
        FPsave
        %saveStep
        exportStep
        storeStep
        saveCounter

        fk
        boundary_nodes

        U_rest

    end
    methods
        function obj=MEA(pg,M,L,T,ionicModelType,U_rest,Di,De,VstoreStep,fk,boundary_nodes)
            obj.pg=pg;
            
            obj.boundary_nodes = boundary_nodes;

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
            obj.Mat2(obj.boundary_nodes, :) = 0;
            obj.Mat2(:, obj.boundary_nodes) = 0;
            obj.Mat2(obj.boundary_nodes, obj.boundary_nodes) = speye(length(obj.boundary_nodes));
            obj.H2=chol(obj.Mat2);
            %[obj.P,obj.Q]=ilu(obj.Mat2);
            
            obj.Vn=zeros(obj.pg.nv,1);
            obj.Vo=zeros(obj.pg.nv,1);
            obj.un=zeros(obj.pg.nv,1);
            obj.uo=zeros(obj.pg.nv,1);
            obj.un=zeros(obj.pg.nv,1);
            obj.fk = fk;
            for k = 1:length(obj.fk)
                obj.Uo{k}=obj.uo'*fk{k};
                obj.Un{k}=obj.un'*fk{k};
                obj.Io{k}=0;%obj.uo'*fk{k};
                obj.In{k}=0;%obj.un'*fk{k};
            end
            obj.Vn(:,1)=U_rest;

            obj.Cel = 1e-10;
            obj.Ri = 1e9;
            Re = 1e6;
            obj.tau = (obj.Ri + Re)*obj.Cel;
            obj.thick = 1e-4;

            obj.storeStep=VstoreStep;
            nSaves = floor(length(T)/obj.storeStep); % + 1;
            obj.Vsave = zeros(obj.pg.nv, nSaves);
            obj.Vsave(:,1) = obj.Vn;
            obj.usave = zeros(obj.pg.nv, nSaves);
            obj.usave(:,1) = obj.un;
            for k = 1:length(obj.fk)
                obj.Usave{k}=0; %obj.uo.*fk{k};
                obj.FPsave{k}=0;
            end
            obj.saveCounter = 1;

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
            if ionicModelType==4
                obj.ionicModel=Botti(obj.V,obj.dt);
            end
            if ionicModelType==5
                obj.ionicModel=Amin(obj.V,obj.dt);
            end
            %obj.exportStep=1;
            
            if obj.ionicModelType==3 || obj.ionicModelType==4 || obj.ionicModelType==5
                exportVTK(obj.Vn*1e3,obj.un*1e3,obj.pg,0,1);
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
                obj.Uo=obj.Un;
                obj.Io=obj.In;

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

                %if obj.factorize
                %    y=obj.H'\rhs;
                %    obj.Vn=obj.H\y;
                %else
                [obj.Vn,~,~]=pcg(obj.Mat,rhs,1e-7,1000,obj.H,obj.H',obj.Vo);
                %obj.Vn=pcg(obj.Mat,rhs,1e-7,1000,obj.H,obj.H',obj.Vo);
                %obj.Vn=pcg(obj.Mat,rhs,1e-7,1000,[],[],obj.Vo);
                %end

                rhs_temp = zeros(obj.pg.nv,1);
                for k = 1:length(obj.fk)
                    rhs_temp = rhs_temp + obj.fk{k}*obj.Io{k}/sum(obj.fk{k})/obj.thick;
                end

                clear rhs

                rhs=-obj.Li*obj.Vn + rhs_temp;
                rhs(obj.boundary_nodes)=0;
                %[obj.un,~,~]=pcg(obj.Mat2,rhs,1e-7,1000,[],[],obj.uo);
                % obj.un=pcg(obj.Mat2,rhs,1e-7,1000,[],[],obj.uo);
                %=pcg(obj.Mat2,rhs,1e-7,1000,obj.H2,obj.H2',obj.uo);
                [obj.un,~,~]=pcg(obj.Mat2,rhs,1e-7,1000,obj.H2,obj.H2',obj.uo);
                clear rhs rhs_temp

                for k = 1:length(obj.fk)
                    obj.Un{k} = obj.uo'*obj.fk{k}/sum(obj.fk{k});
                    dudt=(obj.Un{k}-obj.Uo{k})/obj.dt;
                    rhs=obj.Io{k}./ obj.dt + obj.Cel ./ (obj.tau) .* dudt;
                    den=(1./obj.dt + 1./obj.tau);
                    obj.In{k} = rhs/den;
                end

                % PROTOCOL FOR SPIRAL WAVES
                %if (i==3334)
                %    [X,Y,Z]=obj.pg.getCoo;
                %    which=find(Y<0.5 | X<0.5);
                %    Vn(which)=obj.U_rest;
                %end

                %obj.V(:,i)=Vn;

                if mod(i,obj.storeStep) == 0
                    obj.saveCounter = obj.saveCounter + 1;
                    obj.Vsave(:, obj.saveCounter) = obj.Vn;
                    obj.usave(:, obj.saveCounter) = obj.un;
                    for k = 1:length(obj.fk) 
                        obj.Usave{k}(obj.saveCounter) = obj.Un{k};
                        obj.FPsave{k}(obj.saveCounter) = obj.In{k}.*obj.Ri;
                    end
                end

                if mod(i,obj.exportStep)==0
                    nameCounter=nameCounter+1;
                    if obj.ionicModelType==3 || obj.ionicModelType==4 || obj.ionicModelType==5
                        exportVTK(obj.Vn*1e3,obj.un*1e3,obj.pg,nameCounter,1);
                    else
                        exportVTK(obj.Vn,obj.un,obj.pg,nameCounter,1);
                    end
                    
                end
            end
        end
    end
end



