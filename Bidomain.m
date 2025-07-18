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
        Vsave
        usave
        %saveStep
        exportStep
        storeStep
        saveCounter

        U_rest

    end
    methods
        function obj=Bidomain(pg,M,L,T,ionicModelType,factorize,U_rest,Di,De,VstoreStep)
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
            obj.storeStep=VstoreStep;
            nSaves = floor(length(T)/obj.storeStep); % + 1;
            obj.Vsave = zeros(obj.pg.nv, nSaves);
            obj.Vsave(:,1) = obj.Vn;
            obj.usave = zeros(obj.pg.nv, nSaves);
            obj.usave(:,1) = obj.un;
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
            obj.exportStep=1;
            
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
                    %obj.Vn=pcg(obj.Mat,rhs,1e-7,1000,obj.H,obj.H',obj.Vo);
                    %obj.Vn=pcg(obj.Mat,rhs,1e-7,1000,[],[],obj.Vo);
                end

                clear rhs

                rhs=-obj.Li*obj.Vn;
                [obj.un,~,~]=pcg(obj.Mat2,rhs,1e-7,1000,[],[],obj.uo);
                %obj.un=pcg(obj.Mat2,rhs,1e-7,1000,[],[],obj.uo);
                clear rhs


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
        function plotAndAnalyze(obj)
            % Plot a single time series for one node (e.g., the center)
            centerNode = round(obj.pg.nv / 2);
            timeVec = obj.T(1:obj.storeStep:end); % Tempo associato a Vsave
        
            if obj.ionicModelType==3 || obj.ionicModelType==4 || obj.ionicModelType==5
                timeVec = timeVec*1e3;
                Vsaved = obj.Vsave*1e3;
                usaved = obj.usave*1e3;
            end

            figure;
            plot(timeVec, Vsaved(centerNode,:), 'LineWidth', 2);
            xlabel('Time [ms]');
            ylabel('V [mV]');
            title('Transmembrane Potential at Center Node');

            figure;
            plot(timeVec, usaved(centerNode,:), 'LineWidth', 2);
            xlabel('Time [ms]');
            ylabel('ue [mV]');
            title('Extracellular Potential at Center Node');
        
            % Activation & Repolarization Time Maps
            nNodes = size(obj.Vsave,1);
            activationTimes = zeros(nNodes,1);
            repolarizationTimes = zeros(nNodes,1);
        
            for j = 1:nNodes
                Vj = Vsaved(j,:);
                [t_act, t_repo] = computeActiRepo(timeVec, Vj);                
                activationTimes(j) = t_act;
                repolarizationTimes(j) = t_repo;
            end

            activationTimes_s = activationTimes * 1e-3;
        
            nv = obj.pg.nlv;
            h = obj.pg.h;
            
            
            % Open VTK file
            fid = fopen('AT_RT.vtk', 'w');
            fprintf(fid, '# vtk DataFile Version 3.0\n');
            fprintf(fid, 'Example structured points dataset\n');
            fprintf(fid, 'ASCII\n\n');
            fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
            
            % Adjust for 1D or 2D geometries
            if nv(2) ~= 1 && nv(3) == 1
                nv(3) = 2;
                activationTimes = repmat(activationTimes, 1, 2);
                repolarizationTimes = repmat(repolarizationTimes, 1, 2);
            elseif nv(2) == 1 && nv(3) == 1
                nv(2:3) = 2;
                activationTimes = repmat(activationTimes, 1, 4);
                repolarizationTimes = repmat(repolarizationTimes, 1, 4);
            end
            
            fprintf(fid, 'DIMENSIONS %d %d %d\n', nv(1), nv(2), nv(3));
            fprintf(fid, 'ORIGIN 0 0 0\n');
            fprintf(fid, 'SPACING %f %f %f\n\n', h(1), h(2), h(3));
            
            npoints = prod(nv);
            fprintf(fid, 'POINT_DATA %d\n', npoints);
            
            % Ativation Times
            fprintf(fid, 'SCALARS AT float 1\n');
            fprintf(fid, 'LOOKUP_TABLE default\n');
            fprintf(fid, '%f\n', activationTimes(:));

            % Repolarization Times
            fprintf(fid, 'SCALARS RT float 1\n');
            fprintf(fid, 'LOOKUP_TABLE default\n');
            fprintf(fid, '%f\n', repolarizationTimes(:));

            % Plot activation and repolarization maps (if 2D or 3D)
            % [X,Y,Z] = obj.pg.getCoo;
            % 
            % sz = size(X); 
            % AtMap = reshape(activationTimes, sz);
            % RtMap = reshape(repolarizationTimes, sz);

            % figure;
            % surf(X, Y, Z, AtMap); %,'EdgeColor','none');
            % colorbar;
            % title('Activation Time Map');
            % xlabel('x'); ylabel('y'); zlabel('z');
            % 
            % figure;
            % surf(X, Y, Z, RtMap); %,'EdgeColor','none');
            % colorbar;
            % title('Repolarization Time Map');
            % xlabel('x'); ylabel('y'); zlabel('z');

            % Metodo 3: Average conduction velocity
            E1 = [0.3 0.3];
            E2 = [0.8 0.8];
            [X,Y,Z]=obj.pg.getCoo;
            dist_grid_1 = sqrt((X - E1(1)).^2 + (Y - E1(2)).^2);
            dist_grid_2 = sqrt((X - E2(1)).^2 + (Y - E2(2)).^2);
            [~, which_1] = min(dist_grid_1(:));
            [~, which_2] = min(dist_grid_2(:));
            ds = sqrt((E1(1) - E2(1)).^2 + (E1(2) - E2(2)).^2);
            dt = activationTimes_s(which_2) - activationTimes_s(which_1);
            cv_average = ds / dt;
            disp(['Average velocity considering V = ' num2str(cv_average)])


        end

    end
end



function [t_act, t_repo] = computeActiRepo(timeVec, Vtrace)
    % Trova il massimo e minimo per threshold
    Vmax = max(Vtrace);
    Vmin = min(Vtrace);
    Vthresh = -45;%Vmin + 0.9*(Vmax - Vmin);  % 90% di attivazione
    Vrepo = -50;%Vmin + 0.1*(Vmax - Vmin);    % 10% di ripolarizzazione

    % Trova indice di attivazione e ripolarizzazione
    idx_act = find(Vtrace > Vthresh, 1, 'first');
    idx_repo = find(Vtrace(idx_act:end) < Vrepo, 1, 'first') + idx_act - 1;

    if isempty(idx_act)
        t_act = NaN;
    else
        t_act = timeVec(idx_act);
    end
    if isempty(idx_repo)
        t_repo = NaN;
    else
        t_repo = timeVec(idx_repo);
    end
end

