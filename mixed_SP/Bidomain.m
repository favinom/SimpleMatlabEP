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
        % ionicModelType

        % ionicModel

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

        model0
        model1
        SCELTA % vettore 0/1 per selezionare modello per nodo

        nsd
        iter
        iter_dd
        iter_dd1
        iter_dd2
        h_piccolo
        H_grande
        filename

    end
    methods
        function obj=Bidomain(pg,M,L,T,factorize,U_rest,Di,De,VstoreStep,scelta,nsd)
            obj.pg=pg;
            obj.factorize=factorize;
            obj.nsd=nsd;

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
            obj.Vn=U_rest;
            %obj.Vn(:,1)=U_rest;

            %obj.V=zeros(size(L,1),length(T));
            %obj.V(:,1)=obj.Vn;
            obj.storeStep=VstoreStep;
            nSaves = floor(length(T)/obj.storeStep); % + 1;
            obj.Vsave = zeros(obj.pg.nv, nSaves);
            obj.Vsave(:,1) = obj.Vn;
            obj.usave = zeros(obj.pg.nv, nSaves);
            obj.usave(:,1) = obj.un;
            obj.saveCounter = 1;

            obj.model0 = Botti(obj.Vn, obj.dt);
            obj.model1 = Paci(obj.Vn, obj.dt);
            obj.SCELTA = scelta;
            % 
            % obj.ionicModelType=ionicModelType;
            % if ionicModelType==1
            %     obj.ionicModel=HodgkinHuxley(obj.Vn,obj.dt);
            % end
            % if ionicModelType==2
            %     obj.ionicModel=TenTusscher(obj.V,obj.dt);
            % end
            % if ionicModelType==3
            %     obj.ionicModel=Paci(obj.V,obj.dt);
            % end
            % if ionicModelType==4
            %     obj.ionicModel=Botti(obj.V,obj.dt);
            % end
            % obj.exportStep=1;
            % 
            % if obj.ionicModelType==3 || obj.ionicModelType==4
            %     exportVTK(obj.Vn*1e3,obj.un*1e3,obj.pg,0,1);
            % else
            %     exportVTK(obj.Vn,obj.un,obj.pg,0,1);
            % end
            
            obj.exportStep=1;
            
            exportVTK(obj.Vn*1e3,[],obj.pg,0,1);

        end
        function run(obj)
            nameCounter=0;
            for i=2:length(obj.T)

                [id_S,id]=buildSubdomainIds(obj.pg,obj.nsd,obj.nsd);

                t=obj.T(i);
                
                obj.Vo=obj.Vn;
                obj.uo=obj.un;

                if mod(i,100)==0
                    disp(['t=', num2str(t),'Vmax=',num2str(max(obj.Vo))])
                end
                
                % Calcolo Iion da entrambi i modelli
                Iion0 = obj.model0.getCurr(obj.Vo, i);
                Iion1 = obj.model1.getCurr(obj.Vo, i);

                % Combinazione lineare basata su SCELTA
                Iion = obj.SCELTA .* Iion1 + (1 - obj.SCELTA) .* Iion0;

                %Iion=obj.ionicModel.getCurr(obj.Vo,i);
                if obj.IappStartTime<t && t<obj.IappStopTime
                    Iapp2=obj.Iapp;
                else
                    Iapp2=0*obj.Iapp;
                end
                
                Itot=Iapp2-Iion;
                
                % Iion=obj.ionicModel.getCurr(obj.Vo,i);
                % if obj.IappStartTime<t && t<obj.IappStopTime
                %     Iapp2=obj.Iapp;
                % else
                %     Iapp2=0*obj.Iapp;
                % end
                % 
                % Itot=Iapp2-Iion;

                if obj.nsd == 1
                    rhs=obj.Vo+obj.dt*Itot;
                    rhs=obj.M*rhs;
    
                    if obj.factorize
                        y=obj.H'\rhs;
                        obj.Vn=obj.H\y;
                    else
                        [obj.Vn,~,~,obj.iter]=pcg(obj.Mat,rhs,[],[],[],[],obj.Vo);
                    end

                    clear rhs

                    rhs=-obj.Li*obj.Vn;
                    [obj.un,~,~]=pcg(obj.Mat2,rhs,1e-7,1000,[],[],obj.uo);
    
                    clear rhs
                else
                    rhs=obj.Vo+obj.dt*Itot;
                    rhs=obj.M*rhs;

                    sp = SteklovPoincare(obj.Mat, rhs, obj.Vo);
                    sp.prepare(id_S, id);
                    sp.prepareRhs(id);
                    sp.prepareApp();
                
                    [~,~,~,obj.iter_dd1] = sp.solve(1e-7, 10000);
    
                    obj.Vn(id_S) = sp.xs;
                
                    for ids = 1:length(id)
                        xi = sp.H{ids} \ (sp.H{ids}' \ (sp.bi{ids} - sp.Ais{ids} * sp.xs));
                        obj.Vn(id{ids}) = xi;
                    end

                    clear rhs

                    rhs=-obj.Li*obj.Vn;
                    
                    sp = SteklovPoincare(obj.Mat2, rhs, obj.uo);
                    sp.prepare(id_S, id);
                    sp.prepareRhs(id);
                    sp.prepareApp();
                
                    [~,~,~,obj.iter_dd2] = sp.solve(1e-7, 10000);
    
                    obj.un(id_S) = sp.xs;
                
                    for ids = 1:length(id)
                        xi = sp.H{ids} \ (sp.H{ids}' \ (sp.bi{ids} - sp.Ais{ids} * sp.xs));
                        obj.un(id{ids}) = xi;
                    end
                          
                    clear rhs
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
                end

                if mod(i,obj.exportStep)==0
                    nameCounter=nameCounter+1;
                    exportVTK(obj.Vn*1e3,obj.un*1e3,obj.pg,nameCounter,1);
                end

                obj.iter_dd = obj.iter_dd1+obj.iter_dd2;

                data_row = [i, obj.iter_dd1, obj.iter_dd2, obj.iter_dd];
                dlmwrite(obj.filename, data_row, '-append');
            end
        end
        function plotAndAnalyze(obj)
            % Plot a single time series for one node (e.g., the center)
            centerNode = round(obj.pg.nv / 2);
            timeVec = obj.T(1:obj.storeStep:end); % Tempo associato a Vsave
        
            timeVec = timeVec*1e3;
            Vsaved = obj.Vsave*1e3;
            usaved = obj.usave*1e3;

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
            % activationTimes = zeros(nNodes,1);
            % repolarizationTimes = zeros(nNodes,1);
            % 
            % for j = 1:nNodes
            %     Vj = Vsaved(j,:);
            %     [t_act, t_repo] = computeActiRepo(timeVec, Vj);                
            %     activationTimes(j) = t_act;
            %     repolarizationTimes(j) = t_repo;
            % end
 


            maxBeats = 1; %3; % max numero di AT/RT che vuoi rilevare
            activationTimes = cell(maxBeats, 1);
            repolarizationTimes = cell(maxBeats, 1);
            %keyboard
        
            for beat = 1:maxBeats
                at = zeros(nNodes,1);
                rt = zeros(nNodes,1);
                for j = 1:nNodes
                    if beat == 1
                        Tj{j} =obj.T(1:obj.storeStep:end)*1e3;
                    end
                    Vjj{j} = Vsaved(j,:);
                    [at(j), rt(j), timeVec_new, Vj_new] = computeActiRepo(Tj{j}, Vjj{j});
                    Tj{j}=timeVec_new;
                    Vjj{j}=Vj_new;
        
                    % if length(t_at_list) >= beat && length(t_rt_list) >= beat
                    %     at(j) = t_at_list(beat);
                    %     rt(j) = t_rt_list(beat);
                    % end
                end
                activationTimes{beat} = at;
                repolarizationTimes{beat} = rt;
                clear at rt
            end
        
            nv = obj.pg.nlv;
            h = obj.pg.h;
            

            labels=obj.SCELTA(:);
            
            % Open VTK file
            fid = fopen('AT_RT_labels.vtk', 'w');
            fprintf(fid, '# vtk DataFile Version 3.0\n');
            fprintf(fid, 'Example structured points dataset\n');
            fprintf(fid, 'ASCII\n\n');
            fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
            
            % Adjust for 1D or 2D geometries
            if nv(2) ~= 1 && nv(3) == 1
                nv(3) = 2;
                activationTimes{1} = repmat(activationTimes{1}, 1, 2);
                repolarizationTimes{1} = repmat(repolarizationTimes{1}, 1, 2);
                % activationTimes{2} = repmat(activationTimes{2}, 1, 2);
                % repolarizationTimes{2} = repmat(repolarizationTimes{2}, 1, 2);
                % activationTimes{3} = repmat(activationTimes{3}, 1, 2);
                % repolarizationTimes{3} = repmat(repolarizationTimes{3}, 1, 2);
                labels = repmat(labels, 1, 2);
            elseif nv(2) == 1 && nv(3) == 1
                nv(2:3) = 2;
                activationTimes = repmat(activationTimes, 1, 4);
                repolarizationTimes = repmat(repolarizationTimes, 1, 4);
                labels = repmat(labels, 1, 4);
            end
            
            fprintf(fid, 'DIMENSIONS %d %d %d\n', nv(1), nv(2), nv(3));
            fprintf(fid, 'ORIGIN 0 0 0\n');
            fprintf(fid, 'SPACING %f %f %f\n\n', h(1), h(2), h(3));
            
            npoints = prod(nv);
            fprintf(fid, 'POINT_DATA %d\n', npoints);

            for k = 1:maxBeats
                fprintf(fid, 'SCALARS AT_%d float 1\n', k);
                fprintf(fid, 'LOOKUP_TABLE default\n');
                %fprintf(fid, '%f\n', cellfun(@(x) x(k), activationTimes));
                fprintf(fid, '%f\n', activationTimes{k}(:));
                
                fprintf(fid, 'SCALARS RT_%d float 1\n', k);
                fprintf(fid, 'LOOKUP_TABLE default\n');
                fprintf(fid, '%f\n', repolarizationTimes{k}(:));
            end
            
            % Ativation Times
            % fprintf(fid, 'SCALARS AT float 1\n');
            % fprintf(fid, 'LOOKUP_TABLE default\n');
            % fprintf(fid, '%f\n', activationTimes(:));
            % 
            % % Repolarization Times
            % fprintf(fid, 'SCALARS RT float 1\n');
            % fprintf(fid, 'LOOKUP_TABLE default\n');
            % fprintf(fid, '%f\n', repolarizationTimes(:));

            % Labels
            fprintf(fid, 'SCALARS Labels float 1\n');
            fprintf(fid, 'LOOKUP_TABLE default\n');
            fprintf(fid, '%f\n', labels(:));

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
        end

    end
end



function [t_act, t_repo, timeVec_new, Vtrace_new] = computeActiRepo(timeVec, Vtrace)
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
        timeVec_new = timeVec(idx_repo:end);
        Vtrace_new = Vtrace(idx_repo:end);
    end
end