function plotAndAnalyze(obj,el)
    [X,Y,~]=obj.pg.getCoo;
    number_points = 9;
    P = cell(1, number_points);
    
    % center = [0.1; 0.1];
    d = obj.pg.nlv(1) * obj.pg.h(1) / 4;
    center = [obj.pg.nlv(1) * obj.pg.h(1) /2; obj.pg.nlv(1) * obj.pg.h(1) /2];
    offsets = [-1 0 1] * d;
        
    k = 1;
    for i = 1:sqrt(number_points)  % righe (y)
        for j = 1:sqrt(number_points)  % colonne (x)
            dx = offsets(j);
            dy = offsets(i);  
            P{k} = center + [dx; dy];
            k = k + 1;
        end
    end
    for i = 1:number_points
        dist_point{i} = sqrt((X - P{i}(1)).^2 + (Y - P{i}(2)).^2);
        [~, selectedNodes{i}] = min(dist_point{i}(:));
    end

    timeVec = obj.T(1:obj.storeStep:end); % Tempo associato a Vsave
    Usaved = cell(size(obj.Usave));
    
    if obj.ionicModelType==3 || obj.ionicModelType==4 || obj.ionicModelType==5
        timeVec = timeVec*1e3;
        Vsaved = obj.Vsave*1e3;
        usaved = obj.usave*1e3;
        for k = 1:length(obj.fk)
            %Usaved{k} = obj.Usave{k}*1e3;
            Usaved{k} = obj.FPsave{k}*1e3;
        end
    end

    %[OS, OSpos] = findpeaks(Vsaved(selectedNodes{1}), 'MinPeakHeight', 0);
    %[~, locs] = findpeaks(Vsaved(selectedNodes{5},:), timeVec, 'MinPeakHeight', threshold, 'MinPeakDistance', 0.5);
    %n_beats = length(locs);

    
    % Inserisci manualmente il numero di battiti
    n_beats = 2; %input('Inserisci il numero di battiti: ');


    % % Trova i picchi (fase di depolarizzazione)
    % [~, locs] = findpeaks(V, t, 'MinPeakHeight', threshold, 'MinPeakDistance', d);
    % 
    % % Conta i battiti
    % n_beats = length(locs);


%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% 
% PLOT IN MATLAB AND CSV SINGLE NODE
%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% 
    plot_in_matlab_0d(obj,timeVec,Vsaved,usaved,Usaved,selectedNodes);
    store_csv = 10;

    time_csv = timeVec(1:store_csv:end);
    %V_csv{1} = Vsaved(selectedNodes{1}, 1:store_csv:end); % Modificato per usare selectedNodes
    %u_csv{1} = usaved(selectedNodes{1}, 1:store_csv:end); % Modificato per usare selectedNodes
    
    T = table(time_csv(:), 'VariableNames', {'Time_ms'});
    %T = table(time_csv(:), V_csv{1}(:), u_csv{1}(:), 'VariableNames', {'Time_ms', 'Voltage_mV', 'ue_mV'});

    for k = 1:length(selectedNodes)
        V_csv = Vsaved(selectedNodes{k},1:store_csv:end);
        varName = sprintf('V%d_mV', k);
        T.(varName) = V_csv(:);
    end

    for k = 1:length(selectedNodes)
        u_csv = usaved(selectedNodes{k},1:store_csv:end);
        varName = sprintf('u%d_mV', k);
        T.(varName) = u_csv(:);
    end

    for k = 1:length(obj.fk)
        U_csv = Usaved{k}(1:store_csv:end);
        varName = sprintf('FP%d_mV', k);
        T.(varName) = U_csv(:);
    end
    
    % Scrive la tabella in un file CSV
    filename = './outputs/mea_traces.csv';
    writetable(T, filename);
    
    disp(['File salvato: ', filename]);
    

%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% 
% Activation & Repolarization Time Maps
%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% 
    nNodes = size(obj.Vsave,1);
    %activationTimes = zeros(nNodes,1);
    %Latency = zeros(nNodes,1);
    %repolarizationTimes = zeros(nNodes,1);
    activationTimes = cell(n_beats, 1);
    repolarizationTimes = cell(n_beats, 1);
    
    % for j = 1:nNodes
    %     Vj = Vsaved(j,:);
    %     [t_act, t_repo] = computeActiRepo(timeVec, Vj);
    %     activationTimes(j) = t_act;
    %     repolarizationTimes(j) = t_repo;
    % end

    for beat = 1:n_beats
        at = zeros(nNodes,1);
        rt = zeros(nNodes,1);
        for j = 1:nNodes
            if beat == 1
                Tj{j} =obj.T(1:obj.storeStep:end)*1e3;
                Vjj{j} = Vsaved(j,:);
            end
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

    E1 = el{1};
    [X,Y,~]=obj.pg.getCoo;
    dist_grid_1 = sqrt((X - E1(1)).^2 + (Y - E1(2)).^2);
    [~, which_1] = min(dist_grid_1(:));
    acti_e1 = activationTimes{end}(which_1);
    Latency = activationTimes{end} - acti_e1;

    plot_in_matlab_AT(obj, activationTimes{end}, repolarizationTimes{end}, Latency);
    
    
    nv = obj.pg.nlv;
    h = obj.pg.h;

    export_vtk_AT_RT(activationTimes,repolarizationTimes,nv,h);


%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Conduction Velocity Maps
%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Metodo 1: Gradiente della mappa di attivazione
    activationMap = reshape(activationTimes{end}, nv); % reshape in 2D/3D
    activationTimes_s = activationTimes{end}*1e-3;
    activationMap_s = activationMap *1e-3;
    
    switch nnz(nv > 1)
        case 1  % 1D
            dTx = gradient(activationMap_s, h(1));
            speedMap1 = 1 ./ abs(dTx);
            %speedMap1(speedMap1 > 5) = NaN;
        case 2  % 2D
            [dTx, dTy] = gradient(activationMap_s, h(1), h(2)); % in cm/s
            speedMap1 = 1 ./ sqrt(dTx.^2 + dTy.^2);
            %speedMap1(speedMap1 > 5) = NaN;
        case 3  % 3D
            [dTx, dTy, dTz] = gradient(activationMap_s, h(1), h(2), h(3));
            speedMap1 = 1 ./ sqrt(dTx.^2 + dTy.^2 + dTz.^2);
    end
    

    
    % Metodo 2: Average conduction velocity
    E1 = el{1};
    E2 = el{end};
    [X,Y,~]=obj.pg.getCoo;
    dist_grid_1 = sqrt((X - E1(1)).^2 + (Y - E1(2)).^2);
    dist_grid_2 = sqrt((X - E2(1)).^2 + (Y - E2(2)).^2);
    [~, which_1] = min(dist_grid_1(:));
    [~, which_2] = min(dist_grid_2(:));
    ds = sqrt((E1(1) - E2(1)).^2 + (E1(2) - E2(2)).^2);
    dt = activationTimes_s(which_2) - activationTimes_s(which_1);
    cv_average1 = ds / dt;
    disp(['Average velocity considering V sugli elettrodi= ' num2str(cv_average1)])

    P1 = selectedNodes{1};
    P2 = selectedNodes{end};
    [X,Y,~]=obj.pg.getCoo;
    ds = sqrt((X(P1)-X(P2)).^2 + (Y(P1)-Y(P2)).^2);
    dt = activationTimes_s(P2) - activationTimes_s(P1);
    cv_average2 = ds / dt;
    disp(['Average velocity considering V su due punti = ' num2str(cv_average2)])
    
    % 
    % 
    % % Metodo 3: Interpolazione a partire dalle Ue e ue
    % nk = length(Usaved);
    % activationTimes_el = zeros(nk, 1);
    % for k = 1:nk
    %     [~, idx_max] = max(Usaved{k});
    %     activationTimes_el(k) = timeVec(idx_max)*1e-3;  % in secondi
    % end
    % 
    % v_el = zeros(nk, 1);
    % for k = 1:nk
    %     if k == 1
    %         v_el(k) = NaN;  % nessuna velocità relativa a sé stesso
    %     else
    %         dx = el{k}(1) - el{1}(1);
    %         dy = el{k}(2) - el{1}(2);
    %         dist = sqrt(dx^2 + dy^2);  % cm
    %         dt = activationTimes_el(k) - activationTimes_el(1);  % s
    %         if dt > 0
    %             v_el(k) = dist / dt;  % cm/s
    %         else
    %             v_el(k) = NaN;  % evitiamo valori negativi o infiniti
    %         end
    %     end
    % end
    % 
    % disp(['Average velocity considering V = ' num2str(v_el(end))])
    
    

    plot_in_matlab_cv(obj, speedMap1{beat})

    export_vtk_CV(speedMap1, nv, h);

end
















function [t_act, t_repo, timeVec_new, Vtrace_new] = computeActiRepo(timeVec_old, Vtrace_old)
    % Trova il massimo e minimo per threshold
    % Vmax = max(Vtrace);
    % Vmin = min(Vtrace);
    Vthresh = -45;%Vmin + 0.9*(Vmax - Vmin);  % 90% di attivazione
    Vrepo = -50;%Vmin + 0.1*(Vmax - Vmin);    % 10% di ripolarizzazione

    % Trova indice di attivazione e ripolarizzazione
    idx_act = find(Vtrace_old > Vthresh, 1, 'first');
    idx_repo = find(Vtrace_old(idx_act:end) < Vrepo, 1, 'first') + idx_act - 1;

    if isempty(idx_act)
        t_act = NaN;
    else
        t_act = timeVec_old(idx_act);
    end
    if isempty(idx_repo)
        t_repo = NaN;
    else
        t_repo = timeVec_old(idx_repo);
        timeVec_new = timeVec_old(idx_repo+300:end);
        Vtrace_new = Vtrace_old(idx_repo+300:end);
    end
end


% 
% 
% 
% 
% function neighbors = getNeighbors(idx, nv)
%     % getNeighbors: restituisce gli indici dei nodi adiacenti in una griglia regolare
%     % idx: indice del nodo centrale (lineare)
%     % nv: dimensioni della griglia [nx, ny, nz]
% 
%     % Calcola le coordinate (i,j,k) del nodo
%     [I, J, K] = ind2sub(nv, idx);
% 
%     neighbors = [];
% 
%     directions = [ ...
%         -1,  0,  0;  % sinistra
%          1,  0,  0;  % destra
%          0, -1,  0;  % giù
%          0,  1,  0;  % su
%          0,  0, -1;  % dietro
%          0,  0,  1]; % davanti
% 
%     for d = 1:size(directions,1)
%         i2 = I + directions(d,1);
%         j2 = J + directions(d,2);
%         k2 = K + directions(d,3);
% 
%         if i2 >= 1 && i2 <= nv(1) && ...
%            j2 >= 1 && j2 <= nv(2) && ...
%            k2 >= 1 && k2 <= nv(3)
%             idx2 = sub2ind(nv, i2, j2, k2);
%             neighbors(end+1) = idx2; %#ok<AGROW>
%         end
%     end
% end
% 



function plot_in_matlab_0d(obj,timeVec,Vsaved,usaved,Usaved,selectedNodes)
  
    figure;
    hold on;
    for node = 1:length(selectedNodes)
        plot(timeVec, Vsaved(selectedNodes{node},:), 'LineWidth', 2);
    end
    xlabel('Time [ms]');
    ylabel('V [mV]');
    title('Transmembrane Potential at Selected Nodes');
    hold off;
    
    figure;
    hold on;
    for node = 1:length(selectedNodes)
        plot(timeVec, usaved(selectedNodes{node},:), 'LineWidth', 2);
    end
    xlabel('Time [ms]');
    ylabel('ue [mV]');
    title('Extracellular Potential at Selected Nodes');
    hold off;

    figure;
    for k = 1:length(obj.fk)
        plot(timeVec, Usaved{k}, 'LineWidth', 2);
        hold on
    end
    xlabel('Time [ms]');
    ylabel('U [mV]');
    title('Field Potential of electrodes')

end



function plot_in_matlab_AT(obj, activationTimes, repolarizationTimes, Latency)
    % Plot activation and repolarization maps (if 2D or 3D)
    [X,Y,Z] = obj.pg.getCoo;

    sz = size(X);
    AtMap = reshape(activationTimes, sz);
    RtMap = reshape(repolarizationTimes, sz);
    LatencyMap = reshape(Latency, sz);
    
    figure;
    surf(X, Y, Z, AtMap); %,'EdgeColor','none');
    colorbar;
    title('Activation Time Map');
    xlabel('x'); ylabel('y'); zlabel('z');
    view(2)

    figure;
    surf(X, Y, Z, RtMap); %,'EdgeColor','none');
    colorbar;
    title('Repolarization Time Map');
    xlabel('x'); ylabel('y'); zlabel('z');
    view(2)

    figure;
    surf(X, Y, Z, LatencyMap); %,'EdgeColor','none');
    colorbar;
    title('Latency Map');
    xlabel('x'); ylabel('y'); zlabel('z');
    view(2)

    %figure
    %contour(AtMap,100)
end


function export_vtk_AT_RT(activationTimes, repolarizationTimes, nv, h)

    fid = fopen('./outputs/mea_AT_RT.vtk', 'w');
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'Example structured points dataset\n');
    fprintf(fid, 'ASCII\n\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');

    % Adjust for 1D or 2D geometries
    if nv(2) ~= 1 && nv(3) == 1
        nv(3) = 2;
        for i = 1:length(activationTimes)
            activationTimes{i} = repmat(activationTimes{i}, 1, 2);
            repolarizationTimes{i} = repmat(repolarizationTimes{i}, 1, 2);
        end
    elseif nv(2) == 1 && nv(3) == 1
        nv(2:3) = 2;
        for i = 1:length(activationTimes)
            activationTimes{i} = repmat(activationTimes{i}, 1, 4);
            repolarizationTimes{i} = repmat(repolarizationTimes{i}, 1, 4);
        end
    end

    fprintf(fid, 'DIMENSIONS %d %d %d\n', nv(1), nv(2), nv(3));
    fprintf(fid, 'ORIGIN 0 0 0\n');
    fprintf(fid, 'SPACING %f %f %f\n\n', h(1), h(2), h(3));

    npoints = prod(nv);
    fprintf(fid, 'POINT_DATA %d\n', npoints);

    % Write each AT_i and RT_i from cell arrays
    for i = 1:length(activationTimes)
        AT = activationTimes{i};
        RT = repolarizationTimes{i};

        % Ensure column vector
        AT = AT(:);
        RT = RT(:);

        % Activation Time i
        fprintf(fid, 'SCALARS AT_%d float 1\n', i);
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fprintf(fid, '%f\n', AT);

        % Repolarization Time i
        fprintf(fid, 'SCALARS RT_%d float 1\n', i);
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fprintf(fid, '%f\n', RT);
    end

    fclose(fid);
end



function plot_in_matlab_cv(obj, speedMap1)
    [X, Y, Z] = obj.pg.getCoo;
    sz = size(X);  % dimensione mesh

    CVmap1 = reshape(speedMap1, sz);

    % Plot metodo 1
    figure;
    surf(X, Y, Z, CVmap1, 'EdgeColor', 'none');
    xlabel('x [cm]'); ylabel('y [cm]'); zlabel('z [cm]');
    colorbar; view(2); axis equal tight;
    clim([10 25]); 

    
    % Plot della mappa interpolata
    % figure;
    % imagesc(linspace(min(X(:)), max(X(:)), 200), ...
    %         linspace(min(Y(:)), max(Y(:)), 200), Vinterp);
    % set(gca, 'YDir', 'normal'); axis equal tight;
    % colorbar;
    % title('Mappa interpolata della velocità di conduzione');
    % xlabel('x [cm]');
    % ylabel('y [cm]');

end



function export_vtk_CV(speedMaps, nv, h)
    % Esporta più conduction velocity map in un file VTK (una per campo)
    % speedMaps: cell array di mappe di velocità

    if ~iscell(speedMaps)
        error('speedMaps deve essere una cella.');
    end

    % Adatta le dimensioni per 1D o 2D
    if nv(2) ~= 1 && nv(3) == 1
        nv(3) = 2;
        for i = 1:length(speedMaps)
            speedMaps{i} = repmat(speedMaps{i}, 1, 2);
        end
    elseif nv(2) == 1 && nv(3) == 1
        nv(2:3) = 2;
        for i = 1:length(speedMaps)
            speedMaps{i} = repmat(speedMaps{i}, 1, 4);
        end
    end

    % Scrittura file
    fid = fopen('./outputs/mea_CV_maps.vtk', 'w');
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'Conduction Velocity Maps\n');
    fprintf(fid, 'ASCII\n\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS %d %d %d\n', nv(1), nv(2), nv(3));
    fprintf(fid, 'ORIGIN 0 0 0\n');
    fprintf(fid, 'SPACING %f %f %f\n\n', h(1), h(2), h(3));
    
    npoints = prod(nv);
    fprintf(fid, 'POINT_DATA %d\n', npoints);

    % Scrive ogni mappa come SCALARS CVi
    for i = 1:length(speedMaps)
        CV = speedMaps{i}(:);
        if length(CV) ~= npoints
            fclose(fid);
            error('La mappa %d ha dimensioni non compatibili con nv.', i);
        end
        fprintf(fid, 'SCALARS CV%d float 1\n', i);
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fprintf(fid, '%f\n', CV);
    end

    fclose(fid);
end
