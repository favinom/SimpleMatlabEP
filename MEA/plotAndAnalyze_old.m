function plotAndAnalyze_old(obj,el)
    % Plot a single time series for one node (e.g., the center)
    centerNode = round(obj.pg.nv / 2);
    timeVec = obj.T(1:obj.storeStep:end); % Tempo associato a Vsave
    Usaved = cell(size(obj.Usave));
    
    if obj.ionicModelType==3 || obj.ionicModelType==4 || obj.ionicModelType==5
        timeVec = timeVec*1e3;
        Vsaved = obj.Vsave*1e3;
        usaved = obj.usave*1e3;
        for k = 1:length(obj.fk)
            Usaved{k} = obj.Usave{k}*1e3;
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% 
% PLOT IN MATLAB AND CSV SINGLE NODE
%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% 
  
    plot_in_matlab_0d(obj,timeVec,Vsaved,usaved,Usaved,centerNode);

    store_csv = 5;

    time_csv = timeVec(1:store_csv:end);
    V_csv = Vsaved(centerNode,1:store_csv:end);
    u_csv = usaved(centerNode,1:store_csv:end);
    
    T = table(time_csv(:), V_csv(:), u_csv(:), 'VariableNames', {'Time_ms', 'Voltage_mV', 'ue_mV'});

    % Aggiungi ogni U_csv{k} come colonna separata
    for k = 1:length(obj.fk)
        U_csv = Usaved{k}(1:store_csv:end);
        varName = sprintf('U%d_mV', k);
        T.(varName) = U_csv(:);
    end
    
    
    filename = './outputs/voltage_traces.csv';
    writetable(T, filename);
    writetable(T, filename);
    
    disp(['File salvato: ', filename]);
    

%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% 
% Activation & Repolarization Time Maps
%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% 
     
    nNodes = size(obj.Vsave,1);
    activationTimes = zeros(nNodes,1);
    repolarizationTimes = zeros(nNodes,1);
    
    for j = 1:nNodes
        Vj = Vsaved(j,:);
        [t_act, t_repo] = computeActiRepo(timeVec, Vj);
        activationTimes(j) = t_act;
        repolarizationTimes(j) = t_repo;
    end

    plot_in_matlab_AT(obj, activationTimes, repolarizationTimes);
    
    
    nv = obj.pg.nlv;
    h = obj.pg.h;

    export_vtk_AT_RT(activationTimes,repolarizationTimes,nv,h);


%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Conduction Velocity Maps
%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Metodo 1: Gradiente della mappa di attivazione
    activationMap = reshape(activationTimes, nv); % reshape in 2D/3D
    activationTimes_s = activationTimes*1e-3;
    activationMap_s = activationMap *1e-3;
    
    switch nnz(nv > 1)
        case 1  % 1D
            dTx = gradient(activationMap, h(1));
            speedMap1 = 1 ./ abs(dTx);
            %speedMap1(speedMap1 > 5) = NaN;
        case 2  % 2D
            [dTx, dTy] = gradient(activationMap_s, h(1), h(2)); % in cm/s
            speedMap1 = 1 ./ sqrt(dTx.^2 + dTy.^2);
            %speedMap1(speedMap1 > 5) = NaN;
        case 3  % 3D
            [dTx, dTy, dTz] = gradient(activationMap, h(1), h(2), h(3));
            speedMap1 = 1 ./ sqrt(dTx.^2 + dTy.^2 + dTz.^2);
    end
    

    % Metodo 2: distanza su tempo
    [X, Y, Z] = obj.pg.getCoo;
    X = X(:); Y = Y(:); Z = Z(:);
    nNodes = length(activationTimes);
    speedMap2 = NaN(nNodes, 1);
    min_dt = 1e-6;
    for i = 1:nNodes
        neighbors = getNeighbors(i, obj.pg.nlv);  
        for j = neighbors
            % Calcola la distanza tra i nodi i e j
            dx = X(i) - X(j);
            dy = Y(i) - Y(j);
            dz = Z(i) - Z(j);
            d = sqrt(dx^2 + dy^2 + dz^2); % in centimetri
            dt = abs(activationTimes_s(i) - activationTimes_s(j)); % in secondi
            if dt > min_dt
                v = d / dt; % in cm/s
                if isnan(speedMap2(i)) || v < speedMap2(i)
                    speedMap2(i) = v;
                end
            end
        end
    end



    % Metodo 3: Average conduction velocity
    E1 = el{1};
    E2 = el{end};
    [X,Y,Z]=obj.pg.getCoo;
    dist_grid_1 = sqrt((X - E1(1)).^2 + (Y - E1(2)).^2);
    dist_grid_2 = sqrt((X - E2(1)).^2 + (Y - E2(2)).^2);
    [~, which_1] = min(dist_grid_1(:));
    [~, which_2] = min(dist_grid_2(:));
    ds = sqrt((E1(1) - E2(1)).^2 + (E1(2) - E2(2)).^2);
    dt = activationTimes_s(which_2) - activationTimes_s(which_1);
    cv_average = ds / dt;
    disp(['Average velocity considering V = ' num2str(cv_average)])
    


    % Metodo 4: Interpolazione a partire dalle Ue
    nk = length(Usaved);
    activationTimes_el = zeros(nk, 1);
    for k = 1:nk
        [~, idx_max] = max(Usaved{k});
        activationTimes_el(k) = timeVec(idx_max)*1e-3;  % in secondi
    end
    
    v_el = zeros(nk, 1);
    for k = 1:nk
        if k == 1
            v_el(k) = NaN;  % nessuna velocità relativa a sé stesso
        else
            dx = el{k}(1) - el{1}(1);
            dy = el{k}(2) - el{1}(2);
            dist = sqrt(dx^2 + dy^2);  % cm
            dt = activationTimes_el(k) - activationTimes_el(1);  % s
            if dt > 0
                v_el(k) = dist / dt;  % cm/s
            else
                v_el(k) = NaN;  % evitiamo valori negativi o infiniti
            end
        end
    end

    disp(['Average velocity considering V = ' num2str(v_el(end))])
    
    % Plot delle velocità per elettrodi
    % figure;
    % plot(1:nk, v_el, 'o-', 'LineWidth', 1.5);
    % xlabel('Elettrodo');
    % ylabel('Velocità [cm/s]');
    % title('Velocità di conduzione relativa all'elettrodo 1');
    % grid on;
    
    
    
    
    % Coordinate centri degli elettrodi
    el_x = cellfun(@(c) c(1), el);
    el_y = cellfun(@(c) c(2), el);

    % Coordinate degli elettrodi validi
    valid_idx = ~isnan(v_el);
    % Coordinate degli elettrodi validi
    vx = el_x(valid_idx);
    vy = el_y(valid_idx);
    vz = v_el(valid_idx);
    
    % Crea griglia sul dominio completo
    [Xq, Yq] = meshgrid(linspace(min(X(:)), max(X(:)), 500), ...
                        linspace(min(Y(:)), max(Y(:)), 500));
    
    % Interpolazione primaria (solo nei punti interni)
    Vinterp_linear = griddata(vx, vy, vz, Xq, Yq, 'linear');
    
    % Seconda interpolazione per riempire i NaN (con nearest neighbor)
    F = scatteredInterpolant(vx(:), vy(:), vz(:), 'linear', 'nearest');  % extrapola col valore più vicino
    Vinterp_extrap = F(Xq, Yq);

    
    % Combina: usa i valori "buoni" della prima interpolazione, e completa con la seconda
    Vinterp_combined = Vinterp_linear;
    Vinterp_combined(isnan(Vinterp_combined)) = Vinterp_extrap(isnan(Vinterp_combined));
    
    % Resize e trasposizione per allineare alla mappa
    Vinterp_resized = imresize(Vinterp_combined, [nv(2), nv(1)]);
    Vinterp_resized = Vinterp_resized';  % perché imresize lavora con righe-colonne
    

    plot_in_matlab_cv(obj, speedMap1, speedMap2, Vinterp_resized)

    export_vtk_CV(speedMap1, speedMap2, nv, h, Vinterp_resized);

end
















function [t_act, t_repo] = computeActiRepo(timeVec, Vtrace)
    Vthresh = -45;
    Vrepo = -50;

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






function neighbors = getNeighbors(idx, nv)
    % getNeighbors: restituisce gli indici dei nodi adiacenti in una griglia regolare
    % idx: indice del nodo centrale (lineare)
    % nv: dimensioni della griglia [nx, ny, nz]

    % Calcola le coordinate (i,j,k) del nodo
    [I, J, K] = ind2sub(nv, idx);

    neighbors = [];

    directions = [ ...
        -1,  0,  0;  % sinistra
         1,  0,  0;  % destra
         0, -1,  0;  % giù
         0,  1,  0;  % su
         0,  0, -1;  % dietro
         0,  0,  1]; % davanti

    for d = 1:size(directions,1)
        i2 = I + directions(d,1);
        j2 = J + directions(d,2);
        k2 = K + directions(d,3);

        if i2 >= 1 && i2 <= nv(1) && ...
           j2 >= 1 && j2 <= nv(2) && ...
           k2 >= 1 && k2 <= nv(3)
            idx2 = sub2ind(nv, i2, j2, k2);
            neighbors(end+1) = idx2; %#ok<AGROW>
        end
    end
end




function plot_in_matlab_0d(obj,timeVec,Vsaved,usaved,Usaved,centerNode)
  
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
    
    figure;
    for k = 1:length(obj.fk)
        plot(timeVec, Usaved{k}, 'LineWidth', 2);
        hold on
    end
    xlabel('Time [ms]');
    ylabel('U [mV]');
    title('Field Potential of electrodes')

end



function plot_in_matlab_AT(obj, activationTimes, repolarizationTimes)
    % Plot activation and repolarization maps (if 2D or 3D)
    [X,Y,Z] = obj.pg.getCoo;

    sz = size(X);
    AtMap = reshape(activationTimes, sz);
    RtMap = reshape(repolarizationTimes, sz);
    
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
end


function export_vtk_AT_RT(activationTimes,repolarizationTimes,nv,h)


    fid = fopen('./outputs/AT_RT.vtk', 'w');
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


end



function plot_in_matlab_cv(obj, speedMap1, speedMap2, Vinterp)
    [X, Y, Z] = obj.pg.getCoo;
    sz = size(X);  % dimensione mesh

    CVmap1 = reshape(speedMap1, sz);
    CVmap2 = reshape(speedMap2, sz);

    % Plot metodo 1
    figure;
    surf(X, Y, Z, CVmap1, 'EdgeColor', 'none');
    xlabel('x [cm]'); ylabel('y [cm]'); zlabel('z [cm]');
    colorbar; view(2); axis equal tight;
    %clim([3 5]); 

    % Plot metodo 2
    figure;
    surf(X, Y, Z, CVmap2, 'EdgeColor', 'none');
    xlabel('x [cm]'); ylabel('y [cm]'); zlabel('z [cm]');
    colorbar; view(2); axis equal tight;
    %clim([0 0.01]); 

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



function export_vtk_CV(speedMap1, speedMap2, nv, h, Vinterp)
    % Esporta le conduction velocity map e la mappa interpolata Vinterp

    if nv(2) ~= 1 && nv(3) == 1
        nv(3) = 2;
        speedMap1 = repmat(speedMap1, 1, 2);
        speedMap2 = repmat(speedMap2, 1, 2);
        Vinterp = repmat(Vinterp, 1, 2);
    elseif nv(2) == 1 && nv(3) == 1
        nv(2:3) = 2;
        speedMap1 = repmat(speedMap1, 1, 4);
        speedMap2 = repmat(speedMap2, 1, 4);
        Vinterp = repmat(Vinterp, 1, 4);
    end

    fid = fopen('./outputs/CV_maps.vtk', 'w');
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'Conduction Velocity Maps\n');
    fprintf(fid, 'ASCII\n\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS %d %d %d\n', nv(1), nv(2), nv(3));
    fprintf(fid, 'ORIGIN 0 0 0\n');
    fprintf(fid, 'SPACING %f %f %f\n\n', h(1), h(2), h(3));
    fprintf(fid, 'POINT_DATA %d\n', prod(nv));

    fprintf(fid, 'SCALARS CV1 float 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%f\n', speedMap1(:));

    fprintf(fid, 'SCALARS CV2 float 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%f\n', speedMap2(:));

    fprintf(fid, 'SCALARS CVinterp float 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%f\n', Vinterp(:));

    fclose(fid);
end

