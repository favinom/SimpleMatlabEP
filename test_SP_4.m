clear; close all; clc;

mesh_sizes = [64, 128, 256, 512, 1024]; % puoi aggiungere 1024 se vuoi
Xf = 1; Yf = 1;

for idx = 1:length(mesh_sizes)
    nex = mesh_sizes(idx);
    ney = mesh_sizes(idx);
    
    fprintf('\n--- MESH SIZE: %d x %d ---\n', nex, ney);

    % Mesh generation
    x = linspace(0,Xf,nex+1)'; hx = diff(x);
    y = linspace(0,Yf,ney+1)'; hy = diff(y);
    [X,Y] = ndgrid(x,y);
    
    % Assembla matrici
    [M,L] = assembleMatrices_old(hx,hy);
    L = L + 0.0001*M; % Aggiungo massa per regolarizzazione

    node_id = 1:(nex+1)*(ney+1);
    id{1} = find((X<Xf/2) & (Y<Yf/2));
    id{2} = find((X>Xf/2) & (Y<Yf/2));
    id{3} = find((X<Xf/2) & (Y>Yf/2));
    id{4} = find((X>Xf/2) & (Y>Yf/2));
    id_S = node_id;

    for j = 1:length(id)
        id_S = setdiff(id_S,id{j});
    end

    % Decomposizione dominio
    L_S = L(id_S,id_S);
    for j = 1:length(id)
        Lii{j} = L(id{j},id{j});
        H{j} = chol(Lii{j});
        LSi{j} = L(id_S,id{j});
        LiS{j} = L(id{j},id_S);
    end

    % Vettore RHS
    V = ones(nex+1,ney+1);
    b = V(id_S)';

    % --------------------
    % PCG con Domain Decomposition
    % --------------------
    fcn_fact = @(x) app_SP_fact(x,L_S,H,LSi,LiS);

    tic;
    [x_dd,flag_dd,relres_dd,iter_dd] = pcg(fcn_fact,b,1e-7,1000);
    time_dd = toc;

    % --------------------
    % PCG globale
    % --------------------
    b_full = V(:);
    tic;
    [x_full,flag_full,relres_full,iter_full] = pcg(L,b_full,1e-8,1000);
    time_full = toc;

    % --------------------
    % Confronto risultati
    % --------------------
    x_dd_full = zeros(size(L,1),1);
    x_dd_full(id_S) = x_dd;
    
    res_full = norm(b_full - L*x_full);
    res_dd = norm(b_full - L*x_dd_full);
    err_rel_dd = norm(x_full(id_S) - x_dd) / norm(x_full(id_S));
    
    fprintf('  Tempo PCG globale:     %.4f s (iter: %d)\n', time_full, iter_full);
    fprintf('  Residuo globale:       %.2e\n', res_full);
    fprintf('  Tempo PCG con DD:      %.4f s (iter: %d)\n', time_dd, iter_dd);
    fprintf('  Residuo con DD:        %.2e\n', res_dd);
    fprintf('  Errore relativo (DD):  %.2e\n', err_rel_dd);

    % Pulizia
    clear id id_S Lii H LSi LiS;
end