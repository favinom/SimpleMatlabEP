clear; close all; clc;

mesh_sizes = [64, 128, 256, 512];  % Puoi aggiungere 512, 1024 se vuoi 
Xf = 1; Yf = 1;

for idx = 1:length(mesh_sizes)
    nex = mesh_sizes(idx);
    ney = mesh_sizes(idx);
    
    fprintf('\n--- MESH SIZE: %d x %d ---\n', nex, ney);
    
    % Griglia
    x = linspace(0, Xf, nex+1)'; hx = diff(x);
    y = linspace(0, Yf, ney+1)'; hy = diff(y);
    [X, Y] = ndgrid(x, y);
    node_id = 1:(nex+1)*(ney+1);

    % Matrici
    [M, L] = assembleMatrices_old(hx, hy);
    L = L + 0.001 * M; % regolarizzazione

    % Forzante
    f = @(x, y) sin(pi*x) .* sin(pi*y);
    F = f(X, Y);
    b = M * F(:);
    
    % Decomposizione dominio
    id{1} = find((X < Xf/2) & (Y < Yf/2));
    id{2} = find((X > Xf/2) & (Y < Yf/2));
    id{3} = find((X < Xf/2) & (Y > Yf/2));
    id{4} = find((X > Xf/2) & (Y > Yf/2));

    id_S = node_id;
    for j = 1:length(id)
        id_S = setdiff(id_S, id{j});
    end

    % Soluzione globale
    fprintf('  [*] PCG Full...\n');
    tic;
    [x_full, flag_full, res_full, iter_full] = pcg(L, b, 1e-6, 1500);
    time_full = toc;
    if flag_full ~= 0
        disp('pcg non converge')
    end

    % Metodo di Steklov-Poincaré
    fprintf('  [*] PCG con Steklov-Poincaré...\n');
    x0 = zeros(size(b));
    sp = SteklovPoincare(L, b, x0);
    sp.prepare(id_S, id);
    sp.prepareRhs(id);
    sp.prepareApp();

    tic;
    [x_dd, flag_dd, res_dd, iter_dd] = sp.solve(1e-7, 1000);
    time_dd = toc;
    if flag_dd ~= 0
        disp('pcg non converge')
    end

    % Ricostruzione soluzione completa da SP
    x_rec = x0;
    x_rec(id_S) = sp.xs;

    for i = 1:length(id)
        xi = sp.H{i} \ (sp.H{i}' \ (sp.bi{i} - sp.Ais{i} * sp.xs));
        x_rec(id{i}) = xi;
    end

    % Errore relativo
    err_rel = norm(x_full - x_rec) / norm(x_full);

    % Risultati
    fprintf('  Tempo PCG globale:     %.4f s (iter: %d)\n', time_full, iter_full);
    fprintf('  Tempo SP (DD):         %.4f s (iter: %d)\n', time_dd, iter_dd);
    fprintf('  Errore relativo (DD):  %.2e\n', err_rel);

    % Pulizia
    clear id id_S sp;
end
