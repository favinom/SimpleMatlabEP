clear; close all; clc;

addpath('../.')

mesh_sizes = [64, 128, 256, 512, 1024];

Xf = 1; Yf = 1;

% Nome file con mesh e dominio
filename = sprintf('results.csv');
fid = fopen(filename, 'w');
fprintf(fid, 'h,H,iter_dd,err_rel\n');
fclose(fid);

for idx = 1:length(mesh_sizes)
    nex = mesh_sizes(idx);
    ney = mesh_sizes(idx);
    % nex = mesh_size;
    % ney = mesh_size;
    
    hx=Xf/nex;
    hy=Xf/ney;
    
    pg=PointGrid([nex ney]+1,[hx hy]);
    
    % Matrici
    [M, L] = assembleMatrices(pg);
    L = L + 0.001 * M; % regolarizzazione
    
    % Forzante
    [X,Y]=pg.getCoo;
    f = @(x, y) sin(pi*x) .* sin(pi*y);
    F = f(X, Y);
    b = M * F(:);
    
    % % Soluzione globale
    % fprintf('  [*] PCG Full...\n');
    % tic;
    [x_full, flag_full, res_full, iter_full] = pcg(L, b, 1e-5, 5000);
    % time_full = toc;
    if flag_full ~= 0
        disp('pcg non converge')
    end
        
    nsd_x=[1,2,3,4,5,6,7,8,9,10];
    
    
    for sd = 1:length(nsd_x)
        nsd=nsd_x(sd);
        
        fprintf('\n--- NUMBER SUBDOMAINS %d ---\n', nsd*nsd);
        
        [id_S,id]=buildSubdomainIds(pg,nsd,nsd);
        % plotSubid(pg,id_S,id);
    
        % Metodo di Steklov-Poincaré
        fprintf('  [*] PCG con Steklov-Poincaré...\n');
        x0 = zeros(size(b));
        sp = SteklovPoincare(L, b, x0);
        sp.prepare(id_S, id);
        sp.prepareRhs(id);
        sp.prepareApp();
    
        tic;
        [x_dd, flag_dd, res_dd, iter_dd] = sp.solve(1e-7, 10000);
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
        h = 1/nex;
        H = 1/nsd;
        fprintf('  h = %.4f , \t H = %.4f \n', h, H);
        % fprintf('  Tempo PCG globale:     %.4f s (iter: %d)\n', time_full, iter_full);
        % fprintf('  Tempo SP (DD):         %.4f s (iter: %d)\n', time_dd, iter_dd);
        fprintf('  Numero iterazioni SP (DD):         iter: %d\n', iter_dd);
        fprintf('  Errore relativo (DD):  %.2e\n', err_rel);
    
        % Visualizzazione confronto
        % figure
        % subplot(1,2,1)
        % surf(X, Y, reshape(x_full, size(X)))
        % shading interp
        % title('Soluzione completa')
        % 
        % subplot(1,2,2)
        % surf(X, Y, reshape(x_rec, size(X)))
        % shading interp
        % title('Soluzione SP')
    
    
        % Pulizia
        clear id id_S sp;
    
       
        % Riga da scrivere: h, H, iter_dd, err_rel
        data_row = [h, H, iter_dd, err_rel];
        
      
        % Aggiunta della riga
        dlmwrite(filename, data_row, '-append');
    
    end
end
