clear all
close all

%% Domini e griglia e matrici
Xf = 1; 
nex = 128;
x = linspace(0, Xf, nex+1)'; 
hx = diff(x);

Yf = 1; 
ney = 128;
y = linspace(0, Yf, ney+1)'; 
hy = diff(y);

[X, Y] = ndgrid(x, y);
node_id = 1:(nex+1)*(ney+1);

[M, L] = assembleMatrices_old(hx, hy);
L = L + 0.001 * M; % sistema regolarizzato

%% Forzante
f = @(x, y) sin(pi*x) .* sin(pi*y); % ad esempio
F = f(X, Y); 
b = M * F(:); % membro di destra

%% Decomposizione in sottodomini
id{1} = find((X < Xf/2) & (Y < Yf/2));
id{2} = find((X > Xf/2) & (Y < Yf/2));
id{3} = find((X < Xf/2) & (Y > Yf/2));
id{4} = find((X > Xf/2) & (Y > Yf/2));

id_S = node_id;
for j = 1:length(id)
    id_S = setdiff(id_S, id{j});
end

%% Soluzione completa (PCG su tutto il dominio)
disp('PCG Full:')
tic
[x_full, flag_full, res_full, iter_full] = pcg(L, b, 1e-6, 1500);
toc

%% Soluzione con Steklov-Poincaré
disp('PCG con Steklov-Poincaré:')
x0 = zeros(size(b));
sp = SteklovPoincare(L, b, x0);

sp.prepare(id_S, id);
sp.prepareRhs(id);
sp.prepareApp();
tic
[x_dd, flag_dd, res_dd, iter_dd] = sp.solve(1e-7, 1000);
toc

%% Ricostruzione soluzione completa
x_rec = x0;
x_rec(id_S) = sp.xs;

for i = 1:length(id)
    % risolvi localmente
    xi = sp.H{i} \ (sp.H{i}' \ (sp.bi{i} - sp.Ais{i} * sp.xs));
    x_rec(id{i}) = xi;
end

%% Visualizzazione confronto
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

%% Errore relativo
err_rel = norm(x_full - x_rec) / norm(x_full);
fprintf('Errore relativo: %.2e\n', err_rel);
