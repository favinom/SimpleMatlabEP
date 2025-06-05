function exportVTK(V, u, pg, i, flag_MB)

nv = pg.nlv;
h = pg.h;

% Format timestep index with leading zeros (e.g., 0001)
istr = sprintf('%04d', i);

% Open VTK file
fid = fopen(['./outputs/sol_', istr, '.vtk'], 'w');
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'Example structured points dataset\n');
fprintf(fid, 'ASCII\n\n');
fprintf(fid, 'DATASET STRUCTURED_POINTS\n');

% Adjust for 1D or 2D geometries
if nv(2) ~= 1 && nv(3) == 1
    nv(3) = 2;
    V = repmat(V, 1, 2);
    u = repmat(u, 1, 2);
elseif nv(2) == 1 && nv(3) == 1
    nv(2:3) = 2;
    V = repmat(V, 1, 4);
    u = repmat(u, 1, 4);
end

fprintf(fid, 'DIMENSIONS %d %d %d\n', nv(1), nv(2), nv(3));
fprintf(fid, 'ORIGIN 0 0 0\n');
fprintf(fid, 'SPACING %f %f %f\n\n', h(1), h(2), h(3));

npoints = prod(nv);
fprintf(fid, 'POINT_DATA %d\n', npoints);

% Transmembrane Potential
fprintf(fid, 'SCALARS Transmembrane_Potential float 1\n');
fprintf(fid, 'LOOKUP_TABLE default\n');
fprintf(fid, '%f\n', V(:));

if flag_MB
    % Extracellular Potential
    fprintf(fid, 'SCALARS Extracellular_Potential float 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%f\n', u(:));
end

fclose(fid);
end

