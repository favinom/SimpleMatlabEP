function exportVTK_MEA(pg, cellData)

nv = pg.nlv;
h = pg.h;

% Open VTK file
fid = fopen(['electrodes.vtk'], 'w');
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'Example structured points dataset\n');
fprintf(fid, 'ASCII\n\n');
fprintf(fid, 'DATASET STRUCTURED_POINTS\n');


fprintf(fid, 'DIMENSIONS %d %d %d\n', nv(1), nv(2), nv(3));
fprintf(fid, 'ORIGIN 0 0 0\n');
fprintf(fid, 'SPACING %f %f %f\n\n', h(1), h(2), h(3));


if nv(2) ~= 1 && nv(3) == 1
    nv(3) = 2;
    cellData = repmat(cellData, 1, 2);
elseif nv(2) == 1 && nv(3) == 1
    nv(2:3) = 2;
    cellData = repmat(cellData, 1, 4);
end

% ===============================
% Export CELL DATA (dummy example)
% ===============================

% Compute number of cells in each dimension
nc = nv - 1;
ncells = prod(nc);
fprintf(fid, '\nCELL_DATA %d\n', ncells);

% Dummy scalar cell data (e.g., average of V over each cell)
% For simplicity, let's use zeros or ones — replace with real data if needed
%cellData = ones(ncells, 1);  % Placeholder

fprintf(fid, 'SCALARS Example_Cell_Scalar float 1\n');
fprintf(fid, 'LOOKUP_TABLE default\n');
fprintf(fid, '%f\n', cellData);

fclose(fid);

end

