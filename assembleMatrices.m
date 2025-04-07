function [M,L]=assembleMatrices(pg)

dim=pg.dim;

if dim==0
    M=1;
    L=0;
end

if dim==1
    [M,L]=assembleMatrices1D(pg);
end

if dim==2
    [M,L]=assembleMatrices2D(pg);
end

if dim==3
    [M,L]=assembleMatrices3D(pg);
end
