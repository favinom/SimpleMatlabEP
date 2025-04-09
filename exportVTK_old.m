function exportVTK_old(V,pg,i)

nv=pg.nlv;
h=pg.h;

istr=num2str(i);
if i<10
    istr=['0',istr];
end
if i<100
    istr=['0',istr];
end
if i<1000
    istr=['0',istr];
end
if i<10000
    istr=['0',istr];
end

fid=fopen(['sol_',istr,'.vtk'],'w');
fprintf(fid,'# vtk DataFile Version 3.0\n');


fprintf(fid,'Example structured points dataset\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'\n');
fprintf(fid,'DATASET STRUCTURED_POINTS\n');
if nv(2)~=1 && nv(3)==1
    nv(3)=2*nv(3);
    V=[V V];
elseif nv(2)==1 && nv(3)==1
    nv(2:3)=2*nv(2:3);
    V=[V V V V];
end
fprintf(fid,'DIMENSIONS %d %d %d\n',nv(1),nv(2),nv(3));

fprintf(fid,'ORIGIN 0 0 0\n');
fprintf(fid,'SPACING %f %f %f\n',h(1),h(2),h(3));
fprintf(fid,'\n');

fprintf(fid,'POINT_DATA %d\n',prod(nv));
fprintf(fid,'SCALARS Potential float\n');
fprintf(fid,'LOOKUP_TABLE default\n');

fprintf(fid,'%f\n',V);
%fprintf(fid,'%f\n',V);

fclose(fid);

return
