function dummy = postprocess(inp)

fid = fopen(inp);

m = fscanf(fid,'%d',1);
u = fscanf(fid,'%g',[m-1,m-1]);

fclose(fid);

%x = (1/m)*[0:m]';

surf(u);

