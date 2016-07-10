function visualise

% Matlab function to visualise the result of a least squares 
% fit calculation. "inp" is the input file, "outp" the output file

fid = fopen('lsdata.txt');
n = fscanf(fid,'%d',1);

t = fscanf(fid,'%g',n);
y = fscanf(fid,'%g',n);
fclose(fid);

plot(t,y,'k.','MarkerSize',20)
hold on
 
fid = fopen('lsoutput.txt');
d = fscanf(fid,'%d',1);

x = fscanf(fid,'%g',d+1);
fclose(fid);

tfine = linspace(t(1),t(n))';

plot(tfine,polyval(x(d+1:-1:1),tfine));
hold off

