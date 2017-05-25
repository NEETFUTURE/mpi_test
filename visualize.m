
fig = figure;
files = dir('data_mpi_test');
num = size(files);
n = 0;

for i = 1:num(1)
    if(files(i).isdir == 1)
        continue
    end
    n = n+1;
    fname = strcat(files(i).folder,"\",files(i).name);
    id = fopen(fname, 'r');
    mydata{i} = fread(id,'double');
    plot(mydata{i});
    ylim([-1 1]);
    drawnow
    F1(n) = getframe(gcf);
    fclose(id);
end

%%
fig = figure;
set(fig,'doublebuffer','on')
movie(fig,F1,1)
