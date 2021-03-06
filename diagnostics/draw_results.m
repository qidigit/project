clear
close all

set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultAxesFontSize', 16)

% model resolutions and integration time
% depends on the problems settings
NY = 256;
m = ceil(NY/2);
M_rsv = 170;
M = 171;
TT = 50;
dt=600;
dta=120*600;
stri=dta/dt;
nstep = (TT*24*60*60)/dt/stri;

nproc = 2;
ind_name = ['../model_data/ind_r', sprintf('%03d', 0), '.bin'];
file = fopen(ind_name, 'rb');
ind  = fread(file, 3, 'int');
sind = fread(file, 3, 'int');
leg_points = fread(file, ceil(NY/2), 'double');
fclose(file);

lon = 180/pi*(0:(pi/NY):2*pi-pi/NY)';
lat = 180/pi*[-asin(leg_points(end:-1:1)); asin(leg_points)];
[X, Y] = meshgrid(lon, lat);

u = zeros(NY, 2*NY);
v = zeros(NY, 2*NY);
vort = zeros(NY, 2*NY);

route = '../model_data';

%open files
for i = 1:nproc
    ug_name = ['../model_data/u_r', sprintf('%03d', i-1), '.bin'];
    vg_name = ['../model_data/v_r', sprintf('%03d', i-1), '.bin'];
    vorg_name = ['../model_data/vorg_r', sprintf('%03d', i-1), '.bin'];
%     vors_name = ['../model_data/vors_r', sprintf('%03d', i-1), '.bin'];

    file1(i) = fopen(ug_name, 'rb');
    file2(i) = fopen(vg_name, 'rb');
    file3(i) = fopen(vorg_name, 'rb');
end

for l = 1:nstep+1
    
    
for i = 1:nproc
    ind_name = ['../model_data/ind_r', sprintf('%03d', i-1), '.bin'];
%     ug_name = ['../model_data/u_r', sprintf('%03d', i-1), '.bin'];
%     vg_name = ['../model_data/v_r', sprintf('%03d', i-1), '.bin'];
%     vorg_name = ['../model_data/vorg_r', sprintf('%03d', i-1), '.bin'];
%     vors_name = ['../model_data/vors_r', sprintf('%03d', i-1), '.bin'];
    
    file = fopen(ind_name, 'rb');
    ind  = fread(file, 3, 'int');
    fclose(file);
    
    for j = 1:ind(3)
        cur_ind = 1+ind(1)+(j-1)*ind(2);
        if cur_ind <= m
            cur_ind = m + cur_ind;
        else
            cur_ind = 2*m+1 -cur_ind;
        end
        u(cur_ind, :) = fread(file1(i), 2*NY, 'double');
    end
    
    for j = 1:ind(3)
        cur_ind = 1+ind(1)+(j-1)*ind(2);
        if cur_ind <= m
            cur_ind = m + cur_ind;
        else
            cur_ind = 2*m+1 -cur_ind;
        end
        v(cur_ind, :) = fread(file2(i), 2*NY, 'double');
    end

    for j = 1:ind(3)
        cur_ind = 1+ind(1)+(j-1)*ind(2);
        if cur_ind <= m
            cur_ind = m + cur_ind;
        else
            cur_ind = 2*m+1 -cur_ind;
        end
        vort(cur_ind, :) = fread(file3(i), 2*NY, 'double');
    end
    
end
    % draw u, v, vort
    figure(1)
    set(gcf,'color','w','position',[25,100,800,400])
    contour(X, Y , u);
    caxis([-10 40]);
    title(['velocity u, time = ', num2str((l-1)/60/60*dta), ' h']);
    colorbar;
    export_fig([route, '/u',num2str(1000+l),'.png'])
    
    figure(2)
    set(gcf,'color','w','position',[725,100,800,400])
    contour(X, Y , v);
    caxis([-10 10]);
    title(['velocity v, time = ', num2str((l-1)/60/60*dta), ' h']);
    colorbar;
    export_fig([route, '/v',num2str(1000+l),'.png'])
    
    figure(3)
    set(gcf,'color','w','position',[25,800,1600,400])
    pcolor(X, Y , vort);
    caxis([-4e-5 4e-5]);
    ylim([-90 90])
    colormap jet
    shading flat
    title(['vorticity, time = ', num2str((l-1)/60/60*dta), ' h']);
    colorbar;
    export_fig([route, '/vort',num2str(1000+l),'.png'])
    
    pause;
    

end

for i=1:nproc
fclose(file1(i));
fclose(file2(i));
fclose(file3(i));
end



