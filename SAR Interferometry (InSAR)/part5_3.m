clc; clear all; close all;
% DECORRELATION (partially decorrelated)
h = 693E3; %m
B = 200; %m
f = 5.4E9; %Hz
c = 3e8;
lambda = c/f;

rho_r = 5;%m SLANT RANGE RESOLUTION
theta = deg2rad(35); %rad

N = 21; %Images - N-1
sat = [[0:B:(N-1)*B]',h*ones(N,1)];% satellites positions

% mountain
dz = 0.005; % sinks 5mm per acquisition

span = 1E3; % 1 km along the ground range
center = h*tan(theta); % center of the mountain to preserve the incidence angle
rho_g = rho_r/sin(theta);%ground range resolution
y = (center-span/2:rho_r/4:center+span/2)';
z = gaussmf(y,[span/10, center]);

Bt = 1;
% shift thee gaussian to zero and make the reasonable height
height = 20; %20 m height 

z = height*(z - min(z))/( max(z) - min(z) );% scale
for i = 1:N
    Z(:,i) = z - dz*(i-1);
end
figure
plot(y,Z,'-')
grid on

% complex reflectivity
Np = length(z);% number of points/scatters;

t = complex(rand(Np, 1), rand(Np,1));

p = [y,Z,t];% points 
%% slant range axis
close all
r_min = inf;
r_max = 0;
for i = 1:N-1
    dist_m = sqrt(sum((sat(1,:) - p(:,[1,i+1])).^2,2));
    dist_s = sqrt(sum((sat(i+1,:)  - p(:,[1,i+1])).^2,2));
    
    mx = max([max(dist_m),max(dist_s)]);
    if mx>r_max
       r_max = mx;
    end
    
    mn = min([min(dist_m),min(dist_s)]);
    if mn<r_min
       r_min = mn;
    end
end


r_min = r_min*0.9999; 
r_max = r_max*1.0001;

range = r_min:rho_r:r_max;
%% acquire images
I = zeros(length(sat), length(range));

for n = 1:N%satelite/ image number
    R_n = sqrt(sum((sat(n,:) - p(:,[1,N+1])).^2,2));
    for r = 1:length(range)% range axis
        I(n,r) = sum(p(:,end).*sinc((range(r)-R_n)/rho_r).*exp(-1j*4*pi/lambda*R_n));
    end
end
R_master = sqrt(sum((sat(1,:) - p(:,1:2)).^2,2))';
%%
close all
figure 

for i = 1:N
    subplot(N,1,i)
    plot(range,abs(I(i,:)))
    grid on
    title(sprintf('Sat %0.1i', i))
end

%% coregistration
geo = true; % false doesn't work!
% GEO, works
%1
if geo == true
    y_ref = p(:,1);
    z_ref = zeros(length(p(:,2)),1);
    p_ref = [y_ref, z_ref];
    %2
    for n = 1:size(sat, 1)
        R_n_ref(n,:)= sqrt(sum((sat(n,:) - p_ref).^2,2));
    end
    
    for n = 1:size(I,1)
        I_n_c_r(n,:) = interp1(range, I(n,:), R_n_ref(n,:));
    end
    
    % plot in ground range

    figure 

    for i = 1:N
        subplot(N,1,i)
        plot(y,abs(I_n_c_r(i,:)))
        grid on
        title(sprintf('Sat %0.1i', i))
    end
end
%% Interferometric phase 
if geo == true
    close all;
    %1
    for i = 2:N
        phi(:,i-1) = unwrap(angle(I_n_c_r(1,:).*conj(I_n_c_r(i,:))));

    end

    for j = 1:Np
%         A1 = ones(N-1,1);
%         A2 = transpose(-4*pi/lambda*Bt*(1:N-1));
%         A3 =transpose(-4*pi/(lambda*sin(theta)*R_master(j))*B*(1:N-1)*cos(theta));
%         A = [A1, A2, A3];
%         param(j,:) = pinv(A)*phi(j,:)';
        A1 = transpose(-4*pi/lambda*Bt*(1:N-1));
        A2 =transpose(-4*pi/(lambda*sin(theta)*R_master(j))*B*(1:N-1)*cos(theta));
        A = [A1, A2];
        param(j,:) = pinv(A)*phi(j,:)';
    end

    mean(param(:,1))%mean v_p
    mean(param(:,2))%mean q <-wrong
    cond(A)
end
