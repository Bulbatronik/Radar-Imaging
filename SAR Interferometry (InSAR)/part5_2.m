clc; clear all; close all;
% SUBSIDIDENCE
h = 693E3; %m
B = 200; %m
f = 5.4E9; %Hz
c = 3e8;
lambda = c/f;

rho_r = 5;%m SLANT RANGE RESOLUTION
theta = deg2rad(35); %rad
sat = [0, h;
       B, h]; % satellites positions

% mountain
span = 1E3; % 1 km along the ground range
center = h*tan(theta); % center of the mountain to preserve the incidence angle
rho_g = rho_r/sin(theta);%ground range resolution
y = (center-span/2:rho_r/4:center+span/2)';
z1 = zeros(size(y));
z = gaussmf(y,[1000*span, center]);

% shift thee gaussian to zero and make the reasonable height
height = -0.03; %3 cm down
z2 = height*(z - min(z))/( max(z) - min(z) );% scale
figure
plot(y,z1,'.')
hold on
plot(y,z2,'.')
grid on
xlabel('y, [m]')
ylabel('z, [m]')
%axis equal
legend('1st acquisition', '2nd acquisition')

% complex reflectivity
Np = length(z);% number of points/scatters;
t = complex(rand(Np, 1), rand(Np,1));

p = [y,z1,z2,t];% points 
%% slant range axis
close all
dist_M = sqrt(sum((sat(1,:) - p(:,[1,2])).^2,2));
dist_S = sqrt(sum((sat(2,:)  - p(:,[1,3])).^2,2));

r_min = min([min(dist_M),min(dist_S)])*0.9999; 
r_max = max([max(dist_M),max(dist_S)])*1.0001;

range = r_min:rho_r:r_max;
%% acquire images
I = zeros(length(sat), length(range));

for n = 1:length(sat)%satelite/ image number
    R_n = sqrt(sum((sat(n,:) - p(:,[1,n+1])).^2,2));
    for r = 1:length(range)% range axis
        I(n,r) = sum(p(:,4).*sinc((range(r)-R_n)/rho_r).*exp(-1j*4*pi/lambda*R_n));
    end
end
%%
close all
figure 
subplot(2,1,1)
plot(range,abs(I(1,:)))
title('Master')
grid on
subplot(2,1,2)
plot(range,abs(I(2,:)))
title('Slave')
grid on
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
    close all
    figure 
    subplot(2,1,1)
    plot(y,abs(I_n_c_r(1,:)))
    title('Master')
    grid on
    subplot(2,1,2)
    plot(y,abs(I_n_c_r(2,:)))
    title('Slave')
    grid on
end
%% Interferometric phase 
if geo == true
    close all;
    %1
    % during th class: a(p)*exp(phi)
    % in the hw: a(p)*exp(-phi)
    % CONJUGATE THE ENTIRE IMAGE
    Image = I_n_c_r(1,:).*conj(I_n_c_r(2,:));
    figure;
    plot(y, angle(Image))% if GEO -> plot in y
    title('Original phase')
    
    % 2+3 compensate
    d_R_n_ref = R_n_ref(1,:)-R_n_ref(2,:);
    
    Image_ref = exp(-1j*4*pi/lambda.*d_R_n_ref);
    figure;
    plot(y, angle(Image_ref))
    title('Ref. phase')
    
    Image_flat = Image.*conj(Image_ref);
    phase_comp = angle(Image_flat);
    figure
    plot(y, phase_comp)
    title('compensated phase')
    
    phase_comp_filt = angle(movmean(exp(1j*(phase_comp)), 80));
    %phase_comp_filt = movmedian(phase_comp_filt, 50);
    figure
    plot(y, phase_comp_filt)
    title('filtered phase')
    
    %5
    %MANUALLY
    phase_unwr = phase_comp_filt;
    for i=2:length(phase_comp_filt)% make a decision
        difference = phase_comp_filt(i)-phase_comp_filt(i-1);
        if difference > pi
            phase_unwr(i:end) = phase_unwr(i:end) - 2*pi;
        elseif difference < -pi
            phase_unwr(i:end) = phase_unwr(i:end) + 2*pi;
        end
    end
    
    %phase_unwr = unwrap(phase_comp_filt);
    figure
    plot(y, phase_unwr)
    title('unwrapped phase')
    %6
    
    d_r = -phase_unwr.*lambda/(4*pi)/cos(theta);
    figure
    plot(y, d_r,'r')
    hold on
    plot(y,z2,'g')
    grid on
    legend('estimated','true')
end