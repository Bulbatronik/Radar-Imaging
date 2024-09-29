clc; clear all; close all;
% DECORRELATION (partially decorrelated)
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
z = gaussmf(y,[span/10, center]);

% shift thee gaussian to zero and make the reasonable height
height = 20; %200 m height (mounain)
z = height*(z - min(z))/( max(z) - min(z) );% scale
figure
plot(y,z,'.')
xlabel('y, [m]')
ylabel('z, [m]')
grid on

% complex reflectivity
Np = length(z);% number of points/scatters;

%Cholesky decomposition
T = rand(Np,2) + 1j*rand(Np,2);
rho = 1;%correlation 
R = [1 rho; rho 1]+10e-16*eye(2);
L = chol(R);
T = T*L;
corr(abs(T(:,1)),abs(T(:,2)))

tm = T(:,1);
ts = T(:,2);
p = [y,z,tm,ts];% points 
%% slant range axis
close all
dist_M = sqrt(sum((sat(1,:) - p(:,1:2)).^2,2));
dist_S = sqrt(sum((sat(2,:)  - p(:,1:2)).^2,2));

r_min = min([min(dist_M),min(dist_S)])*0.9999; 
r_max = max([max(dist_M),max(dist_S)])*1.0001;

range = r_min:rho_r:r_max;
%% acquire images
I = zeros(length(sat), length(range));

for n = 1:length(sat)%satelite/ image number
    R_n = sqrt(sum((sat(n,:) - p(:,1:2)).^2,2));
    for r = 1:length(range)% range axis
        I(n,r) = sum(p(:,n+2).*sinc((range(r)-R_n)/rho_r).*exp(-1j*4*pi/lambda*R_n));
    end
end
R_master = sqrt(sum((sat(1,:) - p(:,1:2)).^2,2))';

close all
figure 
subplot(2,1,1)
plot(range,abs(I(1,:)))
xlabel('r, [m]')
ylabel('Magnitude')
title('Master')
grid on
subplot(2,1,2)
plot(range,abs(I(2,:)))
xlabel('r, [m]')
ylabel('Magnitude')
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
    xlabel('y, [m]')
    ylabel('Magnitude')
    title('Master')
    grid on
    subplot(2,1,2)
    plot(y,abs(I_n_c_r(2,:)))
    xlabel('y, [m]')
    ylabel('Magnitude')
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
    xlabel('y, [m]')
    ylabel('Rad')
    
    % 2+3 compensate
    d_R_n_ref = R_n_ref(1,:)-R_n_ref(2,:);
    
    Image_ref = exp(-1j*4*pi/lambda.*d_R_n_ref);
    figure;
    plot(y, angle(Image_ref))
    title('Ref. phase')
    xlabel('y, [m]')
    ylabel('Rad')
    
    
    Image_flat = Image.*conj(Image_ref);
    phase_comp = angle(Image_flat);
    figure
    plot(y, phase_comp)
    title('compensated phase')
    xlabel('y, [m]')
    ylabel('Rad')
    
    
    %phase_comp_filt = movmedian(phase_comp, 20);% median?
    phase_comp_filt = angle(movmean(exp(1j*(phase_comp)), 80));
    %phase_comp_filt = movmedian(phase_comp_filt, 50);
    figure
    plot(y, phase_comp_filt)
    title('filtered phase')
    xlabel('y, [m]')
    ylabel('Rad')
    
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
    xlabel('y, [m]')
    ylabel('z, [m]')
    %6
    k_z = -4*pi/lambda*B*cot(theta)./R_master;
    z_est = phase_unwr./k_z;
    figure
    plot(y, z_est,'r')
    hold on
    plot(y,z,'g')
    grid on
    legend('estimated','true')
    xlabel('y, [m]')
    ylabel('z, [m]')
    
end
%% Slant range
if geo == false
    % 1
    y_ref = p(:,1);
    z_ref = zeros(length(p(:,2)),1);
    p_ref = [y_ref, z_ref];
    % 2
    R_m_ref= sqrt(sum((sat(1,:) - p_ref).^2,2))'; % 1   801 <- #points
    % 3 WHY NAN IF Z = 0?
    z_r = interp1(R_m_ref, p_ref(:,2)', range,'linear','extrap');
    z_r(isnan(z_r)) = 0;
    y_r = interp1(R_m_ref, p_ref(:,1)', range,'linear','extrap');
    y_r(isnan(y_r)) = 0;
    sum(~isnan(z_r))
    sum(~isnan(y_r))
    % 4
    R_n_ref(1,:) = sqrt((sat(1,1) - y_r).^2 + (sat(1,2) + z_r).^2);
    
    for n = 2:size(sat, 1)
        for r = 1:length(y_r)
            R_n_ref(n,r) = sqrt((sat(n,1) - y_r(r))^2 + (sat(n,2) + z_r(r))^2);
        end
        d_r(n-1,:) = R_n_ref(n,:) - R_n_ref(1,:);
    end
    
    I_n_c_r(1,:) = I(1,:);
    for i = 2:size(I,1)
        I_n_c_r(i,:) = interp1(range, I(i,:), range + d_r(i-1,:));% NaNs!!!!
    end
    I_n_c_r(isnan(I_n_c_r)) = 0;

    close all
    figure 
    subplot(2,1,1)
    plot(range,abs(I_n_c_r(1,:)))
    title('Master')
    axis([min(range) max(range) 0 8])
    grid on
    subplot(2,1,2)
    plot(range,abs(I_n_c_r(2,:)))% not 2, because skip master
    title('Slave')
    axis([min(range) max(range) 0 8])
    grid on
end
%%
if geo == false
    %1
    % during th class: a(p)*exp(phi)
    % in the hw: a(p)*exp(-phi)
    % CONJUGATE THE ENTIRE IMAGE
    Image = I_n_c_r(1,:).*conj(I_n_c_r(2,:));
    figure;
    plot(range, angle(Image))% if GEO -> plot in y
    title('Original phase')
    
    % 2+3 compensate
    d_R_n_ref = R_n_ref(1,:)-R_n_ref(2,:);
    
    Image_ref = exp(-1j*4*pi/lambda.*d_R_n_ref);
    figure;
    plot(range, angle(Image_ref))
    title('Ref. phase')
    
    Image_flat = Image.*conj(Image_ref);
    phase_comp = angle(Image_flat);
    figure
    plot(range, phase_comp)
    title('compensated phase')
    
    %phase_comp_filt = movmedian(phase_comp, 20);% median?
    phase_comp_filt = movmean(phase_comp, 120);
    %phase_comp_filt = movmedian(phase_comp_filt, 50);
    figure
    plot(range, phase_comp_filt)
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
    plot(range, phase_unwr)
    title('unwrapped phase')
    %6
    k_z = -4*pi/lambda*B*cot(theta)./R_n_ref(1,:);
    z_est = phase_unwr./k_z;
    figure
    plot(range, z_est,'r')
    hold on
    
    z_1 = interp1(R_m_ref, p(:,2)', range, 'linear','extrap');
    z_1(isnan(z_1)) = 0;
    y_1 = interp1(R_m_ref, p(:,1)', range, 'linear','extrap');
    y_1(isnan(y_1)) = 0;

    plot(y_1,z_1,'g')
    grid on
    legend('estimated','true')
end