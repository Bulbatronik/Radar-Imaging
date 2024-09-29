% PART 1 The monostatic array
clear,close all,clc
%% parameters
f = 77e9; %77GHz
c = 3e8;
lambda = c/f;
rho_theta = deg2rad(3.6); %angular resolution
B = 2e9; %2GHz
T = 1/B; % Duration

rho_r = c/(2*B); % range resolution
d = lambda/4;% antialiasing F_spatial=1/d>=2*2/lambda

L = lambda/(2*rho_theta);% length
Ntx = ceil(L/d+1);

% FAR FIELD: R0>=2*L^2/lambda
% sensors' position
x_s = (0:Ntx-1).'*d;
pos_a = [x_s- mean(x_s) zeros(Ntx,1)];

% target position
xt = -10;
yt = 16;

figure;
plot(pos_a(:,1), pos_a(:,2), 'y-o');
grid on
hold on
disp(sqrt(xt^2+yt^2))
axis equal
plot(xt, yt,'r*')
legend('array', 'target')
%%
%data simulation
r_ax = 0:d:40;% axis of parameter r

data = zeros(length(r_ax), Ntx);% range-compress matrix
R_n = zeros(Ntx,1);

for ii = 1:Ntx% forward problem
    R_n(ii) = sqrt((pos_a(ii,1)-xt)^2+(pos_a(ii,2)-yt)^2);
    data(:, ii) = sinc((r_ax-R_n(ii))/rho_r)*exp(-1j*4*pi*R_n(ii)/lambda);   
end 

 close all
 figure;
 imagesc(-floor(Ntx/2):1:floor(Ntx/2),r_ax, abs(data))
 xlabel('Sensor')
 ylabel('Distance, [m]')
%% range estimation
[~, pos] = max(data); % peaks of the sincs
R_n_est = r_ax(pos);
%%
Nfft = 1024;

data_ft = zeros(length(r_ax), Nfft);
for j = 1:length(r_ax)
    data_ft(j,:) = fftshift(fft(data(j,:), Nfft));
end
%%
% horizontal axis
f_s = 1/d;
df = f_s/Nfft;
frq = (-(Nfft)/2:(Nfft)/2-1)*df;


close all
figure;
imagesc(frq, r_ax, abs(data_ft))
xlabel('frequency, [m^-1]')
ylabel('Distance, [m]')
%%
[~,I] = max(abs(data_ft(:)));
[~, I_col] = ind2sub(size(data_ft),I);

% frequency
f_est = frq(I_col);

theta_est = rad2deg(asin(f_est*lambda/2))
theta_true = rad2deg(atan((xt/yt)))
error = theta_true - theta_est
%% multipath
clc;
D = 5;% position of the wall on the x-axis

% pos of the image
xi = -xt+2*D;%%
yi = yt;

disp(sqrt(xt^2+yt^2))
disp(sqrt(xi^2+yi^2))

figure;
plot(pos_a(:,1), pos_a(:,2), 'y-o');
grid on
hold on
axis equal
plot(xt, yt,'r*')
xline(D)
plot(xi, yi,'g^')
legend('array', 'target', 'refr. surf.', 'image')
%%
%data simulation
r_ax = 0:d:40;% axis of parameter r

data_mult = zeros(length(r_ax), Ntx);% range-compress matrix
R_n_m = zeros(Ntx,1);

for ii = 1:Ntx% forward problem
    R_n_m(ii) = sqrt((pos_a(ii,1)-xt)^2+(pos_a(ii,2)-yt)^2)+...
                sqrt((pos_a(ii,1)-xi)^2+(pos_a(ii,2)-yi)^2);
    data_mult(:, ii) = sinc((r_ax-1/2*R_n_m(ii))/rho_r)*exp(-1j*2*pi/lambda*R_n_m(ii));
end 
%%
close all
figure;
imagesc(-floor(Ntx/2):1:floor(Ntx/2),r_ax,abs(data_mult))
xlabel('Sensor')
ylabel('Distance, [m]')
%% range estimation
[~, pos_mult] = max(data_mult); % peaks of the sincs
R_n__m_est = r_ax(pos_mult);
%%
Nfft = 1024;

data_mult_ft = zeros(length(r_ax), Nfft);
for j = 1:length(r_ax)
    data_mult_ft(j,:) = fftshift(fft(data_mult(j,:), Nfft));
end
%%
close all
figure;
imagesc(frq, r_ax,abs(data_mult_ft))
xlabel('frequency, [m^-1]')
ylabel('Distance, [m]')
%%
clc;
[M,I] = max(abs(data_mult_ft(:)));
[I_row, I_col] = ind2sub(size(data_mult_ft),I);

f_est = frq(I_col);

theta_est = rad2deg(asin(f_est*lambda/2))
theta_true_ghost = rad2deg(asin(1/2*(sin(atan(xt/yt))+sin(atan(xi/yi)))))
theta_true_image = rad2deg(atan(xi/yi))
error = theta_true_ghost - theta_est

close all
figure;
plot(pos_a(:,1), pos_a(:,2), 'y-o');
grid on
hold on
axis equal
plot(xt, yt,'r*')
xline(D)
plot(xi, yi,'g^')
plot(mean(R_n__m_est)*sin(deg2rad(theta_est)),...
     mean(R_n__m_est)*cos(deg2rad(theta_est)),'bx')
legend('array', 'target','refr. surf.', 'image', 'ghost')
title('Monostatic array')
%% Plot all on one figure (WARNING: IT TAKES A LOT OF TIME!)
close all
figure;
full = abs(data_mult_ft+data_ft);

angles = asin((-Nfft/2:Nfft/2-1).*df*lambda/2);
angles = rad2deg(angles);

%I1 = uimage(angles, r_ax, full)
uimage(angles, r_ax, full)
xlabel('Theta, [deg]')
ylabel('Distance, [m]')
hold on;
% Plot cross at target
plot(rad2deg(atan(xt/yt)),sqrt(xt^2+yt^2),'r+', 'MarkerSize', 10, 'LineWidth', .2);
% Plot cross at image
plot(rad2deg(asin(1/2*(sin(atan(xt/yt))+sin(atan(xi/yi))))),0.5*(sqrt(xt^2+yt^2)+sqrt(xi^2+yi^2)),'y+', 'MarkerSize', 10, 'LineWidth', .2);
legend('Target pred.','Image pred.')
%save('MONO.mat', 'I1')
%%
monost_pos = [mean(R_n__m_est)*sin(deg2rad(theta_est)), mean(R_n__m_est)*cos(deg2rad(theta_est))];
save('monostatic_pos.mat','monost_pos');
