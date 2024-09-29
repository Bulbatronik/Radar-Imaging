% PART 3 The MIMO radar
clear,close all,clc
%% parameters
f = 77e9; %77GHz
c = 3e8;
lambda = c/f;
rho_theta = deg2rad(3.6); %angular resolution
B = 2e9; %2GHz
T = 1/B; % Duration

rho_r = c/(2*B); % range resolution

% RX
d_RX = lambda/2;% preserved
N_RX = 8;

% TX
d_TX = N_RX*(lambda/2);
N_TX = 4;

% virtual channels
d = lambda/4;
N = N_TX*N_RX;

% FAR FIELD: R0>=2*L^2/lambda
% sensors' position
x_TX = (0:N_TX-1).'*d_TX;
pos_TX = [x_TX- mean(x_TX) zeros(N_TX,1)];

x_RX = (0:N_RX-1).'*d_RX;
pos_RX = [x_RX- mean(x_RX) zeros(N_RX,1)];

% target position
xt = -10;
yt = 16;

figure;
plot(pos_RX(:,1), pos_RX(:,2), 'g-o');
grid on
hold on
plot(pos_TX(:,1), pos_TX(:,2), 'r-o');
axis equal
plot(xt, yt,'b*')
legend('RX','TX', 'target')
%%
%data simulation
r_ax = 0:d:40;% axis of parameter r

data_mimo = zeros(length(r_ax), N_RX, N_TX);% range-compress matrix
R_m = zeros(N_TX,1);
R_n = zeros(N_RX,1);

for jj = 1:N_TX% forward problem
    R_m(jj) = sqrt((pos_TX(jj,1)-xt)^2+(pos_TX(jj,2)-yt)^2);
    for ii = 1:N_RX
        R_n(ii) = sqrt((pos_RX(ii,1)-xt)^2+(pos_RX(ii,2)-yt)^2);

        data_mimo(:, ii,jj) = sinc((r_ax-0.5*(R_n(ii)+R_m(jj)))/rho_r)...
                        *exp(-1j*2*pi*(R_n(ii)+R_m(jj))/lambda);   
    end
end 
%% 
% Data was collected for the MIMO system. Now data_mimo matrix will be
% converted to an data_virt m x n -> m*n virt channels

data_virt = data_mimo(:,:,1);
for j = 2:N_TX
    data_virt = [data_virt , data_mimo(:,:,j)];%length(r_ax) x N_RX*N_TX
end

close all
figure;
imagesc(-floor(N/2):1:floor(N/2),r_ax,abs(data_virt))
xlabel('Virtual sensors')
ylabel('Distance, [m]')
%% range estimation
%R_n_m_est = zeros
[~, pos] = max(data_virt); % peaks of the sincs
R_n_m_est = r_ax(pos);
%%
Nfft = 1024;
data_ft = zeros(length(r_ax), Nfft);

for j = 1:length(r_ax)
    data_ft(j,:) = fftshift(fft(data_virt(j,:), Nfft));
end
%%
% horizontal axis
f_s = 1/d;
df = f_s/Nfft;
frq = (-(Nfft)/2:(Nfft)/2-1)*df;

close all
figure;
imagesc(frq, r_ax,abs(data_ft))
xlabel('frequency, [m^-1]')
ylabel('Distance, [m]')
%%
[~,I] = max(abs(data_ft(:)));
[~, I_col] = ind2sub(size(data_ft),I);

% fraquency
f_est = frq(I_col);

theta_est = rad2deg(asin(f_est*lambda/2))
theta_true = rad2deg(atan((xt/yt)))
error = theta_true - theta_est
%% multipath
clc;
D = 5;% position of the wall on the x-axis

% pos of the image
xi = -(xt-2*D);
yi = yt;

disp(sqrt(xt^2+yt^2))
disp(sqrt(xi^2+yi^2))

close all
figure;
plot(pos_RX(:,1), pos_RX(:,2), 'g-o');
grid on
hold on
plot(pos_TX(:,1), pos_TX(:,2), 'r-o');
axis equal
plot(xt, yt,'b*')
xline(D)
plot(xi, yi,'c^')
legend('RX','TX', 'target', 'refr. surf.', 'image')
%%
%data simulation
r_ax = 0:d:40;% axis of parameter r

data_mimo_mult = zeros(length(r_ax), N_RX, N_TX);% range-compress matrix
R_m = zeros(N_TX,1);
R_n = zeros(N_RX,1);

for jj = 1:N_TX% forward problem
    R_m(jj) = sqrt((pos_TX(jj,1)-xt)^2+(pos_TX(jj,2)-yt)^2);
    for ii = 1:N_RX
        R_n(ii) = sqrt((pos_RX(ii,1)-xi)^2+(pos_RX(ii,2)-yi)^2);

        data_mimo_mult(:, ii,jj) = sinc((r_ax-0.5*(R_n(ii)+R_m(jj)))/rho_r)...
                        *exp(-1j*2*pi*(R_n(ii)+R_m(jj))/lambda);   
    end
end 
%% 
% Data was collected for the MIMO system. Now data_mimo matrix will be
% converted to an data_virt m x n -> m*n virt channels

data_virt_mult = data_mimo_mult(:,:,1);
for j = 2:N_TX
    data_virt_mult = [data_virt_mult , data_mimo_mult(:,:,j)];%length(r_ax) x N_RX*N_TX
end

close all
figure;
imagesc(-floor(N/2):1:floor(N/2),r_ax,abs(data_virt_mult))
xlabel('Sensor')
ylabel('Distance, [m]')
%% range estimation
[~, pos_mult] = max(data_virt_mult); % peaks of the sincs
R_n__m_est = r_ax(pos_mult);
%%
Nfft = 1024;

data_mult_ft = zeros(length(r_ax), Nfft);
for j = 1:length(r_ax)
    data_mult_ft(j,:) = fftshift(fft(data_virt_mult(j,:), Nfft));
end
%%
close all
figure;
imagesc(frq, r_ax, abs(data_mult_ft))
xlabel('frequency, [m^-1]')
ylabel('Distance, [m]')
%%
clc;
[M,I] = max(abs(data_mult_ft(:)));
[I_row, I_col] = ind2sub(size(data_mult_ft),I);

f_est = frq(I_col);

theta_est = rad2deg(asin(f_est*lambda/2))%!!!!!

%what you would expect from monostatic
theta_true_ghost_m = rad2deg(asin(1/2*(sin(atan(xt/yt))+sin(atan(xi/yi)))))

%what you would expect from bistatic
theta_true_ghost_b = rad2deg(asin(sin(atan(xi/yi))))

theta_true_image = rad2deg(atan(xi/yi))
error_m = theta_true_ghost_m - theta_est
error_b = theta_true_ghost_b - theta_est

close all
figure;
plot(pos_RX(:,1), pos_RX(:,2), 'g-o');
grid on
hold on
plot(pos_TX(:,1), pos_TX(:,2), 'r-o');
axis equal
plot(xt, yt,'b*')
xline(D)
plot(xi, yi,'c^')
plot(mean(R_n__m_est)*sin(deg2rad(theta_est)),...
     mean(R_n__m_est)*cos(deg2rad(theta_est)),'yx')
legend('RX','TX', 'target', 'refr. surf.', 'image', 'ghost')
title('MIMO array')
%% Plot all on one figure (WARNING: IT TAKES A LOT OF TIME!)
close all
figure;
full = abs(data_mult_ft+data_ft);

angles = asin((-Nfft/2:Nfft/2-1).*df*lambda/2);
angles = rad2deg(angles);

%I3 = uimage(angles, r_ax, full)
uimage(angles, r_ax, full)
xlabel('Theta, [deg]')
ylabel('Distance, [m]')
hold on;
% Plot cross at true
plot(rad2deg(atan(xt/yt)),sqrt(xt^2+yt^2),'r+', 'MarkerSize', 10, 'LineWidth', .2);
% closed form solution is comlicated ti find
legend('Target pred.')
%save('MIMO.mat', 'I3')
%%
MIMO_pos = [mean(R_n__m_est)*sin(deg2rad(theta_est)),mean(R_n__m_est)*cos(deg2rad(theta_est))];
save('mimo_pos.mat','MIMO_pos');