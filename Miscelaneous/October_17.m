%MARCO.MANZONI@POLIMI.IT
clear;
clc;
%% 
dt  = 0.1;%PRF, 1/Fs
t = (-250:1:250)*dt;
N = length(t);

%% signals
x = my_rectpuls(t);% rectpuls

figure;
grid on
plot(t, x);
xlabel('Time, [s]')
ylabel('Amplitude')
%% shifting and stretching signals
Tx = 2;
tau_x = 1;

x = rectpuls((t-tau_x)/Tx);

Ty = 3;
tau_y = -1;
y = rectpuls((t-tau_y)/Ty);

%%
s = x+y;
z = x.*y;

figure;
subplot(4,1,1)
plot(t, x);
xlabel('Time, [s]')
ylabel('Amplitude')
title('x(t)')

subplot(4,1,2)
plot(t, y);
xlabel('Time, [s]')
ylabel('Amplitude')
title('y(t)')

subplot(4,1,3)
plot(t, s);
xlabel('Time, [s]')
ylabel('Amplitude')
title('s(t)')

subplot(4,1,4)
plot(t, z);
xlabel('Time, [s]')
ylabel('Amplitude')
title('z(t)')

%% coherent summation(interference)
y = y.*exp(1j*2*pi*3*t);%add a frequency

figure;
plot(t, abs(x+y));
xlabel('Time, [s]')
ylabel('Amplitude')
title('Interference')
% in-phae or out-phase -> distructive or constructive interference
% ALWAYS CONSIDER THAT SIGNALS ARE COMPLEX

%% 
clear;
close
clc;

dt = 1e-2;
t_obs = 2;

t = -t_obs/2:dt:t_obs/2;
Nt = length(t);

f0 = 0.1/dt;
x = exp(1j*2*pi*f0*t);

Nf = 4*Nt;%any
fs = 1/dt;%sampling freq

df = fs/Nf;
f = (-Nf/2:Nf/2-1)*df;
%f = f0;

W = exp(-1j*2*pi*f'*t);

X = W*x.';%only transpose
X = X*dt;

figure;
plot(f, abs(X))%, '*')
xlabel('Frequency')
ylabel('Amplitude')
title('abs(X)')

%% inverse
%df = f(2)-f(1);
y = df*W'*X;%Hermissian

%% INCOHERENT BACKPROJECTION
%Sn = rect(r-Rn)
clear;
close;
clc;

%sensors
N = 5;
dx = 1;

% positions
x_s = (0:N-1)*dx;
x_s = x_s - mean(x_s);
y_s = zeros(size(x_s));

%target
x_t = 0;
y_t = 3;%m

%signal
wid = 0.1;

%data simulation
r_ax = 0:wid/10:20;

data = zeros(length(r_ax), N);

for ii = 1:N% forward problem
    R = sqrt((x_s(ii)-x_t)^2+(y_s(ii)-y_t)^2);
    data(:, ii) = rectpuls((r_ax-R)/wid);%sinc also possible
end

figure;
imagesc(data)
xlabel('Sensor')
ylabel('Distance, [m]')

%% find isorange curves (backprojection)

% define a grid
x_grid = -3:wid/10:3;
y_grid = -0:wid/10:5;

% read the doc of meshgrid !!!!!!!!!!!!!!!!!!!!!!!!
[X,Y] = meshgrid(x_grid, y_grid);

I = zeros(size(X));

for ii = 1:N%backward problem
    R = sqrt((x_s(ii)-X).^2+(y_s(ii)-Y).^2);% NOT 2 NESTED LOOPS
    I = I+interp1(r_ax, data(:, ii), R, 'linear', NaN); % to draw curves on the grid
end

figure;
imagesc(x_grid, y_grid, I); axis xy;
xlabel('x')
ylabel('y')