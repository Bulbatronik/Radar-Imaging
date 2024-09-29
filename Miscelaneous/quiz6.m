lambda = 1;%blue
phiax = -pi/2:0.01:pi/2;
fx = sin(phiax)/lambda;
N = 5;
dx = 3;
apat = sin(pi*N*fx*dx)./sin(pi*fx*dx);
polarplot(phiax, abs(apat), 'b')
rlim([0 5])


hold on

lambda = 1.2;%green
phiax = -pi/2:0.01:pi/2;
fx = sin(phiax)/lambda;
N = 5;
dx = 3;
apat = sin(pi*N*fx*dx)./sin(pi*fx*dx);
polarplot(phiax, abs(apat), 'g')
rlim([0 5])


lambda = 1.4;%red
phiax = -pi/2:0.01:pi/2;
fx = sin(phiax)/lambda;
N = 5;
dx = 3;
apat = sin(pi*N*fx*dx)./sin(pi*fx*dx);
polarplot(phiax, abs(apat), 'r')
rlim([0 5])
legend('1', '1.2', '1.4')