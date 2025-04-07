%---------------------------------------------------------------------
%  COE-603  Controle adaptativo
%
%  Script para simular o exemplo 4
%                                             Leonardo S. da C. Tanaka 
%                                                       30/mar/25, Rio
%---------------------------------------------------------------------
clear all;
clc;

%---------------------------------------------------------------------
disp('-------------------------------')
disp('Script para simular o exemplo 1')
disp(' ')
disp('Caso: Planta ............. n = 1')
disp('      Grau relativo ..... n* = 1')
disp('      Parâmetros ........ np = 2')
disp(' ')
disp('Algoritmo: MRAC Direto')
disp(' ')
disp('-------------------------------')

%------------------------------------------------ Initialization -----
tfinal = 20;    %Simulation interval
st = 0.01;      %Sample time to workspace

s = tf('s');    %trick!

%--------------------------------------------------------- Plant -----
ap = 2;
kp = 1;

P = kp/(s-ap)
P = ss(P);

%----------------------------------------------- Reference model -----
am = 1;
km = 1;

M = 1/(s+am)
M = ss(M);

%--------------------------------------------- Initial condition -----
yp0  = 3
x0   = yp0;

ym0  = 0;
xm0  = ym0;

%----------------------------------- Reference signal parameters -----
DC = 2   %Constant

As = 1   %Sine wave amplitude
ws = 5  %Frequency

%------------------------------------------------- Matching gain -----
thetas = [-(ap+am)/kp; km/kp];   %theta*

%----------------------------------------- Adaptation parameters -----
gamma1 = 2*[1 0;0 1];       %Adaptation gains
gamma2 = 100*[1 0;0 1];
theta0 = [0;0];       %Adaptation inicial condition

%---------------------------------------------------- Simulation -----
gamma = gamma1
sim('MRAC_direto',tfinal);

yp1 = yp;   %Save results
e01 = e0;
theta1 = theta;
u1 = u;

%---------------------------------------------------- Simulation -----
gamma = gamma2
sim('MRAC_direto',tfinal);

yp2 = yp;   %Save results
e02 = e0;
theta2 = theta;
u2 = u;

%----------------------------------------------- Print eps plots -----
figure(1)
clf
subplot(211)
plot(t,e01,t,e02,'Linew',0.5);
grid on
title({'$e_0$'},'FontSize',10,'Interpreter','latex')
par1 = strcat('$e_0\;(\gamma=',strrep(mat2str(gamma1), ' ', '\ '),')$');
par2 = strcat('$e_0\;(\gamma=',strrep(mat2str(gamma2), ' ', '\ '),')$');
legend({par1,par2},...
    'FontSize',8,'Interpreter','latex','Location','NorthEast')
sublaby("   ");
print -dpng images\fig04a.png

Theta1 = thetas(1)*ones(size(t));
Theta2 = thetas(2)*ones(size(t));

figure(2)
clf
subplot(211)
plot(t,theta1,t,theta2,t,Theta1,t,Theta2,'Linew',0.5);
grid on
title({'$\theta, \theta^*$'},'FontSize',10,'Interpreter','latex')
par1 = strcat('$\theta_{1}\;(\gamma=',strrep(mat2str(gamma1), ' ', '\ '),')$');
par2 = strcat('$\theta_{2}\;(\gamma=',strrep(mat2str(gamma1), ' ', '\ '),')$');
par3 = strcat('$\theta_{1}\;(\gamma=',strrep(mat2str(gamma2), ' ', '\ '),')$');
par4 = strcat('$\theta_{2}\;(\gamma=',strrep(mat2str(gamma2), ' ', '\ '),')$');
legend({par1,par2,par3,par4,'$\theta_1^*$','$\theta_2^*$'},...
    'FontSize',8,'Interpreter','latex','Location','NorthEast')
sublaby("   ");
print -dpng images\fig04b.png

figure(3)
clf
subplot(211)
hold on
plot(t,yp1)
plot(t,yp2,t,r,t,ym,'Linew',0.5)
grid on
title({'$r, y_m, y_p$'},'FontSize',10,'Interpreter','latex')
par1 = strcat('$y\;(\gamma=',strrep(mat2str(gamma1), ' ', '\ '),')$');
par2 = strcat('$y\;(\gamma=',strrep(mat2str(gamma2), ' ', '\ '),')$');
legend({par1,par2,'$r$','$y_m$'},...
    'FontSize',8,'Interpreter','latex','Location','NorthEast')
sublaby("   ");
V = axis;
axis([V(1) V(2) 0 2.5 ]);
print -dpng images\fig03c.png

dims = size(t);
thetas_matrix = ones([dims(1),2])*[thetas(1) 0;0 thetas(2)];

ttheta1 = theta1 - thetas_matrix;
ttheta2 = theta2 - thetas_matrix;

figure(4)
clf
hold on
plot(e01,ttheta1)
plot(e02,ttheta2)
grid on
%axis equal
title({'$e_0 \times \tilde\theta$'},'FontSize',10,'Interpreter','latex')
xlabel({'$e_0$'},'FontSize',10,'Interpreter','latex')
ylabel({'$\tilde\theta$'},'FontSize',10,'Interpreter','latex')
par1 = strcat('$e_0 \times \tilde\theta_1\;(\gamma=',strrep(mat2str(gamma1), ' ', '\ '),')$');
par2 = strcat('$e_0 \times \tilde\theta_2\;(\gamma=',strrep(mat2str(gamma1), ' ', '\ '),')$');
par3 = strcat('$e_0 \times \tilde\theta_1\;(\gamma=',strrep(mat2str(gamma2), ' ', '\ '),')$');
par4 = strcat('$e_0 \times \tilde\theta_2\;(\gamma=',strrep(mat2str(gamma2), ' ', '\ '),')$');
legend({par1,par2,par3,par4},...
    'FontSize',8,'Interpreter','latex','Location','SouthEast')
sublaby("   ");
print -dpng images\fig04d.png

figure(5)
clf
subplot(211)
hold on
plot(t,u1)
plot(t,u2,'Linew',0.5)
grid on
title({'$u$'},'FontSize',10,'Interpreter','latex')
par1 = strcat('$u\;(\gamma=',strrep(mat2str(gamma1), ' ', '\ '),')$');
par2 = strcat('$u\;(\gamma=',strrep(mat2str(gamma2), ' ', '\ '),')$');
legend({par1,par2},...
    'FontSize',8,'Interpreter','latex','Location','SouthEast')
sublaby("   ");
print -dpng images\fig04e.png

%------------------------------------------------- Display plots -----
figure(6)
clf

subplot(221)
hold on
plot(t,e01)
plot(t,e02,'Linew',0.5);
grid on
title({'$e_0$'},'FontSize',10,'Interpreter','latex')
par1 = strcat('$\gamma=',strrep(mat2str(gamma1), ' ', '\ '),'$');
par2 = strcat('$\gamma=',strrep(mat2str(gamma2), ' ', '\ '),'$');
legend({par1,par2},...
    'FontSize',10,'Interpreter','latex','Location','SouthEast')

subplot(222)
hold on
plot(t,theta1,t,Theta1)
plot(t,theta2,t,Theta2,'r','Linew',0.5);
grid on; 
title({'$\theta, \theta^*$'},'FontSize',10,'Interpreter','latex')
par1c = strcat('$\theta_{1}\;(\gamma=',strrep(mat2str(gamma1), ' ', '\ '),')$');
par2c = strcat('$\theta_{2}\;(\gamma=',strrep(mat2str(gamma1), ' ', '\ '),')$');
par3c = strcat('$\theta_{1}\;(\gamma=',strrep(mat2str(gamma2), ' ', '\ '),')$');
par4c = strcat('$\theta_{2}\;(\gamma=',strrep(mat2str(gamma2), ' ', '\ '),')$');
legend({par1c,par2c,'$\theta_1^*$',par3c,par4c,'$\theta_1^*$'},...
    'FontSize',10,'Interpreter','latex','Location','NorthEast')

subplot(223)
hold on
plot(t,yp1);
plot(t,yp2,t,r,t,ym,'Linew',0.5);
grid on
title({'$r, y_m, y_p$'},'FontSize',10,'Interpreter','latex')
legend({par1,par2,'$r$','$y_m$'},...
    'FontSize',10,'Interpreter','latex','Location','SouthEast')

subplot(224)
hold on
plot(t,u1)
plot(t,u2,'Linew',0.5);grid;
grid on
title({'$u$'},'FontSize',10,'Interpreter','latex')
legend({par1,par2},...
    'FontSize',10,'Interpreter','latex','Location','SouthEast')

%--------------------------------------- Impressão dos diagramas -----
% open_system('MRAC_111');
% print -depsc2 -sMRAC_111 MRAC-111.eps
% 
% open_system('MRAC_111/Plant');
% print -depsc2 -sMRAC_111/Plant plant.eps
% 
% open_system('MRAC_111/Reference model');
% print -depsc2 '-sMRAC_111/Reference model' reference-model.eps
% 
% open_system('MRAC_111/Adaptation');
% print -depsc2 -sMRAC_111/Adaptation adaptation.eps
% 
% open_system('MRAC_111/Reference signal');
% print -depsc2 '-sMRAC_111/Reference signal' reference-signal.eps
% 
% close_system('MRAC_111');
%---------------------------------------------------------------------


