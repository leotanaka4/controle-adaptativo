%---------------------------------------------------------------------
%  COE-835  Controle adaptativo
%
%  Script para simular o exemplo 1
%                                                       Ramon R. Costa
%                                                       01/out/09, Rio
%---------------------------------------------------------------------
clear all;
clc;

%---------------------------------------------------------------------
disp('-------------------------------')
disp('Script para simular o exemplo 1')
disp(' ')
disp('Caso: Planta ............. n = 1')
disp('      Grau relativo ..... n* = 1')
disp('      Par�metros ........ np = 1')
disp(' ')
disp('Algoritmo: Standard MRAC')
disp(' ')
disp('-------------------------------')

%------------------------------------------------ Initialization -----
tfinal = 10;    %Simulation interval
st = 0.05;      %Sample time to workspace

s = tf('s');    %trick!

%--------------------------------------------------------- Plant -----
ap = 2;

P = 1/(s-ap)
P = ss(P);

%----------------------------------------------- Reference model -----
am = 1;

M = 1/(s+am)
M = ss(M);

%--------------------------------------------- Initial condition -----
yp0  = 0
x0   = yp0;

ym0  = 0;
xm0  = ym0;

%----------------------------------- Reference signal parameters -----
DC = 1   %Constant

As = 0   %Sine wave amplitude
ws = 10  %Frequency

%------------------------------------------------- Matching gain -----
thetas = -ap - am   %theta*

%----------------------------------------- Adaptation parameters -----
gamma1 = 2;       %Adaptation gains
gamma2 = 100;
theta0 = 0;       %Adaptation inicial condition

%---------------------------------------------------- Simulation -----
gamma = gamma1
sim('MRAC_111',tfinal);

yp1 = yp;   %Save results
e01 = e0;
theta1 = theta;
u1 = u;

%---------------------------------------------------- Simulation -----
gamma = gamma2
sim('MRAC_111',tfinal);

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
par1 = strcat('$e_0\;(\gamma=',num2str(gamma1),')$');
par2 = strcat('$e_0\;(\gamma=',num2str(gamma2),')$');
legend({par1,par2},...
    'FontSize',8,'Interpreter','latex','Location','NorthEast')
sublaby("   ");
print -depsc2 fig01a.eps

Thetas = thetas*ones(size(t));

figure(2)
clf
subplot(211)
plot(t,theta1,t,theta2,t,Thetas,'Linew',0.5);
grid on
title({'$\theta, \theta^*$'},'FontSize',10,'Interpreter','latex')
par1 = strcat('$\theta\;(\gamma=',num2str(gamma1),')$');
par2 = strcat('$\theta\;(\gamma=',num2str(gamma2),')$');
legend({par1,par2,'$\theta^*$'},...
    'FontSize',8,'Interpreter','latex','Location','NorthEast')
sublaby("   ");
print -depsc2 fig01b.eps

figure(3)
clf
subplot(211)
hold on
plot(t,yp1)
plot(t,yp2,t,r,t,ym,'Linew',0.5)
grid on
title({'$r, y_m, y_p$'},'FontSize',10,'Interpreter','latex')
par1 = strcat('$y\;(\gamma=',num2str(gamma1),')$');
par2 = strcat('$y\;(\gamma=',num2str(gamma2),')$');
legend({par1,par2,'$r$','$y_m$'},...
    'FontSize',8,'Interpreter','latex','Location','NorthEast')
sublaby("   ");
V = axis;
axis([V(1) V(2) 0 2.5 ]);
print -depsc2 fig01c.eps

ttheta1 = theta1 - thetas;
ttheta2 = theta2 - thetas;

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
par1 = strcat('$e_0 \times \tilde\theta\;(\gamma=',num2str(gamma1),')$');
par2 = strcat('$e_0 \times \tilde\theta\;(\gamma=',num2str(gamma2),')$');
legend({par1,par2},...
    'FontSize',8,'Interpreter','latex','Location','SouthEast')
%sublaby("   ");
print -depsc2 fig01d.eps

figure(5)
clf
subplot(211)
hold on
plot(t,u1)
plot(t,u2,'Linew',0.5)
grid on
title({'$u$'},'FontSize',10,'Interpreter','latex')
par1 = strcat('$u\;(\gamma=',num2str(gamma1),')$');
par2 = strcat('$u\;(\gamma=',num2str(gamma2),')$');
legend({par1,par2},...
    'FontSize',8,'Interpreter','latex','Location','SouthEast')
sublaby("   ");
print -depsc2 fig01e.eps

%------------------------------------------------- Display plots -----
figure(6)
clf

subplot(221)
hold on
plot(t,e01)
plot(t,e02,'Linew',0.5);
grid on
title({'$e_0$'},'FontSize',10,'Interpreter','latex')
par1 = strcat('$\gamma=',num2str(gamma1),'$');
par2 = strcat('$\gamma=',num2str(gamma2),'$');
legend({par1,par2},...
    'FontSize',10,'Interpreter','latex','Location','SouthEast')

subplot(222)
hold on
plot(t,theta1)
plot(t,theta2,t,Thetas,'r','Linew',0.5);
grid on; 
title({'$\theta, \theta^*$'},'FontSize',10,'Interpreter','latex')
legend({par1,par2,'$\theta^*$'},...
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

%--------------------------------------- Impress�o dos diagramas -----
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


