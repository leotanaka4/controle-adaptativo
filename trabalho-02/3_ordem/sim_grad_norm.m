%---------------------------------------------------------------------
%  COE-603  Controle adaptativo
%
%  Script para simular o exemplo 1
%                                             Leonardo S. da C. Tanaka 
%                                                       30/mar/25, Rio
%---------------------------------------------------------------------
clear all;
clc;

%---------------------------------------------------------------------
disp('-------------------------------')
disp('Script para simular o exemplo 1')
disp(' ')
disp('Caso: Planta ............. n = 3')
disp('      Grau relativo ..... n* = 1')
disp('      Parâmetros ........ np = 6')
disp(' ')
disp('Algoritmo: Gradiente Normalizado')
disp(' ')
disp('-------------------------------')

%------------------------------------------------ Initialization -----
tfinal = 100;    %Simulation interval
st = 0.01;      %Sample time to workspace

s = tf('s');    %trick!
%--------------------------------------------------------- Filter ----
lambda2 = 3;
lambda1 = 2;
lambda0 = 1;
%--------------------------------------------------------- Plant -----
a2 = 2;
a1 = 2.5;
a0 = 1.5;
b2 = 3;
b1 = 2;
b0 = 0.8;

P = (b2*s^2+b1*s+b0)/(s^3+a2*s^2+a1*s+a0)
P = ss(P);

%--------------------------------------------- Initial condition -----
yp0  = [0; 0; 0]
x0   = yp0;

%----------------------------------- Reference signal parameters -----
DC = 1   %Constant
As = 2   %Sine wave amplitude
ws = 0.1*pi  %Frequency

%------------------------------------------------- Matching gain -----
theta_star1 = [b0; b1; b2; a0 - lambda0; a1 - lambda1; a2 - lambda2]   %theta* parametrização 1
theta_star2 = [b0; b1; b2; a0; a1;a2]                                  %theta* parametrização 2
theta_star3 = [1/b0; b1/b0; b2/b0; a0/b0; a1/b0; a2/b0]                %theta* parametrização 3
%----------------------------------------- Adaptation parameters -----
gamma = 1*eye(6);             %Adaptation gains
theta0 = [0;0;0;0;0;0];       %Adaptation inicial condition
kappa = 1e-9;
P0 = eye(6);
%---------------------------------------------------- Simulation -----
theta_star = theta_star1;
sim('gradiente_normalizado',tfinal);

yp1 = yp;   %Save results
e01 = e0;
theta1 = theta;
u1 = u;
%---------------------------------------------------- Simulation 2 ---
theta_star = theta_star2;
sim('gradiente_normalizado',tfinal);

yp2 = yp;   %Save results
e02 = e0;
theta2 = theta;
u2 = u;
%---------------------------------------------------- Simulation 3 ---
theta_star = theta_star3;
sim('gradiente_normalizado_p3',tfinal);

yp3 = yp;   %Save results
e03 = e0;
theta3 = theta;
u3 = u;
%---------------------------------------------------- Simulation 4---
theta_star = theta_star1;
sim('rascunho_LS',tfinal);

yp4 = yp;   %Save results
e04 = e0;
theta4 = theta;
u4 = u;

%---------------------------------------------------- Simulation 5---
theta_star = theta_star2;
sim('rascunho_LS',tfinal);

yp5 = yp;   %Save results
e05 = e0;
theta5 = theta;
u5 = u;

%---------------------------------------------------- Simulation 6---
theta_star = theta_star3;
sim('rascunho_LS',tfinal);

yp6 = yp;   %Save results
e06 = e0;
theta6 = theta;
u6 = u;


%----------------------------------------------- Print eps plots -----
figure(1)
clf

% 1º gráfico: saída yp e referência r
subplot(3,2,1)
hold on
plot(t, yp1, 'LineWidth', 1.2, 'DisplayName', '$y$')
plot(t, r, 'LineWidth', 1.2, 'DisplayName', '$r$')
grid on
title('$y$ e $r$', 'Interpreter','latex')
legend('Interpreter','latex','Location','Best')

% 2º gráfico: erros das 3 simulações
subplot(3,2,3)
hold on
plot(t, e01, 'LineWidth', 1.2, 'DisplayName', '$\varepsilon_1$')
plot(t, e02, 'LineWidth', 1.2, 'DisplayName', '$\varepsilon_2$')
plot(t, e03, 'LineWidth', 1.2, 'DisplayName', '$\varepsilon_3$')
grid on
title('Erros $\varepsilon_i$', 'Interpreter','latex')
legend('Interpreter','latex','Location','Best')

% 3º gráfico: theta1 (4 curvas), com theta_star1 como linha pontilhada
subplot(3,2,5)
hold on
colors = lines(4);
for i = 1:4
    plot(t, theta1(:,i), 'Color', colors(i,:), 'LineWidth', 1.2, ...
        'DisplayName', strcat('$\theta_', num2str(i), '$'))
    h = plot(t, theta_star1(i)*ones(size(t)), '--', 'Color', colors(i,:), 'LineWidth', 1.2);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
grid on
title('Parametrização 1', 'Interpreter','latex')
legend('Interpreter','latex','Location','Best')

% 4º gráfico: theta2 (4 curvas), com theta_star2 como linha pontilhada
subplot(3,2,2)
hold on
for i = 1:4
    plot(t, theta2(:,i), 'Color', colors(i,:), 'LineWidth', 1.2, ...
        'DisplayName', strcat('$\theta_', num2str(i), '$'))
    h = plot(t, theta_star2(i)*ones(size(t)), '--', 'Color', colors(i,:), 'LineWidth', 1.2);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
grid on
title('Parametrização 2', 'Interpreter','latex')
legend('Interpreter','latex','Location','Best')

% 5º gráfico: theta3 (4 curvas), com theta_star3 como linha pontilhada
subplot(3,2,4)
hold on
for i = 1:4
    plot(t, theta3(:,i), 'Color', colors(i,:), 'LineWidth', 1.2, ...
        'DisplayName', strcat('$\theta_', num2str(i), '$'))
    h = plot(t, theta_star3(i)*ones(size(t)), '--', 'Color', colors(i,:), 'LineWidth', 1.2);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
grid on
title('Parametrização 3', 'Interpreter','latex')
legend('Interpreter','latex','Location','Best')

% 6º gráfico: controle u1
subplot(3,2,6)
hold on

% Normas dos thetas
norm_theta1 = vecnorm(theta1')';  % ||theta_1(t)||
norm_theta2 = vecnorm(theta2')';  % ||theta_2(t)||
norm_theta3 = vecnorm(theta3')';  % ||theta_3(t)||

% Normas dos theta_star (constantes)
norm_theta_star1 = norm(theta_star1);
norm_theta_star2 = norm(theta_star2);
norm_theta_star3 = norm(theta_star3);

% Usar cores padrão do MATLAB
h1 = plot(t, norm_theta1, 'LineWidth', 1.2, 'DisplayName', '$\|\theta_1\|$');
h2 = plot(t, norm_theta2, 'LineWidth', 1.2, 'DisplayName', '$\|\theta_2\|$');
h3 = plot(t, norm_theta3, 'LineWidth', 1.2, 'DisplayName', '$\|\theta_3\|$');

% Obter as cores automaticamente usadas
c1 = h1.Color;
c2 = h2.Color;
c3 = h3.Color;

% Plot theta*_i normas (pontilhado, mesma cor, sem legenda)
h1s = plot(t, norm_theta_star1 * ones(size(t)), '--', 'Color', c1, 'LineWidth', 1.2);
h1s.Annotation.LegendInformation.IconDisplayStyle = 'off';

h2s = plot(t, norm_theta_star2 * ones(size(t)), '--', 'Color', c2, 'LineWidth', 1.2);
h2s.Annotation.LegendInformation.IconDisplayStyle = 'off';

h3s = plot(t, norm_theta_star3 * ones(size(t)), '--', 'Color', c3, 'LineWidth', 1.2);
h3s.Annotation.LegendInformation.IconDisplayStyle = 'off';

grid on
title('Parametrizações 1, 2 e 3', 'Interpreter','latex')
legend('Interpreter','latex','Location','Best')

%grafico do LS - checar?

figure(2)
clf
% 1º gráfico: saída yp e referência r
subplot(3,2,1)
hold on
plot(t, yp4, 'LineWidth', 1.2, 'DisplayName', '$y$')
plot(t, r, 'LineWidth', 1.2, 'DisplayName', '$r$')
grid on
title('$y$ e $r$', 'Interpreter','latex')
legend('Interpreter','latex','Location','Best')

% 2º gráfico: erros das 3 simulações
subplot(3,2,3)
hold on
plot(t, e04, 'LineWidth', 1.2, 'DisplayName', '$\varepsilon_1$')
plot(t, e05, 'LineWidth', 1.2, 'DisplayName', '$\varepsilon_2$')
plot(t, e06, 'LineWidth', 1.2, 'DisplayName', '$\varepsilon_3$')
grid on
title('Erros $\varepsilon_i$', 'Interpreter','latex')
legend('Interpreter','latex','Location','Best')

% 3º gráfico: theta4 (4 curvas), com theta_star1 como linha pontilhada
subplot(3,2,5)
hold on
colors = lines(4);
for i = 1:4
    plot(t, theta4(:,i), 'Color', colors(i,:), 'LineWidth', 1.2, ...
        'DisplayName', strcat('$\theta_', num2str(i), '$'))
    h = plot(t, theta_star1(i)*ones(size(t)), '--', 'Color', colors(i,:), 'LineWidth', 1.2);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
grid on
title('Parametrização 1 (LS)', 'Interpreter','latex')
legend('Interpreter','latex','Location','Best')

% 4º gráfico: theta5 (4 curvas), com theta_star2 como linha pontilhada
subplot(3,2,2)
hold on
for i = 1:4
    plot(t, theta5(:,i), 'Color', colors(i,:), 'LineWidth', 1.2, ...
        'DisplayName', strcat('$\theta_', num2str(i), '$'))
    h = plot(t, theta_star2(i)*ones(size(t)), '--', 'Color', colors(i,:), 'LineWidth', 1.2);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
grid on
title('Parametrização 2 (LS)', 'Interpreter','latex')
legend('Interpreter','latex','Location','Best')

% 5º gráfico: theta6 (4 curvas), com theta_star3 como linha pontilhada
subplot(3,2,4)
hold on
for i = 1:4
    plot(t, theta6(:,i), 'Color', colors(i,:), 'LineWidth', 1.2, ...
        'DisplayName', strcat('$\theta_', num2str(i), '$'))
    h = plot(t, theta_star3(i)*ones(size(t)), '--', 'Color', colors(i,:), 'LineWidth', 1.2);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
grid on
title('Parametrização 3 (LS)', 'Interpreter','latex')
legend('Interpreter','latex','Location','Best')

% 6º gráfico: normas dos thetas
subplot(3,2,6)
hold on

norm_theta4 = vecnorm(theta4')';
norm_theta5 = vecnorm(theta5')';
norm_theta6 = vecnorm(theta6')';

norm_theta_star1 = norm(theta_star1);
norm_theta_star2 = norm(theta_star2);
norm_theta_star3 = norm(theta_star3);

h1 = plot(t, norm_theta4, 'LineWidth', 1.2, 'DisplayName', '$\|\theta_1\|$');
h2 = plot(t, norm_theta5, 'LineWidth', 1.2, 'DisplayName', '$\|\theta_2\|$');
h3 = plot(t, norm_theta6, 'LineWidth', 1.2, 'DisplayName', '$\|\theta_3\|$');

c1 = h1.Color;
c2 = h2.Color;
c3 = h3.Color;

plot(t, norm_theta_star1 * ones(size(t)), '--', 'Color', c1, 'LineWidth', 1.2, ...
    'HandleVisibility','off');
plot(t, norm_theta_star2 * ones(size(t)), '--', 'Color', c2, 'LineWidth', 1.2, ...
    'HandleVisibility','off');
plot(t, norm_theta_star3 * ones(size(t)), '--', 'Color', c3, 'LineWidth', 1.2, ...
    'HandleVisibility','off');

grid on
title('Parametrizações 1, 2 e 3 (LS)', 'Interpreter','latex')
legend('Interpreter','latex','Location','Best')


%{
figure(2)
clf
subplot(211)
plot(t, theta1, t, Theta1, t, Theta2, t, Theta3, t, Theta4, 'LineWidth', 0.5);
grid on
title({'$\theta, \theta^*$'}, 'FontSize', 10, 'Interpreter', 'latex')

par1 = strcat('$\theta_{1}\;(\gamma=', strrep(mat2str(gamma1), ' ', '\ '), ')$');
par2 = strcat('$\theta_{2}\;(\gamma=', strrep(mat2str(gamma1), ' ', '\ '), ')$');
par3 = strcat('$\theta_{3}\;(\gamma=', strrep(mat2str(gamma1), ' ', '\ '), ')$');
par4 = strcat('$\theta_{4}\;(\gamma=', strrep(mat2str(gamma1), ' ', '\ '), ')$');

legend({'$\theta$', par1, par2, par3, par4}, ...
    'FontSize', 8, 'Interpreter', 'latex', 'Location', 'NorthEast')

sublaby("   ");  % Verifique se sublaby está corretamente definido

print -dpng images/fig01b.png

figure(3)
clf
subplot(211)
hold on
plot(t,yp1,t,r,t,ym,'Linew',0.5)
grid on
title({'$r, y_m, y_p$'},'FontSize',10,'Interpreter','latex')
par1 = strcat('$y\;(\gamma=',strrep(mat2str(gamma1), ' ', '\ '),')$');
legend({par1,par2,'$r$','$y_m$'},...
    'FontSize',8,'Interpreter','latex','Location','NorthEast')
sublaby("   ");
V = axis;
axis([V(1) V(2) 0 2.5 ]);
print -dpng images\fig01c.png

dims = size(t);
thetas_matrix = ones([dims(1),4])*[thetas(1) 0 0 0;0 thetas(2) 0 0;0 0 thetas(3) 0;0 0 0 thetas(4)];

ttheta1 = theta1 - thetas_matrix;

figure(4)
clf
hold on
plot(e01,ttheta1)
grid on
%axis equal
title({'$e_0 \times \tilde\theta$'},'FontSize',10,'Interpreter','latex')
xlabel({'$e_0$'},'FontSize',10,'Interpreter','latex')
ylabel({'$\tilde\theta$'},'FontSize',10,'Interpreter','latex')
par1 = strcat('$e_0 \times \tilde\theta_1\;(\gamma=',strrep(mat2str(gamma1), ' ', '\ '),')$');
par2 = strcat('$e_0 \times \tilde\theta_1\;(\gamma=',strrep(mat2str(gamma2), ' ', '\ '),')$');
legend({par1,par2},...
    'FontSize',8,'Interpreter','latex','Location','SouthEast')
sublaby("   ");
print -dpng images\fig03d.png

figure(5)
clf
subplot(211)
hold on
plot(t,u1,'Linew',0.5)
grid on
title({'$u$'},'FontSize',10,'Interpreter','latex')
par1 = strcat('$u\;(\gamma=',strrep(mat2str(gamma1), ' ', '\ '),')$');
legend({par1},...
    'FontSize',8,'Interpreter','latex','Location','SouthEast')
sublaby("   ");
print -dpng images\fig03e.png

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
%}
