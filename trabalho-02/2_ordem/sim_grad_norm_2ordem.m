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
disp('Caso: Planta ............. n = 2')
disp('      Grau relativo ..... n* = 1')
disp('      Parâmetros ........ np = 4')
disp(' ')
disp('Algoritmo: Gradiente Normalizado')
disp(' ')
disp('-------------------------------')

%------------------------------------------------ Initialization -----
tfinal = 100;    %Simulation interval
st = 0.01;      %Sample time to workspace

s = tf('s');    %trick!

%--------------------------------------------------------- Filter ----
lambda1 = 2;
lambda0 = 1;

%--------------------------------------------------------- Plant -----
a1 = 2.5;
a0 = 1.5;
b1 = 2;
b0 = 0.8;

P = (b1*s+b0)/(s^2+a1*s+a0)
P = ss(P);

%--------------------------------------------- Initial condition -----
yp0  = [0; 0]
x0   = yp0;

%----------------------------------- Reference signal parameters -----
DC = 1   %Constant
As = 2   %Sine wave amplitude
ws = 0.1*pi  %Frequency

%------------------------------------------------- Matching gain -----
theta_star1 = [b0; b1; a0 - lambda0; a1 - lambda1]   %theta* parametrização 1
theta_star2 = [b0; b1; a0; a1]                       %theta* parametrização 2
theta_star3 = [1/b0; b1/b0; a0/b0; a1/b0]            %theta* parametrização 3
%----------------------------------------- Adaptation parameters -----
gamma = 1*eye(4);         %Adaptation gains
theta0 = [0;0;0;0];       %Adaptation inicial condition
kappa = 0;
P0 = 1*eye(4);
%-------------------------------------------------- GNSimulation -----
theta_star = theta_star1;
sim('gradiente_normalizado_2ordem',tfinal);

yp1 = yp;   %Save results
e01 = e0;
theta1 = theta;
u1 = u;
%-------------------------------------------------- GNSimulation 2 ---
theta_star = theta_star2;
sim('gradiente_normalizado_2ordem',tfinal);

yp2 = yp;   %Save results
e02 = e0;
theta2 = theta;
u2 = u;
%-------------------------------------------------- GNSimulation 3 ---
theta_star = theta_star3;
sim('gradiente_normalizado_2ordem_p3',tfinal);

yp3 = yp;   %Save results
e03 = e0;
theta3 = theta;
u3 = u;
%-------------------------------------------------- LSSimulation -----
theta_star = theta_star1;
sim('least_square_2ordem',tfinal);

yp4 = yp;   %Save results
e04 = e0;
theta4 = theta;
u4 = u;
%-------------------------------------------------- LSSimulation 2 ---
theta_star = theta_star2;
sim('least_square_2ordem',tfinal);

yp5 = yp;   %Save results
e05 = e0;
theta5 = theta;
u5 = u;
%-------------------------------------------------- LSSimulation 3 ---
theta_star = theta_star3;
sim('least_square_2ordem_p3',tfinal);

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

% Salvar Figura 1 (após a subplot com o método GN)
saveas(gcf, '../images/figura2_gn.png')

%---------------------------------------------------------------------------------------------------------
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

% Salvar Figura 2 (após a subplot com o método LS)
saveas(gcf, '../images/figura2_ls.png')

