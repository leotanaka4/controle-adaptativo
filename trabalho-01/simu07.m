%---------------------------------------------------------------------
%  COE-603  Controle adaptativo
%
%  Script para simular o exemplo 7
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
disp('      Par�metros ........ np = 2')
disp(' ')
disp('Algoritmo: MRAC Direto')
disp(' ')
disp('-------------------------------')

%------------------------------------------------ Initialization -----
tfinal = 20;    %Simulation interval
st = 0.01;      %Sample time to workspace

s = tf('s');    %trick!

%--------------------------------------------------------- Plant -----
ap = 10;
kp = 1;

P = kp/(s-ap)
P = ss(P);

%----------------------------------------------- Reference model -----
am = 1;
km = 10;

M = km/(s+am)
M = ss(M);

%--------------------------------------------- Initial condition -----
yp0  = 0
x0   = yp0;

ym0  = 0;
xm0  = ym0;

%----------------------------------- Reference signal parameters -----
DC = 1   %Constant

As = 0   %Sine wave amplitude
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

subplot(2,1,1)

% Curvas com cores padrão (azul e laranja) e linha sólida
plot(t, e01, 'LineWidth', 1); hold on;
plot(t, e02, 'LineWidth', 1);

grid on

% Ajuste automático do eixo y com margem
all_errors = [e01(:); e02(:)];
ymin = min(all_errors);
ymax = max(all_errors);
padding = 0.1 * (ymax - ymin);
ylim([ymin - padding, ymax + padding])

title('$e_0$', 'FontSize', 12, 'Interpreter', 'latex')
xlabel('Tempo (s)', 'FontSize', 10, 'Interpreter', 'latex')
ylabel('$e_0$', 'FontSize', 10, 'Interpreter', 'latex')

par1 = ['$e_0\;(\gamma=', strrep(mat2str(gamma1), ' ', '\ '), ')$'];
par2 = ['$e_0\;(\gamma=', strrep(mat2str(gamma2), ' ', '\ '), ')$'];

legend(par1, par2, 'FontSize', 9, 'Interpreter', 'latex', 'Location', 'NorthEast')

sublaby("   "); 
print -dpng images\07_km=10_degrau\fig07a.png

Theta1 = thetas(1) * ones(size(t));
Theta2 = thetas(2) * ones(size(t));

figure(2)
clf
subplot(2,1,1)

% Curvas com linhas sólidas
plot(t, theta1, 'LineWidth', 1); hold on;
plot(t, theta2, 'LineWidth', 1);
plot(t, Theta1, 'k--', 'LineWidth', 1); % Referência theta1* (linha preta tracejada)
plot(t, Theta2, 'k-.', 'LineWidth', 1); % Referência theta2* (linha preta ponto-traço)

grid on

% Título e rótulos
title({'$\theta, \theta^*$'}, 'FontSize', 12, 'Interpreter', 'latex')
xlabel('Tempo (s)', 'FontSize', 10, 'Interpreter', 'latex')
ylabel('$\theta$', 'FontSize', 10, 'Interpreter', 'latex')

% Legenda
par1 = ['$\theta_{1}\;(\gamma=', strrep(mat2str(gamma1), ' ', '\ '), ')$'];
par2 = ['$\theta_{2}\;(\gamma=', strrep(mat2str(gamma1), ' ', '\ '), ')$'];
par3 = ['$\theta_{1}\;(\gamma=', strrep(mat2str(gamma2), ' ', '\ '), ')$'];
par4 = ['$\theta_{2}\;(\gamma=', strrep(mat2str(gamma2), ' ', '\ '), ')$'];

legend({par1, par2, par3, par4, '$\theta_1^*$', '$\theta_2^*$'}, ...
    'FontSize', 9, 'Interpreter', 'latex', 'Location', 'NorthEast')

% Ajuste automático dos limites do eixo y
all_thetas = [theta1(:); theta2(:); Theta1(:); Theta2(:)];
ymin = min(all_thetas);
ymax = max(all_thetas);
padding = 0.1 * (ymax - ymin);
ylim([ymin - padding, ymax + padding])

sublaby("   ");  % Mantido conforme seu código original
print -dpng images\07_km=10_degrau\fig07b.png

figure(3)
clf
subplot(2,1,1)
hold on

% Curvas com linewidth mais visível
plot(t, yp1, 'LineWidth', 1)
plot(t, yp2, 'LineWidth', 1.5)
plot(t, r,   'k--', 'LineWidth', 1)    % Referência r (tracejada preta)
plot(t, ym, 'LineWidth', 1)    % Modelo ym

grid on

% Título e eixos
title({'$r,\ y_m,\ y_p$'}, 'FontSize', 12, 'Interpreter', 'latex')
xlabel('Tempo (s)', 'FontSize', 10, 'Interpreter', 'latex')
ylabel('$y$', 'FontSize', 10, 'Interpreter', 'latex')

% Legenda melhorada
par1 = ['$y\;(\gamma=', strrep(mat2str(gamma1), ' ', '\ '), ')$'];
par2 = ['$y\;(\gamma=', strrep(mat2str(gamma2), ' ', '\ '), ')$'];

legend({par1, par2, '$r$', '$y_m$'}, ...
    'FontSize', 9, 'Interpreter', 'latex', 'Location', 'NorthEast')

% Ajuste automático do eixo y com margem
all_outputs = [yp1(:); yp2(:); r(:); ym(:)];
ymin = min(all_outputs);
ymax = max(all_outputs);
padding = 0.1 * (ymax - ymin);
ylim([ymin - padding, ymax + padding])

sublaby("   ");  % Mantido conforme o original
print -dpng images\07_km=10_degrau\fig07c.png

dims = size(t);
thetas_matrix = ones([dims(1), 2]) * [thetas(1), 0; 0, thetas(2)];

ttheta1 = theta1 - thetas_matrix;
ttheta2 = theta2 - thetas_matrix;

figure(4)
clf
hold on

% Plot com linhas mais espessas
plot(e01, ttheta1(:,1), 'LineWidth', 1)
plot(e01, ttheta1(:,2), 'LineWidth', 1)
plot(e02, ttheta2(:,1), 'LineWidth', 1)
plot(e02, ttheta2(:,2), 'LineWidth', 1)

grid on

% Título e rótulos com LaTeX
title({'$e_0 \times \tilde\theta$'}, 'FontSize', 12, 'Interpreter', 'latex')
xlabel('$e_0$', 'FontSize', 10, 'Interpreter', 'latex')
ylabel('$\tilde\theta$', 'FontSize', 10, 'Interpreter', 'latex')

% Legenda com descrição clara de cada curva
par1 = ['$e_0 \times \tilde\theta_1\;(\gamma=', strrep(mat2str(gamma1), ' ', '\ '), ')$'];
par2 = ['$e_0 \times \tilde\theta_2\;(\gamma=', strrep(mat2str(gamma1), ' ', '\ '), ')$'];
par3 = ['$e_0 \times \tilde\theta_1\;(\gamma=', strrep(mat2str(gamma2), ' ', '\ '), ')$'];
par4 = ['$e_0 \times \tilde\theta_2\;(\gamma=', strrep(mat2str(gamma2), ' ', '\ '), ')$'];

legend({par1, par2, par3, par4}, ...
    'FontSize', 9, 'Interpreter', 'latex', 'Location', 'SouthEast')

% Ajuste dos limites com margem
all_etheta = [ttheta1(:); ttheta2(:)];
ymin = min(all_etheta);
ymax = max(all_etheta);
padding_y = 0.1 * (ymax - ymin);
ylim([ymin - padding_y, ymax + padding_y])

xvals = [e01(:); e02(:)];
xmin = min(xvals);
xmax = max(xvals);
padding_x = 0.05 * (xmax - xmin);
xlim([xmin - padding_x, xmax + padding_x])

% Ativa o sub-rótulo se necessário
sublaby("   ");
print -dpng images\07_km=10_degrau\fig07d.png

figure(5)
clf
subplot(2,1,1)
hold on

% Plotagem com linhas mais visíveis
plot(t, u1, 'LineWidth', 1)
plot(t, u2, 'LineWidth', 1)

grid on

% Título e rótulos
title({'$u$'}, 'FontSize', 12, 'Interpreter', 'latex')
xlabel('Tempo (s)', 'FontSize', 10, 'Interpreter', 'latex')
ylabel('$u$', 'FontSize', 10, 'Interpreter', 'latex')

% Legenda clara e interpretável
par1 = ['$u\;(\gamma=', strrep(mat2str(gamma1), ' ', '\ '), ')$'];
par2 = ['$u\;(\gamma=', strrep(mat2str(gamma2), ' ', '\ '), ')$'];

legend({par1, par2}, ...
    'FontSize', 9, 'Interpreter', 'latex', 'Location', 'SouthEast')

% Ajuste automático do eixo y com margem
all_u = [u1(:); u2(:)];
umin = min(all_u);
umax = max(all_u);
padding = 0.1 * (umax - umin);
ylim([umin - padding, umax + padding])

sublaby("   ");
print -dpng images\07_km=10_degrau\fig07e.png

%------------------------------------------------- Display plots -----
figure(6)
clf

% --- Subplot 1: e0 ---
subplot(2,2,1)
hold on
plot(t, e01, 'LineWidth', 1)
plot(t, e02, 'LineWidth', 1)
grid on

title('$e_0$', 'FontSize', 12, 'Interpreter', 'latex')
xlabel('Tempo (s)', 'FontSize', 10, 'Interpreter', 'latex')
ylabel('$e_0$', 'FontSize', 10, 'Interpreter', 'latex')

par1 = ['$\gamma=', strrep(mat2str(gamma1), ' ', '\ '), '$'];
par2 = ['$\gamma=', strrep(mat2str(gamma2), ' ', '\ '), '$'];
legend({par1, par2}, 'FontSize', 9, 'Interpreter', 'latex', 'Location', 'SouthEast')

% Ajuste automático do eixo y com margem
y_all = [e01(:); e02(:)];
pad = 0.1 * (max(y_all)-min(y_all));
ylim([min(y_all)-pad, max(y_all)+pad])
hold off


% --- Subplot 2: theta e theta* ---
subplot(2,2,2)
hold on
plot(t, theta1, 'LineWidth', 1)
plot(t, theta2, 'LineWidth', 1)
plot(t, Theta1, 'k--', 'LineWidth', 1)
plot(t, Theta2, 'k-.', 'LineWidth', 1)
grid on

title('$\theta,\ \theta^*$', 'FontSize', 12, 'Interpreter', 'latex')
xlabel('Tempo (s)', 'FontSize', 10, 'Interpreter', 'latex')
ylabel('$\theta$', 'FontSize', 10, 'Interpreter', 'latex')

par1c = ['$\theta_{1}\;(\gamma=', strrep(mat2str(gamma1), ' ', '\ '), ')$'];
par2c = ['$\theta_{2}\;(\gamma=', strrep(mat2str(gamma1), ' ', '\ '), ')$'];
par3c = ['$\theta_{1}\;(\gamma=', strrep(mat2str(gamma2), ' ', '\ '), ')$'];
par4c = ['$\theta_{2}\;(\gamma=', strrep(mat2str(gamma2), ' ', '\ '), ')$'];
legend({par1c, par2c, '$\theta_1^*$', par3c, par4c, '$\theta_2^*$'}, ...
    'FontSize', 9, 'Interpreter', 'latex', 'Location', 'NorthEast')

% Ajuste y-limits
yt = [theta1(:); theta2(:); Theta1(:); Theta2(:)];
pad = 0.1 * (max(yt)-min(yt));
ylim([min(yt)-pad, max(yt)+pad])
hold off


% --- Subplot 3: r, y_m, y_p ---
subplot(2,2,3)
hold on
plot(t, yp1, 'LineWidth', 1)
plot(t, yp2, 'LineWidth', 1)
plot(t, r,   'k--', 'LineWidth', 1)
plot(t, ym, 'LineWidth', 1)
grid on

title('$r,\ y_m,\ y_p$', 'FontSize', 12, 'Interpreter', 'latex')
xlabel('Tempo (s)', 'FontSize', 10, 'Interpreter', 'latex')
ylabel('$y$', 'FontSize', 10, 'Interpreter', 'latex')

legend({par1, par2, '$r$', '$y_m$'}, ...
    'FontSize', 9, 'Interpreter', 'latex', 'Location', 'SouthEast')

% Ajuste y-limits
yo = [yp1(:); yp2(:); r(:); ym(:)];
pad = 0.1 * (max(yo)-min(yo));
ylim([min(yo)-pad, max(yo)+pad])
hold off


% --- Subplot 4: u ---
subplot(2,2,4)
hold on
plot(t, u1, 'LineWidth', 1)
plot(t, u2, 'LineWidth', 1)
grid on

title('$u$', 'FontSize', 12, 'Interpreter', 'latex')
xlabel('Tempo (s)', 'FontSize', 10, 'Interpreter', 'latex')
ylabel('$u$', 'FontSize', 10, 'Interpreter', 'latex')

legend({par1, par2}, 'FontSize', 9, 'Interpreter', 'latex', 'Location', 'SouthEast')

% Ajuste y-limits
uu = [u1(:); u2(:)];
pad = 0.1 * (max(uu)-min(uu));
ylim([min(uu)-pad, max(uu)+pad])

% Se você utiliza sublabels:
sublaby("   ");
hold off
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


