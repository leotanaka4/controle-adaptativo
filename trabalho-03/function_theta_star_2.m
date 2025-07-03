clc; clear; close all;

%% Definir Planta
kp = 0.3;
P = tf(kp*[1 4 4], [1 0 0 0 0]);

%% Definir Modelo
km = 1;
M = tf(km, [1 2 1]);

%% Polinômio do Observador
A0 = [1 10]; 

%% Rodar cálculo dos parâmetros ideais
[theta1, theta_n, theta2, theta_2n] = controle2DOF(P, M, A0);

%% Código do Ramon pra conferir - slide 331
% syms s
% P = (1/3)*(s+2)^2/(s^4)
% S = [1 ; s ; s^2]
% theta1 = [-19 -21 -6]
% thetan = -30
% theta2 = [27 75 60]
% theta2n = 3
% %Filters
% Lambda = (s+1)^3
% F = (theta1*S)
% G = (theta2*S) + thetan*Lambda
% %Feedforward & feedback
% Hfb = G/(Lambda - F)
% Hff = theta2n*Lambda/(Lambda - F)
% %Closed-loop transfer function
% M = P*Hff/(1 - P*Hfb)
% M = simplify(M)
% pretty(M)