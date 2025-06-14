clc; clear; close all;

%% Definir Planta
kp = 2;
P = tf(kp, [1 4 3]);

%% Definir Modelo
km = 1;
M = tf(km, [1 5 6]);

%% Polinômio do Observador
A0 = [1 10]; 

%% Rodar cálculo dos parâmetros ideais
[theta1, theta_n, theta2, theta_2n] = controle2DOF(P, M, A0);
