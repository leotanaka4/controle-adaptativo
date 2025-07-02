clc; clear; close all;

%% Definir Planta
kp = 0.1;
b = 1;
P = tf(kp*[1 6 12 8], [1 0 0 0 0]);

%% Definir Modelo
km = 1;
bm = 2;
M = tf(km*[1], [1 1]);

%% Polinômio do Observador
A0 = 1; 

%% Rodar cálculo dos parâmetros ideais
[theta1, theta_n, theta2, theta_2n,Lambda] = controle2DOF(P, M, A0);

n=4;
Ap = [0 0 0 0];
Bp = [6 12 8];
Am = [5 8.75 6.25 1.5];
Bm = [4 4.75 3/2];
ks = kp/km;

for i=1:n-1
k = n - i;
thetas(i) = Bm(k) - Bp(k);
end
thetas(n) = (Ap(1) - Am(1))/kp;
for i=1:n-1
k = n + 1 - i;
thetas(n+i) = (Ap(k) - Am(k))/kp - thetas(n)*Bm(k-1);
end
thetas(2*n) = km/kp;
thetas
