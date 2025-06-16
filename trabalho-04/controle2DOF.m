function [theta1, theta_n, theta2, theta_2n, L] = controle2DOF(P, M, A0)
    % ===============================
    % Função para cálculo dos vetores θ* (controle 2DOF)
    % Entrada:
    %  - P : função de transferência da planta
    %  - M : função de transferência do modelo de referênciaAdd commentMore actions
    %  - A0: polinômio característico do observador (como vetor)
    % Saída:
    %  - theta1 : vetor do filtro na entrada
    %  - theta_n: ganho direto na saída
    %  - theta2 : vetor do filtro na saída
    %  - theta_2n: ganho no feedforward (referência)
    %  - L: Filtro Lambda  
    % ===============================
    
    syms s

    % ===============================
    % Obter dados
    % ===============================
    [N_p, D_p] = tfdata(P, 'v');
    [N_m, D_m] = tfdata(M, 'v');

    % ===============================
    % Normalização
    % ===============================
    kp = N_p(find(N_p ~= 0, 1, 'first'));
    km = N_m(find(N_m ~= 0, 1, 'first'));
    Np = poly2sym(N_p, s)/kp;
    Dp = poly2sym(D_p, s);
    Nm = poly2sym(N_m, s)/km;
    Dm = poly2sym(D_m, s);
    A0_sym = poly2sym(A0, s);

    % ===============================
    % Graus dos polinômios
    % ===============================
    n = length(D_p) - 1;
    m = length(N_p(find(N_p~=0, 1, 'first'):end)) - 1;
    n_star = n - m;

    %% ==============================
    % Cálculo de Λ(s)
    % ===============================
    Lambda = Nm * A0_sym;

    % Verificar grau de Lambda
    Lambda_coeffs = sym2poly(expand(Lambda));
    grau_Lambda = length(Lambda_coeffs) - 1;

    % Se grau de Lambda for menor que n - 1, multiplicar por (s + λ_i) até atingir o grau desejado
    k = 1;  % índice para gerar lambdas distintos
    while grau_Lambda < n - 1
        %lambda_k = sym(['lambda', num2str(k)]);
        lambda_k = 2*k;
        fator = s + lambda_k;
        Nm = expand(Nm * fator);
        Dm = expand(Dm * fator);
        Lambda = Nm * A0_sym;
        grau_Lambda = grau_Lambda + 1;
        k = k + 1;
    end

    %% ==============================
    % Resolver equação Diofantina:
    % H(s)*D_p + G(s) = Dm * A0
    % ===============================
    deg_H = n_star - 1;
    deg_G = n - 1;

    H = sym('h', [1 deg_H]);
    G = sym('g', [1 deg_G + 1]);

    H_poly = poly2sym([1 H], s);
    G_poly = poly2sym(G, s);

    % Montagem da equação Diofantina
    Eq = expand(H_poly * Dp - kp*G_poly - Dm * A0_sym);

    % Coletar os coeficientes em potências de s
    coeffs_Eq = fliplr(coeffs(Eq, s));

    % Sistema de equações
    eqns = coeffs_Eq == zeros(1, length(coeffs_Eq));

    % Resolver equações

    vars = [H G];
    sol = solve(eqns, vars);

    % Obter nomes dos campos da solução
    fields = fieldnames(sol);

    %% Processar H
    h_fields = fields(startsWith(fields, 'h'));

    if isempty(h_fields)
        H_sol = [];  % Caso não exista H (grau zero)
        H_final = 1; % Constante 1, como em H(s) = 1
    else
        h_fields = sort(h_fields); % Ordena lexicograficamente
        H_cells = cellfun(@(c) sol.(c), h_fields, 'UniformOutput', false);
        H_sol = [H_cells{:}]; % Vetor linha com os coeficientes

        % Construir H(s) = 1 + H_sol(1)*s^1 + H_sol(2)*s^2 + ...
        grau_H = length(H_sol);
        H_final = 1; % Começa com o termo constante 1
        for i = 1:grau_H
            H_final = H_final + H_sol(i) * s^i;
        end
        H_final = expand(H_final);
    end

    %% Processar G
    g_fields = fields(startsWith(fields, 'g'));

    if isempty(g_fields)
        G_sol = [];  % Caso não exista G
        G_final = 0; % Polinômio nulo
    else
        g_fields = sort(g_fields); % Ordena lexicograficamente
        G_cells = cellfun(@(c) sol.(c), g_fields, 'UniformOutput', false);
        G_sol = [G_cells{:}]; % Vetor linha com os coeficientes

        % Construir G(s) = G_sol(1)*s^n + G_sol(2)*s^(n-1) + ... + G_sol(end)
        grau_G = length(G_sol) - 1;
        G_final = 0;
        for i = 0:grau_G
            G_final = G_final + G_sol(i+1) * s^(grau_G - i);
        end
        G_final = expand(G_final);
    end

    %% ==============================
    % Cálculo de F(s)
    % F = Λ - N_p * H
    % ===============================
    F = expand(Lambda - Np * H_final);

    %% Cálculo dos thetas
    theta_2n = km / kp;

    % Representações simbólicas
    G_coeffs = G_cells;
    G_coeffs = G_coeffs(:).'; % garantir vetor linha

    [Lambda_coeffs, ~] = coeffs(Lambda, s);
    Lambda_coeffs = Lambda_coeffs(:).';

    % Graus
    grau_G = length(G_coeffs) - 1;
    grau_L = length(Lambda_coeffs) - 1;
    len_theta = n - 1;

    % Variáveis simbólicas
    theta_n_sym = sym('theta_n_scalar');
    theta2_sym = sym('t', [1 len_theta]);

    % θ_2(s) + θ_n * Λ(s)
    theta2_poly = poly2sym(theta2_sym, s);
    theta_n_term = theta_n_sym * Lambda;
    G_est = expand(theta2_poly + theta_n_term);

    % Igualar coeficientes
    G_est_coeffs = fliplr(coeffs(G_est, s));
    G_coeffs_target = [zeros(1, length(G_est_coeffs) - length(G_coeffs)), G_coeffs];

    % Resolver sistema
    eqns = G_est_coeffs == G_coeffs_target;
    vars = [theta2_sym, theta_n_sym];
    sol = solve(eqns, vars);

    % Extrair soluções
    theta2 = double(arrayfun(@(i) sol.(sprintf('t%d', i)), 1:len_theta));
    theta_n = double(sol.theta_n_scalar);
    
    % θ1 -> Coeficientes de F(s) normalizados por Λ
    theta1 = sym2poly(F);
    theta1 = double(theta1);
    
    L = double(fliplr(coeffs(Lambda,s)));
    %% Exibir resultados
    fprintf('\n=== Resultados ===\n');
    fprintf('theta_1  = %s\n', num2str(theta1));
    fprintf('theta_n  = %s\n', num2str(theta_n));
    fprintf('theta_2  = %s\n', num2str(theta2));
    fprintf('theta_2n = %s\n', num2str(theta_2n));
    fprintf('Lambda = %s\n', Lambda);
end
