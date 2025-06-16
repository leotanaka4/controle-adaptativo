function [theta1, theta_n, theta2, theta_2n] = controle2DOF(P, M, A0)
    % ===============================
    % Função para cálculo dos vetores θ* (controle 2DOF)
    % ===============================
    
    syms s

    %% Obter dados da planta e modelo
    [N_p, D_p] = tfdata(P, 'v');
    [N_m, D_m] = tfdata(M, 'v');

    % Normalização
    kp = N_p(find(N_p ~= 0, 1, 'first'));
    km = N_m(find(N_m ~= 0, 1, 'first'));
    Np = poly2sym(N_p, s)/kp;
    Dp = poly2sym(D_p, s);
    Nm = poly2sym(N_m, s)/km;
    Dm = poly2sym(D_m, s);
    A0_sym = poly2sym(A0, s);

    %% Grau dos polinômios
    n = length(D_p) - 1;
    m = length(N_p(find(N_p~=0, 1, 'first'):end)) - 1;
    n_star = n - m;

    %% Cálculo de Λ(s)
    Lambda = Nm * A0_sym;

    %% Equação Diofantina
    deg_H = max(n_star - 1, 0);
    deg_G = n - 1;

    H = sym('h', [1 deg_H]);
    G = sym('g', [1 deg_G + 1]);
    H_poly = poly2sym([1 H], s); % H(s) = 1 + h1*s + h2*s^2 ...
    G_poly = poly2sym(G, s);

    Eq = expand(H_poly * Dp - kp * G_poly - Dm * A0_sym);
    coeffs_Eq = fliplr(coeffs(Eq, s));

    % Ignorar o termo de maior grau SE E SOMENTE SE a planta for de grau maior
    grau_Eq = feval(symengine, 'degree', Eq, s);
    grau_LHS = feval(symengine, 'degree', expand(H_poly * Dp), s);
    grau_RHS = feval(symengine, 'degree', expand(Dm * A0_sym), s);

    if grau_LHS > grau_RHS
        coeffs_Eq = coeffs_Eq(2:end); % Remove o termo de grau mais alto
    end

    % Monta sistema
    eqns = coeffs_Eq == zeros(1, length(coeffs_Eq));
    vars = [H G];
    sol = solve(eqns, vars);
    fields = fieldnames(sol);

     %% Extrair coeficientes de H
    h_fields = fields(startsWith(fields, 'h'));
    if isempty(h_fields)
        H_sol = [];
    else
        h_fields = sort(h_fields);
        H_cells = cellfun(@(c) sol.(c), h_fields, 'UniformOutput', false);
        H_sol = double([H_cells{:}]); % Vetor linha com os coeficientes
    end

    %% Extrair coeficientes de G
    g_fields = fields(startsWith(fields, 'g'));
    if isempty(g_fields)
        G_sol = [];
    else
        g_fields = sort(g_fields); % Ordena lexicograficamente
        G_cells = cellfun(@(c) sol.(c), g_fields, 'UniformOutput', false);
        G_sol = double([G_cells{:}]); % Vetor linha com os coeficientes
    end

    %% Polinômios finais
    H_final = poly2sym([1 H_sol], s);
    G_final = poly2sym(G_sol, s);

    %% Cálculo de F(s)
    F = expand(Lambda - Np * H_final);

      %% Cálculo dos thetas
    theta_2n = km / kp;

    % Representações simbólicas
    G_coeffs = sym2poly(G_final);
    G_coeffs = G_coeffs(:).'; % garantir vetor linha

    Lambda_coeffs = sym2poly(Lambda);
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
    theta2 = double( arrayfun(@(i) sol.(sprintf('t%d', i)), 1:len_theta) );
    theta_n = double(sol.theta_n_scalar);

    % Corrigir theta1 tamanho
    theta1_raw = sym2poly(F);
    if length(theta1_raw) >= len_theta
        theta1 = theta1_raw(end - len_theta + 1:end);
    else
        theta1 = [zeros(1, len_theta - length(theta1_raw)), theta1_raw];
    end

    %% Exibir resultados
    fprintf('\n=== Resultados ===\n');
    fprintf('theta_1  = %s\n', num2str(theta1));
    fprintf('theta_n  = %s\n', num2str(theta_n));
    fprintf('theta_2  = %s\n', num2str(theta2));
    fprintf('theta_2n = %s\n', num2str(theta_2n));
end
