syms s
[num, den] = tfdata(P, 'v');
Ps = poly2sym(num, s) / poly2sym(den, s);
S = [1 ; s ]%; s^2]
%Filters
Lambda = (s+1)*(s+2)
F = (theta1*S)
G = (theta2*S) + theta_n*Lambda
%Feedforward & feedback
Hfb = G/(Lambda - F)
Hff = theta_2n*Lambda/(Lambda - F)
%Closed-loop transfer function
Mod = Ps*Hff/(1 - Ps*Hfb)
Mod = simplify(Mod)
pretty(Mod)
M