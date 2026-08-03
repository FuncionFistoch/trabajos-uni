function y = miu_ss(c0, m0, a0, teta, rho, n, delta, alfa, g, pi_, x, tol)
syms c m a;
ff = -c/teta*(rho + n + delta - alfa*(1-g)*(a-m)^(alfa-1));
gg = (a-m)^(alfa)*(1-g) - c - (pi_-delta)*m - (n+delta)*a;
hh = c^teta - m/x*((1-g)*alfa*(a-m)^(alfa-1) + (pi_-delta));
J = [diff(ff,c), diff(ff,m), diff(ff,a);
    diff(gg,c), diff(gg,m), diff(gg,a);
    diff(hh,c), diff(hh,m), diff(hh,a)];
F = [ff;
    gg;
    hh];
J_eval = matlabFunction(J, 'Vars', {c, m, a});
F_eval = matlabFunction(F, 'Vars', {c, m, a});
error = 1;
X = [c0;
    m0;
    a0];
w = 0;
while error > tol
    w = w + 1;
    X = X - J_eval(X(1), X(2), X(3))^-1 * F_eval(X(1), X(2), X(3));
    E = F_eval(X(1), X(2), X(3));
    error = sqrt(E(1)^2 + E(2)^2 + E(3)^2);
end
display(J_eval(1.35, 8.79, 16.44));
E = [1 0 0;
    0 1 0;
    0 0 0];
lambda = eig(J_eval(1.35, 8.79, 16.44), E);
display(lambda);
y = [X; w];
display(y);
end
