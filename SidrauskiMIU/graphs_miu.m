%%Parametros y Estado Inicial
c0=0.96; m0=2.59; a0=7.64;
teta=1.5; rho=5.129*10^-2; n=0.069*10^-2; delta=1.227*10^-2;
alfa=0.3; g=0.11; pi_=0.02; x=0.15; tol=10^-8;

resultado = miu_ss(c0,m0,a0,teta,rho,n,delta,alfa,g,pi_,x,tol);
css=resultado(1); mss=resultado(2); ass=resultado(3);
fprintf('c*=%.4f, m*=%.4f, a*=%.4f\n', css, mss, ass);

%$Datos del Peru 2010 - 2024
anios = (2010:2024)';
m_data = [2.59 2.81 3.22 3.36 3.50 3.43 3.35 3.62 3.93 4.05 6.42 5.31 4.55 4.41 4.86]';
k_data = [5.05 5.45 5.86 6.18 6.78 6.90 6.96 6.97 6.80 6.85 6.82 6.77 6.77 6.91 6.61]';
c_data = [0.96 0.97 1.02 1.03 1.08 1.08 1.09 1.10 1.08 1.10 1.08 1.09 1.09 1.11 1.09]';
a_data = k_data + m_data;   % activos totales = capital + dinero real

%%
function root = local_bracket(fun, lo, hi, npts)
    ks = linspace(lo, hi, npts);
    vals = nan(size(ks));
    for i = 1:length(ks)
        try
            v = fun(ks(i));
            if isreal(v), vals(i) = v; end
        catch
            vals(i) = NaN;
        end
    end
    idx = find(sign(vals(1:end-1)).*sign(vals(2:end)) < 0, 1);
    if isempty(idx)
        root = NaN;
    else
        try
            root = fzero(fun, [ks(idx), ks(idx+1)]);
        catch
            root = NaN;
        end
    end
end

%%Activos Totales vs Activos Monetarios (Curvas Estacionarias)
a_vals = linspace(mss+0.1, 1.8*ass, 300);
k_gg = nan(size(a_vals));
k_hh = nan(size(a_vals));
for i = 1:length(a_vals)
    a = a_vals(i);
    kmin = 1e-4; kmax = a - 1e-4;
    Fgg = @(k) k.^alfa*(1-g) - (pi_-delta)*(a-k) - (n+delta)*a - css;
    Fhh = @(k) css^teta - (a-k)/x.*((1-g)*alfa*k.^(alfa-1) + (pi_-delta));
    k_gg(i) = local_bracket(Fgg, kmin, kmax, 200);
    k_hh(i) = local_bracket(Fhh, kmin, kmax, 200);
end
m_gg = a_vals - k_gg;   % locus a-dot=0
m_hh = a_vals - k_hh;   % locus hh=0 (condición óptima)

%%Activos Totales vs Consumo Per Capita (Curvas Estacionarias)
a_vals2 = linspace(mss+0.1, 1.8*ass, 300);
c_gg = nan(size(a_vals2));
for i = 1:length(a_vals2)
    a = a_vals2(i);
    if a <= mss, continue; end
    k = a - mss;
    Fgg = @(c) k.^alfa*(1-g) - c - (pi_-delta)*mss - (n+delta)*a;
    c_gg(i) = local_bracket(Fgg, 1e-4, 10, 200);
end

%%Recorrido Real Peruano Activos Monetarios y Totales
figure; hold on;
plot(a_vals, m_gg, 'k-', 'LineWidth', 2);
plot(a_vals, m_hh, 'Color',[0 0.6 0.6], 'LineWidth', 2);
plot(a_data, m_data, 'm-o', 'LineWidth', 1.8, 'MarkerFaceColor','m');
text(a_data(1), m_data(1), ' 2010', 'FontSize', 9);
text(a_data(end), m_data(end), ' 2024', 'FontSize', 9);
plot(ass, mss, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 7);
xlabel('a (activos totales)'); ylabel('m (activo monetario)');
title('Diagrama de fase: m vs a (con datos reales 2010-2024)');
legend('\Delta a=0','condición óptima (hh=0)','Perú 2010-2024','Estado Estacionario','Location','best');
grid on;

%%Recorrido Real Peruano Consumo Per Capita y Activos Totales
figure; hold on;
xline(ass, 'k-', 'LineWidth', 2);
plot(a_vals2, c_gg, 'Color',[0 0.6 0.6], 'LineWidth', 2);
plot(a_data, c_data, 'm-o', 'LineWidth', 1.8, 'MarkerFaceColor','m');
text(a_data(1), c_data(1), ' 2010', 'FontSize', 9);
text(a_data(end), c_data(end), ' 2024', 'FontSize', 9);
plot(ass, css, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 7);
xlabel('a (activos totales)'); ylabel('c (consumo per cápita)');
title('Diagrama de fase: c vs a (con datos reales 2010-2024)');
legend('\dot c=0','\dot a=0','Perú 2010-2024','Estado Estacionario','Location','best');
grid on;
