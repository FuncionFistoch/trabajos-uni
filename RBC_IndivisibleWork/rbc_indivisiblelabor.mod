% Preambulo: Definicion de variables endogenas, exogenas y los parametros
var c, i, y, k, n, r, w, A;
varexo e;

parameters
beta, delta, alfa, rho, teta, n0, B, sigma_ee,
c_ss, i_ss, y_ss, k_ss, n_ss, r_ss, w_ss, A_ss;

% Calibracion del modelo: Parametros y Estados Estacionarios
beta     = 0.984;
delta    = 0.0122723;
alfa     = 0.333;
rho      = 0.979;
teta     = 1.5;
n0       = 0.29;
sigma_ee = 0.072;
B        = teta/(n0*log(1/(1-n0)));
r_ss     = 1/beta + delta - 1;
n_ss     = r_ss*(1-alfa)/(B*(r_ss-alfa*delta));
k_ss     = (alfa/r_ss)^(1/(1-alfa))*n_ss;
i_ss     = delta*(alfa/r_ss)^(1/(1-alfa))*n_ss;
y_ss     = (alfa/r_ss)^(alfa/(1-alfa))*n_ss;
c_ss     = y_ss - i_ss;
w_ss     = B*c_ss;
A_ss     = 1;

% Modelo RBC con trabajo indivisible
model;
1/c = beta*(1/c(+1))*(1+r(+1)-delta);
y = A*(k(-1))^(alfa)*(n)^(1-alfa);
ln(A) = rho*ln(A(-1)) + e;
y = c + i;
k = i + (1-delta)*k(-1);
w = (1-alfa)*y/n;
r = alfa*y/k(-1);
w = B*c;
end;

% Inicializando con los valores del Estado Estacionario
initval;
n = n_ss;
r = r_ss;
A = A_ss;
k = k_ss;
i = i_ss;
y = y_ss;
c = c_ss;
w = w_ss;
end;

% Verificacion del Estado Estacionario
resid;
steady;

% Shock de productividad
shocks;
var e = (sigma_ee)^2;
end;

% Verificacion del modelo
check;

% Solucion del modelo
stoch_simul(order=1, irf=120);
