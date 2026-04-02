%%
% Nume si prenume: Firsov Denis
%

clearvars
clc

%% Magic numbers (replace with received numbers)
m = 3;
n = 13;

%% Process data and experiment setup (fixed, do not modify)
u_star = 0.15+n*0.045; % trapeze + PRBS amplitudes
delta = 0.02;
delta_spab = 0.015;

E = 12;  % converter source voltage

umin = 0; umax = 0.98; % input saturation
assert(u_star < umax-0.1)
ymin = 0; ymax = 1/(1-u_star)*E*2; % output saturation

% circuit components + parasitic terms
R = 15;
rL = 10e-3;
rC = 0.2;
rDS1 = 0.01;
rDS2 = 0.01;
Cv = 600e-6/3*m;
Lv = 40e-3*3/m;

% (iL0,uC0)
rng(m+10*n)
x0_slx = [(-1)^(n+1)*E/R,E/3/(1-u_star)];

Ts = 1e-5*(1+2*(u_star-0.15)/u_star); % fundamental step size
Ts = round(Ts*1e6)/1e6;

% input white noise power and sampling time
whtn_pow_in = 1e-11*(Ts*1e4)/2; 
whtn_Ts_in = Ts*2;
whtn_seed_in = 23341+m+2*n;
q_in = (umax-umin)/pow2(11); % input quantizer (DAC)

% output white noise power and sampling time
whtn_pow_out = 1e-7*E*(Ts*1e4/50)*(1+(50*u_star)*(u_star-0.15))/3; 
whtn_Ts_out = Ts*2;
whtn_seed_out = 23342-m-2*n;
q_out = (ymax-ymin)/pow2(11); % output quantizer (ADC)

meas_rep = 13+ceil(n/2); % data acquisition hardware sampling limitation

%% Input setup (can be changed/replaced/deleted)
% Tfin = 2;
% t0 = Tfin/3;
% t1 = Tfin/2;

t1 = 1.2;
N = 5; % la alegere intre 2-11 dar minim 4 cel mai ok
tr = 0.11; % *2 sau *3 daca ai mai multe oscilatii la urcare, daca nu poti merge conform teoriei intre 10%-90% timpul de urcare

p = round(tr/N/Ts);

DeltaT = 2*p*(2^N-1)*Ts; %putem inmulti cu 2 3 ca sa avem perioade de Spab

[input_LUT_dSPACE, Tfin] = generate_input_signal(Ts,t1,DeltaT,N,p,u_star,delta,delta_spab);

%% Data acquisition (use t, u, y to perform system identification)
out = sim("convertor_R2022b.slx");

t = out.tout;
u = out.u;
y = out.y;

subplot(211) 
plot(t,u)
subplot(212)
plot(t,y)

%% System identification
i1 = 74676;
i2 = 122441;
i3 = 154468;
i4 = 203501;

N = 19;

t_id = t(i1:N:i2);
u_id = u(i1:N:i2);
u_id = u_id - mean(u_id);
y_id = y(i1:N:i2);
y_id = y_id - mean(y_id);


t_vd = t(i3:N:i4);
u_vd = u(i3:N:i4);
u_vd = u_vd - mean(u_vd);
y_vd = y(i3:N:i4);
y_vd = y_vd - mean(y_vd);

figure
subplot(221)
plot(t_id,u_id)
subplot(223)
plot(t_id,y_id)

subplot(222)
plot(t_vd, u_vd)
subplot(224)
plot(t_vd, y_vd)


dat_id = iddata(y_id,u_id, t_id(2)-t_id(1));
dat_vd = iddata(y_vd,u_vd, t_vd(2)-t_vd(1));
%% Validare Metoda Armax + Rafinare Model(PEM)
model_armax = armax(dat_id, [6,3,5,2])
figure,
resid(model_armax,dat_vd)
figure,
compare(model_armax,dat_vd)

armax_pem = pem(dat_vd, model_armax);

figure; resid(dat_vd, armax_pem);
figure; compare(dat_vd, armax_pem);
%% Validare Metoda Oe
model_oe = oe(dat_id,[6,4,3])
figure,
resid(model_oe,dat_vd)
figure,
compare(model_oe,dat_vd)
%% Validare Metoda N4sid + Rafinare Model(PEM)
model_n4 = n4sid(dat_id, 1:15)
figure,
resid(model_n4,dat_vd)
figure,
compare(model_n4,dat_vd)

n4_pem = pem(dat_vd, model_n4);

figure; resid(dat_vd, n4_pem);
figure; compare(dat_vd, n4_pem);
