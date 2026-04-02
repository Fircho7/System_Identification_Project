%%
% Nume si prenume: Firsov Denis
%

clearvars
clc

%% Magic numbers (replace with received numbers)
m = 3; 
n = 13; 

%% Process data (fixed, do not modify)
a1 = 2*(0.15+(m+n/20)/30)*(1000+n*300);
a2 = (1000+n*300);
b0 = (2.2+m+n)/5.5;

rng(m+10*n)
x0_slx = [(-1)^n*(-m/10-rand(1)*m/5); (-1)^m*(n/20+rand(1)*n/100)];

%% Experiment setup (fixed, do not modify)
Ts = 20/a1/1e4; % fundamental step size
Tfin = 36/a1; % simulation duration

gain = 15;
umin = -gain; umax = gain; % input saturation
ymin = -b0*gain/1.8; ymax = b0*gain/1.8; % output saturation

whtn_pow_in = 1e-9*5*(((m-1)*5+n)/5)/2; % input white noise power and sampling time
whtn_Ts_in = Ts*3;
whtn_seed_in = 23341+m+2*n;
q_in = (umax-umin)/pow2(9); % input quantizer (DAC)

whtn_pow_out = 1e-8*5*(((m-1)*8+n)/5)/2; % output white noise power and sampling time
whtn_Ts_out = Ts*5;
whtn_seed_out = 23342-m-2*n;
q_out = (ymax-ymin)/pow2(9); % output quantizer (ADC)

u_op_region = -(m+n/5)/2; % operating point

%% Input setup (can be changed/replaced/deleted)
u0 = 0;     % fixed
ust = 3;  % must be modified (saturation)
t1 = 12/a1; % recommended 

%% Data acquisition (use t, u, y to perform system identification)
out = sim("circuit_electric_R2022b.slx");

t = out.tout;
u = out.u;
y = out.y;

plot(t,u,t,y)
shg

% Savitzky-Golay aplicat
uf = sgolayfilt(u, 1, 19);
yf = sgolayfilt(y, 1, 19);

figure;
subplot(212);
plot(t,uf,t,yf);
subplot(221);
plot(t,u,t,uf);
subplot(222);
plot(t,y,t,yf);
%% System identification

i1 = 4915;
i2 = 5970;
i3 = 11081;
i4 = 11965;

u0 = mean(u(i1:i2));
ust = mean(u(i3:i4));

y0 = mean(y(i1:i2));
yst = mean(y(i3:i4));

K = (yst-y0)/(ust-u0)
%% Partea reala a polilor

i5 = 6036;
i6 = 12045;

t_aux = t(i5:i6);
y_aux = abs(y(i5:i6) - yst);

% Tukey53H, Median_filt aplicat
y_aux = y_aux(:);
y_f = tukey53H(y_aux);
y_f2 = median_filt(y_f,7);

figure;
subplot(212);
plot(t_aux, y_aux);
title('Raspunsul Indicial Nefiltrat');

subplot(221);
plot(t_aux, y_f, 'LineWidth', 1);
title('Raspunsul Indicial Aplicat Tukey53H');

subplot(222);
plot(t_aux(4:end-3), y_f2)
title('Raspunsul Indicial Aplicat Median filt');
grid on
%%
i7 = 12
i8 = 847
i9 = 1722
i10 = 2603
i11 = 3401

t_reg = t_aux([i7,i8,i9,i10,i11]);
y_reg = log(y_aux([i7,i8,i9,i10,i11]));

figure
plot(t_reg,y_reg)

A_reg = [sum(t_reg.^2), sum(t_reg); sum(t_reg), length(t_reg)];
B_reg = [sum(t_reg.*y_reg); sum(y_reg)];

theta = A_reg\B_reg;

Re = theta(1);
%% Partra Imaginara a polilor

i12 = 5995;
i13 = 7746;

Tosc = (t(i12) - t(i13));
Im = 2*pi/Tosc;

%% zeta(tita) si wn 
wn = sqrt(Re^2+Im^2);
zeta = -Re/wn;

% Validare

A = [0,1; -(wn.^2), -2.*zeta*wn];
B = [0; K*(wn^2)];
C = [1 0];
D = 0;

sys = ss(A,B,C,D);

ysim2 = lsim(sys,uf,t,[y(1), 10000]);
figure
plot(t,uf,t,yf,t,ysim2)

J = 1/sqrt(length(t)*norm(yf-ysim2)) *100
er = norm(yf-ysim2)/norm(yf-mean(y))*100
