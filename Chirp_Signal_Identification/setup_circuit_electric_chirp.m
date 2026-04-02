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
Tfin = 36/a1*10; % simulation duration

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
wf = 4900; %wf = wn wosc,usor de citit aprox din rasp la treapta(step)
fmin = wf/2/pi/5;
fmax = wf/2/pi*2;
Ain = 0.85; % ales arbitrat destul de mare incat sa nu confundam cu zgomotul

%% Data acquisition (use t, u, y to perform system identification)
out = sim("circuit_electric_R2022bchirp.slx");

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
Ay = (-6.46307 + 12.0644) / 2
Au = 0.85;
K = Ay/Au

Ayr = (-4.00824 + 14.5079) / 2
Aur = 0.85;
Mr = Ayr/Aur

wr = pi/(0.054256-0.0535605)

r = roots([4*Mr^2,0,-4*Mr^2,0,K^2])

zeta = 0.2777

wn = wr/sqrt(1-2*zeta^2)

A = [0,1; -wn^2, -2*zeta*wn];
B = [0; K*wn^2];
C = [1 0];
D = 0;
sys = ss(A,B,C,D);

ysim2 = lsim(sys,uf,t,[y(1), 10000]);
figure
plot(t,uf,t,yf,t,ysim2)

J = 1/sqrt(length(t)*norm(yf-ysim2)) *100
er = norm(yf-ysim2)/norm(yf-mean(y))*100 
