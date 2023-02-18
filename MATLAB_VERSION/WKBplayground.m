%% WKB approximations of many kinds
% In a 2-D model, k must satisfy k*tanh(kh) = -2*rho*j*omega*Y
% The complex function f(k) = k*tanh(kh) + 2*rho*j*omega*Y has inf roots
% Need the root in the 4th quadrant with the smallest imaginary part
% Here are the parameters we will use, from Steele 1979:

rho = 1e-6; h = 1; L = 35; % density, height, length

M0 = 1.5e-6; R0 = 2e-3; S0 = 10e3; % value at base

M1 = 0; R1 = 0; S1 = -0.2;% exponential constants, e.g. S = S0 e^(S1 x)

X = 5000;
x = linspace(0,L,X);

M = M0*exp(M1*x); R = R0*exp(R1*x); S = S0*exp(S1*x); % spatial maps

F = 500;
freq = linspace(1000,10000,F); % frequency in Hz
omega = 2*pi*freq; % radian frequency

Z = (S./(1j*omega.')) + R + M.*(1j*omega.'); % F x X array, OC impedance
Y = 1./Z; % F x X, OC admittance
%% 1-D k approximation

[k_1D_0, P_1D_0, V_1D_0] = WKB_1D(x,omega,Z,rho,h,0);
[k_1D_1, P_1D_1, V_1D_1] = WKB_1D(x,omega,Z,rho,h,1);
[k_1D_2, P_1D_2, V_1D_2] = WKB_1D(x,omega,Z,rho,h,2);
%%
figure
subplot(2,1,1)
plot(x,10*log10(abs(V_1D_0(1,:))),x,10*log10(abs(V_1D_1(1,:))),...
    x,10*log10(abs(V_1D_2(1,:))))
hold on
plot(x,10*log10(abs(V_1D_0(F/2,:))),x,10*log10(abs(V_1D_1(F/2,:))),...
    x,10*log10(abs(V_1D_2(F/2,:))))
xlim([0,L])
ylim([-60,20])
xlabel("x (mm)")
ylabel("Velocity (dB, unscaled)")
legend("0-order, 1 kHz","1-order, 1 kHz","2-order, 1 kHz",...
    "0-order, 5kHz","1-order, 5kHz","2-order, 5kHz","location",'sw')
title("Cochlear Partition Velocity 1-D WKB Model)")
grid on

subplot(2,1,2)
plot(x,unwrap(angle(V_1D_0(1,:))),x,unwrap(angle(V_1D_1(1,:))),x,unwrap(angle(V_1D_2(1,:))))
hold on
plot(x,unwrap(angle(V_1D_0(F/2,:))),x,unwrap(angle(V_1D_1(F/2,:))),x,unwrap(angle(V_1D_2(F/2,:))))
ylim([-pi*8,pi])
xlim([0,L])
xlabel("x (mm)")
ylabel("Velocity Phase re Stapes (rad)")
yticks(-8*pi:pi:pi)
yticklabels(["-8\pi","-7\pi","-6\pi","-5\pi","-4\pi","-3\pi","-2\pi","-1\pi","0","\pi"])
grid on
%% 2-D long-wave k approximation 1st order

[k_2D_LW, P_2D_LW, V_2D_LW] = WKB_2D_LW(x,omega,Z,rho);
%%
figure
subplot(2,1,1)
plot(x,10*log10(abs(V_2D_LW(1,:))),x,10*log10(abs(V_1D_2(1,:))))
hold on
plot(x,10*log10(abs(V_2D_LW(F/2,:))),x,10*log10(abs(V_1D_2(F/2,:))))
xlim([0,L])
ylim([-60,20])
xlabel("x (mm)")
ylabel("Velocity (dB, unscaled)")
legend("2 kHz 2-D", "2 kHz 1-D (order 2)","5kHz 2-D","5kHz 1-D (order 2)","location",'sw')
title("Cochlear Partition Velocity 2-D Long-Wave WKB Model)")
grid on

subplot(2,1,2)
plot(x,unwrap(angle(V_2D_LW(1,:))),x,unwrap(angle(V_1D_2(1,:))))
hold on
plot(x,unwrap(angle(V_2D_LW(F/2,:))),x,unwrap(angle(V_1D_2(F/2,:))))
ylim([-pi*8,pi])
xlim([0,L])
xlabel("x (mm)")
ylabel("Velocity Phase re Stapes (rad)")
yticks(-8*pi:pi:pi)
yticklabels(["-8\pi","-7\pi","-6\pi","-5\pi","-4\pi","-3\pi","-2\pi","-1\pi","0","\pi"])
grid on

%% 2-D short-wave

[k_2D_SW, P_2D_SW, V_2D_SW] = WKB_2D_SW(x,omega,Z,rho,h);
%%
figure
subplot(2,1,1)
plot(x,10*log10(abs(V_2D_LW(1,:))),x,10*log10(abs(V_2D_SW(1,:))))
hold on
plot(x,10*log10(abs(V_2D_LW(F/2,:))),x,10*log10(abs(V_2D_SW(F/2,:))))
xlim([0,L])
ylim([-60,20])
xlabel("x (mm)")
ylabel("Velocity (dB, unscaled)")
legend("1 kHz 2-D Long-Wave", "1 kHz 2-D Short-Wave",...
    "5kHz 2-D Long-Wave","5kHz 2-D Short-Wave","location",'sw')
title("Cochlear Partition Velocity 2-D Long-Wave WKB Model)")
grid on

subplot(2,1,2)
plot(x,unwrap(angle(V_2D_LW(1,:))),x,unwrap(angle(V_2D_SW(1,:))))
hold on
plot(x,unwrap(angle(V_2D_LW(F/2,:))),x,unwrap(angle(V_2D_SW(F/2,:))))
ylim([-pi*8,pi])
xlim([0,L])
xlabel("x (mm)")
ylabel("Velocity Phase re Stapes (rad)")
yticks(-8*pi:pi:pi)
yticklabels(["-8\pi","-7\pi","-6\pi","-5\pi","-4\pi","-3\pi","-2\pi","-1\pi","0","\pi"])
grid on

%% Iterative Solution via Viergever's method
[k_NR, P_NR, V_NR] = WKB_walkingRF(x,omega,Z,rho,h,20,1);

%% Iterative solution using Alessandro's method
[k_a, P_a, V_a] = WKB_stablePoint(x,omega,Z,rho,h,20,0);

%%
figure
subplot(2,1,1)
loglog(freq/1000,abs(V_NR(:,X/2)),freq/1000,abs(V_NR(:,X/4)),...
    "LineWidth",2)
hold on
loglog(freq/1000,abs(V_2D_LW(:,X/2)),freq/1000,abs(V_2D_LW(:,X/4)))
loglog(freq/1000,abs(V_a(:,X/2)),freq/1000,abs(V_a(:,X/4)))
xlim([freq(1),freq(end)]/1000)
xlabel("Frequency (kHz)")
ylabel("Velocity (unscaled)")
legend("Apex RF","Base RF","Apex LW","Base LW","Apex A","Base A","location",'nw')
title("Cochlear Partition Velocity Comparisons")

ylim([0.001,2000])
subplot(2,1,2)
semilogx(freq/1000,unwrap(angle(V_NR(:,X/2))),freq/1000,...
    unwrap(angle(V_NR(:,X/4))),"LineWidth",2)
hold on
semilogx(freq/1000,unwrap(angle(V_2D_LW(:,X/2))),freq/1000,...
    unwrap(angle(V_2D_LW(:,X/4))))
semilogx(freq/1000,unwrap(angle(V_a(:,X/2))),freq/1000,...
    unwrap(angle(V_a(:,X/4))))
ylim([-12*pi,pi])
yticks(-12*pi:2*pi:2*pi)
yticklabels(["-12\pi","-10\pi","-8\pi","-6\pi","-4\pi","-2\pi","0","2\pi"])
xlim([freq(1),freq(end)]/1000)
grid on
%%
figure
subplot(2,1,1)
plot(x,10*log10(abs(V_NR(1,:))),x,10*log10(abs(V_NR(F/2,:))),...
    "LineWidth",2)
hold on
plot(x,10*log10(abs(V_2D_LW(1,:))),x,10*log10(abs(V_2D_LW(F/2,:))))
plot(x,10*log10(abs(V_a(1,:))),x,10*log10(abs(V_a(F/2,:))))

xlabel("x (mm)")
ylabel("Velocity (unscaled)")
legend("Apex RF","Base RF","Apex LW","Base LW","Apex A","Base A",...
    "location",'sw')
title("Cochlear Partition Velocity 2-D WKB Model")

ylim([-100,30])

subplot(2,1,2)

plot(x,unwrap(angle(V_NR(1,:))),x,...
    unwrap(angle(V_NR(200,:))),"LineWidth",2)
hold on
plot(x,unwrap(angle(V_2D_LW(1,:))),x,...
    unwrap(angle(V_2D_LW(F/2,:))))
plot(x,unwrap(angle(V_a(1,:))),x,...
    unwrap(angle(V_a(F/2,:))))
ylim([-8*pi,pi])
yticks(-8*pi:pi:pi)
yticklabels(["-8\pi","-7\pi","-6\pi","-5\pi","-4\pi","-3\pi","-2\pi","-1\pi","0","\pi"])
