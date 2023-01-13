clear; close all; clc;
%% Recreating figures from Steele 1979
% Here are the parameters we will use, from Steele 1979:

rho = 1e-6; h = 1; L = 35; % density, height, length

M0 = 1.5e-6; R0 = 2e-3; S0 = 10e3; % value at base

M1 = 0; R1 = 0; S1 = -0.2;% exponential constants, e.g. S = S0 e^(S1 x)

X = 5000;
x = linspace(0,L,X);
M = M0*exp(M1*x); R = R0*exp(R1*x); S = S0*exp(S1*x); % spatial maps

%% Compute velocity for Fig 2 via walking method
freq = 400*sqrt(2).^(0:9); % frequency array for Fig 2 in Hz
omega = 2*pi*freq; % radian frequency
Z = (S./(1j*omega.')) + R + M.*(1j*omega.'); % F x X array, OC impedance

[~, ~, V] = WKB_walkingRF(x,omega,Z,rho,h,20,0.1);
%% Figure 2
figure(2)
subplot(2,1,1)
plot(x,10*log10(abs(V)),...
    "LineWidth",2);%,"Color",'k')
ylim([-70,30])
xlim([0,35])
xlabel("DISTANCE FROM STAPES, mm")
ylabel("MAGNITUDE OF VELOCITY (dB)")
yticks(-70:10:30)
yticklabels(["","-60","","-40","","-20","","0","","20",""])

subplot(2,1,2)
hold on
for ff = 1:length(freq)
    plot(x,unwrap(angle(V(ff,:))),"LineWidth",2);%,"Color",'k')
end
hold off
xlim([0,35])
ylim([-9*pi,pi])
yticks(-8*pi:pi:pi)
yticklabels(["-8\pi","-7\pi","-6\pi","-5\pi","-4\pi","-3\pi","-2\pi","-1\pi","0","\pi"])
xlabel("DISTANCE FROM STAPES, mm")
ylabel("PHASE OF VELOCITY (radians)")

%% Compute velocity for Fig 3 via walking method
freq = linspace(400,12800,500); % frequency array for Fig 2 in Hz
omega = 2*pi*freq; % radian frequency
Z = (S./(1j*omega.')) + R + M.*(1j*omega.'); % F x X array, OC impedance

[~, ~, Vf] = WKB_walkingRF(x,omega,Z,rho,h,20,0.1);
V = Vf(:,715*(1:6));
%% Figure 3
figure(3)

subplot(2,1,1)
semilogx(freq,10*log10(abs(V)'),...
    "LineWidth",2);%,"Color",'k')
ylim([-70,30])
xlim([0,35])
xlabel("FREQUNCY (kHz)")
ylabel("MAGNITUDE OF VELOCITY (dB)")
yticks(-70:10:30)
yticklabels(["","-60","","-40","","-20","","0","","20",""])
xlim([freq(1),freq(end)])

subplot(2,1,2,'XScale','log')
hold on
for xx = 1:length(x)
    plot(freq,unwrap(angle(V(:,xx))),"LineWidth",2);%,"Color",'k')
end
hold off

xlim([freq(1),freq(end)])
ylim([-9*pi,pi])
yticks(-8*pi:pi:pi)
yticklabels(["-8\pi","-7\pi","-6\pi","-5\pi","-4\pi","-3\pi","-2\pi","-1\pi","0","\pi"])