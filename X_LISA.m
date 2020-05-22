%% TQ MBHB X signal F-Statistic Testing  2020/03/20

%% generate X signal (waveform model)
clear
clc
tic;
c=3.0*10^8;
G=6.67*10^(-11);
yr=365*86400;
fm=1/yr;
fsc=1/(86400*3.65);
m0=2*10^30;
pc=3.08*10^16;
AU=1.5*10^11;

%% parameters for SMBHB
m1=10^6*m0;
m2=2*10^5*m0;
theta=1.6398;
phi=4.1226;
tc=0.76*yr;
iota=0.7647;
psi=1.373;
z=1.0;
DL=6.828*10^9*pc;
phic=0.9849;


%% waveform model
dt=50;
t=0:dt:tc-3000*dt;   % before meger 1500 hours
m = m1+m2;
eta=m1*m2/m^2; 
Mc=m*eta^(3/5)/m0;
%m = (1+z)*m;  % redshift mass m_z
tau=c^3*eta/(5*G*m)*(tc-t);

wt=c^3/(8*G*m)*(tau.^(-3/8)+(743/2688+11/32*eta)*tau.^(-5/8)-3*pi/10*tau.^(-3/4)...
    +(1855099/14450688+56975/258048*eta+371/2048*eta^2)*tau.^(-7/8));
PHI=-2/eta*(tau.^(5/8)+(3715/8064+55/96*eta).*tau.^(3/8)-3*pi/4*tau.^(1/4)...
    +(9275495/14450688+284875/258048*eta+1855/2048*eta^2)*tau.^(1/8))...
    +wt*AU/c*sin(theta).*cos(2*pi*fm*t-phi);  %phi_Doppler

xx=(G*m*wt/c^3).^(2/3);
hp=2*G*m*eta/(c^2*DL)*(1+cos(iota)^2).*xx.*cos(PHI+phic);
hc=-4*G*m*eta/(c^2*DL)*cos(iota).*xx.*sin(PHI+phic);

k=1000;
xmax=0.0286;
taper = 1/2*(1-tanh(k*(xx-xmax)));

%% Spacecraft to Earth center in Earth Baycentric System 
eta0=0;
k=zeros(3,1);
kai=zeros(1,3);
q1=zeros(length(t),3);
q2=zeros(length(t),3);
q3=zeros(length(t),3);

for i=1:3
kai(i)=2*(i-1)*pi/3;
end
omega=2*pi/yr;
L=5*10^9/c;
R1=L/3^(1/2);
%% PHYSICAL REVIEW D 81, 063008 (2010) Mock LISA data challenge for GWD
% q1(:,1)=L/3^(1/2)*1/4*(cos(2*omega*t.'-kai(1))-3*cos(kai(1))); % x
% q1(:,2)=L/3^(1/2)*1/4*(sin(2*omega*t.'-kai(1))-3*sin(kai(1))); % y
% q1(:,3)=L/3^(1/2)*1/4*(-12^(1/2)*cos(omega*t.'-kai(1)));  %% spacecraft 1 z
% 
% q2(:,1)=L/3^(1/2)*1/4*(cos(2*omega*t.'-kai(2))-3*cos(kai(2)));  % x
% q2(:,2)=L/3^(1/2)*1/4*(sin(2*omega*t.'-kai(2))-3*sin(kai(2))); % y
% q2(:,3)=L/3^(1/2)*1/4*(-12^(1/2)*cos(omega*t.'-kai(2)));  %% spacecraft 2 z
% 
% q3(:,1)=L/3^(1/2)*1/4*(cos(2*omega*t.'-kai(3))-3*cos(kai(3))); % x
% q3(:,2)=L/3^(1/2)*1/4*(sin(2*omega*t.'-kai(3))-3*sin(kai(3))); % y
% q3(:,3)=L/3^(1/2)*1/4*(-12^(1/2)*cos(omega*t-kai(3)));  %% spacecraft 3 z


%% PHYSICAL REVIEW D 70, 022003 (2004) Optimal fittering of LISA data
sig=zeros(1,3);
for i=1:3
sig(i)=3/2*pi-2*(i-1)*pi/3;
end
q1(:,1)=L/3^(1/2)*1/4*(-sin(2*omega*t.'+2*sig(1))-3*sin(2*sig(1))); % x
q1(:,2)=L/3^(1/2)*1/4*(cos(2*omega*t.'+2*sig(1))-3*cos(2*sig(1))); % y
q1(:,3)=L/3^(1/2)*1/4*(-12^(1/2)*cos(omega*t.'-2*sig(1)));  %% spacecraft 1 z

q2(:,1)=L/3^(1/2)*1/4*(-sin(2*omega*t.'+2*sig(2))-3*sin(2*sig(2))); % x
q2(:,2)=L/3^(1/2)*1/4*(cos(2*omega*t.'+2*sig(2))-3*cos(2*sig(2))); % y
q2(:,3)=L/3^(1/2)*1/4*(-12^(1/2)*cos(omega*t.'-2*sig(2)));  %% spacecraft 2 z

q3(:,1)=L/3^(1/2)*1/4*(-sin(2*omega*t.'+2*sig(3))-3*sin(2*sig(3))); % x
q3(:,2)=L/3^(1/2)*1/4*(cos(2*omega*t.'+2*sig(3))-3*cos(2*sig(3))); % y
q3(:,3)=L/3^(1/2)*1/4*(-12^(1/2)*cos(omega*t.'-2*sig(3)));  %% spacecraft 3 z
n1=(q2-q3)/(L);  n2=(q3-q2)/(L);  n3=(q1-q2)/(L);  % unit vector
% n1=(q2-q3);  n2=(q3-q2);  n3=(q1-q2);  % unit vector

% figure
% plot(t,q1(:,1))

%% GW basis 
beta=1.6398;  lambda=4.1226;   %beta=theta   lambda=phi
kk=-[cos(beta)*cos(lambda) cos(beta)*sin(lambda) sin(beta)];
uu=[sin(beta)*cos(lambda) sin(beta)*sin(lambda) -cos(beta)];
vv=[sin(lambda) -cos(lambda) 0];

u(1,:)=-1/2*((uu*n1.').^2-(vv*n1.').^2);
u(2,:)=-1/2*((uu*n2.').^2-(vv*n2.').^2);
u(3,:)=-1/2*((uu*n3.').^2-(vv*n3.').^2);
v(1,:)=(uu*n1.').*(vv*n1.');
v(2,:)=(uu*n2.').*(vv*n2.');
v(3,:)=(uu*n3.').*(vv*n3.');

kn2=kk*n2.';   kn3=kk*n3.';
kq2=kk*q2.'/L;   kq3=kk*q3.'/L;

x=wt*L;
% X(1,:)=4*G*m*eta/c^3*xx.*(u(2,:).*(sin((1+kn2).*x/2)/((1+kn2).*x/2)...
%     .*cos(PHI+x/2.*kq2-3*x/2)+sin((1-kn2).*x/2)/((1-kn2).*x/2)...
%     .*cos(PHI+x/2.*kq2-5*x/2))-u(3,:).*(sin((1+kn3).*x/2)/((1+kn3).*x/2)...
%     .*cos(PHI+x/2.*kq3-5*x/2)+sin((1-kn3).*x/2)/((1-kn3).*x/2)...
%     .*cos(PHI+x/2.*kq3-3*x/2)));
% X(2,:)=4*G*m*eta/c^3*xx.*(v(2,:).*(sin((1+kn2).*x/2)/((1+kn2).*x/2)...
%     .*cos(PHI+x/2.*kq2-3*x/2)+sin((1-kn2).*x/2)/((1-kn2).*x/2)...
%     .*cos(PHI+x/2.*kq2-5*x/2))-v(3,:).*(sin((1+kn3).*x/2)/((1+kn3).*x/2)...
%     .*cos(PHI+x/2.*kq3-5*x/2)+sin((1-kn3).*x/2)/((1-kn3).*x/2)...
%     .*cos(PHI+x/2.*kq3-3*x/2)));
% 
% X(3,:)=4*G*m*eta/c^3*xx.*(u(2,:).*(sin((1+kn2).*x/2)/((1+kn2).*x/2)...
%     .*sin(PHI+x/2.*kq2-3*x/2)+sin((1-kn2).*x/2)/((1-kn2).*x/2)...
%     .*sin(PHI+x/2.*kq2-5*x/2))-u(3,:).*(sin((1+kn3).*x/2)/((1+kn3).*x/2)...
%     .*sin(PHI+x/2.*kq3-5*x/2)+sin((1-kn3).*x/2)/((1-kn3).*x/2)...
%     .*sin(PHI+x/2.*kq3-3*x/2)));
% X(4,:)=4*G*m*eta/c^3*xx.*(v(2,:).*(sin((1+kn2).*x/2)/((1+kn2).*x/2)...
%     .*sin(PHI+x/2.*kq2-3*x/2)+sin((1-kn2).*x/2)/((1-kn2).*x/2)...
%     .*sin(PHI+x/2.*kq2-5*x/2))-v(3,:).*(sin((1+kn3).*x/2)/((1+kn3).*x/2)...
%     .*sin(PHI+x/2.*kq3-5*x/2)+sin((1-kn3).*x/2)/((1-kn3).*x/2)...
%     .*sin(PHI+x/2.*kq3-3*x/2)));

sinckn2=sin((1+kn2).*x/2)./((1+kn2).*x/2);
sinckn22=sin((1-kn2).*x/2)./((1-kn2).*x/2);
sinckn3=sin((1+kn3).*x/2)./((1+kn3).*x/2);
sinckn32=sin((1-kn3).*x/2)./((1-kn3).*x/2);

coskq23=cos(PHI+x/2.*kq2-3*x/2);
coskq25=cos(PHI+x/2.*kq2-5*x/2);
sinkq23=sin(PHI+x/2.*kq2-3*x/2);
sinkq25=sin(PHI+x/2.*kq2-5*x/2);
coskq35=cos(PHI+x/2.*kq3-5*x/2);
coskq33=cos(PHI+x/2.*kq3-3*x/2);
sinkq35=sin(PHI+x/2.*kq3-5*x/2);
sinkq33=sin(PHI+x/2.*kq3-3*x/2);

X(1,:)=4*G*m*eta/c^3*xx.*(u(2,:).*(sinckn2...
    .*coskq23+sinckn22...
    .*coskq25)-u(3,:).*(sinckn3...
    .*coskq35+sinckn32...
    .*coskq33));
X(2,:)=4*G*m*eta/c^3*xx.*(v(2,:).*(sinckn2...
    .*coskq23+sinckn22...
    .*coskq25)-v(3,:).*(sinckn3...
    .*coskq35+sinckn32...
    .*coskq33));

X(3,:)=4*G*m*eta/c^3*xx.*(u(2,:).*(sinckn2...
    .*sinkq23+sinckn22...
    .*sinkq25)-u(3,:).*(sinckn3...
    .*sinkq35+sinckn32...
    .*sinkq33));
X(4,:)=4*G*m*eta/c^3*xx.*(v(2,:).*(sinckn2...
    .*sinkq23+sinckn22...
    .*sinkq25)-v(3,:).*(sinckn3...
    .*sinkq35+sinckn32...
    .*sinkq33));

a(1,:)=c/DL*(1/2*(1+cos(iota)^2)*cos(phic)*cos(2*psi)-cos(iota)*sin(phic)*sin(2*psi));
a(2,:)=c/DL*(1/2*(1+cos(iota)^2)*cos(phic)*sin(2*psi)+cos(iota)*sin(phic)*cos(2*psi));
a(3,:)=-c/DL*(1/2*(1+cos(iota)^2)*sin(phic)*cos(2*psi)+cos(iota)*cos(phic)*sin(2*psi));
a(4,:)=-c/DL*(1/2*(1+cos(iota)^2)*sin(phic)*sin(2*psi)-cos(iota)*cos(phic)*cos(2*psi));

sigma=10^(-27);
n=sigma*randn(1,length(t));
X_signal=2*x.*sin(x).*(a(1,:)*X(1,:)+a(2,:)*X(2,:)+a(3,:)*X(3,:)+a(4,:)*X(4,:));

toc

figure
plot(t,X_signal,'.')
xlabel('t(s)')
ylabel('MBHB X signal') 
title('LISA MBHB X signal')

figure
Fs=1/dt;
T=max(t);
N=length(t)-1;
f=Fs*(1:N/2)/N;
psd=abs(fft(X_signal)).^2/T*dt^2;
plot(log10(f),1/2*log10(2*psd(1:N/2)))
xlabel('log(f) Hz')
ylabel('hf Hz^{-1/2}')
title('LISA MBHB X signal ASD')

% %% LISA X_noise
% fs=1/(2*pi*L);
% 
% S_ps=1.0*10^(-22);
% S_acc=9.0*10^(-30);
% 
% T=365*86400;
% dt=50;
% Fs=1/dt;
% t=0:dt:T;
% Nbin=2^ceil(log2(length(t)));
% %f = Fs*(1:(Nbin)/2)/Nbin;
% f=linspace(1/T,1/(2*dt),Nbin);
% hf_M=zeros(1,length(f));
% for i=1:length(f)
% hf_M(i)=1/(2*L*c)*(4*S_ps+8*(1+cos(f(i)/fs).^2)*S_acc/(2*pi*f(i)).^4).^(1/2);
% end
% hf_X=2*abs(sin(f/fs)).*hf_M;
% figure
% plot(log10(f),2*log10(hf_X))
% title('LISA X noise PSD')
% xlabel('log (f) Hz')
% ylabel('Sn(f) Hz^{-1/2}')

%% LISA TDI Sx(f)
Nbin=2^ceil(log2(length(t)));
N=length(t);
f = Fs*(1:(length(t))/2)/length(t);
% f=linspace(1/T,1/(2*dt),N);
% S_pm=2.54*10^(-48)*f.^(-2);
S_pm=2.8*10^(-41)*(1+(f/3*10^(-3)).^2).^2.*(f/10^(-4)).^(-8/3);
S_shot=5.3*10^(-38)*f.^2;
Sx_LISA=(4*(sin(4*pi*f*L)).^2+32*(sin(2*pi*f*L)).^2).*S_pm+16*(sin(2*pi*f*L)).^2.*S_shot;
hold on
plot(log10(f),1/2*log10(Sx_LISA))
title('LISA MBHB X signal')
xlabel('log (f) Hz')
ylabel('log(hf) Hz^{-1/2}')
set(gca,'XLim',[-5,-2])
%% SNR calculate
% SNR=(4*real(sum((abs(fft(X_signal))).^2./Sx_LISA)/T))^(1/2);
hf=fft(X_signal);
SNR=(4*real(sum((2*abs(hf(1:(N)/2)).^2*dt^2)./Sx_LISA)/T))^(1/2);
