close all
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
% tc=0.76*yr;
iota=0.7647;
psi=1.373;
z=0.5;
DL=2.941*10^9*pc;
% z=1.0;
% DL=6.828*10^9*pc;
phic=0.9849;


%% waveform model
dt=15;
T=15*2^21;
tc=T+1500;
t=1:dt:T;
m = m1+m2;
eta=m1*m2/m^2; 
%  m = (1+z)*m;  % redshift mass m_z
tau=c^3*eta/(5*G*m)*(tc-t);

wt=c^3/(8*G*m)*(tau.^(-3/8)+(743/2688+11/32*eta)*tau.^(-5/8)-3*pi/10*tau.^(-3/4)...
    +(1855099/14450688+56975/258048*eta+371/2048*eta^2)*tau.^(-7/8));
PHI=-2/eta*(tau.^(5/8)+(3715/8064+55/96*eta).*tau.^(3/8)-3*pi/4*tau.^(1/4)...
    +(9275495/14450688+284875/258048*eta+1855/2048*eta^2)*tau.^(1/8))...
    +wt*AU/c*sin(theta).*cos(2*pi*fm*t-phi);  %phi_Doppler

xx=(G*m*wt/c^3).^(2/3);
hp=2*G*m*eta/(c^2*DL)*(1+cos(iota)^2).*xx.*cos(PHI+phic);
hc=-4*G*m*eta/(c^2*DL)*cos(iota).*xx.*sin(PHI+phic);

k=150;
xmax=1/7;
taper = 1/2*(1-tanh(k*(xx-xmax)));

%% Spacecraft to Earth center in Earth Baycentric System 
thetas=-4.7/180*pi;
phis=120.5/180*pi;
lambda0=0;
R1=10^8/c;  % unit of m
L=R1*3^(1/2);
k=zeros(3,1);
alpha=zeros(length(t),3);
q1=zeros(length(t),3);
q2=zeros(length(t),3);
q3=zeros(length(t),3);

for i=1:3
k(i)=2*pi/3*(i-1)+lambda0;
alpha(:,i)=2*pi*fsc*t+k(i);
end

q1(:,1)=R1*(cos(phis)*sin(thetas)*sin(alpha(:,1))+cos(alpha(:,1))*sin(phis)); % x
q1(:,2)=R1*(sin(phis)*sin(thetas)*sin(alpha(:,1))-cos(alpha(:,1))*cos(phis)); % y
q1(:,3)=R1*(-sin(alpha(:,1))*cos(thetas));  %% spacecraft 1 z

q2(:,1)=R1*(cos(phis)*sin(thetas)*sin(alpha(:,2))+cos(alpha(:,2))*sin(phis)); % x
q2(:,2)=R1*(sin(phis)*sin(thetas)*sin(alpha(:,2))-cos(alpha(:,2))*cos(phis)); % y
q2(:,3)=R1*(-sin(alpha(:,2))*cos(thetas));  %% spacecraft 2 z

q3(:,1)=R1*(cos(phis)*sin(thetas)*sin(alpha(:,3))+cos(alpha(:,3))*sin(phis)); % x
q3(:,2)=R1*(sin(phis)*sin(thetas)*sin(alpha(:,3))-cos(alpha(:,3))*cos(phis)); % y
q3(:,3)=R1*(-sin(alpha(:,3))*cos(thetas));  %% spacecraft 3 z
n1=(q2-q3)/(3^(1/2)*R1);  n2=(q3-q2)/(3^(1/2)*R1);  n3=(q1-q2)/(3^(1/2)*R1);  % unit vector

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

kn1=kk*n1.';  kn2=kk*n2.';   kn3=kk*n3.';
kq1=kk*q1.'/L;  kq2=kk*q2.'/L;   kq3=kk*q3.'/L;

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

sinckn2=sin((1+kn2).*x/2)/((1+kn2).*x/2);
sinckn22=sin((1-kn2).*x/2)/((1-kn2).*x/2);
sinckn3=sin((1+kn3).*x/2)/((1+kn3).*x/2);
sinckn32=sin((1-kn3).*x/2)/((1-kn3).*x/2);

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

sinckn1=sin((1+kn1).*x/2)/((1+kn1).*x/2);
sinckn12=sin((1-kn1).*x/2)/((1-kn1).*x/2);

coskq33=cos(PHI+x/2.*kq3-3*x/2);
coskq35=cos(PHI+x/2.*kq3-5*x/2);
coskq13=cos(PHI+x/2.*kq1-3*x/2);
coskq15=cos(PHI+x/2.*kq1-5*x/2);
sinkq33=sin(PHI+x/2.*kq3-3*x/2);
sinkq35=sin(PHI+x/2.*kq3-5*x/2);
sinkq13=sin(PHI+x/2.*kq1-3*x/2);
sinkq15=sin(PHI+x/2.*kq1-5*x/2);

Y(1,:)=4*G*m*eta/c^3*xx.*(u(3,:).*(sinckn3...
    .*coskq33+sinckn32...
    .*coskq35)-u(1,:).*(sinckn1...
    .*coskq15+sinckn12...
    .*coskq13));
Y(2,:)=4*G*m*eta/c^3*xx.*(v(3,:).*(sinckn3...
    .*coskq33+sinckn32...
    .*coskq35)-v(1,:).*(sinckn1...
    .*coskq15+sinckn12...
    .*coskq13));
Y(3,:)=4*G*m*eta/c^3*xx.*(u(3,:).*(sinckn3...
    .*sinkq33+sinckn32...
    .*sinkq35)-u(1,:).*(sinckn1...
    .*sinkq15+sinckn12...
    .*sinkq13));
Y(4,:)=4*G*m*eta/c^3*xx.*(v(3,:).*(sinckn3...
    .*sinkq33+sinckn32...
    .*sinkq35)-v(1,:).*(sinckn1...
    .*sinkq15+sinckn12...
    .*sinkq13));

a(1,:)=c/DL*(1/2*(1+cos(iota)^2)*cos(phic)*cos(2*psi)-cos(iota)*sin(phic)*sin(2*psi));
a(2,:)=c/DL*(1/2*(1+cos(iota)^2)*cos(phic)*sin(2*psi)+cos(iota)*sin(phic)*cos(2*psi));
a(3,:)=-c/DL*(1/2*(1+cos(iota)^2)*sin(phic)*cos(2*psi)+cos(iota)*cos(phic)*sin(2*psi));
a(4,:)=-c/DL*(1/2*(1+cos(iota)^2)*sin(phic)*sin(2*psi)-cos(iota)*cos(phic)*cos(2*psi));

load('TQXnoise.txt')   %% length(n)=2^21=2097152
n=TQXnoise(:,2);       %% TianQin X time domain noise
X_signal=2*x.*sin(x).*(a(1,:)*X(1,:)+a(2,:)*X(2,:)+a(3,:)*X(3,:)+a(4,:)*X(4,:));
Y_signal=2*x.*sin(x).*(a(1,:)*Y(1,:)+a(2,:)*Y(2,:)+a(3,:)*Y(3,:)+a(4,:)*Y(4,:));
% dat=X_signal+n;
figure
plot(t,X_signal)
X_signal=X_signal.*taper;
Y_signal=Y_signal.*taper;
hold on
plot(t,X_signal,'r-')
xlabel('t (s)')
ylabel('TQ MBHB h(t)')
legend('X_{signal}','X_{signal} with taper')
%% TQ TDI Sx(f)
Nbin=2^ceil(log2(length(t)));
N=length(t);  Fs=1/dt;  T=max(t);
f = Fs*(1:(length(t))/2)/length(t);
% f=linspace(1/T,1/(2*dt),Nbin);
S_pm=2.8*10^(-41)*(1+(f/3*10^(-3)).^2).^2.*(f/10^(-4)).^(-8/3);
S_shot=4.3865*10^(-40)*f.^2;
Sx_TQ=(4*sin(4*pi*f*L).^2+32*sin(2*pi*f*L).^2).*S_pm+16*sin(2*pi*f*L).^2.*S_shot;
figure
plot(log10(f),log10(Sx_TQ),'m-')
title('TQ X MBHB')
xlabel('log (f) Hz')
ylabel('log(hf) Hz^{-1/2}')

%% SNR calculate
hf_X=fft(X_signal);
SNR_X=(4*real(sum((2*abs(hf_X(1:(N)/2)).^2*dt^2)./Sx_TQ)/T))^(1/2);

hf_Y=fft(Y_signal);
SNR_Y=(4*real(sum((2*abs(hf_Y(1:(N)/2)).^2*dt^2)./Sx_TQ)/T))^(1/2);

SNR=(SNR_X^2+SNR_Y^2)^(1/2);
%% cut off waveform for PSO search
% N=2^ceil(log2(length(t)));
% ht=zeros(1,N);
% ht(1:length(t))=X_signal;
ht(1,:)=X_signal;
ht(2,:)=Y_signal;
NFFT=length(t);
df=1/max(t);
hf(1,:)=fft(ht(1,:),NFFT);
hf(2,:)=fft(ht(2,:),NFFT);
%% F-Statistic for MBHB
sinx=sin(x);
A1f(1,:)=fft(2*x.*sinx.*X(1,:).*taper,NFFT);
A1f(2,:)=fft(2*x.*sinx.*X(2,:).*taper,NFFT);
A1f(3,:)=fft(2*x.*sinx.*X(3,:).*taper,NFFT);
A1f(4,:)=fft(2*x.*sinx.*X(4,:).*taper,NFFT);

A2f(1,:)=fft(2*x.*sinx.*Y(1,:).*taper,NFFT);
A2f(2,:)=fft(2*x.*sinx.*Y(2,:).*taper,NFFT);
A2f(3,:)=fft(2*x.*sinx.*Y(3,:).*taper,NFFT);
A2f(4,:)=fft(2*x.*sinx.*Y(4,:).*taper,NFFT);

M=zeros(4,4);  N=zeros(4,1);
for k=1:4
    for j=1:4
        tmp_X=(A1f(k,:).*conj(A1f(j,:)));
        tmp_Y=(A2f(k,:).*conj(A2f(j,:)));
        M(k,j)=4*real(2*sum(tmp_X(1:NFFT/2)./Sx_TQ)*df)+4*real(2*sum(tmp_Y(1:NFFT/2)./Sx_TQ)*df);
    end
end

for l=1:4
    tmpX=A1f(l,:).*conj(hf(1,:));
    tmpY=A2f(l,:).*conj(hf(2,:));
    N(l,1)=4*real(2*sum(tmpX(1:NFFT/2)./Sx_TQ)*df)+4*real(2*sum(tmpY(1:NFFT/2)./Sx_TQ)*df);
end
a=M\N;

a1=a(1);  a2=a(2);  a3=a(3);  a4=a(4);
Ap = ((a1+a4)^2+(a2-a3)^2)^(1/2)+((a1-a4)^2+(a2+a3)^2)^(1/2);
Ac = ((a1+a4)^2+(a2-a3)^2)^(1/2)-((a1-a4)^2+(a2+a3)^2)^(1/2);
AA = (Ap+(Ap^2-Ac^2)^(1/2));

%% Extrinsic parameter
iota_est = acos(Ac/(AA)); 
psi_est = -1/2*atan2((Ap*a4-Ac*a1),-(Ac*a2+Ap*a3));
cc = sign(sin(2*psi));
DL_est = 2*c/AA;
phic_est = -atan2(cc*(Ap*a4-Ac*a1),cc*(Ap*a2+Ac*a3));

toc
LLR1 = 1/2*N'*(inv(M))*N;
LLR2 = 1/2*N'*(M\N);

