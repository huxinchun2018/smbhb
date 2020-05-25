function [fitVal,varargout] = TDI_LISA_searchfunc(xVec,params)
%A benchmark test function for CRCBPSO
%F = CRCBPSOTESTFUNC(X,P)
%Compute the Rastrigin fitness function for each row of X.  The fitness
%values are returned in F. X is standardized, that is 0<=X(i,j)<=1. P has
%two arrays P.rmin and P.rmax that are used to convert X(i,j) internally to
%actual coordinate values before computing fitness: X(:,j) ->
%X(:,j)*(rmax(j)-rmin(j))+rmin(j).
%
%For standardized coordinates, F = infty if a point X(i,:) falls
%outside the hypercube defined by 0<=X(i,j)<=1.
%
%[F,R] =  CRCBPSOTESTFUNC(X,P)
%returns the real coordinates in R.
%
%[F,R,Xp] = CRCBPSOTESTFUNC(X,P)
%Returns the standardized coordinates in Xp. This option is to be used when
%there are special boundary conditions (such as wrapping of angular
%coordinates) that are better handled by the fitness function itself.

%Soumya D. Mohanty, Aug 2015
%Just a renamed version of the rastrigin benchmark function.

%Soumya D. Mohanty
%June, 2011
%April 2012: Modified to switch between standardized and real coordinates.

%Shihan Weerathunga
%April 2012: Modified to add the function rastrigin.

%Soumya D. Mohanty
%May 2016: New optional output argument introduced in connection with
%handling of special boundary conditions.

%Soumya D. Mohanty
%Dec 2017: Modified PTAPSOTESTFUNC (just renaming) for the LDAC school.

%Soumya D. Mohanty
%Dec 2018: Changed name
%==========================================================================

%rows: points
%columns: coordinates of a point
[nrows,~]=size(xVec);

%storage for fitness values
fitVal = zeros(nrows,1);
validPts = ones(nrows,1);

%Check for out of bound coordinates and flag them
validPts = crcbchkstdsrchrng(xVec);
%Set fitness for invalid points to infty
fitVal(~validPts)=inf;
%Convert valid points to actual locations
xVec(validPts,:) = s2rv(xVec(validPts,:),params);


for lpc = 1:nrows
    if validPts(lpc)
        % Only the body of this block should be replaced for different fitness
        % functions
        x = xVec(lpc,:);
        fitVal(lpc) = ssrqc(x, params);
    end
end

%Return real coordinates if requested
if nargout > 1
    varargout{1}=xVec;
    if nargout > 2
        varargout{2} = r2sv(xVec,params);
    end
end


function ssrVal = ssrqc(x,params)
%params.c^3   x(1)=tc  x(2)=m  x(3)=eta  x(4)=theta  x(5)=phi

tau=params.c^3*x(3)/(5*params.G*10^6*x(2)*params.m0)*(x(1)*params.yr-params.dataX);

wt=params.c^3/(8*params.G*10^6*x(2)*params.m0)*(tau.^(-3/8)+(743/2688+11/32*x(3))*tau.^(-5/8)-3*pi/10*tau.^(-3/4)...
    +(1855099/14450688+56975/258048*x(3)+371/2048*x(3)^2)*tau.^(-7/8));
PHI=-2/x(3)*(tau.^(5/8)+(3715/8064+55/96*x(3)).*tau.^(3/8)-3*pi/4*tau.^(1/4)...
    +(9275495/14450688+284875/258048*x(3)+1855/2048*x(3)^2)*tau.^(1/8))...
    +wt*params.AU/params.c*sin(x(4)).*cos(2*pi*params.fm*params.dataX-x(5));

xx=(params.G*10^6*x(2)*params.m0*wt/(params.c^3)).^(2/3);
k=150;
xmax=1/7;
taper = 1/2*(1-tanh(k*(xx-xmax)));

q1=zeros(length(params.dataX),3);
q2=zeros(length(params.dataX),3);
q3=zeros(length(params.dataX),3);

sig=zeros(1,3);
for i=1:3
sig(i)=3/2*pi-2*(i-1)*pi/3;
end
L=params.R1*3^(1/2);
omega=2*pi/params.yr;
q1(:,1)=L/3^(1/2)*1/4*(-sin(2*omega*params.dataX.'+2*sig(1))-3*sin(2*sig(1))); % x
q1(:,2)=L/3^(1/2)*1/4*(cos(2*omega*params.dataX.'+2*sig(1))-3*cos(2*sig(1))); % y
q1(:,3)=L/3^(1/2)*1/4*(-12^(1/2)*cos(omega*params.dataX.'-2*sig(1)));  %% spacecraft 1 z

q2(:,1)=L/3^(1/2)*1/4*(-sin(2*omega*params.dataX.'+2*sig(2))-3*sin(2*sig(2))); % x
q2(:,2)=L/3^(1/2)*1/4*(cos(2*omega*params.dataX.'+2*sig(2))-3*cos(2*sig(2))); % y
q2(:,3)=L/3^(1/2)*1/4*(-12^(1/2)*cos(omega*params.dataX.'-2*sig(2)));  %% spacecraft 2 z

q3(:,1)=L/3^(1/2)*1/4*(-sin(2*omega*params.dataX.'+2*sig(3))-3*sin(2*sig(3))); % x
q3(:,2)=L/3^(1/2)*1/4*(cos(2*omega*params.dataX.'+2*sig(3))-3*cos(2*sig(3))); % y
q3(:,3)=L/3^(1/2)*1/4*(-12^(1/2)*cos(omega*params.dataX.'-2*sig(3)));  %% spacecraft 3 z
n1=(q2-q3)/(L);  n2=(q3-q2)/(L);  n3=(q1-q2)/(L);  % unit vector


% figure
% plot(params.dataX,q1(:,1))

%% GW basis 
%beta=theta   lambda=phi
kk=-[cos(x(4))*cos(x(5)) cos(x(4))*sin(x(5)) sin(x(4))];
uu=[sin(x(4))*cos(x(5)) sin(x(4))*sin(x(5)) -cos(x(4))];
vv=[sin(x(5)) -cos(x(5)) 0];

u(1,:)=-1/2*((uu*n1.').^2-(vv*n1.').^2);
u(2,:)=-1/2*((uu*n2.').^2-(vv*n2.').^2);
u(3,:)=-1/2*((uu*n3.').^2-(vv*n3.').^2);
v(1,:)=(uu*n1.').*(vv*n1.');
v(2,:)=(uu*n2.').*(vv*n2.');
v(3,:)=(uu*n3.').*(vv*n3.');

kn2=kk*n2.';   kn3=kk*n3.';
kq2=kk*q2.'/(3^(1/2)*params.R1);   kq3=kk*q3.'/(3^(1/2)*params.R1);

L=params.R1*3^(1/2);
x=wt*L;

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

X(1,:)=4*params.G*10^6*x(2)*params.m0*x(3)/params.c^3*xx.*(u(2,:).*(sinckn2...
    .*coskq23+sinckn22...
    .*coskq25)-u(3,:).*(sinckn3...
    .*coskq35+sinckn32...
    .*coskq33));
X(2,:)=4*params.G*10^6*x(2)*params.m0*x(3)/params.c^3*xx.*(v(2,:).*(sinckn2...
    .*coskq23+sinckn22...
    .*coskq25)-v(3,:).*(sinckn3...
    .*coskq35+sinckn32...
    .*coskq33));

X(3,:)=4*params.G*10^6*x(2)*params.m0*x(3)/params.c^3*xx.*(u(2,:).*(sinckn2...
    .*sinkq23+sinckn22...
    .*sinkq25)-u(3,:).*(sinckn3...
    .*sinkq35+sinckn32...
    .*sinkq33));
X(4,:)=4*params.G*10^6*x(2)*params.m0*x(3)/params.c^3*xx.*(v(2,:).*(sinckn2...
    .*sinkq23+sinckn22...
    .*sinkq25)-v(3,:).*(sinckn3...
    .*sinkq35+sinckn32...
    .*sinkq33));



%% F-Statistic for MBHB
sinx=sin(x);
A(1,:)=2*x.*sinx.*X(1,:).*taper;
A(2,:)=2*x.*sinx.*X(2,:).*taper;
A(3,:)=2*x.*sinx.*X(3,:).*taper;
A(4,:)=2*x.*sinx.*X(4,:).*taper;

NFFT=2^ceil(log2(length(params.dataX)));
A1=zeros(4,NFFT);
A1(:,1:length(params.dataX))=A;

A1f(1,:) = fft(A1(1,:));
A1f(2,:) = fft(A1(2,:));
A1f(3,:) = fft(A1(3,:));
A1f(4,:) = fft(A1(4,:));

df=1/max(params.dataX);
M=zeros(4,4);N=zeros(4,1);
for k=1:4
    for j=1:4
        tmp=(A1f(k,:).*conj(A1f(j,:)));
        M(k,j)=4*real(2*sum(tmp(1:NFFT/2)./params.Sx_LISA)*df);
    end
end

for l=1:4
    tmp2=A1f(l,:).*conj(params.dataY);
    N(l,1)=4*real(2*sum(tmp2(1:NFFT/2)./params.Sx_LISA)*df);
end


%[normSigVec, ~] = normSig(qc, params.sampFreq, params.psd, 1);

%Compute fitness
ssrVal = -1/2*N'*(M\N);

