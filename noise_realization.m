%% start from Time Domain . 2020/5/22

paras=[2^21,1/15];  %2^21 is the number of time series
N=paras(1)/2+1;
fs=paras(2);   %sampling frequency
c=3.0*10^8;
L=1.7*10^8;
fb=(1:N)/N*fs/2;  %the frequency bin
% S_pm=2.5*10^-48*fb.^-2;%the one-side noise power spectral density of proof mass;
% S_op=1.8*10^-37*fb.^2;%the one-side noise power spectral density of optical path;
%the noise' transfer function of AET.
S_pm=2.8*10^(-41)*(1+(fb/3*10^(-3)).^2).^2.*(fb/10^(-4)).^(-8/3);
S_shot=4.3865*10^(-40)*fb.^2; 
% S_shot(i)=5.3*10^(-38)*f(i).^2;   %LISA
S_X=(4*sin(4*pi*fb*L/c).^2+32*sin(2*pi*fb*L/c).^2).*S_pm...
   +16*sin(2*pi*fb*L/c).^2.*S_shot;
%generate the time-series random noise satifying the TDI X combination noise PSD
%spectrum under the choosed sampling frequency.
h_fX1=sqrt(1/2)*S_X.^(1/2);
X_PHI=unifrnd(0,2*pi,1,length(fb)); %generate the uniform distributions random phase.
h_fX2=h_fX1.*exp(1i*X_PHI);

h_fX=[h_fX2,conj(fliplr(h_fX2(2:end)))];%let hf conjugate symmetry,so the ifft time series is real
X_noise1=sqrt(fs*2*N)*ifft(h_fX);%factor sqrt(fs*2*N) is the relationship between continuous and discrete spectrum
X_noise=real(X_noise1(1:end-1));

t=1:1/fs:2*(N-1)/fs;
figure
plot(t,X_noise)
xlabel('t (s)')
ylabel('X noise for TQ')

%Gaussian noise test
figure
histfit(X_noise*10^23)
xlabel('noise amplitude')
ylabel('number')

% PSD test
psd_X=abs(fft(X_noise)).^2/fs/(2*(N-1));
figure
plot(log10(fb),log10(2*psd_X(1:N)))
xlabel('log(f) Hz')
ylabel('TQ TDI X PSD')
hold on
plot(log10(fb),log10(S_X))
%% start from frequncy domain 2020/5/22

% M=2^21; % Number of points in the PSD
% fmax=1/30; % Maximal frequency in the PSD
% f=linspace(10^(-7),fmax,M); %frequency axis for evaluating the PSD
% 
% %PSD function evaluated at the discreet points
% 
% c=3.0*10^8;
% L=1.7*10^8;
% 
% S_ps=1.0*10^(-24);
% S_acc=1.0*10^(-30);
% 
% S_pm=2.8*10^(-41)*(1+(f/3*10^(-3)).^2).^2.*(f/10^(-4)).^(-8/3);
% S_shot=4.3865*10^(-40)*f.^2; 
% % S_shot(i)=5.3*10^(-38)*f(i).^2;   %LISA
% S_X=(4*sin(4*pi*f*L/c).^2+32*sin(2*pi*f*L/c).^2).*S_pm...
%    +16*sin(2*pi*f*L/c).^2.*S_shot;
% 
% PSDfkt=S_X;
% %Add the symmetric part to the PSD function needed for getting real value
% PSDfkt2=[PSDfkt fliplr(PSDfkt(2:end))];
% 
% rn=randn(1,length(PSDfkt2));  %generate time series based on the method of Mina Abdallah
% sig_t=ifft(fft(rn).*sqrt(PSDfkt2))*fmax^(1/2);
% sig=sig_t(1:M);
% %corresponding sample frequency in the time domain
% fsamp=2*fmax;
% t=1:1/fsamp:M/fsamp;

% %use matlab periodogram to estimate PSD
% [pest,fest]=periodogram(sig,[],2^16,fsamp);
% %manual FFT based estimate
% % see https://en.wikipedia.org/wiki/Spectral_density#Power_spectral_density
% pest2=fft(sig).*conj(fft(sig));        %% hxc add pest2=abs(fft(sig)).^2;
% pest2=2*pest2(1:M)/fsamp/length(sig);  %% sig psd calculate
% fest2=fsamp*(0:(M-1))/M;               %% time correponding to frequency . max(fest2)=1/15


