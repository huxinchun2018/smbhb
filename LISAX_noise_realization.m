function OUT=LISAX_noise_realization(paras)
N=paras(1)/2+1;
fs=paras(2);
c=3.0*10^8;
L=5.0*10^9;
fb=(1:N)/N*fs/2;   %the frequency bin
% PSD need to satisfy
S_pm=2.8*10^(-41)*(1+(fb/3*10^(-3)).^2).^2.*(fb/10^(-4)).^(-8/3);
% S_shot=4.3865*10^(-40)*fb.^2;   %TQ
S_shot=5.3*10^(-38)*fb.^2;   %LISA
S_X=(4*sin(4*pi*fb*L/c).^2+32*sin(2*pi*fb*L/c).^2).*S_pm...
   +16*sin(2*pi*fb*L/c).^2.*S_shot;
%generate the time-series random noise satifying the X noise power 
%spectrum under the choosed sampling frequency.
h_fX1=sqrt(1/2)*S_X.^(1/2);
X_PHI=unifrnd(0,2*pi,1,length(fb)-1);%generate the uniform distributions random phase.
h_fX2=h_fX1.*[1,exp(1i*X_PHI)];
h_fX=[h_fX2,conj(fliplr(h_fX2(2:end)))];%let hf conjugate symmetry,so the ifft time series is real
X_noise1=sqrt(fs*2*N)*ifft(h_fX);%factor sqrt(fs*2*N) is the relationship between continuous and discrete spectrum
X_noise=X_noise1(1:end-1);
OUT=real(X_noise);
end