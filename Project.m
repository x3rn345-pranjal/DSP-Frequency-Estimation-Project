fo = 2400;
% Frequency of Signal
Fs = 2*fo;
T = 1/Fs;
L = 100000;
t = (0:L-1)*T;
A = 100;
Nsnr = 20;
%%%%%%%%%%%%%% SIGNAL GENERATION %%%%%%%%%%%%%%
phi = pi;
% Angles are in radians
hs = A*sin( ((2*pi*fo*t)/Fs) + phi );
hc = A*cos( ((2*pi*fo*t)/Fs) + phi );
h  = hs+hc/1i;
h = awgn(h,Nsnr);
%figure(1);
%plot (t,hs);

%%%%%%%%%%%%%% TSCW %%%%%%%%%%%%%%

M = 25000;
%TRIANGUAL WINDOW LENGTH

p = 4;
%FOURTH ORDER SYSTEM

N = p*M;
%Sequence Length

x = triang(M);
w = conv(x,conv(x,conv(x,x))); 
z = zeros(p-1,1); 
w = [w;z];
%%%%%%%%%%%%% MULTIPLICATION OF SIGNALS %%%%%%%%

w = transpose(w);
h = h.*w;

%%%%%%%%%%%%%%%% DFT %%%%%%%%%%%%%%%%%%%

H = fft(h);
H = abs(real(H));
figure(1);
plot((0:L-1),H/L);
title('Frequency Estimation by Non-even Item interpolation FFT');
xlabel('k','FontWeight','bold');
ylabel('X(k)','FontWeight','bold');

%P2 = abs(H/L);
%P1 = P2(1: L/2 + 1);
%P1(2:end-1) = 2*P1(2:end-1);
%f = Fs*(0:L/2)/L;
%figure(2);
%plot(f,P1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[max1, ind1] = max(H);
H(ind1)      = -Inf;
[max2, ind2] = max(H);
H(ind2)      = -Inf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kreq = min(ind1,ind2)-1;
%if ind1 == kreq
if ind2 == kreq
    temp = max2;
    max2 = max1;
    max1 = temp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = (max2 - max1)/(max2 + max1);
alpha = 0.154*(beta^7) + 0.238*(beta^5) + 0.431*(beta^3) + 4.852*beta;
zeta = alpha + 0.5;
f_estimated = Fs*(kreq + zeta)*Fs/N;
f_estimated = round(f_estimated,4,'significant')