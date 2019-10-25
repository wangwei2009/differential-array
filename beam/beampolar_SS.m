function [beamOut] = beampolar_SS(Cf,Cb,N,spacing,tao0)
%UNTITLED9 此处显示有关此函数的摘要
%   此处显示详细说明
f = 1:1:20000;
theta = linspace(0,2*pi,360);   % scaning angle
d = spacing;
c = 340;
% tao0 = d/c;
tao = tao0*cos(theta);

omega = 2*pi*f;

half_bin = size(Cf,2);
fs = 16000;
N_FFT = (half_bin-1)*2;

Nele = size(Cf,1);
H = zeros(Nele,length(omega));

beamOut = zeros(length(theta),length(half_bin)); % beamformer output
B12 = zeros(length(theta),length(half_bin)); % beamformer output
B21 = zeros(length(theta),length(half_bin)); % beamformer output
N12 = zeros(length(theta),length(half_bin)); % beamformer output
M12 = zeros(length(theta),length(half_bin)); % beamformer output
Y12_2 = zeros(length(theta),length(half_bin)); % beamformer output
Y12 = zeros(length(theta),length(half_bin)); % beamformer output

for ang = 1:length(theta)
    for freIndex = 1:half_bin
        omega_k = 2*pi*(freIndex-1)*fs/N_FFT; % normalized digital angular frequency
%         omega_k = omega(freIndex);          % analog angular frequency
        a = [1,exp(-1j*omega_k*tao(ang))];  % signal model,steering vector
        B12(ang,freIndex) = a*(Cf(:,freIndex));
        B21(ang,freIndex) = a*(Cb(:,freIndex));
        N12(ang,freIndex) = a*(N(:,freIndex));
        M12(ang,freIndex) = min(abs(B12(ang,freIndex)),abs(B21(ang,freIndex)));
        Y12_2(ang,freIndex) = max(abs(M12(ang,freIndex))^2-abs(N12(ang,freIndex))^2,0);
        Y12(ang,freIndex) = sqrt(Y12_2(ang,freIndex))*1/sqrt(2*(1-cos(omega_k*tao0)));
%         beamOut(ang,freIndex) = min(abs(a*(Cf(:,freIndex))),abs(a*(Cb(:,freIndex)))); % y = w'*a;
    end   
end
k = 16;
figure,polarplot(linspace(0,2*pi,360),real(B12(:,k)));%rlim([-1 1])
hold on,polarplot(linspace(0,2*pi,360),real(B21(:,k)));
hold on,polarplot(linspace(0,2*pi,360),real(N12(:,k)));
hold on,polarplot(linspace(0,2*pi,360),real(M12(:,k)));
hold on,polarplot(linspace(0,2*pi,360),real(Y12(:,k)));
legend('B12','B21','N12','M12','Y12');
beamOut = Y12;
end

