function [beamOut] = beampolar(weights,spacing,tao0)
%UNTITLED9 此处显示有关此函数的摘要
%   此处显示详细说明
f = 1:1:20000;
theta = linspace(0,2*pi,360);   % scaning angle
d = spacing;
c = 340;
% tao0 = d/c;
tao = tao0*cos(theta);

omega = 2*pi*f;

half_bin = size(weights,2);
fs = 16000;
N_FFT = (half_bin-1)*2;

Nele = size(weights,1);
H = zeros(Nele,length(omega));

beamOut = zeros(length(theta),length(half_bin)); % beamformer output
for ang = 1:length(theta)
    for freIndex = 1:half_bin
        omega_k = 2*pi*(freIndex-1)*fs/N_FFT; % normalized digital angular frequency
%         omega_k = omega(freIndex);          % analog angular frequency
        a = [1,exp(-1j*omega_k*tao(ang))];  % signal model,steering vector
        beamOut(ang,freIndex) = a*(weights(:,freIndex)); % y = w'*a;
    end   
end

if(nargout == 0)
    figure,polarplot(linspace(0,2*pi,360),real(beamOut(:,96)));
%     hold on,polarplot(real(beamOut(:,96)));
end
end

