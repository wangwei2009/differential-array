function [ y ] = ADMA( x,spacing)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
fs = 16000;
X = stft(x);
% output spectral
Y = squeeze(X(:,1,:));

half_bin = size(Y,2);
N_FFT = (half_bin-1)*2;
frameNum = size(X,1);
omega = zeros(1,half_bin);
HL = zeros(1,half_bin);
Hf = zeros(2,half_bin);
Hb = zeros(2,half_bin);

c = 340;
tao0 = spacing/c;

beta = 1;

for k = 2:half_bin
    omega(k) = 2*pi*(k-1)*fs/N_FFT;    

    Hf(:,k) = [1;
               -1*exp(-1j*omega(k)*tao0)];
    Hb(:,k) = [-1*exp(-1j*omega(k)*tao0);
               1;];
    % compensation filter
    if k<=66
        HL(k) = 1/(2*sin(omega(k)*tao0));
    else
        HL(k) = 0.5;
    end
end
        
half_bin = size(X,3);
    for frameIndex = 1:frameNum

        d = squeeze(X(frameIndex,:,1:half_bin));
         Cf = Hf.*HL;
         Cb = Hb.*HL;
         C = (Cf - beta*Cb);
         Yout = sum(d.*(C));
         Y(frameIndex,:) = Yout;
    end
    y = istft(Y);
    y = real(y);

end


