function [ y ] = DMA1( x,omega,Hb,Hf,HL,fs,N,tao0,alpha,beta)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
X = stft(x);
% output spectral
Y = squeeze(X(:,1,:));

frameNum = size(X,1);
half_bin = size(X,3);
    for frameIndex = 1:frameNum
        for k = 2:N/2+1
            omega(k) = 2*pi*(k-1)*fs/N;    
%             if k<16
%                 H(:,k) = 1j/(omega(k)*tao0*(alpha-1))*[1-beta*exp(-1j*omega(k)*tao0*(1-alpha));...
%                                                             -(1-beta)*exp(1j*omega(k)*tao0*alpha)];
%             else
            %         H(:,k) = 1j/((alpha-1)*tao0*omega(k))*[1;-exp(1j*omega(k)*tao0*alpha)];
%             H(:,k) = 1/(1-exp(1j*omega(k)*tao0*(alpha-1)))*[1;-exp(1j*omega(k)*tao0*alpha)];
            Hf(:,k) = [1;
                       -1*exp(-1j*omega(k)*tao0)];
            Hb(:,k) = [-1*exp(-1j*omega(k)*tao0);
                       1;];
            if k<=66
                HL(k) = 1/(2*sin(omega(k)*tao0));
            else
                HL(k) = 0.5;
            end
%             HL(k) = 1/(1-exp(1j*omega(k)*tao0*(alpha-1)));
%             end
        end
        d = squeeze(X(frameIndex,:,1:half_bin));
%          d = fft(x(frameIndex:frameIndex+frameLength-1,:).*hann(frameLength));
         Cf = sum(d.*conj(Hf));
         Cb = sum(d.*conj(Hb));
         C = (Cf - beta*Cb).*conj(HL);
         Y(frameIndex,:) = C;
    end
    y = istft(Y);
    y = real(y);

end


