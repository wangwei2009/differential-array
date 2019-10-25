function [ y ] = DMA1_b2b( x,y,frameLength,inc,omega,Hb,Hf,HL,fs,N,tao0,alpha,beta)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
win = sqrt(hann(frameLength+1));
win = win(1:end-1);
    for i = 1:inc:length(x(:,1))-(frameLength-inc)
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
         d = fft(x(i:i+frameLength-1,:).*win);
         Cf = sum(d(1:N/2+1,:).*Hf',2);
         Cb = sum(d(1:N/2+1,:).*Hb',2);
         C = (Hf - beta*Hb).*HL;
         yd = sum(d(1:129,:).*C.',2);
%          yd = d(1:129,1);
         fftd = [yd;conj(flipud(yd(2:N/2)))];
         y(i:i+frameLength-1) = y(i:i+frameLength-1)+(ifft(fftd).*win);
    end

end

