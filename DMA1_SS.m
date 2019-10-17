function [ y,beam] = DMA1_SS( x,spacing,omega,Hb,Hf,HL,fs,N,tao0,alpha,beta)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
X = stft(x);
% output spectral
Y = squeeze(X(:,1,:));

frameNum = size(X,1);

half_bin = size(X,3);

     B12 = zeros(2,half_bin);
     B21 = zeros(2,half_bin);
     Hn = zeros(2,half_bin);
     HL = zeros(1,half_bin);
     M12 = zeros(frameNum,half_bin);
     
    for frameIndex = 1:frameNum
        d = squeeze(X(frameIndex,:,1:half_bin));
        a = ones(2,half_bin); % [1,exp(-1j*omega(k)*sin0)],input signal broadside
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
            Hn(:,k) = [1;
                        -1];
            Hf(:,k) = [exp(-1j*omega(k)*tao0);
                        -1];
            Hb(:,k) = [1;
                        -exp(-1j*omega(k)*tao0)];
            HL(k) = 1/sqrt(2*(1-cos(omega(k)*tao0)));
            HL_2(k) = 1j/(omega(k)*tao0*(alpha-1));
            

%             if k<=66
%                 HL(k) = 1/(2*sin(omega(k)*tao0));
%             else
%                 HL(k) = 0.5;
%             end
%             HL(k) = 1/(1-exp(1j*omega(k)*tao0*(alpha-1)));
%             end
        end
        Cf = Hf;%.*HL_2;
        Cb = Hb;%.*HL_2;
        
        B12 = B12.*HL;
        B21 = B21.*HL;
        
        B12 = sum(Hf.*d);
        B21 = sum(Hb.*d);
        N12 = sum(Hn.*d);
        M12 = min(abs(B12),abs(B21));
        phase = angle(B12);
        Y12_2 = max(abs(M12).^2-abs(N12).^2,1e-6);
        Y(frameIndex,:) = sqrt(Y12_2).*HL.*(cos(phase)+1j*(sin(phase)));
        
%         M12(frameIndex,:) = min(abs(sum(d.*(Cf))),abs(sum(d.*(Cb))));
%         M12 = min(abs(Cf),abs(Cb));
%         phase = angle(sum(d.*conj(B12)));
%         img = sqrt(-1);
            
%         Y_2 = max(abs(M12(frameIndex,k)).^2-abs(sum(a.*conj(B21))).^2,0.01);
%         Y(frameIndex,:) = sqrt(Y_2).*(cos(phase)+img*(sin(phase))).*HL;
% %         Y(frameIndex,:) = sum(d.*conj(B12));
%         Y(frameIndex,:) = sum(d.*Hf);
        
%          d = fft(x(frameIndex:frameIndex+frameLength-1,:).*hann(frameLength));
%          Cf = sum(d.*conj(Hf));
%          Cb = sum(d.*conj(Hb));
%          C = (Cf - beta*Cb).*conj(HL);
%          Y(frameIndex,:) = C;
    end
    y = istft(Y);
    y = real(y);
    beam = beampolar_SS(Cf,Cb,Hn,spacing,tao0);
%     beam = beampolar(Hb,spacing,tao0);

end


