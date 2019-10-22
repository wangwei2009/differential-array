function [ y] = DMA1_SS( x,spacing,omega,Hb,Hf,HL,fs,N,tao0,alpha,beta)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
X = stft(x);
% output spectral
Y = squeeze(X(:,1,:));

frameNum = size(X,1);

N_FFT = 256;
half_bin = N_FFT/2+1;

 
 theta = linspace(0,2*pi,360);   % scaning angle

 Hn = zeros(2,half_bin);
 HL = zeros(1,half_bin);

 H_1st_f = zeros(2,half_bin);
 H_1st_b = zeros(2,half_bin);
 H_1st_n = zeros(2,half_bin);
 HL_f = zeros(1,half_bin);
 HL_b = zeros(1,half_bin);
 HL_n = zeros(1,half_bin);

 B12 = zeros(length(theta),half_bin); % beamformer output
 B21 = zeros(length(theta),length(half_bin)); % beamformer output
 N12 = zeros(length(theta),length(half_bin)); % beamformer output
M12 = zeros(length(theta),length(half_bin)); % beamformer output
Y12_2 = zeros(length(theta),length(half_bin)); % beamformer output
Y12 = zeros(length(theta),length(half_bin)); % beamformer output

 for ang = 1:length(theta)
    for k = 2:half_bin
        omega_k = 2*pi*(k-1)*fs/N_FFT; % normalized digital angular frequency
        a = [1,exp(-1j*omega_k*tao0*cos(theta(ang)))];  % signal model,steering vector

        theta_target = 0*pi/180;
        theta_null = 180*pi/180;
        HL_f(k) = 1/(1-exp(1j*omega_k*tao0*(cos(theta_null)-cos(theta_target))));

        H_1st_f(:,k) = [1;
                    -exp(1j*omega_k*tao0*cos(theta_null))];                
        Cf = H_1st_f(:,k);%.*HL_f(k);
        
        theta_target = 180*pi/180;
        theta_null = 0*pi/180;
        HL_b(k) = 1/(1-exp(1j*omega_k*tao0*(cos(theta_null)-cos(theta_target))));

        H_1st_b(:,k) = [1;
                    -exp(1j*omega_k*tao0*cos(theta_null))];              
        Cb = H_1st_b(:,k);%.*HL_b(k);
        
        theta_target = 0*pi/180;
        theta_null = 90*pi/180;
        HL_n(k) = 1/(1-exp(1j*omega_k*tao0*(cos(theta_null)-cos(theta_target))));

        H_1st_n(:,k) = [1;
                    -exp(1j*omega_k*tao0*cos(theta_null))];
        Cn = H_1st_n(:,k);%.*HL_n(k);
                
        
        HL(k) = 1/sqrt(2*(1-cos(omega_k*tao0)));
                    
        B12(ang,k) = a*Cf;
        B21(ang,k) = a*Cb;
        N12(ang,k) = a*Cn;
        M12(ang,k) = min(abs(B12(ang,k)),abs(B21(ang,k)));
        Y12_2(ang,k) = max(abs(M12(ang,k))^2-abs(N12(ang,k))^2,0);
        Y12(ang,k) = (sqrt(Y12_2(ang,k))*1/sqrt(2*(1-cos(omega_k*tao0))))^(1/48);
        Y12(ang,k) = (sqrt(Y12_2(ang,k)))^(1/48);
%         beamOut(ang,freIndex) = min(abs(a*(Cf(:,freIndex))),abs(a*(Cb(:,freIndex)))); % y = w'*a;
    end   
end
k = 48;
figure,polarplot(linspace(0,2*pi,360),abs(B12(:,k)));%rlim([-1 1])
hold on,polarplot(linspace(0,2*pi,360),abs(B21(:,k)));
hold on,polarplot(linspace(0,2*pi,360),abs(N12(:,k)));
hold on,polarplot(linspace(0,2*pi,360),abs(M12(:,k)));
hold on,polarplot(linspace(0,2*pi,360),abs(Y12(:,k)));
legend('B12','B21','N12','M12','Y12');

     
     
     
for frameIndex = 1:frameNum
    d = squeeze(X(frameIndex,:,1:half_bin));
    
    for k = 2:N/2+1
        a = d(:,k); % [1,exp(-1j*omega(k)*sin0)],input signal broadside
        omega_k = 2*pi*(k-1)*fs/N;
        
        % backward
        theta_target = 0*pi/180;
        theta_null = 180*pi/180;
        HL_f(k) = 1/(1-exp(1j*omega_k*tao0*(cos(theta_null)-cos(theta_target))));

        H_1st_f(:,k) = [1;
                    -exp(1j*omega_k*tao0*cos(theta_null))];                
        Cf = H_1st_f(:,k);%.*HL_f(k);
        
        
        % forward
        theta_target = 180*pi/180;
        theta_null = 0*pi/180;
        HL_b(k) = 1/(1-exp(1j*omega_k*tao0*(cos(theta_null)-cos(theta_target))));
        H_1st_b(:,k) = [1;
                    -exp(1j*omega_k*tao0*cos(theta_null))];      
%         H_1st_b(:,k) = [exp(1j*omega_k*tao0*cos(theta_target));
%                         -1];
        Cb = H_1st_b(:,k);%.*HL_b(k);
        
        % noise
        theta_target = 0*pi/180;
        theta_null = 90*pi/180;
        HL_n(k) = 1/(1-exp(1j*omega_k*tao0*(cos(theta_null)-cos(theta_target))));
        H_1st_n(:,k) = [1;
                    -exp(1j*omega_k*tao0*cos(theta_null))];
        Cn = H_1st_n(:,k).*HL_n(k);
        
        % fixed beamformer
        B12(frameIndex,k) = a.'*Cf;
        B21(frameIndex,k) = a.'*Cb;
        N12(frameIndex,k) = a.'*Cn;
        M12(frameIndex,k) = min(abs(B12(frameIndex,k)),abs(B21(frameIndex,k)));
        % spectral-subtraction
        Y12_2(frameIndex,k) = max(abs(M12(frameIndex,k))^2-abs(N12(frameIndex,k))^2,0);
        phase = angle(B12(frameIndex,k));
        Y(frameIndex,k) = (sqrt(Y12_2(frameIndex,k))*1/sqrt(2*(1-cos(omega_k*tao0))))*(cos(phase)+1j*(sin(phase)));
%         Y(frameIndex,k) = M12(frameIndex,k);   
%         Y(frameIndex,k) = Y(frameIndex,k)*(cos(phase)+1j*(sin(phase)));   
%             Hf(:,k) = [1;
%                        -1*exp(-1j*omega(k)*tao0)];
%             Hb(:,k) = [-1*exp(-1j*omega(k)*tao0);
%                        1;];
%             Hn(:,k) = [1;
%                         -1];
%             Hf(:,k) = [exp(-1j*omega(k)*tao0);
%                         -1];
%             Hb(:,k) = [1;
%                         -exp(-1j*omega(k)*tao0)];
%             HL(k) = 1/sqrt(2*(1-cos(omega(k)*tao0)));
%             HL_2(k) = 1j/(omega(k)*tao0*(alpha-1));
%             
%             theta_target = 0*pi/180;
%             theta_null = 180*pi/180;
%             HL_f(k) = 1/(1-exp(1j*omega(k)*tao0*(cos(theta_null)-cos(theta_target))));
%             H_1st_f(:,k) = [1;
%                         -exp(1j*omega(k)*tao0*cos(theta_null))];
%             
%             theta_target = 180*pi/180;
%             theta_null = 0*pi/180;
%             HL_f(k) = 1/(1-exp(1j*omega(k)*tao0*(cos(theta_null)-cos(theta_target))));
%             
%             H_1st_b(:,k) = [1;
%                         -exp(1j*omega(k)*tao0*cos(theta_null))];
            
            
            

%             if k<=66
%                 HL(k) = 1/(2*sin(omega(k)*tao0));
%             else
%                 HL(k) = 0.5;
%             end
%             HL(k) = 1/(1-exp(1j*omega(k)*tao0*(alpha-1)));
%             end
     end
%         Cf = Hf;%.*HL_2;
%         Cb = Hb;%.*HL_2;
%         
%         Cf = H_1st_f.*HL_f;
%         Cb = H_1st_b.*HL_b;
%         
%         B12 = B12.*HL;
%         B21 = B21.*HL;
%         
%         B12 = sum(Hf.*d);
%         B21 = sum(Hb.*d);
%         N12 = sum(Hn.*d);
%         M12 = min(abs(B12),abs(B21));
%         phase = angle(B12);
%         Y12_2 = max(abs(M12).^2-abs(N12).^2,1e-6);
%         Y(frameIndex,:) = sqrt(Y12_2).*HL.*(cos(phase)+1j*(sin(phase)));
        
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
%     beam = beampolar_SS(Cf,Cb,Hn,spacing,tao0);
%     beam = beampolar(Hb,spacing,tao0);

end


