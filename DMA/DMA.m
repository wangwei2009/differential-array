function [ y ] = DMA( x,spacing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first-order DMAs
% refer to
% [1]."differential microphone arrays"
%
% example Usage:
%   y = DMA( x,0.02)
%
% Inputs:
%   x        dual-mic input data,[samples,channel]
%   spacing  mic spacing     
%
% Outputs:
%   y            processed data
%
% Created by Wang wei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = stft(x);
% output spectral
Y = squeeze(X(:,1,:));
fs = 16000;

alpha_1_1 = cos(180*pi/180);
theta_target = 0*pi/180;
d = spacing;
frameLength = 256;
inc = frameLength/2;

 t = 27;
 c = (331.3+0.606*t);
 tao0 = d/c;

 N = frameLength;
 omega = zeros(N/2+1,1);
 HL = zeros(1,N/2+1);
 H = zeros(2,N/2+1);
         
frameNum = size(X,1);
half_bin = size(X,3);

for k = 2:N/2+1
    omega(k) = 2*pi*(k-1)*fs/N;
%     HL(k) = 1/(1-exp(1j*omega(k)*tao0*(alpha_1_1-cos(theta_target))));        
    HL(k) = 1j/(omega(k)*tao0*(alpha_1_1-cos(theta_target)));  % approximating e^x with 1+x
    H(:,k) = HL(k)*[1;
                    -exp(1j*omega(k)*tao0*alpha_1_1)];
end
    
for frameIndex = 1:frameNum

    d = squeeze(X(frameIndex,:,1:half_bin));

     Yout = sum(d.*(H));
     Y(frameIndex,:) = Yout;
end
    y = istft(Y);
    y = real(y);
end


