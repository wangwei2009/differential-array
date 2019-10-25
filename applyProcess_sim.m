%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% endfire
% refer to "A Dual-Microphone Speech Enhancement Algorithm
% Based on the Coherence Function"
%
% broadside
% refer to "A coherence-based noise reduction algorithm for binaural
% hearing aids"
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
% clear all;
addpath(genpath('lib'));
addpath(genpath('DMA_SS'));
c = 340; % speed of sound

%%
%% load recorded office noise audio

fs = 16000;

angle = [0,0]/180*pi;
% array spacing
d = 0.025;
r = d/2; 

switch 1
    case 1
        slice = [1,3]; % extract speaker-1
        disp('speaker-1 is in front of mic1')
    case 2
        slice = [2,4]; % extract speaker-2
        disp('speaker-2 is in front of mic2')
    otherwise
        disp('other value')
end

[ sig ] = sim.signal_simulation( r,slice );
x = sig.x;

x1 = x;

frameLength = 256;
inc = frameLength/2;
overlap = frameLength - inc;
     t = 27;
     c = (331.3+0.606*t);
     tao0 = d/c;
     theta0 = 180;
     alpha = cos(theta0/180*pi);
     beta = 1;
     N_FFT = frameLength;
     omega = zeros(N_FFT/2+1,1);
     omega_c = pi/(2*tao0);
     Hf = zeros(2,N_FFT/2+1);
     Hb = zeros(2,N_FFT/2+1);
     HL = zeros(1,N_FFT/2+1);

%% process
y = zeros(size(x,1),1);
% [ out ] = DMA1( x,omega,Hb,Hf,HL,fs,N_FFT,tao0,alpha,beta);
x = [x(:,1),x(:,2)];
[ y] = DMA1_SS( x,d);

%% evaluate
speech = sig.speech;
% [pesq_mos]= pesq_vec(speech, out,fs)
rmpath(genpath('lib'));
rmpath(genpath('DMA_SS'));
visual( x(:,1),y );
% util.fig(out, fs);


