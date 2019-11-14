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
% close all
% clear all;
%addpath(genpath('lib'));
addpath(genpath('DMA_SS'));
c = 340; % speed of sound

%%
%% load recorded office noise audio

fs = 16000;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = size(x,2);
%% 
frameLength = 256;
overlap = 128;
inc = frameLength - overlap;
N_FFT = 256;
% test xmos 4-mic circular array recordings
x = loadwav('wav/xmos/rec/');
d = 0.064;
% x = loadwav('wav/4mic_r0.005/target_2mic_ganrao_180/');
% d = 0.005;
switch 2
    case 1
        x = x(:,[1,3]); % extract speaker-1
        disp('speaker-1 is in front of mic1')
    case 2
        x = x(:,[4,2]); % extract speaker-2
        disp('speaker-2 is in front of mic4')
    otherwise
        disp('other value')
end

x1 = x;

frameLength = 256;
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
[ y] = DMA1_SS( x,d);

%% evaluate
speech = sig.speech;
% [pesq_mos]= pesq_vec(speech, out,fs)
%rmpath(genpath('lib'));
rmpath(genpath('DMA_SS'));
visual( x(:,1),y );
% util.fig(out, fs);


