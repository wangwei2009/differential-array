%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% test DMA_SS
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% clear all;

c = 340; % speed of sound

%%
%% load recorded office noise audio

fs = 16000;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
frameLength = 256;
overlap = 128;
inc = frameLength - overlap;
N_FFT = 256;
% test xmos 4-mic circular array recordings
x = loadwav('wav/4mic_r0.005/target_2mic_ganrao_90/');
d = 0.005*2;
% x = loadwav('wav/xmos/meetingroom_2/');
% d = 0.064;

cfgs = 'config.m';
switch 2
    case 1
        x = x(:,[1,3]); % extract speaker-1
        disp('speaker-1 is in front of mic1')
    case 2
        x = x(:,[2,4]); % extract speaker-2
        disp('speaker-2 is in front of mic4')
    otherwise
        disp('other value')
end
% x = pcmread('wav/631/32dB/STEREO_0107.pcm',2)'*10;
% d = 0.025;
% M = size(x,2);
% x1 = x;

% x = loadwav('wav/cmu/');
% d = 0.03;

frameLength = 256;
overlap = frameLength - inc;
     t = 27;
     c = 340;%(331.3+0.606*t);
     null = 90*pi/180;
     tao0 = sin(null)*d/c;
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
% x = pcmread('wav/STEREO_0024.pcm',2)';
y = zeros(size(x,1),1);
% WPE_out_2ch = fdndlp(x, cfgs);
[ y] = DMA1_SS( x,d);

%% evaluate
%speech = sig.speech;
% [pesq_mos]= pesq_vec(speech, out,fs)

visual( x(:,1),y*2 );
% util.fig(out, fs);


