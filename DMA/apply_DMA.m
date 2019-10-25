%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% test first-order DMA
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% clear all;
%addpath(genpath('lib'));
c = 340; % speed of sound

%%
%% load recorded office noise audio

fs = 16000;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
x = loadwav('../wav/4mic_r0.005/target_2mic_ganrao_90/');
d = 0.005*2;
% x = loadwav('wav/xmos/meetingroom_2/');
% d = 0.064;

cfgs = 'config.m';
switch 1
    case 1
        x = x(:,[3,1]); % extract speaker-1
        disp('speaker-1 is in front of mic1')
    case 2
        x = x(:,[2,4]); % extract speaker-2
        disp('speaker-2 is in front of mic4')
    otherwise
        disp('other value')
end
%% process
% x = pcmread('../wav/STEREO_0024.pcm',2)';
% d = 0.025;
y = DMA( x,d);

%% evaluate
%speech = sig.speech;
% [pesq_mos]= pesq_vec(speech, out,fs)
%rmpath(genpath('lib'));
visual( x(:,1),y);
% util.fig(out, fs);


