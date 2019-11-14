% # COLA check:
% # from scipy import signal
% # print(signal.check_COLA(signal.windows.hamming(400,sym=False),400,300))
% #
frameLength = 256;
overlap = 128;
inc = frameLength - overlap;
N_FFT = 512;

window = sqrt(hamming(frameLength+1));
window = window(1:frameLength);
% window = KaiserBesselDerived(1.5,256);

win = zeros(1,inc*9+overlap)';
frameNum = fix((length(win)-overlap)/inc);
win_ind = zeros(frameNum,length(win));
figure,
for n = 1:frameNum
    win_ind(n,(n-1)*inc+1:(n-1)*inc+frameLength) = window.^2;
    hold on,plot(win_ind(n,:))
end
win_sum = sum(win_ind);  % overlap-add window
hold on,plot(win_sum);
