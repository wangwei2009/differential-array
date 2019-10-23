close all;
theta = [0,90,135]'*pi/180;%0:2*pi/360:2*pi;  %注视方向
theta = linspace(0,2*pi,360);

f = 1:1:8000;

d = 0.020;
c = 340;
tao0 = d/c;
omega = 2*pi*f;

N_FFT = 256;
fs = 16000;
half_bin = N_FFT/2+1;
omega = [1:half_bin]*fs/N_FFT;

H = zeros(2,length(omega));
theta_target = 0*pi/180;             % target
alpha_1_1 = cos(180*pi/180);          % 零点方向
B = zeros(length(theta),length(omega));
for i = 1:length(omega)
    omega_k = 2*pi*omega(i);
    B(:,i) = 1j/(omega_k*tao0*(alpha_1_1-1))*(1-exp(-1*1j*omega_k*tao0*(cos(theta)-alpha_1_1)));
    B(:,i) = (1-exp(-1*1j*omega_k*tao0*(cos(theta)-cos(alpha))));
    H(:,i) = 1j/(omega_k*tao0*(alpha_1_1-cos(theta_target)))*[1;
                                           -exp(1j*omega_k*tao0*alpha_1_1)];
    H(:,i) = 1/(1-exp(1j*omega_k*tao0*(alpha_1_1-cos(theta_target))))*[1;
                                           -exp(1j*omega_k*tao0*alpha_1_1)];
    B(:,i) = 1/(1-exp(1j*omega_k*tao0*(alpha_1_1-cos(theta_target))))*...
                                  (1-exp(-1*1j*omega_k*tao0*(cos(theta)-alpha_1_1)));
%     B(:,i) = 1/(1-alpha_1_1)*...
%                                   (cos(theta)-alpha_1_1);
end
% figure,plot(pow2db(abs(B(1,:))),'b'),ylim([-20,5]),
% set(gca,'XScale','log'),grid on
% hold on
% plot(pow2db(abs(B(90,:))),'r'),ylim([-20,5]),
% hold on,
% plot(pow2db(abs(B(135,:))),'g'),ylim([-20,5]),
% legend('theta = 0','theta = 90','theta = 135');

figure,polarplot(linspace(0,2*pi,360),abs(B(:,16)));
figure, mesh(abs(B(:,:))),title('beampattren');

beamOut = zeros(length(theta),length(f)); % beamformer output
for ang = 1:length(theta)
    for freIndex = 1:length(f)
        omega_k = 2*pi*f(freIndex); % normalized digital angular frequency
%         omega_k = omega(freIndex);          % analog angular frequency
        a = [1,exp(-1j*omega_k*tao0*cos(theta(ang)))];  % signal model,steering vector
        beamOut(ang,freIndex) = a*(H(:,freIndex)); % y = w'*a;
    end   
end

figure,polarplot(linspace(0,2*pi,360),real(beamOut(:,1000)));%rlim([0 1])

% [beamOut] = beampolar(H,d,tao0);

