%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  webrtc window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function window = KaiserBesselDerived(alpha,length)

    half = fix((length + 1) / 2);
    sum = 0.0;
    window = zeros(1,length);

    for i = 0: half
    
        r = (4.0 * i) / length - 1.0;
        sum = sum + I0(pi * alpha * sqrt(1.0 - r * r));
        window(i+1) = sum;
    end
    i = length - 1:-1:half;
%     for i = length - 1; i >= half; --i)
    
        window(length - i) = sqrt(window(length - i) / sum);
        window(i+1) = window(length - i);
%     end
    if (mod(length,2) == 1)
    
        window(half - 1) = sqrt(window(half - 1) / sum);
    end
    
end
function i0 = I0(x)

    y = x / 3.75;

    y = y*y;
    i0 = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
end