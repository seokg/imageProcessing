function [ H, S, I ] = rgb2HSI( R,G,B )

I = 1/3 *(R+G+B);

S = 1 - min(min(R,G),B)./I;

H = 1/(2*pi) * acos (double(2*R-G-B)./(2*sqrt( double((R-G).^2 + (R-B).*(G-B)) )));

end

