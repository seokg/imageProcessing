function [ H, S, I ] = rgb2HSI( R,G,B )

    I = 1/3 .*(R+G+B);

    S=1- (3./(R+G+B)+0.000001).*(min(min(R,G),B));

    H = acosd ((2*R-G-B) ./(2*sqrt( ( (R-G).^2 + (R-B).*(G-B) ) )+0.000001)); % 1/(2*pi)

    H(B>G)=360-H(B>G);
    H=H/360;
end

