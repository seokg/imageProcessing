function [ newR newG newB ] = HSItoRGB( H, S, I )



    Chroma = S .* I;
    Hdash = H*360.0 / 60.0;
    X = Chroma .* (ones(size(Hdash)) - abs(rem(Hdash,2.0) -ones(size(Hdash))));
    Min = I - Chroma;

    newR = zeros(size(H));
    newG = zeros(size(H));
    newB = zeros(size(H));

    newR = newR +(Hdash < 1.0).*Chroma ; 
    newG = newG +(Hdash < 1.0).*X ;

    newR = newR +(Hdash>= 1.0 & Hdash < 2.0).*X ; 
    newG = newG +(Hdash>= 1.0 & Hdash< 2.0).*Chroma ;

    newG = newG +(Hdash >= 2.0 & Hdash < 3.0).*Chroma ; 
    newB = newB +(Hdash >= 2.0 & Hdash < 3.0).*X ;

    newG = newG +(Hdash >= 3.0 & Hdash < 4.0).*X ; 
    newB = newB +(Hdash >= 3.0 & Hdash < 4.0).*Chroma ;

    newR = newR +(Hdash >= 4.0 & Hdash < 5.0).*X ; 
    newB = newB +(Hdash >= 4.0 & Hdash < 5.0).*Chroma ;

    newR = newR+(Hdash >= 5.0 & Hdash < 6.0).*Chroma ; 
    newB = newB +(Hdash >= 5.0 & Hdash < 6.0).*X ;

    newR = newR + Min;
    newG = newG + Min;
    newB = newB + Min;

end

