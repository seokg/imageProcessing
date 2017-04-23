function [ R, G, B ] = RGBchan( rawImg )

    R = double(rawImg(:,1:3:end));
    G = double(rawImg(:,2:3:end));
    B = double(rawImg(:,3:3:end));

end

