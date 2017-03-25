function [ R, G, B ] = RGBchan( rawImg )

    R = double(rawImg(:,1:3:end))/256;
    G = double(rawImg(:,2:3:end))/256;
    B = double(rawImg(:,3:3:end))/256;

end

