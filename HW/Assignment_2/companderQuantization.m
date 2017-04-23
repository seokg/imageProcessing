function [ finalImage ] = companderQuantization( nbitQuant, R1, G1, B1 ,fx)


    level = 2^nbitQuant;
    stepSize = 256 / level;
    start = stepSize/2;
    lowerBound = 0;
    upperBound = stepSize;

    map = zeros(size(fx));
    for i = 1:level
        idx = lowerBound<= fx & upperBound> fx;
        map = map + idx * start;
        
        start = start +stepSize;
        lowerBound = lowerBound + stepSize;
        upperBound = upperBound + stepSize;
        
    end
    
    
    R1 = map(R1);
    G1 = map(G1);
    B1 = map(B1);
    

    
    finalImage =  zeros(size(R1,1),size(R1,2),3);
    finalImage(:,:,1) = R1;
    finalImage(:,:,2) = G1;
    finalImage(:,:,3) = B1;
end

