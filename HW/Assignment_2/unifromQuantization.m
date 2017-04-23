function [ finalImage ] = unifromQuantization( nbitQuant, R1, G1, B1)


    level = 2^nbitQuant;
    stepSize = 256 / level;
    start = stepSize/2;
    lowerBound = 0;
    upperBound = stepSize;

    for i= 1:level   
        
            boundedPixel = lowerBound< R1 & upperBound> R1;
            idx = find(boundedPixel);
            R1(idx) = start;
            
            boundedPixel = lowerBound< G1 & upperBound> G1;
            idx = find(boundedPixel);
            G1(idx) = start;

            boundedPixel = lowerBound< B1 & upperBound> B1;
            idx = find(boundedPixel);
            B1(idx) =start;



            start = start +stepSize;
            lowerBound = lowerBound + stepSize;
            upperBound = upperBound + stepSize;

    end
    
    finalImage =  zeros(size(R1,1),size(R1,2),3);
    finalImage(:,:,1) = R1;
    finalImage(:,:,2) = G1;
    finalImage(:,:,3) = B1;
end

