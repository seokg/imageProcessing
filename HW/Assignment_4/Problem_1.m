function  [output_image, output_image2] = Problem_1()
    % template for ee535 DIP
    % Insert the cod in the designated area below
    
    %% loading directory for image files
    imgdir= uigetdir('Image Directory');
    file = fopen(fullfile(imgdir,'\Gray_baby_512x512.raw'),'rb');
    gray_image = fread(file,fliplr([512,512]),'*uint8')';
    fclose(file);
    
    %% Sample image for test
    centerImg = uint8(ones(32,32)*255);
    padRight = uint8(zeros(256,112));
    padTop = uint8(zeros(112,32));
    midImg = [padTop;centerImg;padTop];
    finalImg = [padRight midImg padRight];    
    
    sampleOutput = fft2(finalImg);
    sampleOutput2 = dft2D(finalImg);
    
    figure; imshow (finalImg,[]); title('orignal');
    figure; imshow (log(1+abs(sampleOutput)),[]); title('matlab');
    figure; imshow (log(1+abs(sampleOutput2)),[]); title('test');
    %% -------------- insert code below --------------%%
    output_image = dft2D(gray_image); % sample code - delete
    output_image2 = fft2(gray_image);
    
    figure; imshow (output_image,[]); title('Problem_1.1');
    figure; imshow (output_image2,[]); title('Problem_1.1');
    

    % abs for magnitude
    % angle for phase
    
    % my code and fft2 is different
    % first do the ifft and compare the reconstructed image for wrong
    % algorithm
    


end

function output = dft1D(input)
    %% functions code
    [~ , col] = size(input);
    output = zeros(1,col);
    if(col==1)
       output = input(col) ;
   
    
    else
        idxOdd = 1;
        idxEven = 1;
        odd = zeros(1, col/2);
        even = zeros(1,col/2);
        for i = 1 : col
            checkOddEven = mod(i,2);
            if (checkOddEven==1)
                odd(idxOdd) = input(i);
                idxOdd = idxOdd + 1;
            else
                even(idxEven) = input(i);
                idxEven = idxEven + 1;
            end 
        end
        evenOut = dft1D(odd);
        oddOut = dft1D(even);
        
        for j = 1: col/2
           expVal = - 1i * 2 * pi * (j-1) / col;
           output(j)= evenOut(j) + oddOut(j) * exp (expVal);
           output(j + col/2) = evenOut(j) - oddOut(j) * exp (expVal);
           
        end
        
    end

end

function  output = dft2D(input)
    [row, col] = size(input);
    output = zeros(row,col);
    
    for i = 1 : row
        output(i,:) = dft1D(input(i,:));
    end
    
    for j = 1: col
       output(:,j) = dft1D(input(:,j)');
    end


end