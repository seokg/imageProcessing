function Problem_1()
%% image Restoration
%    blur lenna image
%    add random gaussian noise -> SNR be 12dB
%    restore the blur image 
%    
%    H:LPF corresponding to Hamming Window in the freq domain
%    implement following restoration filters
%    calc PSNR

    close all;
    clear all;
    
    imgdir = uigetdir('Image Directory');
    file = fopen(fullfile(imgdir,'\lenna_gray_256x256.raw'),'rb');
    grayLenna = fread(file,fliplr([256,256]),'*uint8')';
    fclose(file);

%%  1.1 

% % hamming window / blur image
    window = hammingWindow(size(grayLenna,1));
    blurImage = window.*dft2D(double(grayLenna));
    magnitude = log(1+abs(blurImage));
    figure; imshow(magnitude,[]); title('blur image');
    
    
% adding gaussian nosie
%     matlabNosie = imnoise(grayLenna,'gaussian');
%     figure; imshow(matlabNosie,[]); title('matlab noise');
    noiseLenna = uint8(gausNoise(double(grayLenna)/255)*255);
    figure; imshow(noiseLenna,[]); title('Problem 1.1');    

end

function output = gausNoise(input)
    mean = 0;     % default mean
    var = 0.01;  % default variance
  
    output = input + sqrt(var)*randn(size(input)) + mean;
end

function output = hammingWindow(inputSize)
    w = size(inputSize,1);
    for idx = 1 : inputSize
        w(idx) = 0.54 - 0.46*cos(2 * pi * (idx-1)/(inputSize-1));
    end
    
    output = w.'* w;

end


%% DFT
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
        evenOut = dft1D(odd); % start with 0 in equation 
        oddOut = dft1D(even); % start with 1 in equation
        
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
       output(:,j) = dft1D(output(:,j).');
    end


end

function output = invDft1D(input)

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
%         odd = input(2:2:col);
%         even = input(1:2:col);
%         evenOut = invDft1D(even);
%         oddOut = invDft1D(odd);
        evenOut = invDft1D(odd);
        oddOut = invDft1D(even);        
        for j = 1: col/2
           expVal = 1i * 2 * pi * (j-1) / col;
           output(j)= evenOut(j) + oddOut(j) * exp (expVal);
           output(j + col/2) = (evenOut(j) - oddOut(j) * exp (expVal));
           
        end
        output = output/2;
    end
    
end

function output = invDft2D(input)

    [row, col] = size(input);
    output = zeros(row,col);
    
    for i = 1 : row
        output(i,:) = invDft1D(input(i,:));
    end
    
    for j = 1: col
       output(:,j) = invDft1D(output(:,j).');
    end
    
%     output = 1/(row*col) * output;
end
