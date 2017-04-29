function Problem_1()
    % template for ee535 DIP
    % Insert the cod in the designated area below
    
    %% loading directory for image files
    imgdir= uigetdir('Image Directory');
    file = fopen(fullfile(imgdir,'\Gray_baby_512x512.raw'),'rb');
    gray_image = fread(file,fliplr([512,512]),'*uint8')';
    fclose(file);
    
%     %% Sample image for test
%     % sample image test
%     centerImg = uint8(ones(32,32)*255);
%     padRight = uint8(zeros(256,112));
%     padTop = uint8(zeros(112,32));
%     midImg = [padTop;centerImg;padTop];
%     finalImg = [padRight midImg padRight];    
%     
%     % test
%     sampleOutput = dft2D(finalImg);
%     dftMatrix = shiftResult(sampleOutput);
%     figure; imshow (finalImg,[]); title('orignal');
%     figure; imshow (log(1+abs(dftMatrix)),[]); title('phase');
%     figure; imshow ((angle(dftMatrix)),[]); title('angle');
    
    %% -------------- insert code below --------------%%
    output_image = dft2D(gray_image);
    tempOutput = shiftResult(output_image);
    
    magnitude = log(1+abs(tempOutput));
    phase = log(1+angle(tempOutput));
    figure; imshow (magnitude,[]); title('Problem 1.1 a - FFT magnitude');
    figure; imshow (phase,[]); title('Problem 1.1 a - FFT phase');
    
    
    reconImage = invDft2D(output_image);
    figure; imshow(reconImage,[]); title('Problem1.1 b - reconstrected image')
    
    
    dctImage = dct2D(gray_image);
%     dctImage2 = ourDCT2D(gray_image);
    magnitude = log(1+abs(dctImage));
%     magnitude2 = log(1+abs(dctImage2));
    figure; imshow(magnitude,[]); title('Problem1.1 b - dctImage')
%     figure; imshow(magnitude2,[]); title('Problem1.1 b - dctImage2')
    
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
        evenOut = invDft1D(odd);
        oddOut = invDft1D(even);
        
        for j = 1: col/2
           expVal =  1i * 2 * pi * (j-1) / col;
           output(j)= evenOut(j) + oddOut(j) * exp (expVal);
           output(j + col/2) = evenOut(j) - oddOut(j) * exp (expVal);
           
        end
        
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
    
    output = 1/(row*col) * output;
end

function output = shiftResult(input)
    [row, col] = size(input);
    temp = zeros(row,col);
    for i = 1 : row
        idxRow = i + row/2 ;
        if (idxRow>row)
            idxRow = i - row/2;
        end
        
        for j = 1: col
            idxCol = j + col/2;
            if (idxCol>col)
                idxCol = j - col/2;
            end
            
            temp(idxRow, idxCol) = input(i,j);
            
            
        end
    end
    
    output = temp;
    

end

function output = dct1D(input)
     %% functions code
    [~ , col] = size(input);
    output = zeros(1,col);
    
    
    
    
    y = [input(1:2:col) input(col:-2:2) ];
    dftValue = dft1D(y);
    for j = 1: col
        expVal = -1i*j*pi/(2*col);
        if (j ==1)
            output(j) = sqrt(1/(4*col))*dftValue(j) * exp(expVal)*2;
        else
            output(j) = sqrt(1/(2*col))*dftValue(j) * exp(expVal)*2;
        end
    end    
end

function output = dct2D(input)
    [row, col] = size(input);
    output = zeros(row,col);
    
    for i = 1 : row
        output(i,:) = dct1D(input(i,:));
%         output(i,:) = dct(double(input(i,:)));
        
    end
    
    for j = 1: col
       output(:,j) = dct1D(output(:,j).');
%        output(:,j) = dct(double(output(:,j).'));

    end
end

function [X]=ourDCT(x)

N=length(x);
k=0:(N-1);
n=0:(N-1);

X=2*double(x)*cos(pi/N*(n'+1/2)*k);
end