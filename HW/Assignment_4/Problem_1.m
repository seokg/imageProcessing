function [blocDctImage, reconBlocDCTImage]=Problem_1()
    % template for ee535 DIP
    % Insert the cod in the designated area below
    
    %% loading directory for image files
    imgdir= uigetdir('Image Directory');
    file = fopen(fullfile(imgdir,'\Gray_baby_512x512.raw'),'rb');
    gray_image = fread(file,fliplr([512,512]),'*uint8')';
    fclose(file);

    
    %% -------------- insert code below --------------%%
%     output_image = dft2D(gray_image);
%     tempOutput = shiftResult(output_image);
%     
%     magnitude = log(1+abs(tempOutput));
%     phase = log(1+angle(tempOutput));
%     figure; imshow (magnitude,[]); title('Problem 1.1 a - FFT magnitude');
%     figure; imshow (phase,[]); title('Problem 1.1 a - FFT phase');
%     
%     
%     reconImage = invDft2D(output_image);
%     reconImage2 = ifft2(output_image);
%     figure; imshow(reconImage,[]); title('Problem1.1 b - reconstrected image')
%     
%     
%     dctImage = dct2D(gray_image);
%     magnitude = log(1+abs(dctImage));
%     figure; imshow(magnitude,[]); title('Problem1.2 a - dctmagnitude')
%      
%     
%     reconDCTImage = invDct2D(dctImage);
%     figure; imshow(reconDCTImage,[]); title('Problem1.2 a - dct reconstructed Image')
%     
%     blocDctImage = blocDct2D(gray_image,16);
%     magnitude = log(1+abs(blocDctImage));
%     figure; imshow(magnitude,[]); title('Problem1.2 b - block dctmagnitude')
%     
%     reconBlocDCTImage = blocInvDct2D(blocDctImage,16);
%     figure; imshow(reconBlocDCTImage,[]); title('Problem1.2 a - block dct reconstructed Image')
    
     dftImage = dht2D(double(gray_image));
     magnitude = log(1+abs(dftImage));
    figure; imshow(magnitude*255,[]); title('Problem1.3  - DHT: magnitude')
    reconDhtImage= invDht2D(dftImage);
    figure; imshow(reconDhtImage,[]); title('Problem1.3  - DHT: reconstructed Image')
    
     [downLow, downHigh] = dwt1D(gray_image, 'horizontal');
     [W1, W2] = dwt1D(downLow, 'vertical');
     [W3, W4] = dwt1D(downHigh, 'vertical');
    figure; imshow([W1 W2; W3 W4],[]); title('Problem1.3  - DHT: magnitude')

    
    
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
%% DCT
function output = dct1D(input)
     %% functions code
     

    [~ , col] = size(input);
    if col == 1
        output = input;
     else   
    
        y = [input(1:2:(col)) input(col:-2:2) ];
        output = dft1D(y);
        
        w = zeros(1,col);
        for idx = 1: col
            expVal = -pi*1i*(idx-1)/ (2*col);
            w(idx) = 2*exp(expVal)/sqrt(2*col);
%             if idx == 1
%                 y(idx) = y(idx) * sqrt(2);
%             end
        end 
%           y = sqrt(2*col) * exp(1i*pi*(0:col-1)/(2*col));
%          y(1) = y(1)/sqrt(2);
        w(1) = w(1) / sqrt(2);
        output = w.*output; 
%         for j = 1: col
%             realPart = real(output(j));
%             imgPart = imag(output(j));
% %             expVal = -1i*(j-1)*pi/(2*col);
%             if (j ==1)
%                 output(j) = sqrt(1/col)*(cos(pi*(j-1)/(2*col))*realPart + sin(pi*(j-1)/(2*col))*imgPart);
% %                 output(j) = sqrt(1/(4*col))*output(j) * exp(expVal)*2;
%             else
%                 output(j) = sqrt(2/col)*(cos(pi*(j-1)/(2*col))*realPart + sin(pi*(j-1)/(2*col))*imgPart);
% %                 output(j) = sqrt(1/(2*col))*output(j) * exp(expVal)*2;
%             end
%         end    
    end
     output = real(output);
end

function output = dct2D(input)
    [row, col] = size(input);
    output = zeros(row,col);
    
    for i = 1 : row
        output(i,:) = dct1D(input(i,:));
    end
    
    for j = 1: col
       output(:,j) = dct1D(output(:,j).');
    end
end

function output = invDct1D(input)
    [~, col] = size(input);
    if (col == 1)
       output = input; 
    
    else
        y = zeros(1,col);
        for idx = 1: col
            expVal = pi*1i*(idx-1)/ (2*col);
            y(idx) = exp(expVal)*sqrt(2*col);
%             if idx == 1
%                 y(idx) = y(idx) * sqrt(2);
%             end
        end 
%           y = sqrt(2*col) * exp(1i*pi*(0:col-1)/(2*col));
%          y(1) = y(1)/sqrt(2);
        y(1) = y(1) / sqrt(2);
        inputVal = y.*input; 
%         inputRev = zeros(1,col-1);
%         for j = 2: col
%             inputRev(col-j+1)=input(j); 
%         end
% 
%         inputTemp = [ input(1), y(2:col).*(input(2:col) - 1i*inputRev) ];
        dftVal = invDft1D(inputVal); % fft not wrong

%         dftVal = invDft1D(inputTemp); % fft not wrong
        output = zeros(1,col);
        output(1:2:col) = dftVal(1:(col/2));
        output(2:2:col) = dftVal(col:-1:(col/2+1));
        
    end
    
    output = real(output);
    
end

function output = invDct2D(input)
    [row, col] = size(input);
    output = zeros(row,col);
    
    for i = 1 : row
        output(i,:) = invDct1D(input(i,:));
%         output(i,:) = IDCTImpl(input(i,:));
%          output(i,:) = idct(input(i,:));       
    end
    
    for j = 1: col
       output(:,j) = invDct1D(output(:,j).');
%        output(:,j) = IDCTImpl(output(:,j).');
%           output(:,j) = idct(output(:,j).');    
    end
end

function output = blocDct2D(input, blkSize)

    [row, col] = size(input);
    output = zeros(row,col);
    nCol = col / blkSize;
    nRow = row/ blkSize;

    for i = 1: nCol
        for j = 1 : nRow
            for a = 1: blkSize
               output(a + blkSize*(j-1),(1+blkSize*(i-1)):(blkSize*i)) = dct1D(input(a + blkSize*(j-1),(1+blkSize*(i-1)):(blkSize*i)));
            end
            for b = 1: blkSize
               output((1+blkSize*(j-1)):(blkSize*j),b + blkSize*(i-1)) = dct1D(output((1+blkSize*(j-1)):(blkSize*j),b + blkSize*(i-1)).');
            end 
        end
    end
end

function output = blocInvDct2D(input, blkSize)

    [row, col] = size(input);
    output = zeros(row,col);
    nCol = col / blkSize;
    nRow = row/ blkSize;

    for i = 1: nCol
        for j = 1 : nRow
            for a = 1: blkSize
               output(a + blkSize*(j-1),(1+blkSize*(i-1)):(blkSize*i)) = invDct1D(input(a + blkSize*(j-1),(1+blkSize*(i-1)):(blkSize*i)));
            end
            for b = 1: blkSize
               output((1+blkSize*(j-1)):(blkSize*j),b + blkSize*(i-1)) = invDct1D(output((1+blkSize*(j-1)):(blkSize*j),b + blkSize*(i-1)).');
            end 
        end
    end
end

%% DHT
function output = dht2D(input)
    [row, col] = size(input);
    output = zeros(row,col);
    
    for i = 1 : row
        output(i,:) = ( HMatrix(row)*(input(i,:)).').' / row;
    end
    
    for j = 1: col
       output(:,j) = HMatrix(row) * output(:,j) / col;
    end
end

function output = invDht2D(input)
    [row, col] = size(input);
    output = zeros(row,col);
    
    for i = 1 : row
        output(i,:) = ( HMatrix(row)*input(i,:).').' * row;
    end
    
    for j = 1: col
       output(:,j) = HMatrix(row) * output(:,j) * col;
    end
end

function output = HMatrix(inputsize)

    H = [1 1 ; 1 -1];
    N = log2(inputsize);
    if N == 2
        output = H;
        return;
    end
    while N > 1
        newMatrix = [H H; H -H];
        H = newMatrix;
        N = N-1;
    end
    output = H;
end
%% Wavelet Transform
function [downLow, downHigh] = dwt1D(input, HV)
    [row, col] = size(input);
    lowOutput = zeros(row, col);
    highOutput = zeros(row, col);
    switch HV
        case 'horizontal'
            lowOutput(1,:) = (input(1,:) + input(row,:))/2;
            highOutput(1,:) = (input(1,:) - input(row,:))/2;
            for idx = 2 : row
                lowOutput(idx,:) = (input(idx,:) + input(idx-1,:))/2;
                highOutput(idx,:) = (input(idx,:) - input(idx-1,:))/2;
            end

            downLow = zeros(row/2,col);
            downHigh = zeros(row/2,col);

            for idx = 1:row/2
               downLow(idx,:) = lowOutput(2*idx,:);
               downHigh(idx,:) = highOutput(2*idx,:);
            end
            
        case 'vertical'
            lowOutput(:,1) = (input(:,1) + input(:,row))/2;
            highOutput(:,1) = (input(:,1) - input(:,row))/2;
            for idx = 2 : row
                lowOutput(:,idx) = (input(:,idx) + input(:,idx-1))/2;
                highOutput(:,idx) = (input(:,idx) - input(:,idx-1))/2;
            end

            downLow = zeros(row,col/2);
            downHigh = zeros(row,col/2);

            for idx = 1:row/2
               downLow(:,idx) = lowOutput(:,2*idx);
               downHigh(:,idx) = highOutput(:,2*idx);
            end            
    end
end

function ouput = invDwt1D(low, high, HV)
    [row, col] = size(low);
    
    switch HV
        case 'horizontal'
            lowOutput = zeros(row*2, col);
            highOutput = zeros(row*2, col);

            for i = 1 : row
                lowOutput = low(i,:);
                highOutput = high(i,:);
            end

            [row, col] = size(lowOutput);
            output = zeros(row,col);

            output(1,:) = (lowOutput(1,:) + lowOutput(row,:)) + (highOutput(1,:) - highOutput(row,:));
            for idx = 2 : row
                output(idx,:) = (lowOutput(idx,:) + lowOutput(idx-1,:)) + (highOutput(idx,:) - highOutput(idx-1,:));

            end
    
        case 'vertical'
            lowOutput = zeros(row, col*2);
            highOutput = zeros(row, col*2);

            for i = 1 : col
                lowOutput = low(:,i);
                highOutput = high(:,i);
            end

            [row, col] = size(lowOutput);
            output = zeros(row,col);

            output(:,1) = (lowOutput(:,1) + lowOutput(:,row)) + (highOutput(:,1) - highOutput(:,row));
            for idx = 2 : col
                output(:,idx) = (lowOutput(:,idx) + lowOutput(:,idx-1)) + (highOutput(:,idx) - highOutput(:,idx-1));

            end
    end
    
end
