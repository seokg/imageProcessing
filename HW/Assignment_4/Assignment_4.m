clear
clc

% % assignment 4
% % 20164332 ¼­±¤±Õ



file = fopen(fullfile('Test_images\Gray_baby_512x512.raw'),'rb');
gray_image = fread(file,fliplr([512,512]),'*uint8')';
figure; imshow(gray_image,[]); title('original image');


dctImage = dct2(gray_image);
fftImage = fft2(gray_image);

magnitude = log(1+abs(dctImage));
magnitude2 = log(1+abs(fftImage));



    [row, col] = size(magnitude2);
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
            
            temp(idxRow, idxCol) = magnitude2(i,j);
            
            
        end
    end
    
    magnitude2 = temp;

figure; imshow (magnitude,[]); title('magnitude');
figure; imshow (magnitude2,[]); title('magnitude');



