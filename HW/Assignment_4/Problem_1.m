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
% % Outer Product:
% % --------------
% w1 = hamming(15);
% w2 = w1(:) * w1(:).';

% % Rotational
% % ----------
% w1 = hamming(15);
% [x,y] = meshgrid(-7:7);
% r = sqrt(x.^2 + y.^2);
% w = zeros(size(r));
% w(r<=7) = interp1(linspace(-7,7,15),w1,r(r<=7));

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
    
    output = w* w.';

end