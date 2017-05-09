function Problem_2()

    %% loading directory for image files
    imgdir= uigetdir('Image Directory');
    file = fopen(fullfile(imgdir,'\Gray_baby_512x512.raw'),'rb');
    gray_image = fread(file,fliplr([512,512]),'*uint8')';
    fclose(file);
   
    file = fopen(fullfile(imgdir,'\Color_baby_512x512.raw'),'rb');
    color_image = fread(file,fliplr([512,512*3]),'*uint8')';
    fclose(file);   
    %% -------------- insert code below --------------%%

    
    grayHistEq = histEqualizer(gray_image);
    figure, imshow(grayHistEq,[]), title('2.1 a Histogram equalization Gray Image')
    colorHistEq = histEqualizerColor(color_image);
    figure, imshow(colorHistEq,[]), title('2.1 a Histogram equalization RGB Image')
%     hsiHistEq = histEqualizerHSI(color_image)
%     figure, imshow(hsiHistEq,[]), title('2.1 a Histogram equalization HSI Image')
% 
%     figure, imshow(colorHistEq,[]), title('2.1 a Histogram equalization CMY Image')

    alpha = 0.5;
    grayHistModi = histModification(gray_image, alpha);
    figure, imshow(grayHistModi, []), title('2.1 b Histogram Modification Gray Image')
    
    gamma = 2;
    gammaCorrection(color_image, gamma);
    figure, imshow(grayHistModi, []), title('2.2 a Gamma Correction Color Image')
    
    grayguasNois = guasNoise(gray_image);
    graysaltNois = saltpepperNoise(gray_image);
    figure, imshow(grayguasNois, []), title('2.2 b image with gaussian noise')
    figure, imshow(graysaltNois, []), title('2.2 b image with salt pepper noise')
    
    meansalt3 = meanfilter(graysaltNois, 3);
    meansalt7 = meanfilter(graysaltNois, 7);
    meangaus3 = meanfilter(grayguasNois, 3);
    meangaus7 = meanfilter(grayguasNois, 7);
    figure, imshow(meansalt3, []), title('2.2 b Mean filtered Salt image (3x3 kernel)')
    figure, imshow(meansalt7, []), title('2.2 b Mean filtered Saltimage (7x7 kernel)')    
    figure, imshow(meangaus3, []), title('2.2 b Mean filtered Guassian image (3x3 kernel)')
    figure, imshow(meangaus7, []), title('2.2 b Mean filtered Guassian image (7x7 kernel)')    
    
    mediansalt3 = medianfilter(graysaltNois, 3);
    mediansalt7 = medianfilter(graysaltNois, 7);
    mediangaus3 = medianfilter(grayguasNois, 3);
    mediangaus7 = medianfilter(grayguasNois, 7);
    figure, imshow(mediansalt3, []), title('2.2 b Median filtered Salt image (3x3 kernel)')
    figure, imshow(mediansalt7, []), title('2.2 b Median filtered Salt image (7x7 kernel)')    
    figure, imshow(mediangaus3, []), title('2.2 b Median filtered Guassian image (3x3 kernel)')
    figure, imshow(mediangaus7, []), title('2.2 b Median filtered Guassian image (7x7 kernel)')        
end
%% histogram equalization 
function output = histEqualizer(input)
    [row, col] = size(input);
    
    pixelTotal = row*col;   
    % compute pdf
    pdf = zeros(256,1);
    for i = 1: row
        for j = 1: col
            val = input(i,j);
            pdf(val+1) = pdf(val+1) + 1;
        end
    end
    pdf = pdf/ pixelTotal;
    
    % compute cdf
    cdf= zeros(256,1);
    cdf(1) = pdf(1);
    for i = 2: size(pdf,1)
        cdf(i) = pdf(i) + cdf(i-1);

    end
    
        normalizedCDF = round ((cdf - min(cdf)) / (1 - min(cdf)) * 255 + 0.5);
    output = size(input);
    for i = 1:row
        for j = 1:col
            output(i,j) = normalizedCDF(input(i,j)+1);
        end  
    end
    

end

function output = histEqualizerColor(input)
    [row, col] = size(input);
    
    % seperate RGB
    R = zeros(row, col/3);
    G = zeros(row, col/3);
    B = zeros(row, col/3);
    
    for i = 1 : col/3
        for j = 1 : row
            R(j,i) = input(j,3*(i-1) + 1);
            G(j,i) = input(j,3*(i-1) + 2);
            B(j,i) = input(j,3*(i-1) + 3);
        end
    end
    
    outputR = histEqualizer(R);
    outputG = histEqualizer(G);
    outputB =histEqualizer(B);
    
    output = zeros(row, col/3 , 3);
    
    output(:,:,1) =outputR;
    output(:,:,2) =outputG;
    output(:,:,3) =outputB;
    
    output = output / 255;
    
    output1 = zeros(row, col/3 , 3);
    
    output1(:,:,1) =R;
    output1(:,:,2) =G;
    output1(:,:,3) =B;
    
    output1 = output1/255;
    figure, imshow(output, [])
    figure, imshow(output1, [])
end

function output = histEqualizerHSI(input)
    [row, col] = size(input);
    
    % seperate RGB
    R = zeros(row, col/3);
    G = zeros(row, col/3);
    B = zeros(row, col/3);
    
    for i = 1 : col/3
        for j = 1 : row
            R(j,i) = input(j,3*(i-1) + 1);
            G(j,i) = input(j,3*(i-1) + 2);
            B(j,i) = input(j,3*(i-1) + 3);
        end
    end
    % convert to HSI
    [H S I]= rgb2HSI( R,G,B )
    outputH = histEqualizer(floor(H*255))/255;
    outputS = histEqualizer(floor((S+0.5)*255))/255-0.5;
    outputI =histEqualizer(floor(I*255))/255;
    
   [newR newG newB ] = HSItoRGB( outputH, outputS, outputI )   
    
    output = zeros(row, col/3 , 3);
    
    output(:,:,1) =floor(newR);
    output(:,:,2) =floor(newG);
    output(:,:,3) =floor(newB);
    
    output = output / 255;
    
    figure, imshow(output, [])
end

% HSI RGB conversion
function [ H, S, I ] = rgb2HSI( R,G,B )

    I = 1/3 .*(R+G+B);

    S=1- (3./(R+G+B)+0.000001).*(min(min(R,G),B));

    H = acosd ((2*R-G-B) ./(2*sqrt( ( (R-G).^2 + (R-B).*(G-B) ) )+0.000001)); % 1/(2*pi)

    H(B>G)=360-H(B>G);
    H=H/360;
end

function [ newR newG newB ] = HSItoRGB( H, S, I )



    Chroma = S .* I;
    Hdash = H*360.0 / 60.0;
    X = Chroma .* (ones(size(Hdash)) - abs(rem(Hdash,2.0) -ones(size(Hdash))));
    Min = I - Chroma;

    newR = zeros(size(H));
    newG = zeros(size(H));
    newB = zeros(size(H));

    newR = newR +(Hdash < 1.0).*Chroma ; 
    newG = newG +(Hdash < 1.0).*X ;

    newR = newR +(Hdash>= 1.0 & Hdash < 2.0).*X ; 
    newG = newG +(Hdash>= 1.0 & Hdash< 2.0).*Chroma ;

    newG = newG +(Hdash >= 2.0 & Hdash < 3.0).*Chroma ; 
    newB = newB +(Hdash >= 2.0 & Hdash < 3.0).*X ;

    newG = newG +(Hdash >= 3.0 & Hdash < 4.0).*X ; 
    newB = newB +(Hdash >= 3.0 & Hdash < 4.0).*Chroma ;

    newR = newR +(Hdash >= 4.0 & Hdash < 5.0).*X ; 
    newB = newB +(Hdash >= 4.0 & Hdash < 5.0).*Chroma ;

    newR = newR+(Hdash >= 5.0 & Hdash < 6.0).*Chroma ; 
    newB = newB +(Hdash >= 5.0 & Hdash < 6.0).*X ;

    newR = newR + Min;
    newG = newG + Min;
    newB = newB + Min;

end


%% histogram modification

function output = histModification(input, alpha)
    [row, col] = size(input);
    
    pixelTotal = row*col;   
    % compute pdf
    pdf = zeros(256,1);
    for i = 1: row
        for j = 1: col
            val = input(i,j);
            pdf(val+1) = pdf(val+1) + 1;
        end
    end
    pdf = pdf/ pixelTotal;
    
    % compute cdf
    cdf= zeros(256,1);
    cdf(1) = pdf(1);
    for i = 2: size(pdf,1)
        cdf(i) = pdf(i) + cdf(i-1);
    end
    
    normalizedCDF =  round((1 - exp(-(cdf.^2 / (2*alpha^2))))*255);
    output = size(input);
    for i = 1:row
        for j = 1:col
            output(i,j) = normalizedCDF(input(i,j)+1);
        end  
    end
    
end


%% gamma correction

function output = gammaCorrection(input, gamma)

    gammaCorrection = 1/ gamma;
    
    [row, col] = size(input);
    
    % seperate RGB
    R = zeros(row, col/3);
    G = zeros(row, col/3);
    B = zeros(row, col/3);
    
    for i = 1 : col/3
        for j = 1 : row
            R(j,i) = input(j,3*(i-1) + 1);
            G(j,i) = input(j,3*(i-1) + 2);
            B(j,i) = input(j,3*(i-1) + 3);
        end
    end
    
    
    newR = 255*(R./255).^gammaCorrection;
    newG = 255*(G./255).^gammaCorrection;
    newB = 255*(B./255).^gammaCorrection;

    
   output = zeros(row, col/3 , 3);
    
    output(:,:,1) =newR;
    output(:,:,2) =newG;
    output(:,:,3) =newB;
    
    output = output./255;
    figure, imshow(output, [])

end



%% Noise

function output = guasNoise(input)
  mean = 0;     % default mean
  variance = 0.01;  % default variance
  output = double(input) + sqrt(variance)*randn(size(input));  
end

function output = saltpepperNoise(input)
    d = 0.01;   % default density
    output = input;
    x = rand(size(input));
    output(x < d/2) = 0; % Minimum value
    output(x >= d/2 & x < d) = 1; % Maximum (saturated) value
    
    
end

%% mean filter
function output = meanfilter(input, kernelSize)
    matrixSize = kernelSize^2;
    kernel = ones(kernelSize, kernelSize) ./ matrixSize;
    output = colvolution(input,kernel);

end


function output = colvolution(input, kernel)
    % invert
    [row, col] = size(kernel);
    temp = zeros(row,col);
    for i = 1: row
        for j = 1: col
           temp(i,j) = kernel(row - i +1, col - j +1);
        end 
    end
    kernel = temp;
    
    [row1, col1] = size(input);
    a  = ceil(row /2);
    output = zeros(row1, col1);
    for i = 1 :row1
        for j = 1: col1
            for k = 1:row
                for l = 1:col
                    if (i-a+k >0) &  (i-a+k <= row1) &  (j-a+l >0) & (j-a+l <= col1)
                        output(i,j) = output(i,j) + input(i-a+k,j-a+l) * kernel(k,l);
                    end
                end
            end
            
        end 
    end
     


end
%% median filter

function output = medianfilter(input, kernelSize)
    matrixSize = kernelSize^2;
    kernel = ones(kernelSize, kernelSize);
    output = colvolution4median(input,kernel);

end

function output = colvolution4median(input, kernel)
    % invert
    [row, col] = size(kernel);
    temp = zeros(row,col);
    for i = 1: row
        for j = 1: col
           temp(i,j) = kernel(row - i +1, col - j +1);
        end 
    end
    kernel = temp;
    
    [row1, col1] = size(input);
    a  = ceil(row /2);
    temp = zeros(row*col,1);
    output = zeros(row1, col1);
    for i = 1 :row1
        for j = 1: col1
            
            for k = 1:row
                for l = 1:col
                    if (i-a+k >0) &  (i-a+k <= row1) &  (j-a+l >0) & (j-a+l <= col1)
                        temp((k-1)*col+l,1) =  input(i-a+k,j-a+l);
                    else
                        notsort = 1;
                    end
                end
            end
            
            
            if notsort == 1
                output(i,j) = input(i,j);
            else
                sortMatrix = sort(temp);
                output(i,j)= sortMatrix(a);
                              
            end
            notsort = 0 ;
            temp = zeros(row*col,1);  

            
        end 
    end
 end
%% directional filite