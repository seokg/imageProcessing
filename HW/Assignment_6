function Problem_1()
%% Parallel Beam imag ereconstruction algorithm
    close all;
    clear all;
    % Range: (-1.0, 1.0)
    % number of samples: 256
    % number of views: 180
    % data type : float (4bytes for each sample) 
    
    imgdir = uigetdir('Image Directory');
    file = fopen(fullfile(imgdir,'\shepp_proj_256x180.raw'),'rb');
    
    shepp = fread(file,fliplr([180, 256]),'*float');
    fclose(file);

    figure; imshow(shepp,[]); title('1.1 sigmoid image');    
%%  1.1 Back Projection
    backprojImg = backProj(shepp)
    figure; imshow((backprojImg),[]),title('1.1 back projected image');
%% 1.2 Filter Back Projection (Ram-Lak and Shepp-Logan Filters)

% Ram-Lak
    ramlakImg = filtBackProj(shepp, 'ramlak')
    figure; imshow((ramlakImg),[]),title('1.2 back projected image with ram-lak');
% Sheppen-Logan
    shepploganImg = filtBackProj(shepp, 'shepplogan')
    figure; imshow((shepploganImg),[]),title('1.2 back porjected image with shepp-logane');

%% 1.3 compare the obtained images
% 1.1 vs 1.2

% Ram-Lak vs Sheppen-Logan



end
%% show the following 
% 1.1 back projected image
% 1.2 back projected image with ram-lak
% 1.2 back porjected image with shepp-logan
% 1.3 cut view of back projected image
% 1.3 cut view of back projected image with ram-lack
% 1.2 cut view of back projected image with shepp-logan

%% main functions
function output = backProj(input)
%% back projection witout filter
% input: sigmoid image

    [samples, views] = size(input);
    
    % angle in radian
    theta = 1:1:views;
    theta = theta * (pi / 180);
    
    % size of reconstructed image
    halfSamples = floor(samples/2);
    reconSize = 2 * floor(samples / (2*sqrt(2)));
    halfReconSize = reconSize /2;

    % init back projection image
    output = zeros(reconSize);

   
    
    % setting X Y coordinate 
    x = (1:reconSize) - halfReconSize;
    y = (1:reconSize) - halfReconSize;
    
    % back project from n views
    for i = 1:1:views
        
        pos = zeros(reconSize);
        for u = 1 : reconSize
           for v = 1 : reconSize
               pos(u,v)  = x(u)*cos(theta(i)) + y(v)*sin(theta(i)) + halfSamples;
           end
        end
        output = output + interpolate(1:samples, input(:,i), pos);
    end


output = output * pi / (2*views);
output = output';
end

function output = filtBackProj(input, type)
%% back projection with filter
% input: sigmoid image

    [samples, views] = size(input);
    
    % angle in radian
    theta = 1:1:views;
    theta = theta * (pi / 180);
    
    
    % making it power of 2 zero padding
    newviews = findpower2(views);
    newinput = [input zeros(samples, newviews-views)];    
    % size of reconstructed image
    halfSamples = floor(samples/2);
    reconSize = 2 * floor(samples / (2*sqrt(2)));
    halfReconSize = reconSize /2;
    % Compute the filter
    freqs = -1:2/(samples-1):1;
    tempmyFilter = abs(freqs');
    myFilter = zeros(samples, newviews);
    for i = 1:newviews
        myFilter(:,i) = tempmyFilter;
    end    
    switch type
        case 'ramlak'
            % do nothing
        case 'shepplogan'
            for i = 1:newviews
                myFilter(:,i) = sin(pi*tempmyFilter/ (2)) / pi;
            end    
        otherwise
            error(message('type is wrong... ramlak or shepplogan'))
    end
    H = myFilter;
    
    %FT domain filtering
    ftS = fftShift(dft2D(newinput));
    filteredProj = ftS .* H;
    filteredProj = ifftshift(filteredProj);
    iftS = real(invDft2D(filteredProj));

    
    
    % init back projection image
    output = zeros(reconSize);
   
    % setting X Y coordinate 
    x = (1:reconSize) - halfReconSize;
    y = (1:reconSize) - halfReconSize;
    
    % back project from n views
    for i = 1:1:views
        
        pos = zeros(reconSize);
        for u = 1 : reconSize
           for v = 1 : reconSize
               pos(u,v)  = x(u)*cos(theta(i)) + y(v)*sin(theta(i)) + halfSamples;
           end
        end
        output = output + interpolate(1:samples, iftS(:,i), pos);
    end


output = output * pi / (2*views);
output = output';
end

function output = interpolate(x, y , xprime)
    [row, col] = size(xprime);
    output = zeros(size(xprime));
    for i = 1: row
        for j = 1: col
            idx = findIdx(x,xprime(i,j));
            if idx > 1
                x1 = x(idx);
                x0 = x(idx-1);
                y1 = y(idx);
                y0 = y(idx-1);
            elseif idx <= 1
                idx = 2;
                x1 = x(idx);
                x0 = x(idx-1);
                y1 = y(idx);
                y0 = y(idx-1); 
            elseif idx>length(x)
                idx = length(x);
                x1 = x(idx);
                x0 = x(idx-1);
                y1 = y(idx);
                y0 = y(idx-1); 
            end
            
            
            m = (y1-y0) / (x1-x0);
            output(i,j) = y0 + m * (xprime(i,j) - x0);
            
            
         end
    end
    
end

function output = findIdx(x, input)
    for i = 1 : length(x);
       if x(i) > input
           output = i;
           break;
       end
    end
end

function output = findpower2(input)
    cmp=bitand(input, input -1);
    if cmp == 0
        output = input;
    else
        output = pow2(ceil(log2(input)));
    end
end

function output = findIdx(x, input)
    for i = 1 : length(x);
       if x(i) > input
           output = i;
           break;
       end
    end
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

function output = fftShift(input)
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

%% Filter
% mean filter
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
