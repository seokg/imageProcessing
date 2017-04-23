function [ PSNR ] = findPNSR( ref, img )

    MSE = 1/(size(ref,2)*size(ref,1)*3)* sum(sum(sum((ref-img).^2)));
    PSNR = 10 * log10(255^2 / MSE);

end

