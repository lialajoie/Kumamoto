function [Hurst,Roughness,Rsq2,slope,x_1,yCalc] = fit_lin2log(X,Y,limit)
       
    clear x_1
    clear y_1
    clear X_1
    
    % limit to wavelengths/frequencies < fault length
    if length(limit) == 1
        xfind = X <= log10(limit);
    elseif length(limit) == 2
        xfind = X >= log10(limit(1)) & X <= log10(limit(2));
    end
    x_1 = X(xfind);
    y_1 = Y(xfind);
    X_1 = [ones(length(x_1),1) x_1'];
    
    size_x = size(X_1);
    size_y = size(y_1);
    if size_x(1) ~= size_y(1)
        y_1 = y_1';
    end

    b_1 = X_1\y_1;
    slope = b_1(2);
    yCalc = X_1*b_1;
    Rsq2 = 1 - sum((y_1 - yCalc).^2)/sum((y_1 - mean(y_1)).^2);
    Roughness = (5-abs(b_1(2)))/2;
    Hurst = (abs(b_1(2))-1)/2;
      
end
