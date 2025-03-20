function x8 = resample_to_8bit(x)
    % will resample a vector of x = single(N,1) or double(N,1)
    % to an int8 vector
    [~,indMax] = max(abs(x));
    x8 = x./max(abs(x));
    
    x8         = x8.*(128-(x8(indMax)>0)); 
    x8(x8>127) = 127;
    x8(x8<-128)= -128;
    x8         = int8(x8);
end


