function y = feature_descaling(x,xmax,xmin)
    y = x.*(xmax-xmin) + xmin; 
end

