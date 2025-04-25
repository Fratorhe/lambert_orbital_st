function S = stumpff_S(z)
    if z == 0
        S = 1/6;
    elseif z > 0
        S = (sqrt(z) - sin(sqrt(z))) / (z^(3/2));
    else
        S = (sinh(sqrt(-z)) - sqrt(-z)) / ((-z)^(3/2));
    end
end