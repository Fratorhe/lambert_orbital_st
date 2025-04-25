function C = stumpff_C(z)
    if z == 0
        C = 1/2;
    elseif z > 0
        C = (1 - cos(sqrt(z))) / z;
    else
        C = (cosh(sqrt(-z)) - 1) / (-z);
    end
end