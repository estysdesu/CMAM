function [m] = BSplineMod(i, n)
    switch i
        case i <n
            m = i;
        case i == n
            m = 0;
        otherwise
            m = mod(i, n);
    end
end