classdef BSpline
    properties
        k uint64
        vertices double
        n uint64
        knots double
    end
    methods
        function obj = BSpline(k, vertices, knots)
            if nargin > 0
                obj.k = k;
                obj.vertices = vertices;
                obj.n = length(vertices) - 1;
                obj.knots = knots;
            end
        end
    end≤
    methods (Access = protected)
        function m = ModBSpline(obj, i)
            switch i
            case i < obj.n
                m = i;
            case i == obj.n
                m = 0;
            otherwise
                m = mod(i, obj.n);
            end
        end
    end
end