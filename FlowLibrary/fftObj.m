classdef fftObj
    properties
        fft_name {string} = "Empty"
        fft_id {int8} = 0
        f {double} = []
        c  = 0 + 0i
        A {double} = []
        phase {double} = []
        numPoints {int32} = 0
        df {double} = 1
        f_max {double} = 0
        originalLength int32 = 0
    end

    methods
        function obj = fftObj(fft_name, fft_id, f, c, L)
            if nargin > 4
                obj.originalLength = L;
                obj.fft_name = fft_name;
                obj.fft_id = fft_id;
                obj.f = f;
                obj.c = c;
                obj.A = abs(c);
                obj.phase = angle(c);
                obj.numPoints = length(f);
                obj.df = f(2)-f(1);
                obj.f_max = max(f);
            elseif nargin > 0
                obj.fft_name = fft_name;
                obj.fft_id = fft_id;
                obj.f = f;
                obj.c = c;
                obj.A = abs(c);
                obj.phase = angle(c);
                obj.numPoints = length(f);
                obj.df = f(2)-f(1);
                obj.f_max = max(f);
                obj.originalLength = obj.numPoints;
            end
            
        end
            
    end   
end