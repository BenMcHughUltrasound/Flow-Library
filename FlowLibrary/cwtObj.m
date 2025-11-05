classdef cwtObj
    properties
    cwt_id {int8} = 0;
    frequencyArray = [];
    timeArray {double} = [];
    complex_cwt = 0 + 0i;
    cwtAmp {double} = [];

    end

    methods

        function obj = cwtObj(f, t, cfs, cwt_id)

            obj.frequencyArray = f;
            obj.timeArray = t;
            obj.complex_cwt = cfs;
            obj.cwtAmp = abs(cfs);
            if nargin > 3
                obj.cwt_id = cwt_id;
            end
        end

    end

end