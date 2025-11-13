classdef WaveformObj 

    properties 
        wave_name {string} = "Empty"
        wave_id {int8} = 0;
        time {double} = [];
        volts {double} = [];
        dt {double} = 0;
        fs {double} = 0;
    end
    
    methods
        function obj = WaveformObj(wave_id, wave_name, t, v)
            if nargin > 0
                obj.wave_name = wave_name;
                obj.wave_id = wave_id;
                obj.time = t;
                obj.volts = v;
                obj.dt = t(2)-t(1);
                obj.fs = 1./(t(2)-t(1));
            end
        end
                
    end     
end