% Checks for connected picoscope and clears all variables
if exist('device') == 1
pico_close(device);
Tlogger = [];
end
clear;
clc;
close all;

%%
%constants for analysis
D = 5.528; %Pipe internal diameter in cm
TR_spacing = 3.388599877; %transducer spacing in cm
P_chordal = 11.54456; %Chordal path length in cm
P_Vpath = 11.544; %V path length in cm
cs = c_water(21); %sound speed in water

%set up device
device = pico_setup('IZ186/0055');

%enable channel A with a +/- maxVolts y range
maxVolts = 1;
pico_enable_channel(device, 0, maxVolts,0);

%set up trigger on channel D
pico_set_trigger(device, 3);

% setup timebase parameters for passing into data collection functions
[interval, maxSamples] = pico_timebase(device);

% set time to capture up to in microseconds
timeToCapture = 150;

%set number of averages (max on picoscope is 64, for higher change
%num64Avg) and create waveforms array to store data
numAvg = 64;
num64Avg = 2;
RecordedWaveforms = WaveformObj.empty();
scope = Wave_Lib();

%enable switcher and connect
dev = serialport("COM3", 115200);
configureTerminator(dev, 'LF');

% Return switcher to known position
writeline(dev, '1')
StartTimer = datetime;
%%
i = 0;
numSmooth = 5;
h = figure(1);
forever = 1; 


while forever 
    i = i+1;
    
    %small component that makes the loop run until a key is pressed
    drawnow
    isKeyPressed = ~isempty(get(h,'CurrentCharacter'));
    if isKeyPressed
        pause(1)
    toStop = input('Program stopped. Press c to continue, and s to stop.');
    if toStop == 'c'
        isKeyPressed = false;
    else
        break
    end
    end
    
    %Collect data from picoscope and concatenate to save array
    for k = 1:num64Avg
        [time, chA(:,k), ~] = pico_averaged_waveform(device, numAvg, timeToCapture*1E-6, interval, maxSamples);
        writeline(dev, '2');
        [~, chB(:,k), ~] = pico_averaged_waveform(device, numAvg, timeToCapture*1E-6, interval, maxSamples);
        writeline(dev, '1');
    end
    chA = mean(chA, 2);
    chB = mean(chB, 2);
    volts = cat(2, chA, chB);

    RecordedWaveforms(i) = WaveformObj(i,i,time, volts);

    CutWaveforms(i) = scope.cutWaveObj(RecordedWaveforms(i), 70e-6, timeToCapture);
    [starttimes(i,:), endtimes(i,:)] = scope.DetectArrival(CutWaveforms(i), 0.05);
    WindowedWaveforms(i) = scope.windowWaveObj(CutWaveforms(i), starttimes, endtimes, 0);
    SmoothedWindowedWaveforms(i) = scope.smoothWaveObj(WindowedWaveforms(i));
    [lags(i), correlation(i)] = scope.CalcLag(SmoothedWindowedWaveforms(i), 0.5e-6, 1000);
    [v(i),q(i)] = scope.FlowRate(lags(i), 21, 1, D, 0);
    Tcurrent(i) = datetime-StartTimer;

    %plot and print results
    figure(1);
    
    subplot(2,1,2, "replace")
    hold on
        plot(Tcurrent,lags/1e-9, linestyle="none", marker = 'x', Color='blue', MarkerSize=7)
        %plot(lags*1e9, color = 'blue', Marker='o')
        if i > numSmooth
            dt_err(i) = std(lags(i-numSmooth:i));      
            avdt = movmedian(lags, [numSmooth,0]);
            errorbar(Tcurrent, avdt/1e-9,dt_err/1e-9, color = 'red', Marker='o', MarkerSize=7)
        end
        grid("on")
        xlabel('Time (hrs:min:sec)')
        
        ylabel('Transit time difference (ns)')
    hold off
    
    subplot(2,1,1, "replace")
    hold on
        plot(Tcurrent, q, linestyle="none", marker = 'x', Color='blue', MarkerSize=7)
        %plot(lags*1e9, color = 'blue', Marker='o')
        if i > numSmooth
        v_err(i) = std(v);
        q_err(i) = std(q(i-numSmooth:i));
        avQ = movmedian(q, [numSmooth,0]);
        errorbar(Tcurrent, avQ,q_err, color = 'red', Marker='o', MarkerSize=7)
        end
        grid("on")
        xlabel('Time (hrs:min:sec)')
        ylabel('Flowrate (L/min)')
    hold off
    scope.plt(SmoothedWindowedWaveforms(i),2, 1)
end
pico_close(device)
clear("device")



%%

scope.SaveWaveforms(RecordedWaveforms, 'testData.csv')
