
classdef Wave_Lib

    properties
        WaveStorage = WaveformObj.empty;
        FFTStorage = fftObj.empty;
        
    end



    methods

        function obj = Wave_Lib()
            %Returns a Wave_Viewer object with an internal WaveStorage variable
            obj.WaveStorage = WaveformObj.empty(1,0);
            obj.FFTStorage = fftObj.empty(1,0);
        end


    end

    methods (Static)

        function help()
        %HELP Displays usage information for the Wave_Lib class
        fprintf('\nWave_Lib Class Help\n');
        fprintf('-------------------\n');
        fprintf('Wave_Lib provides utilities for handling waveform and FFT objects.\n\n');

        fprintf('Constructor:\n');
        fprintf('  Wave_Viewer() - Initializes WaveStorage and FFTStorage as empty arrays.\n\n');

        fprintf('Static Methods:\n');
        fprintf('  plt(WaveObj, figNum) - Plots time-domain waveforms from a WaveformObj.\n');
        fprintf('  pltFFT(fftObj, figNum) - Plots frequency-domain FFT from an fftObj.\n');
        fprintf('  LoadWaveforms(name_or_matrix) - Loads waveform data from a file or matrix.\n');
        fprintf('  SaveWaveforms(WaveArrayToSave, fileToSaveTo) - Saves waveform data to a file (not yet implemented).\n');
        fprintf('  WaveformObjToMat(waveObj) - Converts a WaveformObj to a matrix [time, ch1, ch2].\n');
        fprintf('  MatToWaveformObj(Mat) - Converts a matrix to a WaveformObj.\n');
        fprintf('  MatToFFTObj(Mat) - Converts a matrix to an fftObj.\n');
        fprintf('  fftObjToMat(fftObj) - Converts an fftObj to a matrix [frequency, complex amplitudes].\n');
        fprintf('  fftObj(waveObj, fft_name, fft_id, n) - Generates an fftObj from a WaveformObj.\n');
        fprintf('  IFFTObj(fftObj) - Reconstructs a WaveformObj from an fftObj using inverse FFT.\n\n');

        fprintf('Usage Example:\n');
        fprintf('  wave = Wave_Lib.LoadWaveforms("data.csv");\n');
        fprintf('  Wave_Lib.plt(wave(1), 1);\n');
        fprintf('  fft = Wave_Lib.fftObj(wave(1), 1, 1, 1024);\n');
        fprintf('  Wave_Lib.pltFFT(fft, 2);\n');
        fprintf('  wave_reconstructed = Wave_Lib.IFFTObj(fft);\n');
        fprintf('  Wave_Lib.plt(wave_reconstructed, 3);\n');
    end
%{
Plotting functions for plotting and formatting WaveObjs and fftObjs
%}
        function plt(WaveObj, figNum) 
            %{This just plots all waves in a WaveObj on figure figNum and formats the axes
            %with the right labels.
            %}

            figure(figNum) %move to figure figNum

            %loop through all of the volts data in a WaveformObj to plot
            hold on
            for i = 1 : length(WaveObj.volts(1,:))
                plot(WaveObj.time./1e-6,WaveObj.volts(:,i))
            end
            hold off

            fig = gcf; % get current figure data to pass to standardization

            %format axes
            grid("on")
            grid("minor")
            xlim([min(WaveObj.time)/1e-6, max(WaveObj.time)/1e-6])
            xlabel('Time (\mus)')
            ylabel('Amplitude (V)')
            title_text = sprintf('Time domain plot of WaveformObj %d', (WaveObj.wave_name));
            title(title_text)

            STANDARDIZE_FIGURE(fig);

        end
        
        function pltFFT(fftObj, figNum)
            %{Plots an fftObj and labels axes correctly
            figure(figNum)
            hold on
            for i = 1 : length(fftObj.c(1,:))
                           
                Amp = fftObj.A(:,i);
                f_to_plot = fftObj.f;
                plot(f_to_plot/1e6, Amp)
            end
            hold off
                xlim([0, 2])
                grid('on')
                grid('Minor')
                xlabel('Frequency (MHz)')
                ylabel('|FFT(Y)| (V/MHz)')
                title_text = sprintf('Single sided amplitude FFT of fftObj %d', fftObj.fft_name);
                title(title_text)
                fig = gcf;
                STANDARDIZE_FIGURE(fig);
            
        end
       
%{
Load and save functions
%}
        function LoadWaveform = LoadWaveforms(name_or_matrix)
            %{Converts a txt or csv file into a WaveformObj array to
            %perform processing on. Can take array inputs also
            %}
            disp('Loading...')

            %Check if input name_or_matrix is a name or a local variable
            %matrix
            if isstring(name_or_matrix) | ischar(name_or_matrix)
                RawLoadingIn = readmatrix(name_or_matrix);
            else
                RawLoadingIn = name_or_matrix;
            end
            
           disp('Files loaded! Allocating WaveformObj...')

            t = ~any(isnan(RawLoadingIn),2); %check for nan values in file or array
        
            time = RawLoadingIn(t, 1);  %seperate time array and filter nans out     
            AllVolts = RawLoadingIn(t, 2:end); % filter nans out of voltage data

            numMeasures = floor(length(AllVolts(1,:))/2); % calculate number of upstream downstream pairs
            LoadWaveform = WaveformObj.empty(numMeasures, 0); %allocate empty WaveformObj array for speed
            dt = time(2)-time(1); 
            fs = 1./dt;

            if numMeasures == 0
                %handle case of single waveform
                LoadWaveform = WaveformObj(1,1,time, AllVolts, dt, fs);    

            else
                %handle all other cases
                for ii = 1:numMeasures
                    fprintf('\n Allocating WaveformObj %d / %d... \n', ii, numMeasures)
     
                    wave_name = ii; %set to idx as default
                    wave_id = ii;
                    LoadWaveform(ii) = WaveformObj(wave_name, wave_id,time, AllVolts(:,ii:ii+1)); %create WaveformObj and allocate
                end
            end
        end

        function SaveWaveform = SaveWaveforms(WaveArrayToSave, fileToSaveTo)
            % figure this out
        end
%{
Helper functions that convert WaveformObjs and fftObjs to matrices and vice
versa.
%}
        function Mat = WaveformObjToMat(waveObj)
            %converts a WaveformObj to an nx3 array with columns time, ch1,
            %ch2
            Mat = cat(2,waveObj.time, waveObj.volts);
        end

        function waveObj = MatToWaveformObj(Mat)
            time = Mat(:, 1);  %seperate time array and filter nans out     
            AllVolts = Mat(:, 2:end); % seperate voltage data

            %create waveObj
            waveObj = WaveformObj(1, 1,time, AllVolts); %create WaveformObj and allocate
        end
        
        function FFTObj = MatToFFTObj(Mat)
            f = Mat(:,1);
            c = Mat(:,2:end);
            FFTObj = fftObj(1,1,f,c);
        end
        
        function Mat = fftObjToMat(fftObj)
            Mat = [fftObj.f, fftObj.c];
        end

%{
Functions to convert WaveformObjs to fftObjs and vice versa
%}
        function SingleSidedFFT = fftObj(waveObj, fft_name, fft_id, n)
            %{generates an fftObj from a WaveformObj.
            % Arguments: waveObj {WaveformObj} - a WaveformObj to
            %            transform into a fftObj.
            %            fft_id {double} - an id to assign to
            %            SingleSidedFFT
            %            n {int} - number of fft points
            %}

            L = length(waveObj.time);

            if nargin < 4
                n = 2^nextpow2(L);
            else
                n = 2^nextpow2(n);
            end %{ sets number of fft points to next power of 2 after input. If no input is given,
                %  n is set to the next powe%}

            Y = fft(waveObj.volts,n,1);

            c = Y(1:n/2+1,:);
            c(2:end-1,:) = 2*c(2:end-1,:);

            f = transpose(waveObj.fs/n*(0:(n/2)));
       
            SingleSidedFFT = fftObj(fft_name, fft_id, f, c, L);

        end

        function waveObj = IFFTObj(fftObj)
            % Performs ifft on a fftObj, preserving the original length of
            % the signal it is calculated from

            fftMat = Wave_Lib.fftObjToMat(fftObj); %convert to matrix

            fftAmp = fftMat(:,2:end); %seperate complex amplitudes

            N = 2*length(fftAmp(:,1)); %calculate length of full twosided fft spectrum

            fftFull = nan(N-2, length(fftAmp(1,:))); %allocates array for full twosided fft (N-2 to exclude DC and Nyquist)

            signal_complex = nan(N, length(fftAmp(1,:))); %allocates array for ifft results

            %loops through ffts in fftObj
            for ii = 1:length(fftAmp(1,:))
            fftFull(:,ii) = [fftAmp(:,ii); conj(fftAmp(end-1:-1:2,ii))]; %Convert to twosided spectrum
            signal_complex(:,ii) = ifft(fftFull(:,ii), N); %perform ifft
            end

            fs = N*fftObj.df; % reconstruct fs (df = fs/N from fft theory)
            
            dt = 1./fs; %compute dt for WaveformObj output

            N_actual = double(fftObj.originalLength); %get original signal length from fftObj
            time = transpose(dt.*(0:(N_actual-1))); %reconstruct time array from dt and orignal length
            signal_complex = signal_complex(1:N_actual,:); %crop full signal to remove zero padding

            % Return real signal (assume original was real-valued)
            waveObj = WaveformObj(fftObj.fft_name, fftObj.fft_id,time, real(signal_complex));
                


        end
%{
Functions to perform CWT and plot it
%}
        function cwtObjs = cwtObj(WaveObj, freqLimLower, freqLimHigher)
            time = WaveObj.time;
            volts = WaveObj.volts;
            fb = cwtfilterbank(SignalLength=length(time), ...
                SamplingFrequency=WaveObj.fs, ...
                VoicesPerOctave=48, ...
                FrequencyLimits=[freqLimLower freqLimHigher], ...
                Wavelet="Morse");
            for ii = 1:length(volts(1,:))
            [cfs(:,:,ii), frq(:,ii)] = cwt(volts(:,ii), FilterBank = fb);
            cwtObjs(ii) = cwtObj(frq(:,ii), time, cfs(:,:,ii), ii);
            end
        end

        function pltcwt(cwtObj, figNum)
            figure(figNum)
            surface(cwtObj.timeArray,cwtObj.frequencyArray,cwtObj.cwtAmp)
            axis tight
            shading flat
            xlabel("Time (s)")
            ylabel("Frequency (Hz)")
            text3 = sprintf('Amplitude Scaleogram of Signal %d', i);
            title(text3)
            colorbar();
            fig = gcf;
            STANDARDIZE_FIGURE(fig);

        end


        function cutWaveObj = cutWaveObj(waveObj, starttime, endtime)
            time = waveObj.time;
            volts = waveObj.volts;
            signalTime = find(time> starttime &  time< endtime);
            cut_time = time(signalTime);
            cut_volts = volts(signalTime,:);
            cutWaveObj = WaveformObj(1,1, cut_time, cut_volts);
        end


        function windowedWaveObj = windowWaveObj(waveObj, starttime, endtime, windowConst)
            if nargin < 4
                windowConst = 0;
            end
            time = waveObj.time;
            volts = waveObj.volts;
            
            signalTime = find(time> starttime &  time< endtime);
            L = length(signalTime);
            window = zeros(1,length(time));
            window(1,signalTime) = tukeywin(L,windowConst);

            for ii = 1:length(volts(1,:))
            windowedVolts = transpose(window).*volts(:,ii);
            end
            windowedWaveObj = WaveformObj(1,1,time, windowedVolts);
        end


        function interpWaveObj = InterpWave(waveObj, numInt)
            %inputs: A waveformObj with m voltage and a number of points to
            %interpolate per point in its corresponding time array. Returns the interpolated time and amplitude
            %packaged identically to the input.
                time = waveObj.time;
                volts = waveObj.volts;

                xq = time(1):waveObj.dt/numInt:time(end);
                for ii = 1:length(volts(1,:)) %loop across voltage columns
                intWave(:,ii) = transpose(interpn(time, volts(:,ii), xq, 'cubic'));
                end
                interpWaveObj = WaveformObj(1,1,transpose(xq), intWave);
            
        end

        
    end
end

