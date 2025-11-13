classdef Wave_Lib
    %WAVE_LIB
    %   Wave_Lib is a MATLAB class that provides a library of static and
    %   instance methods for working with waveform, FFT, and CWT objects.
    %
    %   It supports:
    %       • Loading, saving, and converting WaveformObj / fftObj data
    %       • Time-domain and frequency-domain plotting
    %       • FFT and inverse FFT operations
    %       • Continuous Wavelet Transforms (CWT)
    %       • Waveform cutting, windowing, smoothing, and interpolation
    %       • Flow and velocity calculations
    %       • Signal arrival detection and thresholding
    %
    %   Dependencies:
    %       WaveformObj, fftObj, cwtObj, and STANDARDIZE_FIGURE must exist
    %       in your MATLAB path.

    %% ====================================================================
    %  Properties
    % =====================================================================
    properties
        % WaveStorage  - container for waveform objects
        WaveStorage = WaveformObj.empty;
        % FFTStorage   - container for FFT objects
        FFTStorage = fftObj.empty;
    end

    %% ====================================================================
    %  Constructor
    % =====================================================================
    methods
        function obj = Wave_Lib()
            %WAVE_LIB Constructor
            %   Initializes the class with empty WaveformObj and fftObj
            %   containers.
            obj.WaveStorage = WaveformObj.empty(1,0);
            obj.FFTStorage  = fftObj.empty(1,0);
        end
    end

    %% ====================================================================
    %  Static Methods
    % =====================================================================
    methods (Static)

        function help()
        %HELP Displays detailed usage information for the Wave_Lib class
        %
        %   Wave_Lib.help()
        %
        %   Shows available methods grouped by functionality category.

            fprintf('\nWave_Lib Class Help\n');
            fprintf('===================\n');
            fprintf(['Wave_Lib provides tools for waveform processing, FFT analysis,\n' ...
                     'plotting, filtering, interpolation, and flow calculations.\n\n']);

            fprintf('Constructor:\n');
            fprintf('  Wave_Lib() - Initializes WaveStorage and FFTStorage as empty arrays.\n\n');

            % --------------------------------------------------------------
            fprintf('------------------------------------------\n');
            fprintf('PLOTTING FUNCTIONS\n');
            fprintf('------------------------------------------\n');
            fprintf('  plt(WaveObj, figNum, [ClearFlag])   - Plots time-domain waveform(s).\n');
            fprintf('  pltFFT(fftObj, figNum)              - Plots FFT amplitude spectra.\n');
            fprintf('  pltcwt(cwtObj, figNum)              - Plots CWT amplitude scalogram.\n\n');

            fprintf('------------------------------------------\n');
            fprintf('LOADING AND SAVING\n');
            fprintf('------------------------------------------\n');
            fprintf('  LoadWaveforms(name_or_matrix)       - Loads waveform data from file or matrix.\n');
            fprintf('  SaveWaveforms(WaveArray, file)      - Saves WaveformObj array to file.\n\n');

            fprintf('------------------------------------------\n');
            fprintf('OBJECT <-> MATRIX CONVERSION\n');
            fprintf('------------------------------------------\n');
            fprintf('  WaveformObjToMat(waveObj)           - Converts WaveformObj to numeric matrix.\n');
            fprintf('  MatToWaveformObj(Mat)               - Converts numeric matrix to WaveformObj.\n');
            fprintf('  MatToFFTObj(Mat)                    - Converts numeric matrix to fftObj.\n');
            fprintf('  fftObjToMat(fftObj)                 - Converts fftObj to numeric matrix.\n\n');

            fprintf('------------------------------------------\n');
            fprintf('FFT AND SIGNAL RECONSTRUCTION\n');
            fprintf('------------------------------------------\n');
            fprintf('  fftObj(waveObj, name, id, n)        - Computes FFT from waveform.\n');
            fprintf('  IFFTObj(fftObj)                     - Reconstructs waveform from FFT.\n\n');

            fprintf('------------------------------------------\n');
            fprintf('CWT (CONTINUOUS WAVELET TRANSFORM)\n');
            fprintf('------------------------------------------\n');
            fprintf('  cwtObj(WaveObj, fLow, fHigh)        - Computes Continuous Wavelet Transform.\n');
            fprintf('  pltcwt(cwtObj, figNum)              - Plots amplitude scaleogram of a CWT object.\n\n');

            fprintf('------------------------------------------\n');
            fprintf('WAVEFORM EDITING AND PROCESSING\n');
            fprintf('------------------------------------------\n');
            fprintf('  cutWaveObj(waveObj, t1, t2)         - Extracts a segment between two times.\n');
            fprintf('  windowWaveObj(waveObj, t1, t2, a)   - Applies a Tukey window to waveform.\n');
            fprintf('  windowfftObj(fftObj, f1, f2, a)     - Applies a Tukey window in frequency domain.\n');
            fprintf('  InterpWave(waveObj, numInt)         - Interpolates waveform using cubic spline.\n');
            fprintf('  smoothWaveObj(WaveObj)              - Smooths waveform with Savitzky-Golay filter.\n');
            fprintf('  normWaveObj(WaveObj)                - Normalizes waveform amplitudes.\n');
            fprintf('  CalcLag(waveObj, maxLag, numInt)    - Computes lag between two waveform channels.\n');
            fprintf('  DetectArrival(waveObj, thresh)      - Detects signal start/stop times.\n');
            fprintf('  FindCrossing(signal, thresh)        - Finds threshold crossing indices.\n\n');

            fprintf('------------------------------------------\n');
            fprintf('FLOW AND PHYSICAL CALCULATIONS\n');
            fprintf('------------------------------------------\n');
            fprintf('  FlowRate(dt, T, FPCF, D_i, angle)   - Calculates water flow velocity & rate.\n');
            fprintf('  calculateV_w(dt, c_s, P, theta)     - Computes water velocity from geometry.\n');
            fprintf('  CalcWaterAngle(c_w)                 - Computes acoustic beam angle in water.\n');
            fprintf('  c_water(T)                          - Calculates sound speed in water (0–60°C).\n\n');

            fprintf('------------------------------------------\n');
            fprintf('USAGE EXAMPLE\n');
            fprintf('------------------------------------------\n');
            fprintf('  wave = Wave_Lib.LoadWaveforms("data.csv");\n');
            fprintf('  Wave_Lib.plt(wave(1), 1);\n');
            fprintf('  fft = Wave_Lib.fftObj(wave(1), 1, 1, 1024);\n');
            fprintf('  Wave_Lib.pltFFT(fft, 2);\n');
            fprintf('  wave_recon = Wave_Lib.IFFTObj(fft);\n');
            fprintf('  Wave_Lib.plt(wave_recon, 3);\n\n');

            fprintf('------------------------------------------\n');
            fprintf('NOTE:\n');
            fprintf(['  This library requires WaveformObj, fftObj, and cwtObj definitions.\n' ...
                     '  STANDARDIZE_FIGURE() is used for consistent figure formatting.\n']);
            fprintf('------------------------------------------\n\n');
        end

        %% ================================================================
        %  PLOTTING FUNCTIONS
        % ================================================================
        function plt(WaveObj, figNum, ClearFlag)
            %PLT Plots WaveformObj in time domain.
            figure(figNum);
            if nargin > 2 && ClearFlag == 1, clf; end

            hold on
            for i = 1:size(WaveObj.volts,2)
                plot(WaveObj.time./1e-6, WaveObj.volts(:,i))
            end
            hold off

            grid("on"); grid("minor");
            xlabel('Time (\mus)'); ylabel('Amplitude (V)');
            title(sprintf('Time domain plot of WaveformObj %d', WaveObj.wave_name));
            STANDARDIZE_FIGURE(gcf);
        end
        
        function pltFFT(fftObj, figNum)
            %PLTFFT Plots FFT magnitude spectrum for an fftObj.
            figure(figNum); hold on
            for i = 1:size(fftObj.c,2)
                plot(fftObj.f/1e6, fftObj.A(:,i));
            end
            hold off
            xlabel('Frequency (MHz)'); ylabel('|FFT(Y)| (V/MHz)');
            title(sprintf('FFT Spectrum of fftObj %d', fftObj.fft_name));
            grid('on'); grid('minor');
            STANDARDIZE_FIGURE(gcf);
        end

        function pltcwt(cwtObj, figNum)
            %PLTCWT Plots a Continuous Wavelet Transform (CWT) amplitude scalogram.
            figure(figNum)
            surface(cwtObj.timeArray, cwtObj.frequencyArray, cwtObj.cwtAmp);
            axis tight; shading flat;
            xlabel("Time (s)"); ylabel("Frequency (Hz)");
            title(sprintf('Amplitude Scalogram of Signal %d', cwtObj.cwt_id));
            colorbar;
            STANDARDIZE_FIGURE(gcf);
        end

        %% ================================================================
        %  LOADING AND SAVING
        % ================================================================
        function LoadWaveform = LoadWaveforms(name_or_matrix)
            %LOADWAVEFORMS Loads waveform data from a file or numeric array.
            disp('Loading...')
            if isstring(name_or_matrix) || ischar(name_or_matrix)
                RawLoadingIn = readmatrix(name_or_matrix);
            else
                RawLoadingIn = name_or_matrix;
            end

            disp('Files loaded! Allocating WaveformObj...');
            t = ~any(isnan(RawLoadingIn),2);
            time = RawLoadingIn(t, 1);
            AllVolts = RawLoadingIn(t, 2:end);

            numMeasures = floor(size(AllVolts,2)/2);
            LoadWaveform = WaveformObj.empty(numMeasures, 0);

            if numMeasures == 0
                LoadWaveform = WaveformObj(1,1,time, AllVolts);
            else
                for ii = 1:numMeasures
                    fprintf('\nAllocating WaveformObj %d / %d...\n', ii, numMeasures);
                    LoadWaveform(ii) = WaveformObj(ii, ii, time, AllVolts(:,2*ii-1:2*ii));
                end
            end
        end

        function SaveWaveforms(WaveformObjArrayToSave, fileToSaveTo)
            %SAVEWAVEFORMS Saves an array of WaveformObjs to CSV file.
            disp('Saving...')
            MatrixToSave(:,1) = WaveformObjArrayToSave(1).time;
            for ii = 1:numel(WaveformObjArrayToSave)
                MatrixToSave = cat(2, MatrixToSave, WaveformObjArrayToSave(ii).volts);
            end
            writematrix(MatrixToSave, fileToSaveTo);
            disp('Done!')
        end

        %% ================================================================
        %  CONVERSION FUNCTIONS
        % ================================================================
        function Mat = WaveformObjToMat(waveObj)
            %WAVEFORMOBJTOMAT Converts a WaveformObj to [time, volts...] matrix.
            Mat = [waveObj.time, waveObj.volts];
        end

        function waveObj = MatToWaveformObj(Mat)
            %MATToWAVEFORMOBJ Converts numeric matrix into a WaveformObj.
            time = Mat(:, 1);
            AllVolts = Mat(:, 2:end);
            waveObj = WaveformObj(1, 1, time, AllVolts);
        end
        
        function FFTObj = MatToFFTObj(Mat)
            %MATTOFFTOBJ Converts numeric matrix to fftObj.
            f = Mat(:,1); c = Mat(:,2:end);
            FFTObj = fftObj(1,1,f,c);
        end
        
        function Mat = fftObjToMat(fftObj)
            %FFTOBJToMAT Converts an fftObj to a numeric matrix.
            Mat = [fftObj.f, fftObj.c];
        end

        %% ================================================================
        %  FFT AND SIGNAL RECONSTRUCTION
        % ================================================================
        function SingleSidedFFT = fftObj(waveObj, fft_name, fft_id, n)
            %FFTOBJ Computes single-sided FFT from a WaveformObj.
            L = length(waveObj.time);
            if nargin < 4, n = 2^nextpow2(L); else, n = 2^nextpow2(n); end

            Y = fft(waveObj.volts, n, 1);
            c = Y(1:n/2+1,:);
            c(2:end-1,:) = 2*c(2:end-1,:);
            f = (waveObj.fs/n)*(0:(n/2))';
            SingleSidedFFT = fftObj(fft_name, fft_id, f, c, L);
        end

        function waveObj = IFFTObj(fftObj)
            %IFFTOBJ Performs inverse FFT to reconstruct WaveformObj.
            fftMat = Wave_Lib.fftObjToMat(fftObj);
            fftAmp = fftMat(:,2:end);
            N = 2*size(fftAmp,1);
            for ii = 1:size(fftAmp,2)
                fftFull = [fftAmp(:,ii); conj(fftAmp(end-1:-1:2,ii))];
                signal_complex(:,ii) = ifft(fftFull, N);
            end
            fs = N*fftObj.df;
            dt = 1/fs;
            N_actual = double(fftObj.originalLength);
            time = (0:(N_actual-1))'*dt;
            waveObj = WaveformObj(fftObj.fft_name, fftObj.fft_id, time, real(signal_complex(1:N_actual,:)));
        end

        %% ================================================================
        %  CWT COMPUTATION
        % ================================================================
        function cwtObjs = cwtObj(WaveObj, freqLimLower, freqLimHigher)
            %C WTOBJ Computes Continuous Wavelet Transform for WaveformObj.
            time = WaveObj.time;
            volts = WaveObj.volts;
            fb = cwtfilterbank('SignalLength', length(time), ...
                               'SamplingFrequency', WaveObj.fs, ...
                               'VoicesPerOctave', 48, ...
                               'FrequencyLimits', [freqLimLower freqLimHigher], ...
                               'Wavelet', "Morse");
            for ii = 1:size(volts,2)
                [cfs(:,:,ii), frq(:,ii)] = cwt(volts(:,ii), FilterBank=fb);
                cwtObjs(ii) = cwtObj(frq(:,ii), time, cfs(:,:,ii), ii);
            end
        end

        %% ================================================================
        %  SIGNAL PROCESSING
        % ================================================================
        function cutWaveObj = cutWaveObj(waveObj, starttime, endtime)
            %CUTWAVEOBJ Extracts waveform segment between two times.
            mask = waveObj.time > starttime & waveObj.time < endtime;
            cutWaveObj = WaveformObj(1,1, waveObj.time(mask), waveObj.volts(mask,:));
        end

        function windowedWaveObj = windowWaveObj(waveObj, starttime, endtime, windowConst)
            %WINDOWWAVEOBJ Applies a Tukey window to a time range.
            if nargin < 4, windowConst = 0; end
            time = waveObj.time; volts = waveObj.volts;
            windowedVolts = zeros(size(volts));
            for ii = 1:size(volts,2)
                range = time > starttime(ii) & time < endtime(ii);
                L = sum(range);
                w = zeros(size(time));
                w(range) = tukeywin(L, windowConst);
                windowedVolts(:,ii) = w .* volts(:,ii);
            end
            windowedWaveObj = WaveformObj(1,1,time, windowedVolts);
        end

        function windowedfftObj = windowfftObj(fftObject, minfreq, maxfreq, windowConst)
            %WINDOWFFTOBJ Applies a Tukey window to an fftObj within frequency limits.
            if nargin < 4, windowConst = 0; end
            freq = fftObject.f; c = fftObject.c;
            windowedComplexAmps = zeros(size(c));
            for ii = 1:size(c,2)
                relfreqs = freq > minfreq & freq < maxfreq;
                L = sum(relfreqs);
                w = zeros(size(freq));
                w(relfreqs) = tukeywin(L, windowConst);
                windowedComplexAmps(:,ii) = w .* c(:,ii);
            end
            windowedfftObj = fftObj(1,1,freq, windowedComplexAmps);
        end

        function interpWaveObj = InterpWave(waveObj, numInt)
            %INTERPWAVE Performs cubic interpolation on WaveformObj data.
            xq = waveObj.time(1):waveObj.dt/numInt:waveObj.time(end);
            for ii = 1:size(waveObj.volts,2)
                intWave(:,ii) = interpn(waveObj.time, waveObj.volts(:,ii), xq, 'cubic')';
            end
            interpWaveObj = WaveformObj(1,1,xq', intWave);
        end

        function [lag, cross_correlation] = CalcLag(waveObj, maxLag, numInt)
            %CALCLAG Calculates lag between two waveform channels.
            maxlagidx = ceil(maxLag/waveObj.dt);
            [xcf, rawlag] = xcorr(waveObj.volts(:,1), waveObj.volts(:,2), maxlagidx);
            cross_correlation = WaveformObj(1,1,rawlag*waveObj.dt, xcf);
            if nargin > 2
                cross_correlation = Wave_Lib.InterpWave(cross_correlation, numInt);
            end
            [~,maxidx] = max(cross_correlation.volts);
            lag = cross_correlation.time(maxidx);
        end

        function smoothed_waveObj = smoothWaveObj(WaveObj)
            %SMOOTHWAVEOBJ Smooths waveform using Savitzky-Golay filter.
            for ii = 1:size(WaveObj.volts,2)
                smoothedVolts(:,ii) = sgolayfilt(WaveObj.volts(:,ii),3,101);
            end
            smoothed_waveObj = WaveformObj(1,1,WaveObj.time, smoothedVolts);
        end

        function normWaveObj = normWaveObj(waveObj)
            %NORMWAVEOBJ Normalizes waveform amplitudes to max = 1.
            for ii = 1:size(waveObj.volts,2)
                normvolts(:,ii) = waveObj.volts(:,ii)./max(waveObj.volts(:,ii));
            end
            normWaveObj = WaveformObj(1,1,waveObj.time, normvolts);
        end

        %% ================================================================
        %  FLOW AND PHYSICAL CALCULATIONS
        % ================================================================
        function [v_w, Q] = FlowRate(dt, T, FPCF, D_i, chordAngle)
            %FLOWRATE Computes flow velocity and volumetric rate.
            c_s = Wave_Lib.c_water(T);
            theta_w = Wave_Lib.CalcWaterAngle(c_s);
            PathLength = 2*D_i/(cos(theta_w)*cosd(chordAngle));
            v_w = Wave_Lib.calculateV_w(dt, c_s, PathLength/100, theta_w) * FPCF * cos(chordAngle);
            Area = 5.528^2*pi/4;
            Q = v_w *100 * Area *60 /1000; % L/min
        end

        function WaterVelocity = calculateV_w(dt, c_s, P_water, theta)
            %CALCULATEV_W Computes water velocity from acoustic geometry.
            WaterVelocity = (dt*c_s^2)/(2*P_water*sin(theta));
        end
            
        function theta_w = CalcWaterAngle(c_w)
            %CALCWATERANGLE Computes acoustic beam angle in water.
            theta_w = asin(sind(38)*c_w/2586);
        end

        function c_l = c_water(T)
            %C_WATER Computes speed of sound in water (empirical 0–60 °C).
            if T < 0 || T > 60
                warning("T outside of 0–60°C range. Result extrapolated.");
            end
            p1 = 0.0001866; p2 = -0.05278; p3 = 4.97; p4 = 1403;
            c_l = p1*T^3 + p2*T^2 + p3*T + p4;
        end

        %% ================================================================
        %  ARRIVAL DETECTION / THRESHOLDING
        % ================================================================
        function [starttimes, stoptimes] = DetectArrival(waveObj, thresh)
            %DETECTARRIVAL Detects signal arrival and end times by threshold.
            time = waveObj.time; volts = waveObj.volts; dt = waveObj.dt;
            for ii = 1:size(volts,2)
                zc = Wave_Lib.FindCrossing(volts(:,ii), 0);
                tc = Wave_Lib.FindCrossing(volts(:,ii), thresh);
                [~, idx] = min(abs(zc - tc(1)));
                stop = zc(idx+4);
                starttimes(ii) = zc(idx)*dt + time(1);
                stoptimes(ii) = stop*dt + time(1);
            end
        end
        
        function CrossingIdx = FindCrossing(signal, thresh)
            %FINDCROSSING Returns indices where signal crosses threshold.
            signal = signal/max(signal);
            sigAbove = signal.*(signal>thresh);
            CrossingIdx = find(diff(sigAbove>0));
        end

    end
end
