classdef Psd
        
    properties (Constant)
        
        % Energy (Power/Freq) units useful for ZEMAX
        A_A_um      = 1e-10*1e-10*1e-6;
        nm_nm_mm    = 1e-9*1e-9*1e-3;
    end


    methods (Static)
        
        
        % Compute the PSD of a signal
        % @param {double 1xm} t - time (s) samples or pos (m)
        % @param {double 1xm} v - amplitude samples (volts) (can also be other unit)    
            
        function [f, psd] = calc(t, v)
            
            n = length(t);
            dt = t(2) - t(1);       % time space delta
            df = 1/dt/n;            % frequency space delta
            period = max(t) - min(t);
            
            debug = false;
            
            if (mod(n, 2) == 0)
                % Even
                f  = -1/dt/2 : df : 1/dt/2 - df;
            else
                % Odd
                f  = -1/dt/2 + df/2 : df : 1/dt/2 - df/2;
            end
            
            % https://en.wikipedia.org/?title=Spectral_density
            %
            % The power spectrum of a time series x(t) describes how the 
            % variance of the data x(t) is distributed over the frequency
            % components into which x(t) may be decomposed.
            %
            % When continuous FT is performed, the units of the FT are
            % the units of the original signal (in this case volts)
            % multiplied by the time/space unit, in this case seconds.  So
            % the units are s*V
            %
            % The MATLAB FFT (DFT) needs to be multiplied by "dt" to get
            % the units and scaling correct.  See my notes "General FT
            % sampling stuff". This is an artifact of the non-symmetric DFT
            % pair that MATLAB defines.

            % One way to show that this needs to be done is to not multiply
            % the fft by the sampling delta and plot abs(v_dft).^2 for
            % different sampling rates (dt) (all better than Nyquist)
            % 
            % Each plot looks different which you know can't be true since
            % the fundamental signal isn't changing; just the way it is
            % sampled

            v_dft = fft(v);
            v_dft = v_dft*dt;
            v_dft = fftshift(v_dft);
            
            % The energy of a signal is the integral of the square of the
            % amplitude.  For a signal with amplitude units of Volts and
            % time/space units of sec, the energy has units of s*V^2
            %
            % PSD coefficients need to have unit of energy, or s*V^2.  This unit
            % can be otbained by dividing v_dft^2 by something with units of 
            % seconds. But which thing with this unit? 
            %
            % By dividing the v_dft^2 vector by period of the time/space signal
            % the energy of the signal at each frequency is obtained. I.E., 
            % v_dft^2/period is the energy of the signal at each frequency.
            %
            % An example of this can be seen by creating a signal that is
            % the sum of a few sinusoids.  psd_example_2_2.m (in this
            % folder) shows this.  Because the energy is evenly distributed
            % between the +1 and -1, the PSD coefficient for +1 and -1 each
            % have half of the energy
            %
            % Note also that power (energy/time) or (energy/space) has units
            % of V^2 and would be obtained by dividing the energy by the 
            % length of the integration period (energy over time).  
            %
            % Power per frequency has units of V^2/Hz.  This is equivalent
            % to s*V^2 or energy.  The units of the PSD coefficients can
            % either be thought of as power per frequency or energy. These
            % units are equivalent.
            %
            % Multiplying the PSD coefficient by df gives the total power
            % within that frequency band.  (Just a property of PSDs)
            %
            % For signals whose mean is zero, the power of the signal is
            % the same as the variance of the signal.  Lets prove it:
            %
            % The variance of a sampled signal is defined as
            %   sum(abs(x-mean(x)).^2)/N,
            %   i.e., it is the average of the abs(x_i - mean(x))^2 values.
            %
            % When mean(x) = 0, the variance is:
            %   sum(abs(x).^2)/N
            %
            % What about the power of the signal?  The energy of a signal
            % is defined as the integral of the square of the signal, i.e.:
            %   sum(x.^2)*dx
            %
            % Power is energy divided by duration or:
            %   sum(x.^2)*dx/period == sum(x.^2)/N
            %
            % Since this is numerically equivalent to the variance, the
            % result holds.  Why did I bother showing that result?
            %
            % Since the area under any part of the PSD curve gives the
            % power of the signal within a frequency band.  IF THE SIGNAL
            % HAS A MEAN OF ZERO (NEED TO SUBTRACT MEAN), the power of the signal is the same as
            % the variance of the signal.  This means the area of the PSD
            % curve for any frequency band gives the variance (square root
            % for RMS) of the surface height within that band. This is a
            % cool result.
            %
            % Note that df == 1/dt/n.  With n == T/dT, df = 1/T.  Thus df
            % doesn't change with how well you sample the signal; it only
            % changes when you sample the signal longer or shorter in time.
            % By multiplying the PSD coefficient by df, you are really
            % dividing by the signal duration in time/space, which converts
            % energy to power.
            %
            % The way to think about frequency resolution is that you
            % convolve the spectrum with a sinc whose width is by the data
            % acquisition period.
            
            v_psd = 1/period*abs(v_dft).^2;

            if (debug)
                figure
                plot(f, v_psd, '.-r')
                ylabel('Energy or Power/Freq');
            end
            
            % Return
            
            % Notes on combining PSDs on the same plot
            % 
            % If you ever need to combine two PSDs on the same plot, for
            % example by measuring the same surface with two different
            % methods with different resolution, the PSD coefficients that are over a
            % shorter period need to be scaled by the ratio of the longer and
            % shorter data windows since the energy scales linearly with
            % the length of the measurement window
            
            psd = v_psd;
           
        end
        
        % Returns the cumulative power up to frequency f vs. f
        % @param {double 1xm} f - frequency samples (units of 1/s or 1/m)
        % @param {double 1xm} psd - psd coefficients (units of energy, i.e., s*V^2 or m * m * m)
        
        function [f, power] = powerCumulative(f, psd)
           
            power = zeros(size(f));
            df = f(2) - f(1);

            for n = 1 : length(f)
               if n == 1
                   power(n) = psd(n) * df;
               else
                   power(n) = power(n - 1) + psd(n) * df;
               end
            end
            
        end
        
        % For optical PSDs, we want to only report the PSD coefficient at
        % frequencies > 0.  In order to conserve energy, we need to sum
        % the +f and -f psd coefficients so the power of the reported
        % PSD (computed using Parseval) is the same as the original PSD.
        % @param {double 1xm} f - frequency samples (units of 1/s or 1/m)
        % @param {double 1xm} psd - psd coefficients (units of energy, i.e., s*V^2 or m * m * m)
            
            
        function [f, psd] = fold(f, psd)
                                   
            
            n = length(f);
            df = f(2) - f(1);
            dt = 1/df/n;
            f_pos = 0: df : 1/dt/2;
            
            debug = false;

            if (mod(n, 2) == 0)

                % Even
                psd_neg_f = fliplr(psd(1:n/2));
                psd_0     = psd(n/2 + 1);
                psd_pos_f = psd(n/2+2:end);

                % neg is 1 longer than pos. assign 0 amplitude to the last
                % element of pos so summing the arrays is simple

                psd_pos_f(end + 1) = 0;
                psd_pos = [psd_0 [psd_neg_f + psd_pos_f]];
                
            else
                % Odd
                psd_neg_f = fliplr(psd(1:floor(n/2)));
                psd_0 = fliplr(psd(ceil(n/2)));
                psd_pos_f = psd(ceil(n/2)+1:end);

                psd_pos = [psd_0 [psd_neg_f + psd_pos_f]];

            end

            if debug
                
                figure
                plot(psd_neg_f, '.-r');
                title('Neg f');

                figure
                plot(psd_pos_f, '.-r');
                title('Pos f');

                figure
                plot(f_pos, psd_pos, '.-r');
            
            end
            
            % Return
            psd = psd_pos;
            f = f_pos;

        end
        
        % Splice a (f, psd) vector set to only periods of a specified band.
        % This is useful when plotting different bands of a PSD with
        % different colors
        % @param {double 1xm} f - frequency samples (units of 1/s or 1/m)
        % @param {double 1xm} psd - psd coefficients (units of energy, i.e., s*V^2 or m * m * m)
        % @param {double 1x1} d_max - the maximum period of the band(units of s or m)
        % @param {double 1x1} d_min - the minimum period of the band (units of s or m)
         
        function [f, psd] = band(f, psd, d_max, d_min)
            
            % Var       Type    Description         Unit
            % f         double  freqency            1/m or 1/s
            % psd       double  PSD                 m^2*m or V^2*s
            % d_max     double  max period          m or s
            % d_min     double  min period          m or s
            
            
            idx = abs(f) >= 1/d_max & abs(f) < 1/d_min;
            f = f(idx);
            psd = psd(idx);
            
        end
        
        
        % Compute the power of a (f, psd) vector set
        % @param {double 1xm} f - frequency samples (units of 1/s or 1/m)
        % @param {double 1xm} psd - psd coefficients (units of energy, i.e., s*V^2 or m * m * m)
        
        function [o] = power(f, psd)
            
            if (length(f) < 2)
                o = 0;
                return;
            end
            
            df = f(2) - f(1);
            o = sum(psd)*df;
            
        end
        
        
        % To compute the rms of the surface height for a particular
        % frequency band just do:         
        %{
        [f_lsf, psd_lsf] = Psd.band(f_psd, psd, d_lsf_max, d_lsf_min);
        [f_msf, psd_msf] = Psd.band(f_psd, psd, d_msf_max, d_msf_min);
        [f_hsf, psd_hsf] = Psd.band(f_psd, psd, d_hsf_max, d_hsf_min);

        rms_lsf = sqrt(Psd.power(f_lsf, psd_lsf));
        rms_msf = sqrt(Psd.power(f_msf, psd_msf));
        rms_hsf = sqrt(Psd.power(f_hsf, psd_hsf));
        
        %}
        
        
        % Generate a random PSD from a nominal slope and y-intercept
        % @param {struct} params - config
        % @param {double 1xm} params.x - log x values (period or freq)
        % @param {double 1x1} params.icept - log y intercept (m*m*m or nm*nm*mm)
        % @param {double 1x1} params.slope - loglog slope.  When x values
        % are period, the loglog slope is positive (left to right - more
        % energy in longer periods).  The slope is negative when the x
        % values are frequency.
        % @param {double 1x1} params.order - order of polynomial fit
        % @returns {double 1xorder} - polynomial coefficients of PSD
        
        function [p] = rand(params)
            
            % Linear portion
            y1 = params.icept + params.slope*params.x;

            % Ripple 
            y2 = zeros(size(params.x));
            for n = 1:5
                f = 6*rand(1, 1)/(max(params.x) - min(params.x));
                phase = rand(1, 1)*2*pi;
                
                a = sin(2*pi*f*params.x + phase);
                y2 = y2 + a;
            end
            
            % Polynomial fit
            y3 = y1 - 0.1*y2;
            p = polyfit(params.x, y3, params.order);
            
            
        end
        
        % Convert a PSD in SI units to log10 SI units space 
        % @param {double 1xm} f - frequency samples (units of 1/s or 1/m)
        % @param {double 1xm} psd - psd coefficients (units of energy, i.e., s*V^2 or m * m * m)
        % @returns {double 1xm} fLog10 - the log10 of frequency samples in
        % SI uints
        % @returns {double 1xm} eLog10 - the log10 of energy samples in
        % SI uints
        
        function [fLog10, eLog10] = log(f, e)
            
            fLog10 = log10(f);
            eLog10 = log10(e);
 
        end
        
        % Convert a PSD in SI units to log10 Zygo units: 
        % freq      == 1/mm
        % energy    == nm*nm*mm
        % @param {double 1xm} f - frequency samples (units of 1/s or 1/m)
        % @param {double 1xm} psd - psd coefficients (units of energy, i.e., s*V^2 or m * m * m)
        % @returns {double 1xm} fLog10 - the log10 of frequency samples in
        % Zygo uints
        % @returns {double 1xm} eLog10 - the log10 of energy samples in
        % Zygo uints
        
        function [fLog10Z, eLog10Z] = logZ(f, e)
            
            [fLog10, eLog10] = Psd.log(f, e);
            
            % In log10 space, metric unit conversion is add/sub
            
            fLog10Z = fLog10 - 3;   % 1/m -> 1/mm
            eLog10Z = eLog10 + 21;  % m*m*m -> nm*nm*mm (9)(9)(3)
            
            % Also, converting from freq to period in log10 space is done
            % by adding a negative sign to the frequency
            
        end
        
        
        
        
        
    end

    
end

