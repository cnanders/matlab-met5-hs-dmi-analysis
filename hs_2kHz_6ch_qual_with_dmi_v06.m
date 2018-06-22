% hs_2kHz


clear
clc
close all

m_per_diff_over_sum = 120e-6/2;  % 200 um diameter pinhole.  As the wafer moves dZ, spot on diode moves 2*dZ verified in ZEMAX

freq = 2000; % Hz 

f_min = 10;
f_max = 500;

% addpath(genpath(fullfile(pwd, '../../Functions')));

[file_name, file_path] = uigetfile( ...
    '*','Select the log file',...
    '~/Documents/Matlab/Data/MET5/HS/2018-06-22-VME' ... // 06-13  // 2014-12-22 //2014-10-23
);

fprintf('%s \n\n', [file_path, file_name]);

% textread help:
% u     integer
% d     signed integer
% f     floating-point


tic

hFile = fopen([file_path, file_name]);

% %d = signed integer, 32-bit
cFormat = [...
    '%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,', ...
    '%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,', ...
    '%d,%d,%d,%d' ...
];

ceData = textscan(...
    hFile, cFormat, -1, ...
    'headerlines', 1 ...
);
fclose(hFile);

% The last line of the text file contains partial data.  Truncate all sample
% arrays to the length of dmi4 sample appay

numSamples = length(ceData{28});
for n = 1 : length(ceData)
    while length(ceData{n}) > numSamples
        ceData{n}(end) = [];
    end
end
    
toc

stMap = struct();

stMap.ch6top1 = 1;
stMap.ch6bot1 = 2;
stMap.ch5top1 = 3;
stMap.ch5bot1 = 4;
stMap.ch4top1 = 5;
stMap.ch4bot1 = 6;
stMap.ch3top1 = 7;
stMap.ch3bot1 = 8;
stMap.ch2top1 = 9;
stMap.ch2bot1 = 10;
stMap.ch1top1 = 11;
stMap.ch1bot1 = 12;

stMap.ch6top2 = 13;
stMap.ch6bot2 = 14;
stMap.ch5top2 = 15;
stMap.ch5bot2 = 16;
stMap.ch4top2 = 17;
stMap.ch4bot2 = 18;
stMap.ch3top2 = 19;
stMap.ch3bot2 = 20;
stMap.ch2top2 = 21;
stMap.ch2bot2 = 22;
stMap.ch1top2 = 23;
stMap.ch1bot2 = 24;

stMap.dmi1 = 25;
stMap.dmi2 = 26;
stMap.dmi3 = 27;
stMap.dmi4 = 28;

% side === indicates to the capacitor used.  Do not confuse this with
% "side" of the diode.  "a", and "b" refer to "top" and "bottom" halves of
% the dual-diode, respectively
% A single measurement returns the "capacitor/side" of the ADC, and the
% counts from the top and bottom of each diode of each channel

% average the values from each cap within the 1 ms acquisition period

% top
top(1, :) = (ceData{stMap.ch1top1} + ceData{stMap.ch1top2}) / 2;
top(2, :) = (ceData{stMap.ch2top1} + ceData{stMap.ch2top2}) / 2;
top(3, :) = (ceData{stMap.ch3top1} + ceData{stMap.ch3top2}) / 2;
top(4, :) = (ceData{stMap.ch4top1} + ceData{stMap.ch4top2}) / 2;
top(5, :) = (ceData{stMap.ch5top1} + ceData{stMap.ch5top2}) / 2;
top(6, :) = (ceData{stMap.ch6top1} + ceData{stMap.ch6top2}) / 2; 

% bottom
bot(1, :) = (ceData{stMap.ch1bot1} + ceData{stMap.ch1bot2}) / 2;
bot(2, :) = (ceData{stMap.ch2bot1} + ceData{stMap.ch2bot2}) / 2;
bot(3, :) = (ceData{stMap.ch3bot1} + ceData{stMap.ch3bot2}) / 2;
bot(4, :) = (ceData{stMap.ch4bot1} + ceData{stMap.ch4bot2}) / 2;
bot(5, :) = (ceData{stMap.ch5bot1} + ceData{stMap.ch5bot2}) / 2;
bot(6, :) = (ceData{stMap.ch6bot1} + ceData{stMap.ch6bot2}) / 2;    

diff = top - bot;
sum = top + bot;
dos = double(diff)./double(sum);
z = dos * m_per_diff_over_sum * 1e9;

% CLOCK TIMES ARE SOURCE HEAD TIMES
cecLabels = {'z 5:30 (1)',  'z 9:30 (2)',  'z 1:30 (3)',  'ang 0:30 (4)', 'ang 4:30 (5)', 'ang 8:30 (6)'};


%{
figure
plot(dos')
title('diff/sum')
legend(cecLabels)

figure
plot(a')
title('top')
legend(cecLabels)

figure
plot(b')
title('bottom')
legend(cecLabels)
%}

%{

figure('Name', [file_path, file_name])

subplot(141)
plot(dos')
title('diff/sum')
legend(cecLabels)

subplot(142)
plot(sum')
title('sum')
legend(cecLabels)

subplot(143)
plot(top')
title('top')
legend(cecLabels)

subplot(144)
plot(bot')
title('bottom')
legend(cecLabels)

%}

%{
figure('Name', sprintf('%s SUM', file_name))
plot(sum')
title('sum')
legend(cecLabels)
%}


dos_z_avg = (dos(1, :) + dos(2, :) + dos(3, :)) / 3;
nm_z_avg = dos_z_avg * m_per_diff_over_sum ./ 1e-9;

lPlotAvg = false;

if lPlotAvg
    figure('Name', sprintf('%s dos z avg', file_name))
    plot(nm_z_avg, '.-b')
    title('z avg (DOS * 120e-6/2) of three central channels')
    ylabel('nm')
    xlabel('sample (@ 1 kHz cap sides averaged)')
end

% Need to average from each capacitor

% Compute and plot PSD

t = [0 : length(dos) - 1] * 1e-3;

channels = [1, 2, 3];
%channels = [4]


hFigOverlap = figure;
dScreenSize = get(0, 'ScreenSize');
dWidthFigure = 1200;
dHeightFigure = 500;
dPos = [...
    (dScreenSize(3) - dWidthFigure) / 2 ...
    (dScreenSize(4) - dHeightFigure) / 2 ...
    dWidthFigure ...
    dHeightFigure ...
];
set(hFigOverlap,'Position', dPos);
set(hFigOverlap, 'Name', file_name)
hAxesPsd = subplot(211);
hAxesCas = subplot(212);
% hold(hAxesPsd, 'on'); Need to turn on hold after first loglog
% hold(hAxesCas, 'on');
ceLabelPsd = {}
ceLabelCas = {}




h = figure;
dScreenSize = get(0, 'ScreenSize');
dWidthFigure = 1200;
dHeightFigure = 500;
dPos = [...
    (dScreenSize(3) - dWidthFigure) / 2 ...
    (dScreenSize(4) - dHeightFigure) / 2 ...
    dWidthFigure ...
    dHeightFigure ...
];
set(h,'Position', dPos);
set(h, 'Name', file_name)



for n = 1 : length(channels)
    
    channel = channels(n);
    pos = z(channel, :);
    
    
    % Remove DC
    pos = pos - mean(pos);
    
    subplot(length(channels), 3, n * 3 - 2)
    plot(t, pos , '.-b');
    xlabel('Time (s)')
    ylabel('Z (nm)');
    cTitle = [...
        'Z vs. Time ', ...
        sprintf('%s ', cecLabels{channel}), ...
        sprintf('PV = %1.2f nm ', max(pos) - min(pos)), ...
        sprintf('RMS = %1.2f nm ', std(pos)) ...
    ];
    title(cTitle)
    
    % PSD
    
    [freq_psd, energy_psd] = Psd.calc(t, pos - mean(pos));
    [freq_psd, energy_psd] = Psd.fold(freq_psd, energy_psd);

    subplot(length(channels), 3, n * 3 - 1)
    loglog(freq_psd, energy_psd, '-b');
    xlabel('Freq (Hz)')
    ylabel('PSD (nm^2/Hz)');
    
    % Push to cumulative plot
    loglog(hAxesPsd, freq_psd, energy_psd);
    ceLabelPsd{end + 1} = cecLabels{channel};
    hold(hAxesPsd, 'on') % need to do after first loglog
        
    
    % Compute the RMS of the height signal within 98 - 102 Hz
    [f_band, energy_band] = Psd.band(freq_psd, energy_psd, 1/f_min, 1/f_max);
    
    z_rms_band = sqrt(Psd.power(f_band, energy_band));
    
    % Compute percentage of signal energy
    band_power_frac = Psd.power(f_band, energy_band)/Psd.power(freq_psd, energy_psd);
    
    cTitle = [ ...
       sprintf('PSD of %s ', cecLabels{channel}), ...
       sprintf(...
        'RMS of [%1.0fHz, %1.0fHz] band = %1.2f nm (%1.1f %%)', ...
        f_min, ...
        f_max, ...
        z_rms_band, ...
        band_power_frac * 100 ...
       ) ...
    ];
    title(cTitle);
    
    
    % Power cumulative
    [f, powerc] = Psd.powerCumulative(f_band, energy_band);
    subplot(length(channels), 3, n * 3);
    semilogx(f, sqrt(powerc), '-b');
    xlabel('Freq (Hz)')
    ylabel('RMS Cumulative (nm)');
    
    % Push to overlap plot
    semilogx(hAxesCas, f, sqrt(powerc));
    ceLabelCas{end + 1} = cecLabels{channel};
    hold(hAxesCas, 'on')


end

%{

figure('Name', [file_path, file_name])


for n = 1 : length(channels)
    
    channel = channels(n);
    pos = z(channel, :);
    
    [freq_psd, energy_psd] = Psd.calc(t, pos - mean(pos));
    [freq_psd, energy_psd] = Psd.fold(freq_psd, energy_psd);

    subplot(length(channels), 1, channel)
    loglog(freq_psd, energy_psd, '-b');
    xlabel('Freq (Hz)')
    ylabel('Power/Freq of Diff/Sum');
    
    % Compute the RMS of the height signal within 98 - 102 Hz
    f_min = 98;
    f_max = 102;
    [f_band, energy_band] = Psd.band(freq_psd, energy_psd, 1/f_min, 1/f_max);
    
    z_rms_band = sqrt(Psd.power(f_band, energy_band))
    
    % Compute percentage of signal energy
    band_power_frac = Psd.power(f_band, energy_band)/Psd.power(freq_psd, energy_psd)
    
    cTitle = [ ...
       sprintf('PSD of %s', cecLabels{channel}), ...
       sprintf(...
        'RMS of [%1.0fHz, %1.0fHz] band = %1.2f nm (%1.1f %%)', ...
        f_min, ...
        f_max, ...
        z_rms_band, ...
        band_power_frac * 100 ...
       ) ...
    ];
    title(cTitle);

end

%}



dDMI_SCALE = 632.9907/4096; % dmi axes come in units of 1.5 angstroms
                                  % Convert to nm
                                 
dErrU_ret = double(ceData{stMap.dmi1});
dErrV_ret = double(ceData{stMap.dmi2});
dErrU_waf = double(ceData{stMap.dmi3});
dErrV_waf = double(ceData{stMap.dmi4});

dXDat_ret = dDMI_SCALE * 1/sqrt(2) * (dErrU_ret + dErrV_ret);
dYDat_ret = -dDMI_SCALE * 1/sqrt(2) * (dErrU_ret - dErrV_ret);
dXDat_waf = -dDMI_SCALE * 1/sqrt(2) * (dErrU_waf + dErrV_waf);
dYDat_waf = dDMI_SCALE * 1/sqrt(2) * (dErrU_waf - dErrV_waf);

dPosDmi(1, :) = dXDat_ret;
dPosDmi(2, :) = dYDat_ret;
dPosDmi(3, :) = dXDat_waf; 
dPosDmi(4, :) = dYDat_waf;



cecLabels = {'X reticle', 'Y reticle', 'X wafer', 'Y wafer'};
channels = [1, 2, 3, 4];


% Remove DC and tilt from DMI data

for n = 1 : length(channels)
    
    channel = channels(n);
    pos = dPosDmi(channel, :);
    
    % Remove DC
    % pos = pos - mean(pos);
    
    % Remove DC and tilt
    p = polyfit(t, pos, 1);
    pos = pos - polyval(p, t);
    
    
    dPosDmi(channel, :) = pos;
end



h = figure;
dScreenSize = get(0, 'ScreenSize');
dWidthFigure = 1200;
dHeightFigure = 600;
dPos = [...
    (dScreenSize(3) - dWidthFigure) / 2 ...
    (dScreenSize(4) - dHeightFigure) / 2 ...
    dWidthFigure ...
    dHeightFigure ...
];
set(h,'Position', dPos);
set(h, 'Name', file_name)

for n = 1 : length(channels)
    
    channel = channels(n);
    pos = dPosDmi(channel, :);
    
    % Vanilla
    
    subplot(length(channels), 3, n * 3 - 2)
    plot(t, pos , '.-b');
    xlabel('Time (s)')
    ylabel(sprintf('%s (nm)', cecLabels{channel}));
    cTitle = [...
        sprintf('%s vs. Time ', cecLabels{channel}), ...
        sprintf('PV = %1.2f nm ', max(pos) - min(pos)), ...
        sprintf('RMS = %1.2f nm ', std(pos)) ...
    ];
    title(cTitle)
    
    % PSD
    
    [freq_psd, energy_psd] = Psd.calc(t, pos - mean(pos));
    [freq_psd, energy_psd] = Psd.fold(freq_psd, energy_psd);

    subplot(length(channels), 3, n * 3 - 1)
    loglog(freq_psd, energy_psd, '-b');
    xlabel('Freq (Hz)')
    ylabel('PSD (nm^2/Hz)');
    
    % Push to cumulative plot
    loglog(hAxesPsd, freq_psd, energy_psd);
    ceLabelPsd{end + 1} = cecLabels{channel};
    
    % Compute the RMS of the height signal within 98 - 102 Hz
    
    [f_band, energy_band] = Psd.band(freq_psd, energy_psd, 1/f_min, 1/f_max);
    
    z_rms_band = sqrt(Psd.power(f_band, energy_band));
    
    % Compute percentage of signal energy
    band_power_frac = Psd.power(f_band, energy_band)/Psd.power(freq_psd, energy_psd);
    
    cTitle = [ ...
       sprintf('PSD of %s', cecLabels{channel}), ...
       sprintf(...
        'RMS of [%1.0fHz, %1.0fHz] band = %1.2f nm (%1.1f %%)', ...
        f_min, ...
        f_max, ...
        z_rms_band, ...
        band_power_frac * 100 ...
       ) ...
    ];
    title(cTitle);
    
    % Power Cumulative of 10Hz onward
    
    

    [f, powerc] = Psd.powerCumulative(f_band, energy_band);
    subplot(length(channels), 3, n * 3);
    semilogx(f, sqrt(powerc), '-b');
    xlabel('Freq (Hz)')
    ylabel('RMS Cumulative (nm)');
    
    % Push to overlap plot
    semilogx(hAxesCas, f, sqrt(powerc));
    ceLabelCas{end + 1} = cecLabels{channel};

end

legend(hAxesPsd, ceLabelPsd);
title(hAxesPsd, 'PSD')
xlabel(hAxesPsd, 'Freq (Hz)');
ylabel(hAxesPsd, 'PSD (nm^2/Hz)');

cTitle = sprintf(...
    'Cumulative Amplitude Spectrum [%1.0fHz, %1.0fHz]', ...
    f_min, ...
    f_max ...
);
title(hAxesCas, cTitle);
legend(hAxesCas, ceLabelCas);
xlabel(hAxesCas, 'Freq (Hz)');
ylabel(hAxesCas, 'Cumulative Amplitude RMS (nm)');
figure(hFigOverlap); % bring to front


%{
timeStart = 1;
timeEnd = 1.3;
figure
subplot(211)
hold on
plot(t, dXDat_ret, '.-b');
plot(t, dXDat_waf, '.-r');
xlim([timeStart timeEnd])
legend({'x reticle', 'x wafer'})
subplot(212)
hold on
plot(t, dYDat_ret, '.-b');
plot(t, dYDat_waf, '.-r');
xlim([timeStart timeEnd])
legend({'y reticle', 'y wafer'})


figure
plot(xcorr(1/1000, dXDat_ret, dXDat_waf))
%}



