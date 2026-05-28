clear; clc; close all;

% Diode
Vx = 0.49; % Diode DC voltage drop (V)
Rd = 80e-3; % Diode resistance (ohms)
diodeCost = 0.78; % ($ total)
diodeCurrentRating = 2; % (A average)
diodePeakCurrentRating = 50; % (A peak)
diodeReverseVoltageRating = 200; % (V peak)

% Inductor
L = 300e-3; % Inductance (H)
RDCR = 6; % Inductor DC resistance (ohms)
inductorCurrentRating = 1; % (A RMS)
inductorCost = 31.04; % ($ total)

% Capacitor
C = 5e-6; % Capacitance (F)
RESR = 0.004; % Capacitor Equivalent Series Resistance (ohms)
capacitorVoltageRating = 400; % (V)
capacitorCurrentRating = 9; % (A rms)
capacitorCost = 2.37; % ($ total)

f = 60; % Hz


% Set the design number for the Right summary Output: (1, 2, 3)
DesignNumber = 1;



% Parameters
Vll_rms = 208;
%f = 60;
w = 2*pi*f;
%L = 0.3;
%C = 5e-6;
%R = 50;

% Non-Ideal Parameters
%Vx = 0.49;
%Rd = 0.08;
R_DCR = RDCR;
R_ESR = RESR;

% Line-Line/Neutral Voltage Calculations
Vln_rms = Vll_rms/sqrt(3);
Vll_pk = Vll_rms * sqrt(2); 
Vln_pk = Vll_pk / sqrt(3);


% Exceed Rating and CCM Switch
% Checks if Ratings at exceeded at some point
CCM=0;
SWdiodeReverseVoltageRating = 0;
SWdiodeCurrentRating = 0;
SWdiodePeakCurrentRating = 0;
SWinductorCurrentRating = 0;
SWcapacitorVoltageRating = 0;
SWcapacitorCurrentRating = 0;

results = [];
RVals = 100:10:500;
for Rindx = 1:length(RVals)
    R = RVals(Rindx);
% DC Calculations
V0 = 3*Vll_pk/pi - 2*Vx;
I0 = V0/(R + R_DCR+ 2*Rd);
Pavg = I0^2*R;

% Samples and Harmonics array
wt = linspace(0, 2*pi/w, 1001);
wt = wt(1:end-1);
harmonics = 6:6:length(wt);
num_harmonics = length(harmonics);

% Initialize Phasor variables
VRectPhasor = 0*harmonics;
IRectPhasor = 0*harmonics;
vLoadPhasor = 0*harmonics;
IloadPhasor = 0*harmonics;
ICPhasor = 0*harmonics;

% AC Calculations
for indx = 1:num_harmonics
    n = harmonics(indx);

    VRectPhasor(n) = (-1)^(n/6) * -6 * Vll_pk/(pi * (n^2 - 1));
    
    % Impedance Calcuations
    ZL = 1j * n*w * L + R_DCR;
    ZC = 1 / (1j * n*w * C) + R_ESR;
    ZeqRC = 1/(1/R+1/ZC);
    
    IRectPhasor(n) = VRectPhasor(n)/(ZeqRC+ZL+ 2*Rd);
    vLoadPhasor(n) = VRectPhasor(n) * (ZeqRC/(ZeqRC+ZL+ 2*Rd));
    IloadPhasor(n) = vLoadPhasor(n)/R;
    ICPhasor(n) = vLoadPhasor(n)/ZC;
    Pavg = Pavg + abs(vLoadPhasor(n))*abs(IloadPhasor(n)) * cos(angle(vLoadPhasor(n))-angle(IloadPhasor(n)))/2;
end



% Build Waveform from phasor
vRectWaveform = V0 + 0*wt; 
iRectWaveform = I0 + 0*wt; 
vLoadWaveform = I0*R + 0*wt;
vRectWaveform_ideal =  V0 - 2*Vx + 0*wt;
iCWaveform = 0*wt;
for indx = 1:length(harmonics)
    vRectWaveform_ideal = vRectWaveform_ideal + abs(VRectPhasor(indx)) * cos(indx*w*wt + angle(VRectPhasor(indx)));
    vRectWaveform = vRectWaveform + abs(VRectPhasor(indx)) * cos(indx*w*wt + angle(VRectPhasor(indx)));
    iRectWaveform = iRectWaveform + abs(IRectPhasor(indx)) * cos(indx*w*wt + angle(IRectPhasor(indx)));
    iCWaveform = iCWaveform + abs(ICPhasor(indx)) * cos(indx*w*wt + angle(ICPhasor(indx)));
    vLoadWaveform = vLoadWaveform + abs(vLoadPhasor(indx)) * cos(indx*w*wt + angle(vLoadPhasor(indx)));
end


% Verify CCM
if (min(iRectWaveform)>0 && CCM==0)
    CCM=1;
end

% Current and Voltage Ratings
%  Diode Voltage/Current Ratings
iDiodePeak = max(iRectWaveform);
if (iDiodePeak > diodePeakCurrentRating && SWdiodePeakCurrentRating == 0)
    SWdiodePeakCurrentRating = 1;
end
vDiodePeak = max(vRectWaveform);
if (vDiodePeak > diodeReverseVoltageRating && SWdiodeReverseVoltageRating == 0)
    SWdiodeReverseVoltageRating = 1;
end
iDiodeAvg = mean(iRectWaveform);
if (iDiodeAvg > diodeCurrentRating && SWdiodeCurrentRating == 0)
    SWdiodeCurrentRating = 1;
end

% Inductor Current Rating
iLrms = rms(iRectWaveform);
if (iLrms > inductorCurrentRating && SWinductorCurrentRating == 0)
    SWinductorCurrentRating = 1;
end

% Capacitor Voltage and Current Rating
vCPeak = rms(vLoadWaveform);
if (vCPeak > capacitorVoltageRating && SWcapacitorVoltageRating == 0)
    SWcapacitorVoltageRating = 1;
end
iCrms = rms(iCWaveform);
if (iCrms > capacitorCurrentRating && SWcapacitorCurrentRating == 0)
    SWcapacitorCurrentRating = 1;
end



vLoadPkPk = max(vLoadWaveform) - min(vLoadWaveform);
pAvg = (sqrt(3)/pi)*I0*Vln_pk;
pf = 3*pAvg/(3*rms(iRectWaveform)*Vln_rms*sqrt(2/3));
iD1Avg = mean(iRectWaveform)/3;
iD1Pk = max(iRectWaveform);
Pin = mean(iRectWaveform.*vRectWaveform_ideal);
Ploss = rms(iRectWaveform)^2*(2*Rd+R_DCR)+rms(iCWaveform)^2*R_ESR+mean(iRectWaveform)*(2*Vx);
efficiency = (Pin-Ploss)/Pin * 100;
LoadPercentRipple = vLoadPkPk/mean(vLoadWaveform)*100;

    results = [results; R, efficiency, LoadPercentRipple, pf];
end

% Calculate Diode Unit Cost
diodeUnitCost = diodeCost/6;

% Part Number and List Array
diodePartNumbers = {
    'SBR10U300CT-ND';   
    'SBR10U300CT-ND';         
    'SBR10U300CT-ND'       
};

diodeLinks = {
    'https://www.digikey.com/en/products/detail/diodes-incorporated/SBR10U300CT/1778834';   
    'https://www.digikey.com/en/products/detail/diodes-incorporated/SBR10U300CT/1778834';         
    'https://www.digikey.com/en/products/detail/diodes-incorporated/SBR10U300CT/1778834'       
};

inductorPartNumbers = {
    'HM4112-ND';   
    '595-1878-ND';         
    'HM4112-ND'       
};

inductorLinks = {
    'https://www.digikey.com/en/products/detail/hammond-manufacturing/159ZE/455205';   
    'https://www.digikey.com/en/products/detail/signal-transformer/CH-4/1984839';         
    'https://www.digikey.com/en/products/detail/hammond-manufacturing/159ZE/455205'       
};

capacitorPartNumbers = {
    '399-ALS30A471DE350-ND';   
    'ALS31A471DE400';         
    'ALS31A471DE400'       
};

capacitorLinks = {
    'https://www.digikey.com/en/products/detail/kemet/ALS30A471DE350/18166859';   
    'https://www.digikey.com/en/products/detail/kemet/ALS31A471DE400/18196440';         
    'https://www.digikey.com/en/products/detail/kemet/ALS31A471DE400/18196440'       
};

efficiencyArray = {98, 98.5, 98.8};
rippleArray = {1, 0.25, 0.1};
pfArray = {0.89, 0.93, 0.945};



% Printing Design Summary
fprintf('Design %d Summary:\n', DesignNumber);
if CCM
    fprintf('Design Operates in Continous Current Mode\n');
else
    fprintf('Design Operates in Discontinous Current Mode\n');
end
fprintf('─ DIODES ────────────────────────────\n');
fprintf('Part Number:             %s  \n', diodePartNumbers{DesignNumber});
fprintf('Link:                    %s  \n', diodeLinks{DesignNumber});
fprintf('Forward Voltage:         %0.2f V  \n', Vx);
fprintf('Diode Resistance:        %0.2f Ω  \n', Rd);
if SWdiodeReverseVoltageRating
    fprintf('Reverse Voltage Rating:  %0.2f V(Rating Exceeded)  \n', diodeReverseVoltageRating);
else
    fprintf('Reverse Voltage Rating:  %0.2f V  \n', diodeReverseVoltageRating);
end
if SWdiodeCurrentRating
    fprintf('Average Current Rating:  %0.2f A(Rating Exceeded)  \n', diodeCurrentRating);
else
    fprintf('Average Current Rating:  %0.2f A  \n', diodeCurrentRating);
end
if SWdiodePeakCurrentRating
    fprintf('Peak Current Rating:     %0.2f A(Rating Exceeded)  \n', diodePeakCurrentRating);
else
    fprintf('Peak Current Rating:     %0.2f A  \n', diodePeakCurrentRating);
end
fprintf('Cost per unit:           $%0.3f  \n', diodeUnitCost);
fprintf('Total Diode Cost:        $%0.2f  \n\n', diodeCost);


fprintf('─ INDUCTOR ───────────────────────────\n');
fprintf('Part Number:             %s  \n', inductorPartNumbers{DesignNumber});
fprintf('Link:                    %s  \n', inductorLinks{DesignNumber});
fprintf('Inductance:              %0.2e H  \n', L);
fprintf('DCR:                     %0.2e Ω  \n', RDCR);
if SWinductorCurrentRating
    fprintf('Current Rating:          %0.2f A(Rating Exceeded)  \n', inductorCurrentRating);
else
    fprintf('Current Rating:          %0.2f A  \n', inductorCurrentRating);
end
fprintf('Total Inductor Cost:     $%0.2f  \n\n', inductorCost);

fprintf('─ Capacitor ───────────────────────────\n');
fprintf('Part Number:             %s  \n', capacitorPartNumbers{DesignNumber});
fprintf('Link:                    %s  \n', capacitorLinks{DesignNumber});
fprintf('Capacitance:             %0.2e F  \n', C);
fprintf('ECR:                     %0.2e Ω  \n', RESR);
if SWcapacitorVoltageRating
    fprintf('Voltage Rating:          %0.2f V(Rating Exceeded)  \n', capacitorVoltageRating);
else
    fprintf('Voltage Rating:          %0.2f V  \n', capacitorVoltageRating);
end
if SWcapacitorCurrentRating
    fprintf('Current Rating:          %0.2f A(Rating Exceeded)  \n', capacitorCurrentRating);
else
    fprintf('Current Rating:          %0.2f A  \n', capacitorCurrentRating);
end
fprintf('Total Capacitor Cost:    $%0.2f  \n\n', capacitorCost);

fprintf('─ General Summary ───────────────────────────\n');
fprintf('Total Component Cost:    $%0.2f  \n', diodeCost + capacitorCost + inductorCost);
if min(results(:,2)) > efficiencyArray{DesignNumber}
    fprintf('Efficiency is > %0.2f%% for all load conditions\n', efficiencyArray{DesignNumber});
else
    fprintf('Efficiency is not > %0.2f%% for all load conditions\n', efficiencyArray{DesignNumber});
end
if min(results(:,3)) < rippleArray{DesignNumber}
    fprintf('Load voltage ripple is < %0.2f%% for all load conditions\n', rippleArray{DesignNumber});
else
    fprintf('Load voltage ripple is not < %0.2f%% for all load conditions\n', rippleArray{DesignNumber});
end
if min(results(:,4)) > pfArray{DesignNumber}
    fprintf('AC source power factor is > %0.2f for all load conditions\n', pfArray{DesignNumber});
else
    fprintf('AC source power factor is not > %0.3f for all load conditions\n', pfArray{DesignNumber});
end




% Plots
subplot(1,3,1);
plot(results(:,1),results(:,2),'linewidth',1)
ylabel('% Efficiency')
xlabel('Resistance(R)')
grid on
axis tight

subplot(1,3,2);
plot(results(:,1),results(:,3),'linewidth',1)
ylabel('% Load Voltage Ripple')
xlabel('Resistance(R)')
grid on
axis tight
title("Design " + DesignNumber, 'FontSize', 14, 'Color','black');


subplot(1,3,3);
plot(results(:,1),results(:,4),'linewidth',1)
ylabel('Power Factor')
xlabel('Resistance(R)')
grid on
axis tight
