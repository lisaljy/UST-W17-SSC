clc;
clear;

snapshots = 100;
sensorsNumber = 10;
array_1 = 0;
array_2 = 10;
SNRMinimun = -15;
SNRMaximun = 30;
x = SNRMinimun:1:SNRMaximun;
MonteNum = 1000;
ProbabilityAIC = [];
ProbabilityMDL = [];

% Signal Model
% X = A * Signal + Noise;

% Construct A
Omega_1 = pi * sind(array_1);
Omega_2 = pi * sind(array_2);

a1 = [];
a2 = [];
for k = 0:1:(sensorsNumber - 1);
    a1 = [a1;exp(-1j * k * Omega_1)];
    a2 = [a2;exp(-1j * k * Omega_2)];
end

A = [a1,a2];

% Construct X under different SNR
for SNR = SNRMinimun:SNRMaximun
% SNR=30;
Right_AIC = 0;
Right_MDL = 0;

    for MonteNumTest = 1:MonteNum

        % Genrate Signal S
    %     SigamaSignalSquare = 10^(SNR/10);    
    %     Signal_1 = sqrt(SigamaSignalSquare/2) * randn(1,snapshots) + j * sqrt(SigamaSignalSquare/2) * rand(1,snapshots);
    %     Signal_2 = sqrt(SigamaSignalSquare/2) * randn(1,snapshots) + j * sqrt(SigamaSignalSquare/2) * rand(1,snapshots);

        % Assume 50HZ and 100HZ
        t = randn(1,snapshots);

        Signal_1 = exp(j*2*pi*50*t);
        Signal_2 = exp(j*2*pi*100*t);

    %     Signal_1 = sin(2*pi*50*t);
    %     Signal_2 = sin(2*pi*100*t);

        Signal = [Signal_1;Signal_2];

        % Genrate Noise N
        ATimesSignal = A * Signal;    
        
    %     Noise = sqrt(1/2) * randn(sensorsNumber,snapshots) + j * sqrt(1/2) * rand(sensorsNumber,snapshots);
    %     X = ATimesSignal + Noise;

        X = awgn(ATimesSignal,SNR,'measured');

    %     X = ATimesSignal;
    
        XH = X';

        % Acquired correlation matrix R from simulation data
        R = (X * XH) / snapshots;

        % SVD
        eigenvalue = svd (R); 

        %AIC
        %MDL
        AIC = [];
        MDL = [];
        for n=0:1:(sensorsNumber - 1)

           %likelihood function
           sumS = sum(eigenvalue(n+1:sensorsNumber,:),1);
           prodS = prod(eigenvalue(n+1:sensorsNumber,:),1);
           likelihoodFunction = ( sumS/(sensorsNumber-n) )/( prodS^( 1/(sensorsNumber-n) ) );

           AIC = [AIC,2*snapshots*(sensorsNumber-n)*log(likelihoodFunction) + 2*n*(2*sensorsNumber-n)];
           MDL = [MDL,snapshots*(sensorsNumber-n)*log(likelihoodFunction) + 0.5*n*log(snapshots)*(2*sensorsNumber-n)];
    %       AIC = [AIC,snapshots*(sensorsNumber-n)*log(likelihoodFunction) + n*(2*sensorsNumber-n)];

        end
        
        %EGM
%         eigenvalueDescend = sort(eigenvalue,1,'descend');
%         eigenvalueDetaAve = eigenvalueDescend(1,1) - eigenvalueDescend(sensorsNumber,1);
%         eigenvalueDetaAve = eigenvalueDetaAve/(sensorsNumber -1);
%         eigenvalueDeta = [];
%         for i=1:sensorsNumber-1;
%             eigenvalueDeta = [eigenvalueDeta;eigenvalueDescend(i,1) - eigenvalueDescend(i+1,1)];
%         end
%         index = find(eigenvalueDeta <= eigenvalueDetaAve);
%         eigenvalueDetaSmallerThanEigenvalueDetaAve = eigenvalueDeta(index);
%         ............??
        
        %Detection Probability
        AICmin = min(AIC);
        n_AIC = find(AIC == AICmin);
        if n_AIC == 3
            Right_AIC = Right_AIC + 1;
        end
        
        MDLmin = min(MDL);
        n_MDL = find(MDL == MDLmin);
        if n_MDL == 3
            Right_MDL = Right_MDL + 1;
        end                 
        
    end
    
    ProbabilityAIC = [ProbabilityAIC,(Right_AIC/MonteNum)*100];
    ProbabilityMDL = [ProbabilityMDL,(Right_MDL/MonteNum)*100];
    
end

hold on;
title('Detection Probability of non-coherent signal versus SNR');
xlabel('SNR(dB)');
ylabel('Detection Probability');
p1 = plot(x,ProbabilityAIC,'b');
p2 = plot(x,ProbabilityMDL,'r');
legend([p1 p2],'AIC','MDL');
hold off;
