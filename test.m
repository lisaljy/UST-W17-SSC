clc;
clear;

snapshots = 100;
sensorsNumber = 9;
overlapSubArrayNum = sensorsNumber/3;
SubArraySize = sensorsNumber - overlapSubArrayNum +1;
array_1 = 0;
array_2 = 10;
array_3 = 20;
rou = [0.7+0.7j ; 0.6+0.5j ; 0.2+0.4j];
SNRMinimun = -10;
SNRMaximun = 10;
x = SNRMinimun:1:SNRMaximun;
MonteNum = 1000;
ProbabilityAIC = [];
ProbabilityMDL = [];
ProbabilityBIC = [];
ProbabilityHQ = [];

% Signal Model
% X = A * rou * Signal_0 + Noise;

% Construct A
Omega_1 = pi * sind(array_1);
Omega_2 = pi * sind(array_2);
Omega_3 = pi * sind(array_3);

a1 = [];
a2 = [];
a3 = [];
for k = 0:1:(sensorsNumber - 1);
    a1 = [a1;exp(-1j * k * Omega_1)];
    a2 = [a2;exp(-1j * k * Omega_2)];
    a3 = [a3;exp(-1j * k * Omega_3)];
end

A = [a1,a2,a3];

% Construct X under different SNR
% for SNR = SNRMinimun:SNRMaximun
SNR=50;
Right_AIC = 0;
Right_MDL = 0;
Right_BIC = 0;
Right_HQ = 0;

    for MonteNumTest = 1:MonteNum

        % Genrate Signal_0 S
%         SigamaSignalSquare = 10^(SNR/10);    
%         Signal_0 = sqrt(SigamaSignalSquare/2) * randn(1,snapshots) + j * sqrt(SigamaSignalSquare/2) * rand(1,snapshots);

        % Assume 50HZ 
        t = randn(1,snapshots);

        Signal_0 = exp(j*2*pi*50*t);

%         Signal_0 = sin(2*pi*50*t);

        ATimesRouTimesSignal_0 = A * rou * Signal_0;
        
        % Genrate Noise N
                   
%         Noise = sqrt(1/2) * randn(sensorsNumber,snapshots) + j * sqrt(1/2) * rand(sensorsNumber,snapshots);
%         X = ATimesRouTimesSignal_0 + Noise;

        X = awgn(ATimesRouTimesSignal_0,SNR,'measured');

%         X = ATimesRouTimesSignal_0;
    
        XH = X';

        % Acquired correlation matrix R from simulation data
        R = (X * XH) / snapshots;
        %Since under this circumstance, overlapSubArrayNum = 3
        Rk_1 = R(1:SubArraySize , 1:SubArraySize);
        Rk_2 = R(2:SubArraySize+1 , 2:SubArraySize+1);
        Rk_3 = R(3:SubArraySize+2 , 3:SubArraySize+2);
        Rf = (Rk_1+Rk_2+Rk_3) / 3;
        Rb = ((Rk_1).'+(Rk_2).'+(Rk_3).') / 3;
        R2 = (Rf + Rb)/2;
        
        % SVD
        eigenvalue = svd (R); 
        eigenvalue2 = svd (R2); 

        %AIC
        %MDL
        AIC = [];
        MDL = [];
        BIC = [];
        HQ = [];
        for n=0:1:(2*sensorsNumber/3)

           %likelihood function
           sumS = sum(eigenvalue(n+1:2*sensorsNumber/3+1,:),1);
           prodS = prod(eigenvalue(n+1:2*sensorsNumber/3+1,:),1);
           likelihoodFunction = ( sumS/(2*sensorsNumber/3-n) )/( prodS^( 1/(2*sensorsNumber/3-n) ) );
%            likelihoodFunction = ( sumS/(2*sensorsNumber/3+1-n) )/( prodS^( 1/(2*sensorsNumber/3+1-n) ) );
           
%            AIC = [AIC,2*snapshots*(2*sensorsNumber/3+1-n)*log(likelihoodFunction) + 2*n*(2*2*sensorsNumber/3+1-n)];
           AIC = [AIC,2*snapshots*(sensorsNumber-n)*log(likelihoodFunction) + 2*n*(2*sensorsNumber-n)];
%            MDL = [MDL,snapshots*(sensorsNumber-n)*log(likelihoodFunction) + 0.5*n*log(snapshots)*(2*sensorsNumber-n)];
%            BIC = [BIC,2*snapshots*(sensorsNumber-n)*log(likelihoodFunction) + n*log(snapshots)*(2*sensorsNumber-n)];
%            HQ = [HQ,snapshots*(sensorsNumber-n)*log(likelihoodFunction) + 0.5*n*log(log(snapshots))*(2*sensorsNumber-n)];

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
        if n_AIC == 4
            Right_AIC = Right_AIC + 1;
        end
        
%         MDLmin = min(MDL);
%         n_MDL = find(MDL == MDLmin);
%         if n_MDL == 4
%             Right_MDL = Right_MDL + 1;
%         end  
%         
%         BICmin = min(BIC);
%         n_BIC = find(BIC == BICmin);
%         if n_BIC == 4
%             Right_BIC = Right_BIC + 1;
%         end  
%         
%         HQmin = min(HQ);
%         n_HQ = find(HQ == HQmin);
%         if n_HQ == 4
%             Right_HQ = Right_HQ + 1;
%         end  
        
    end
    
    ProbabilityAIC = [ProbabilityAIC,(Right_AIC/MonteNum)*100];
%     ProbabilityMDL = [ProbabilityMDL,(Right_MDL/MonteNum)*100];
%     ProbabilityBIC = [ProbabilityBIC,(Right_BIC/MonteNum)*100];
%     ProbabilityHQ = [ProbabilityHQ,(Right_HQ/MonteNum)*100];
    
% end

% hold on;
% title('Detection Probability of non-coherent signal versus SNR');
% xlabel('SNR(dB)');
% ylabel('Detection Probability');
% p1 = plot(x,ProbabilityAIC,'b');
% p2 = plot(x,ProbabilityMDL,'*r');
% p3 = plot(x,ProbabilityBIC,'g');
% p4 = plot(x,ProbabilityHQ,'k');
% legend([p1 p2 p3 p4],'AIC','MDL','BIC','HQ');
% hold off;
