% Aydin Roozbeh - 9923037 - Coding2 project
% Main code 
close all;
clear;
clc;

% Defining The parity check matrices
H1= [1 0 0 0 1 0 0 0 0 0 0 1;
     0 0 0 1 0 1 0 0 0 0 0 1;
     0 1 0 0 1 0 0 0 0 0 1 0;
     0 0 1 0 0 0 1 0 0 0 1 0;
     0 1 0 0 0 0 0 1 0 1 0 0;
     0 0 1 0 0 1 0 0 0 1 0 0;
     0 0 0 1 0 0 0 1 1 0 0 0;
     1 0 0 0 0 0 1 0 1 0 0 0];

H2= [1 0 0 0 1 0 1 0 0 0 0 1;
     0 0 0 1 0 1 0 0 0 0 0 1;
     0 1 0 0 1 0 0 0 0 0 1 0;
     0 0 1 0 0 0 1 0 0 0 1 0;
     1 1 0 0 0 0 0 1 0 1 0 0;
     0 0 1 0 0 1 0 0 0 1 0 0;
     0 0 0 1 0 0 0 1 1 0 0 0;
     1 0 0 0 0 0 1 0 1 0 0 0];

H3= [1 0 0 0 1 0 0 0 0 0 0 1;
     0 0 0 1 0 1 0 0 0 0 0 1;
     0 1 0 0 1 0 1 0 0 0 1 0;
     0 0 1 0 0 0 1 0 0 0 1 0;
     1 1 0 0 1 0 0 1 0 1 0 0;
     0 0 1 0 0 1 0 0 0 1 0 0;
     0 0 0 1 0 0 0 1 1 0 0 0;
     1 0 0 0 0 0 1 0 1 0 0 0];

% Number of simulations to run
sim = 10;

% Dimensions
n = size(H1,2);
k = n-size(H1,1);

% Rate
Rate = k/n;

% SNR
SNR_db = 1:0.5:6;          % SNR in Decibels
SNR_abs = 10.^(SNR_db/10); % SNR in absolute value
N0 = 1./(Rate.*SNR_abs);   % N0
var = N0/2;                % Noise sample variance based on N0
std = sqrt(var);           % Standard Deviation based on variance

% Bit Mapping: (D=Data , S=Symbol)
% 0 => +1
% 1 => -1
% Data : D = 0/1 , Symbol : S = 1/-1   |   S=-2D+1 , D = (1-S)/2

% Data :
MI = [1 5 10 30 100];
word = 10000;                       % Number of words to be sent
Y_D = zeros(word,12);               % Data bits
Y_S = -2*Y_D + 1;                   % Data to Symbol
Y_dec = zeros(word,12);             % Decoded bits
total = word * 12;

% Average number of bit errors for each SNR value, for each PCM
H1_def = zeros(5,2,11);
H2_def = zeros(5,2,11);
H3_def = zeros(5,2,11);

% ---------------- for H1 ---------------------
% For all SNRs
for m=1:1:size(MI,2)
    for s=1:1:size(SNR_db,2)
        % Averaging over 10 executions
        for a=1:1:sim
    
            % Changing the seed everytime the code runs
            reset(RandStream.getGlobalStream,sum(100*clock));
            noise = std(s)*randn(word,12);      % noise
            Yn = Y_S + noise;                   % Symbol + Noise
    
            for i=1:1:word
               Y_dec(i,:) = decoder_awgn(MI(m) , H3 , Yn(i,:) , SNR_db(s));
            end
    
            % Calculating the bit and frame error rate
            def_bit = sum(sum(double(xor(Y_D , Y_dec))));
            def_frame = cal_frame_error(Y_D , Y_dec);
            
            H1_def(m,1,s)=H1_def(m,1,s)+def_bit;
            H1_def(m,2,s)=H1_def(m,2,s)+def_frame;
        end
    end
end
H1_def = H1_def / sim;

% ---------------- for H2 ---------------------
% For all SNRs
for m=1:1:size(MI,2)
    for s=1:1:size(SNR_db,2)
        % Averaging over 10 executions
        for a=1:1:sim
    
            % Changing the seed everytime the code runs
            reset(RandStream.getGlobalStream,sum(100*clock));
            noise = std(s)*randn(word,12);      % noise
            Yn = Y_S + noise;                   % Symbol + Noise
    
            for i=1:1:word
               Y_dec(i,:) = decoder_awgn(MI(m) , H2 , Yn(i,:) , SNR_db(s));
            end
    
            % Calculating the bit and frame error rate
            def_bit = sum(sum(double(xor(Y_D , Y_dec))));
            def_frame = cal_frame_error(Y_D , Y_dec);
            
            H2_def(m,1,s)=H2_def(m,1,s)+def_bit;
            H2_def(m,2,s)=H2_def(m,2,s)+def_frame;
        end
    end
end
H2_def = H2_def / sim;

% ----------------- for H3 ---------------------
% For all SNRs
for m=1:1:size(MI,2)
    for s=1:1:size(SNR_db,2)
        % Averaging over 10 executions
        for a=1:1:sim
    
            % Changing the seed everytime the code runs
            reset(RandStream.getGlobalStream,sum(100*clock));
            noise = std(s)*randn(word,12);      % noise
            Yn = Y_S + noise;                   % Symbol + Noise
    
            for i=1:1:word
               Y_dec(i,:) = decoder_awgn(MI(m) , H3 , Yn(i,:) , SNR_db(s));
            end
    
            % Calculating the bit and frame error rate
            def_bit = sum(sum(double(xor(Y_D , Y_dec))));
            def_frame = cal_frame_error(Y_D , Y_dec);
            
            H3_def(m,1,s)=H3_def(m,1,s)+def_bit;
            H3_def(m,2,s)=H3_def(m,2,s)+def_frame;
        end
    end
end
H3_def = H3_def / sim;

% --------- plotting --------------
figure(1)
subplot(1,2,1);
title("BER for H1");
hold on;
plot(SNR_db , squeeze(H1_def(1,1,:))/(word*12) , 'o-');
plot(SNR_db , squeeze(H1_def(2,1,:))/(word*12) , 'b--o');
plot(SNR_db , squeeze(H1_def(3,1,:))/(word*12) , 'c-*' );
plot(SNR_db , squeeze(H1_def(4,1,:))/(word*12) , 'g--');
plot(SNR_db , squeeze(H1_def(5,1,:))/(word*12) , 'r-.');
legend("MI=1","MI=5","MI=10","MI=30","MI=100")
grid on;

subplot(1,2,2);
title("FER for H1");
hold on;
plot(SNR_db , squeeze(H1_def(1,1,:))/word , 'o-');
plot(SNR_db , squeeze(H1_def(2,1,:))/word , 'b--o');
plot(SNR_db , squeeze(H1_def(3,1,:))/word , 'c-*' );
plot(SNR_db , squeeze(H1_def(4,1,:))/word , 'g--');
plot(SNR_db , squeeze(H1_def(5,1,:))/word , 'r-.');
legend("MI=1","MI=5","MI=10","MI=30","MI=100")
grid on;

%---------------------

figure(2)
subplot(1,2,1);
title("BER for H2");
hold on;
plot(SNR_db , squeeze(H2_def(1,1,:))/(word*12) , 'o-');
plot(SNR_db , squeeze(H2_def(2,1,:))/(word*12) , 'b--o');
plot(SNR_db , squeeze(H2_def(3,1,:))/(word*12) , 'c-*' );
plot(SNR_db , squeeze(H2_def(4,1,:))/(word*12) , 'g--');
plot(SNR_db , squeeze(H2_def(5,1,:))/(word*12) , 'r-.');
legend("MI=1","MI=5","MI=10","MI=30","MI=100")
grid on;

subplot(1,2,2);
title("FER for H2");
hold on;
plot(SNR_db , squeeze(H2_def(1,1,:))/word , 'o-');
plot(SNR_db , squeeze(H2_def(2,1,:))/word , 'b--o');
plot(SNR_db , squeeze(H2_def(3,1,:))/word , 'c-*' );
plot(SNR_db , squeeze(H2_def(4,1,:))/word , 'g--');
plot(SNR_db , squeeze(H2_def(5,1,:))/word , 'r-.');
legend("MI=1","MI=5","MI=10","MI=30","MI=100")
grid on;

% ------------------

figure(3)
subplot(1,2,1);
title("BER for H3")
hold on;
plot(SNR_db , squeeze(H3_def(1,1,:))/(word*12) , 'o-');
plot(SNR_db , squeeze(H3_def(2,1,:))/(word*12) , 'b--o');
plot(SNR_db , squeeze(H3_def(3,1,:))/(word*12) , 'c-*' );
plot(SNR_db , squeeze(H3_def(4,1,:))/(word*12) , 'g--');
plot(SNR_db , squeeze(H3_def(5,1,:))/(word*12) , 'r-.');
legend("MI=1","MI=5","MI=10","MI=30","MI=100")
grid on;

subplot(1,2,2);
title("FER for H3");
hold on;
plot(SNR_db , squeeze(H3_def(1,1,:))/word , 'o-');
plot(SNR_db , squeeze(H3_def(2,1,:))/word , 'b--o');
plot(SNR_db , squeeze(H3_def(3,1,:))/word , 'c-*' );
plot(SNR_db , squeeze(H3_def(4,1,:))/word , 'g--');
plot(SNR_db , squeeze(H3_def(5,1,:))/word , 'r-.');
legend("MI=1","MI=5","MI=10","MI=30","MI=100")
grid on;
%---------
