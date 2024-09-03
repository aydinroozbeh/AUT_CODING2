% Aydin Roozbeh - 9923037 - Coding2 project
% Main code 
close all;
clear;
clc;
tic;
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

% Maximum Number of iterations
MI = 5;

% Number of simulations to run
sim = 20;

% Dimensions
n = size(H1,2);
k = n-size(H1,1);

% Rate
Rate = k/n;

% SNR
SNR_db = 1:0.5:6;          % SNR in Decibels
SNR_abs = 10.^(SNR_db/10); % SNR in absolute value
SNR_num = size(SNR_db,2);   % Number of different SNRs to be simluated
N0 = 1./(Rate.*SNR_abs);   % N0
var = N0/2;                % Noise sample variance based on N0
std = sqrt(var);           % Standard Deviation based on variance

% Bit Mapping: (D=Data , S=Symbol)
% 0 => +1
% 1 => -1
% Data : D = 0/1 , Symbol : S = 1/-1   |   S=-2D+1 , D = (1-S)/2

% Data :

word = 1000;                       % Number of words to be sent
Y_D = zeros(word,12);               % Data bits
Y_S = -2*Y_D + 1;                   % Data to Symbol
Y_dec = zeros(word,12);             % Decoded bits
total = word * 12;

% Average number of bit and frame errors for each SNR value, for each PCM
H1_def = zeros(2,SNR_num);
H2_def = zeros(2,SNR_num);
H3_def = zeros(2,SNR_num);

% ---------------- for H1 ---------------------

% For all SNRs
for s=1:1:SNR_num
    % Averaging over 10 executions
    for a=1:1:sim

        % Changing the seed everytime the code runs
        reset(RandStream.getGlobalStream,sum(100*clock));
        noise = std(s)*randn(word,12);      % noise
        Yn = Y_S + noise;                   % Symbol + Noise

        for i=1:1:word
           Y_dec(i,:) = decoder_awgn(MI , H1 , Yn(i,:) , SNR_db(s));
        end

        % Calculating the bit and frame error rate
        def_bit = sum(sum(double(xor(Y_D , Y_dec))));
        def_frame = cal_frame_error(Y_D , Y_dec);
        
        H1_def(1,s)=H1_def(1,s)+def_bit;
        H1_def(2,s)=H1_def(2,s)+def_frame;
    end
end
H1_def = H1_def / sim;

% ---------------- for H2 ---------------------

% For all SNRs
for s=1:1:SNR_num
    % Averaging over 10 executions
    for a=1:1:sim

        % Changing the seed everytime the code runs
            
        noise = std(s)*randn(word,12);      % noise
        Yn = Y_S + noise;                   % Symbol + Noise

        for i=1:1:word
           Y_dec(i,:) = decoder_awgn(MI , H2 , Yn(i,:) , SNR_db(s));
        end

        % Calculating the bit and frame error rate
        def_bit = sum(sum(double(xor(Y_D , Y_dec))));
        def_frame = cal_frame_error(Y_D , Y_dec);
        
        H2_def(1,s)=H2_def(1,s)+def_bit;
        H2_def(2,s)=H2_def(2,s)+def_frame;
    end
end
H2_def = H2_def / sim;

% ---------------- for H3 ---------------------

% For all SNRs
for s=1:1:SNR_num
    % Averaging over 10 executions
    for a=1:1:sim

        % Changing the seed everytime the code runs
        reset(RandStream.getGlobalStream,sum(100*clock));
        noise = std(s)*randn(word,12);      % noise
        Yn = Y_S + noise;                   % Symbol + Noise

        for i=1:1:word
           Y_dec(i,:) = decoder_awgn(MI , H1 , Yn(i,:) , SNR_db(s));
        end

        % Calculating the bit and frame error rate
        def_bit = sum(sum(double(xor(Y_D , Y_dec))));
        def_frame = cal_frame_error(Y_D , Y_dec);
        
        H3_def(1,s)=H3_def(1,s)+def_bit;
        H3_def(2,s)=H3_def(2,s)+def_frame;
    end
end
H3_def = H3_def / sim;

% Plotting the results;

figure(1)
title("Bit Error Rate");
hold on;
plot(SNR_db , H1_def(1,:)/total , "LineStyle","-" , Color='green');
plot(SNR_db , H2_def(1,:)/total , "LineStyle","--" , Color='blue');
plot(SNR_db , H3_def(1,:)/total , "LineStyle","-." , Color='red');
legend("H1" , "H2" , "H3");
grid on;    

figure(2)
title("Frame Error Rate");
hold on;
plot(SNR_db , H1_def(2,:)/word , "LineStyle","-" , Color='green');
plot(SNR_db , H2_def(2,:)/word , "LineStyle","--" , Color='blue');
plot(SNR_db , H3_def(2,:)/word , "LineStyle","-." , Color='red');
legend("H1" , "H2" , "H3");
grid on;

toc