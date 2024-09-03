% Aydin Roozbeh - 9923037 - Coding 2 - BPSK AWGN Channel

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

% Maximum Number of iterations
MI = 5;

% Number of simulations to run
sim = 5;

% Dimensions
n = size(H1,2);
k = n-size(H1,1);

% Rate
Rate = k/n;

% SNR
SNR_db = 1:0.5:6;          % SNR in Decibels
SNR_abs = 10.^(SNR_db/10); % SNR in absolute value
SNR_num = size(SNR_db,2);   % Number of different SNRs to be simluated

% Average number of bit and frame errors for each SNR value, for each PCM
H1_def = zeros(3,SNR_num);
H2_def = zeros(3,SNR_num);
H3_def = zeros(3,SNR_num);

% Zero code word
data = zeros(1,12);
data_demod = data;

% Running multiple simulations
for k=1:1:sim
    % ------------- H1 Parity Chech matrix --------------
    for i=1:1:SNR_num
    
        % Frame and Bit errors
        eb = 0;
        ef = 0;
        count = 0;
    
        while(ef<50)
            count=count+1;
            TX = pskmod(data,2); % Baseband BPSK Modulation
            RX = awgn(TX , SNR_db(i)); % Adding AWGN 
            data_demod = pskdemod(RX,2,0); % Baseband BPSK Demodulation
            data_dec = decoder_awgn(MI , H1 , data_demod , SNR_db(i)); % Decoding
           
        end
    
        H1_def(1,i)=H1_def(1,i)+eb;
        H1_def(2,i)=H1_def(2,i)+ef;
        H1_def(3,i)=H1_def(3,i)+count;
    
    end
    
    % ------------- H2 Parity Chech matrix --------------
    for i=1:1:SNR_num
    
        % Frame and Bit errors
        eb = 0;
        ef = 0;
        count = 0;
    
        while(ef<50)
            count=count+1;
            TX = pskmod(data,2,0); % Baseband BPSK Modulation
            RX = awgn(TX , SNR_db(i)); % Adding AWGN 
            data_demod = pskdemod(RX,2,0); % Baseband BPSK Demodulation
            data_dec = decoder_awgn(MI , H2 , data_demod , SNR_db(i)); % Decoding
            
            if(data_dec ~= data)    % Counting the number defect bits and frames
                eb = eb + sum(xor(data,data_dec));
                ef = ef + 1;
            end
        end
    
        H2_def(1,i)=H2_def(1,i)+eb;
        H2_def(2,i)=H2_def(2,i)+ef;
        H2_def(3,i)=H2_def(3,i)+count;
    
    end
    
    % ------------- H3 Parity Chech matrix --------------
    for i=1:1:SNR_num
    
        % Frame and Bit errors
        eb = 0;
        ef = 0;
        count = 0;
    
        while(ef<50)
            count=count+1;
            TX = pskmod(data,2,0); % Baseband BPSK Modulation
            RX = awgn(TX , SNR_db(i)); % Adding AWGN 
            data_demod = pskdemod(RX,2,0); % Baseband BPSK Demodulation
            data_dec = decoder_awgn(MI , H3 , data_demod , SNR_db(i)); % Decoding
            
            if(data_dec ~= data)    % Counting the number defect bits and frames
                eb = eb + sum(xor(data,data_dec));
                ef = ef + 1;
            end
        end
    
        H3_def(1,i)=H3_def(1,i)+eb;
        H3_def(2,i)=H3_def(2,i)+ef;
        H3_def(3,i)=H3_def(3,i)+count;
    
    end
end

% Averaging over all simulations
H1_def = H1_def / sim;
H2_def = H2_def / sim;
H3_def = H3_def / sim;


% Plotting the results;

figure(1)
title("Bit Error Rate");
hold on;
plot(SNR_db , H1_def(1,:)./(H1_def(3,:)*12) , "LineStyle","-" , Color='green');
plot(SNR_db , H2_def(1,:)./(H2_def(3,:)*12) , "LineStyle","--" , Color='blue');
plot(SNR_db , H3_def(1,:)./(H3_def(3,:)*12) , "LineStyle","-." , Color='red');
legend("H1" , "H2" , "H3");
grid on;    

figure(2)
title("Frame Error Rate");
hold on;
plot(SNR_db , H1_def(2,:)./H1_def(3,:) , "LineStyle","-" , Color='green');
plot(SNR_db , H2_def(2,:)./H2_def(3,:) , "LineStyle","--" , Color='blue');
plot(SNR_db , H3_def(2,:)./H3_def(3,:) , "LineStyle","-." , Color='red');
legend("H1" , "H2" , "H3");
grid on;