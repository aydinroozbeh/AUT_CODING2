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
sim = 30;

% Dimensions
n = size(H1,2);
k = n-size(H1,1);

% Rate
Rate = k/n;

% SNR
SNR_db = 1:0.5:6;          % SNR in Decibels
SNR_abs = 10.^(SNR_db/10); % SNR in absolute value
SNR_num = size(SNR_db,2);   % Number of different SNRs to be simluated

MI = [1 5 10 30 100];
MI_num = size(MI,2);

% Average number of bit errors for each SNR value, for each PCM
H1_def = zeros(5,3,11);
H2_def = zeros(5,3,11);
H3_def = zeros(5,3,11);

% Zero code word
data = zeros(1,12);
data_demod = data;

% ---------------- for H1 ---------------------
% For All MI's
for m=1:1:MI_num
    % For all SNRs
    for s=1:1:SNR_num
        % Averaging over 10 executions
        for a=1:1:sim
            % Frame and Bit errors
            eb = 0;
            ef = 0;
            count = 0;
        
            while(ef<50)
                reset(RandStream.getGlobalStream,sum(100*clock));
                S=randi(10000);

                count=count+1;
                TX = pskmod(data,2,0); % Baseband BPSK Modulation
                RX = awgn(TX , SNR_db(s),1,S); % Adding AWGN 
                data_demod = pskdemod(RX,2,0); % Baseband BPSK Demodulation
                data_dec = decoder_awgn(MI(m) , H1 , data_demod , SNR_db(s)); % Decoding
                
                if(data_dec ~= data)    % Counting the number defect bits and frames
                    eb = eb + sum(xor(data,data_dec));
                    ef = ef + 1;
                end
            end
            H1_def(m,1,s)=H1_def(m,1,s)+eb;
            H1_def(m,2,s)=H1_def(m,2,s)+ef;
            H1_def(m,3,s)=H1_def(m,3,s)+count;
        end
    end
end


% ---------------- for H2 ---------------------
% For All MI's
for m=1:1:MI_num
    % For all SNRs
    for s=1:1:SNR_num
        % Averaging over 10 executions
        for a=1:1:sim
            % Frame and Bit errors
            eb = 0;
            ef = 0;
            count = 0;
        
            while(ef<50)
                reset(RandStream.getGlobalStream,sum(100*clock));
                S=randi(10000);

                count=count+1;
                TX = pskmod(data,2,0); % Baseband BPSK Modulation
                RX = awgn(TX , SNR_db(s),1,S); % Adding AWGN 
                data_demod = pskdemod(RX,2,0); % Baseband BPSK Demodulation
                data_dec = decoder_awgn(MI(m) , H2 , data_demod , SNR_db(s)); % Decoding
                
                if(data_dec ~= data)    % Counting the number defect bits and frames
                    eb = eb + sum(xor(data,data_dec));
                    ef = ef + 1;
                end
            end
            H2_def(m,1,s)=H2_def(m,1,s)+eb;
            H2_def(m,2,s)=H2_def(m,2,s)+ef;
            H2_def(m,3,s)=H2_def(m,3,s)+count;
        end
    end
end

% ---------------- for H1 ---------------------
% For All MI's
for m=1:1:MI_num
    % For all SNRs
    for s=1:1:SNR_num
        % Averaging over 10 executions
        for a=1:1:sim
            % Frame and Bit errors
            eb = 0;
            ef = 0;
            count = 0;
        
            while(ef<50)
                reset(RandStream.getGlobalStream,sum(100*clock));
                S=randi(10000);

                count=count+1;
                TX = pskmod(data,2,0); % Baseband BPSK Modulation
                RX = awgn(TX , SNR_db(s),1,S); % Adding AWGN 
                data_demod = pskdemod(RX,2,0); % Baseband BPSK Demodulation
                data_dec = decoder_awgn(MI(m) , H3 , data_demod , SNR_db(s)); % Decoding
                
                if(data_dec ~= data)    % Counting the number defect bits and frames
                    eb = eb + sum(xor(data,data_dec));
                    ef = ef + 1;
                end
            end
            H3_def(m,1,s)=H3_def(m,1,s)+eb;
            H3_def(m,2,s)=H3_def(m,2,s)+ef;
            H3_def(m,3,s)=H3_def(m,3,s)+count;
        end
    end
end


% --------- plotting --------------
figure(1)
subplot(1,2,1);
title("BER for H1");
hold on;
plot(SNR_db , squeeze(H1_def(1,1,:))./ squeeze(H1_def(1,3,:)*12) , 'o-');
plot(SNR_db , squeeze(H1_def(2,1,:))./ squeeze(H1_def(2,3,:)*12) , 'b--o');
plot(SNR_db , squeeze(H1_def(3,1,:))./ squeeze(H1_def(3,3,:)*12) , 'c-*' );
plot(SNR_db , squeeze(H1_def(4,1,:))./ squeeze(H1_def(4,3,:)*12) , 'g--');
plot(SNR_db , squeeze(H1_def(5,1,:))./ squeeze(H1_def(5,3,:)*12) , 'r-.');
legend("MI=1","MI=5","MI=10","MI=30","MI=100")
grid on;

subplot(1,2,2);
title("FER for H1");
hold on;
plot(SNR_db , squeeze(H1_def(1,2,:))./squeeze(H1_def(1,3,:)) , 'o-');
plot(SNR_db , squeeze(H1_def(2,2,:))./squeeze(H1_def(2,3,:)) , 'b--o');
plot(SNR_db , squeeze(H1_def(3,2,:))./squeeze(H1_def(3,3,:)) , 'c-*' );
plot(SNR_db , squeeze(H1_def(4,2,:))./squeeze(H1_def(4,3,:)) , 'g--');
plot(SNR_db , squeeze(H1_def(5,2,:))./squeeze(H1_def(5,3,:)) , 'r-.');
legend("MI=1","MI=5","MI=10","MI=30","MI=100")
grid on;

% --------- plotting --------------
figure(2)
subplot(1,2,1);
title("BER for H2");
hold on;
plot(SNR_db , squeeze(H2_def(1,1,:))./ squeeze(H2_def(1,3,:)*12) , 'o-');
plot(SNR_db , squeeze(H2_def(2,1,:))./ squeeze(H2_def(2,3,:)*12) , 'b--o');
plot(SNR_db , squeeze(H2_def(3,1,:))./ squeeze(H2_def(3,3,:)*12) , 'c-*' );
plot(SNR_db , squeeze(H2_def(4,1,:))./ squeeze(H2_def(4,3,:)*12) , 'g--');
plot(SNR_db , squeeze(H2_def(5,1,:))./ squeeze(H2_def(5,3,:)*12) , 'r-.');
legend("MI=1","MI=5","MI=10","MI=30","MI=100")
grid on;

subplot(1,2,2);
title("FER for H2");
hold on;
plot(SNR_db , squeeze(H2_def(1,2,:))./squeeze(H2_def(1,3,:)) , 'o-');
plot(SNR_db , squeeze(H2_def(2,2,:))./squeeze(H2_def(2,3,:)) , 'b--o');
plot(SNR_db , squeeze(H2_def(3,2,:))./squeeze(H2_def(3,3,:)) , 'c-*' );
plot(SNR_db , squeeze(H2_def(4,2,:))./squeeze(H2_def(4,3,:)) , 'g--');
plot(SNR_db , squeeze(H2_def(5,2,:))./squeeze(H2_def(5,3,:)) , 'r-.');
legend("MI=1","MI=5","MI=10","MI=30","MI=100")
grid on;


% --------- plotting --------------
figure(3)
subplot(1,2,1);
title("BER for H3");
hold on;
plot(SNR_db , squeeze(H3_def(1,1,:))./ squeeze(H3_def(1,3,:)*12) , 'o-');
plot(SNR_db , squeeze(H3_def(2,1,:))./ squeeze(H3_def(2,3,:)*12) , 'b--o');
plot(SNR_db , squeeze(H3_def(3,1,:))./ squeeze(H3_def(3,3,:)*12) , 'c-*' );
plot(SNR_db , squeeze(H3_def(4,1,:))./ squeeze(H3_def(4,3,:)*12) , 'g--');
plot(SNR_db , squeeze(H3_def(5,1,:))./ squeeze(H3_def(5,3,:)*12) , 'r-.');
legend("MI=1","MI=5","MI=10","MI=30","MI=100")
grid on;

subplot(1,2,2);
title("FER for H3");
hold on;
plot(SNR_db , squeeze(H3_def(1,2,:))./squeeze(H3_def(1,3,:)) , 'o-');
plot(SNR_db , squeeze(H3_def(2,2,:))./squeeze(H3_def(2,3,:)) , 'b--o');
plot(SNR_db , squeeze(H3_def(3,2,:))./squeeze(H3_def(3,3,:)) , 'c-*' );
plot(SNR_db , squeeze(H3_def(4,2,:))./squeeze(H3_def(4,3,:)) , 'g--');
plot(SNR_db , squeeze(H3_def(5,2,:))./squeeze(H3_def(5,3,:)) , 'r-.');
legend("MI=1","MI=5","MI=10","MI=30","MI=100")
grid on;


