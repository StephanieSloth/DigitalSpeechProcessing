%% Pitch Detection
% Estimate the speech signal's pitch

% Speech file: [`MaoYiSheng.wav'; ' ...
% Time: 0.9s;
% Sampling rate: 8kHz;
% Quantization precision: 16bits/sample;
% Format: LSB, MSB;
% Frame length: 20 ms (160 samples);
% Total frames：45

clc;clear;
file = 'MaoYiSheng.wav';
[x,Fs] = audioread(file);
time = 0.9;
L = 160; % window length
FrameNum = length(x)/L;

%% SHORT-TIME AUTOCORRELATION FUNCTION
% Rn的行数是第几个窗，列数代表lag index k
for n = 1:FrameNum % analysis position or window shift(n^)
    for k = 1:L % lag
        Rn(n,k) = 0;
        
        for m = 1:L-k % each point
            Rn(n,k) = Rn(n,k)+x(m+(n-1)*L)*x(m+k+(n-1)*L);
        end
  
    end
end
% get the pitch

maxk = zeros(1,FrameNum);
pitch = zeros(1,FrameNum);
% next max magnitude is the pitch 
for n = 1:FrameNum
    maxk(n) = find(Rn(n,15:end) == max(Rn(n,15:end)));
end
maxk = maxk+15;
pitch = maxk/8;% 20ms = 160 points
%% hamming 
% Rn的行数是第几个窗，列数代表lag index k
for n = 1:FrameNum % analysis position or window shift(n^)
    w = hamming(L);
    for k = 1:L % lag
        Rn_ham(n,k) = 0;
        for m = 1:L-k % each point
            Rn_ham(n,k) = Rn_ham(n,k)+x(m+(n-1)*L)*w(m)*x(m+k+(n-1)*L)*w(m+k);
        end
    end
end
% get the pitch

maxk_ham = zeros(1,FrameNum);
pitch_ham = zeros(1,FrameNum);
% next max magnitude is the pitch 
for n = 1:FrameNum
    maxk_ham(n) = find(Rn_ham(n,15:end) == max(Rn_ham(n,15:end)));
end
maxk_ham = maxk_ham+15;
pitch_ham = maxk_ham/8;% 20ms = 160 points
%% THE MODIFIED SHORT-TIME AUTOCORRELATION FUNCTION
for n = 1:FrameNum-1 % analysis position or window shift(n^)
    for k = 1:L % lag
        RnM(n,k) = 0;
        for m = 1:L % all point in window(length = L)
            RnM(n,k) = RnM(n,k)+x(m+(n-1)*L)*x(m+k+(n-1)*L);
        end
  
    end
end
% get the pitch

maxkM = zeros(1,FrameNum);
pitchM = zeros(1,FrameNum);
% next max magnitude is the pitch 
for n = 1:FrameNum-1
    maxkM(n) = find(RnM(n,15:end) == max(RnM(n,15:end)));
end
maxkM = maxkM+15;
pitchM = maxkM/8;% 20ms = 160 points

%% draw
% time
t=linspace(0,0.9,7200);
figure(1);
stem(t,x,'.');
xlabel('Time(s)');title('Time Domain');

% Autocorrelation functions
figure(2);
subplot(221)
k = 1:L;
stem(k,Rn(5,:),'.');grid on;
xlabel('Lag k');ylabel('R(x)');title('Frame 5');
subplot(222)
stem(k,Rn(20,:),'.');grid on;
xlabel('Lag k');ylabel('R(x)');title('Frame 20');
subplot(223)
stem(k,Rn(35,:),'.');grid on;
xlabel('Lag k');ylabel('R(x)');title('Frame 30');
subplot(224)
stem(k,Rn(40,:),'.');grid on;
xlabel('Lag k');ylabel('R(x)');title('Frame 45');
sgtitle('STAF with Rectangular Window');

% Hamming - destory periodicity
figure(3);
subplot(221)
k = 1:L;
stem(k,Rn_ham(5,:),'.');grid on;
xlabel('Lag k');ylabel('R(x)');title('Frame 5');
subplot(222)
stem(k,Rn_ham(20,:),'.');grid on;
xlabel('Lag k');ylabel('R(x)');title('Frame 20');
subplot(223)
stem(k,Rn_ham(35,:),'.');grid on;
xlabel('Lag k');ylabel('R(x)');title('Frame 30');
subplot(224)
stem(k,Rn_ham(40,:),'.');grid on;
xlabel('Lag k');ylabel('R(x)');title('Frame 45');
sgtitle('STAF with Hamming Window');

% Modified
figure(4);
subplot(221)
k = 1:L;
stem(k,RnM(5,:),'.');grid on;
xlabel('Lag k');ylabel('R(x)');title('Frame 5');
subplot(222)
stem(k,RnM(20,:),'.');grid on;
xlabel('Lag k');ylabel('R(x)');title('Frame 20');
subplot(223)
stem(k,RnM(35,:),'.');grid on;
xlabel('Lag k');ylabel('R(x)');title('Frame 30');
subplot(224)
stem(k,RnM(40,:),'.');grid on;
xlabel('Lag k');ylabel('R(x)');title('Frame 45');
sgtitle('Modified STAF');

% pitch
figure(5);
subplot(211);
stem(pitch,'.');
xlabel('frame(n)');ylabel('Pitch(ms)');title('STAF with Rectangular Window ');
grid on;
subplot(212);
stem(pitchM,'.');
xlabel('frame(n)');ylabel('Pitch(ms)');title('Modified STAF');
grid on;
sgtitle('Pitch for each frame');