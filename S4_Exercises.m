% S4 Exercises
% Ch 6 - Fourier Transform

%6.3 Plotting a spectrum

n = 0:99; % number of points
fs = 200; % sampling frequency
Ts = 1/fs; % sampling period

% x is our example signal
x = cos(2*pi*20*n*Ts + pi/4) + ...
    3*cos(2*pi*40*n*Ts - 2*pi/5) + ...
    2*cos(2*pi*60*n*Ts + pi/8);

% calculate X, being the fft of x
X = fft(x);
m = 0:length(X)-1;

% Print frequency resolution for interest
disp(sprintf('Freq resolution is every %5.2f Hz',...
    fs/length(X)));


% Plot frequency magnitude response and phase angles
figure(1)
plot(x)

figure(2)
% Plot magnitudes
subplot(2,1,1);
% stem(m*fs/length(y),abs(y),'b')
half_m = 0:ceil(length(X)/2);
stem(half_m*fs/length(X),abs(X(half_m+1)), 'b');
ylabel('magnitude');
xlabel('frequency(Hz)');
title('Frequency magnitude response');

% Plot phase angles
subplot(2,1,2);
stem(m*fs/length(X),angle(X), 'b');
ylabel('phase angle');
xlabel('frequency(Hz)');
title('Phase angle plot');


% set tolerance to set zero phase anglex
tolerance = 0.00001;
X2 = ceil(abs(X) - tolerance);
X3 = round(X2 ./ (X2+1)); % X3 is a vector of 0s and 1s: are they > tolerance?

% plot phase angles that matter

subplot(2,1,2);
stem(half_m*fs/length(X),angle(X(half_m+1)).*X3(half_m+1), 'b');

%%

% 6.9

%
% Show how the DFT function can represent a triangle wave

clear x

% Make a saw-tooth wave
for i=1:40
    x(i) = i;
end
for i=41:80
    x(i) = i-40;
end
figure(691)
plot(x)


% Find the DFT of it, and scale the result.
% In other words, represent signal x as a sum of sinusoids
[y] = fft(x);

% Scale the result, so that it is the same size as original
y= y/ length(x);
% Convert to polar coordinates (magnitudes and phases)
mag = abs(y);
phi = angle(y);

% plot frequency-magnitude response
figure(692)
stem(abs(y(1:length(y)/2)))

% Now, reconstruct it.
% This shows what the "sum of sinusoids" version looks like
t=0:(1/length(mag)):1;
f = 1; % Our fundamental frequency is 1, since time t=n*m/N
a = 0;
% Show it, adding another sinusoid each time
figure(693)
for k=0:length(mag)-1
    a=a+mag(k+1)*(cos(2*pi*f*k*t+phi(k+1)) ...
        + j*sin(2*pi*f*k*t+phi(k+1)));
    plot(1:length(x), x, 'r', 1:length(a), real(a), 'b-.')
    pause(0.1);
end


%% 6.11 Exercise

clear x

% create impulse function x
x = zeros(200,1);
x(50) = 1;

figure(6111)
plot(x)

% Frequency response is found by performing DFT on a system's output when
% impulse function is given as input.
for n = 2:length(x)
    w(n) = x(n) + x(n-1);
    y(n) = x(n) - x(n-1);
end

% take DFT
W = fft(w);
Y = fft(y);

% plot freq-mag response
figure(6112)
subplot(211)
stem(abs(W(1:length(W)/2)))
subplot(212)
stem(abs(Y(1:length(Y)/2)))

% So, w is a lowpass filter (and not a very good one)
% y is a ... what?


%% Ch 6 Project: Given a signal, find an approximation to it using the fewest sinusoids possible, with an error of less than 10%.

 clear all 



% % Make signal x, saw-tooth wave
% for i=1:40
%     x(i) = i;
% end
% for i=41:80
%     x(i) = 40; %i-40;
% end
% % repeat
% for i=81:120
%     x(i) = i-80;
% end
% for i=121:160
%     x(i) = 40; %40;
% end

filename = 'Kolohe_Nov102012_fasted_none';
load(filename)

subx = O2(ind(1):10:ind(1)+1500)';
x = subx;

figure(61)
plot(x)

% Perform FFT
X = fft(x);
X = X/length(X); % scale result
half_m = 0:ceil(length(X)/2);
fs = 200; 

% find magnitudes and phases
mag = abs(X);
phi = angle(X);

% plot frequency-magnitude response
figure(62)
subplot(211)
stem(half_m*fs/length(X),abs(X(half_m+1)));
xlabel('Frequency, Hz'); ylabel('Magnitude')
subplot(212)
stem(half_m*fs/length(X),angle(X(half_m+1)));
xlabel('Frequency, Hz'); ylabel('Phase Angle');

% approximate
f = 1;    % fundamental frequency
a = 0;    % initial approximation
ten_perc = 0.001*sum(abs(x));   % goal is to have error be within 10% 
err = sum(abs(x - a));      % error is difference between original signal and approximate
count = 0;  % number of sinusoids we have used

step = 1/length(mag);
t = 0:step:1-step;


while ((err > ten_perc) && count < length(mag)) % while our error is still greater than we want, and while we still have sinusoids
    % find next biggest magnitude sinusoid to add
    [v,p] = max(mag);
    % add this sinusoid to our approximation
a=a+mag(p)*(cos(2*pi*f*(p-1)*t+phi(p)) ...
     + j*sin(2*pi*f*(p-1)*t+phi(p)));
    
    % mark magnitude as having already been incorporated
    mag(p) = 0;
    
    % update variables
    err = sum(abs(x-a)); % recompute error
    count = count+1; % tally number of sinusoids used
    
    figure(63)
     plot(1:length(x), x, 'r', 1:length(a), real(a), 'b-.')
     pause(0.1);
end


