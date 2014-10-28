% S3_Exercises
% Example section 3.5

% Get bandstop filter coefficients
B2 = fir1(100,[0.3, 0.5],'stop'); % band-stop filter with cutoff 
                                  % frequencies of 30% and 50% fs/2
x = zeros(1,1000);
x(50) = 1; % have created an impulse here 
Y2 = fft(conv(x,B2));
half = 1:ceil(length(Y2)/2);
plot(half/max(half),abs(Y2(half)));

% Make bandpass filter with cutoff frequencies of 0.2, 0.4, 0.6, 0.8
B5 = fir1(100,[0.2, 0.4, 0.6, 0.8],'bandpass');
Y5 = fft(conv(x,B5));
figure(2)
plot(half/max(half),abs(Y5(half)));

%% Examples section 3.7

% IIR
N = 20;
y(1) = 10;
a = 0.5; % or other examples in text
for n = 2:N+1
    y(n) = a*y(n-1);
end
figure(3)
plot(y)

%% Examples section 3.8.2

x = [0 0 1 5 1 -2 -3 -2 0 0];
y = [1 5 1 -2 -3 -2 0 0 0 0];

conv(x,y(length(y):-1:1))/length(x);
ind = max(abs(ans))

%% Example 3.83
x = [0 0 1 5 1 -2 -3 -2 0 0];
y = [1 5 1 -2 -3 -2 0 0 0 0];
yn = -y;

% compute cross-covariance
% only works when signal average is zero ****
conv(x,yn(length(yn):-1:1))/length(x)

%% Example 3.84

% signals
x = [0 0 1 5 1 -2 -3 -2 0 0];
y = [0 0 1 5 1 -2 -3 -2 0 0];
yn = -y;
y10 = y*10;
yimpulse = [0 0 0 0 400 0 0 0 0 0];
   

N = length(x);
sxx = x*x.' - sum(x)*sum(x)/N; % autocovariance x
syy = y*y.' - sum(y)*sum(y)/N; % autocovariance y
sxy = x*y.' - sum(x)*sum(y)/N; % cross-covariance x,y

% calculate cross-correlation value, rho
rho = sxy / sqrt(sxx*syy);

%% 3.8.6

x = [ 0 0 1 5 1 -2 -3 -2 0 0 ];
y = [ 1 5 1 -2 -3 -2 0 0 0 0 ];
yn = -y;
y10 = y*10;
yimpulse = [ 0 0 0 0 400 0 0 0 0 0 ];

figure(4)
plot(x); hold on; plot(y,'r')

% find offset
[val, x_peak] = max(x);
[val, y_peak] = max(y);

k = abs(x_peak-y_peak);

% **** k = 2 BUT BELOW THEY USE 8 IN THE BOOK -- WHY?
rho = correlate(8, x, y)

%% Exercise 3.8

% Make up a test signal of at least 10 values, and filter it with the 
% following coefficients. 

clear all

% a) b[k] = {0.4,0.4}
f = 100;
t = 0:0.0001:(2*1/f);
x = 3*cos(2*pi*100*t + (pi/6)) + rand(1,length(t));

% plot original
figure(8); clf; hold on
plot(x,'b')

% convolve with filter coefficients
ba = [0.4, 0.4];
Ya = conv(x,ba,'same');

% plot filtered
plot(Ya,'r')

% b) b[k] = {0.2, 0.2, 0.2, 0.2};
bb = [0.2, 0.2, 0.2, 0.2];
Yb = conv(x,bb,'same');
plot(Yb,'k')

% c) b[k] = {.1,.1, .1, .1, .1, .1, .1, .1}
bc = [.1,.1, .1, .1, .1, .1, .1, .1];
Yc = conv(x,bc,'same');
plot(Yc,'g')

% Best at smoothing = c -- averages over more values
% Can look at frequency response of filter
x = zeros(1,1000);
x(50) = 1; % have created an impulse here 
Y2 = fft(conv(x,bc));
half = 1:ceil(length(Y2)/2);
figure(2)
plot(half/max(half),abs(Y2(half)));


%% Exercise 3.14

% For the input signal x[n] = [6, 1, 3, 5, 1, 4], and a filter 
% h[k] = [-1, 3, -1], find y[n], where y[n] = h[k] * x[n].
% Use the algorithm given in Section 3.2.

clear all

x = [6, 1, 3, 5, 1, 4];
h = [-1, 3, -1];

y = conv(x,h);
