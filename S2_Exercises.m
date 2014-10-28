% S2 Exercises
% 4.1, 4.3

%% 4.1
% What are the minimum and maximum values for x(t)?
% x(t) = 3cos(2pi 100t + pi / 6)

f = 100;
t = 0:0.0001:(2*1/f);
x = 3*cos(2*pi*100*t + (pi/6));
figure(41)
plot(t,x)

% minimum and maximum values should be amplitude, +/-3. 
% check:
disp(min(x))
disp(max(x))

%% 4.3
% The analog signal, x(t), is:
  % x(t) = 3*cos(2*pi*2000*t + (pi/4)) + 2*cos(2*pi*5000*t) + cos(2*pi*11000*t - (pi/7)) 
 
% Plot each frequency component (in the time-domain) separately.
    % find fo
fo = gcd(2000,gcd(5000,11000));
    % set time
t = 0:0.000001:(2*1/fo);
    % separate sinusoid into components
x1 = 3*cos(2*pi*2000*t + (pi/4));
x2 = 2*cos(2*pi*5000*t);
x3 = cos(2*pi*11000*t - (pi/7)); 
    % entire sinusoid
x = 3*cos(2*pi*2000*t + (pi/4)) + 2*cos(2*pi*5000*t) + cos(2*pi*11000*t - (pi/7));

    % plot components and total signal
subplot(411)
plot(t,x1);
subplot(412)
plot(t,x2)
subplot(413)
plot(t,x3)


% Graph this function (in the time-domain).
subplot(414)
plot(t,x)
xlabel('Time (s)')
ylabel('Amplitude')

% Represent x(t) in terms of a fundamental frequency, amplitudes, and phases.
% fundamental frequency = gcd(2000,5000,11000) = 1000;
% amplitudes = 3,2,1;
% phase angles = pi/4, 0, pi/7;

x_new = 3*cos(2*pi*2*fo*t + (pi/4)) + 2*cos(2*pi*5*fo*t) + cos(2*pi*11*fo*t - (pi/7));
hold on
plot(t,x_new,'r')

%% 5.15
% Assume the signal x(t) = 3*cos(2*pi*2000*t + pi/4) + 2*cos(2*pi*5000*t) +
% cos(2*pi*11000*t - pi/7); is sampled at 10 kHz
% a. What does the x(n) equation become?
n = 1:100;
x(n) = 3*cos(2*pi*n*0.2+pi/4) + 2*cos(2*pi*n*0.5) + cos(2*pi*n*1.1 - pi/7);

% b. Plot x[n].
plot(n,x,'r')


%% Ch 5 Project

% Simulate sampling a signal made of frequencies: 200 Hz, 400 Hz, 1400 Hz, 
% and 1600 Hz at 1200 samples/second. Show that aliasing combines two 
% frequencies into one

freq = [200, 400, 1400, 1600];     % example frequencies
phase = [0, 0, 0, 0];              % example phases
mag = [1, 1, 1, 1];                % example magnitudes

fs = 1200; % sampling frequency
Ts = 1/fs; % sampling period

num_points = 200;   % How many points to use
num_samples = 11;   % How many samples to simulate reading

step = 1/(min(freq)*num_points);  % Calculate step size
t = 0:step:2*(1/min(freq));       % Create "time" index
n = 0:num_samples-1;              % Create sample index
                 
   for m=1:length(freq)
       x(m, :) = mag(m)*cos(2*pi*freq(m)*t + phase(m));
   end
   
plot(t,x); hold on

    for m=1 : length(freq)
        x2(m, n+1) = mag(m)*cos(2*pi*freq(m)*n*Ts + phase(m));
    end

figure(1)
 subplot(2,1,1);
   plot(t, x(1,:), 'b', t, x(3,:), 'g');
   subplot(2,1,2);
   plot(t, x(2,:), 'r', t, x(4,:), 'y');

      figure(2)
   subplot(2,1,1);
   plot(n,x2(1,:),'bx-', n, x2(3,:), 'go-');
   subplot(2,1,2);
   plot(n,x2(2,:),'rx-', n, x2(4,:), 'yo-');
