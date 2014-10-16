% MAS500 Session 1
% Ch 2.12 Exercises


%%
% 1. Create a signal x1, where x1 = 2cos(2?15t).
% Make t a range large enough to show 2 repetitions of x1.
% Plot x1.
% What are the minimum and maximum values of x1 from the graph?

% 2 repetitions = 2*f = 2* 15 = 30
t = 0:0.01:(2*1/15);
x1 = 2*cos(2*pi*15*t);
figure(1)
plot(t,x1)

% create a more smooth plot, by increasing resolution of t
t = 0:0.0001:(2*1/15);
x1 = 2*cos(2*pi*15*t);
hold on; plot(t,x1,'r')

minimum = min(x1)
maximum = max(x1)

%%
% 2. Let x2 = x1 + 2.5cos(2?30t), and adjust t as needed.

common = gcd(15,30)
t = 0:0.01:(2*1/common);
x2 = 2*cos(2*pi*15*t) + 2.5*cos(2*pi*30*t);
figure(2)
plot(t,x2)

minimum = min(x2)
maximum = max(x2)

%%
% 4. Find y1[n] = x1[n] ? .9x1[n ? 1], using x1 from Exercise 3.
% Plot y1 (adjusting t as needed).

for n = 2:length(x1)
    y1(n) = x1(n) - 0.9*x1(n-1);
end
figure(4)
plot(y1)

%%
% 5. Create y2[n] to be the average of the current and previous two values 
% of x2[n].
% Plot y2.

for n = 3:length(x2)
    y2(n) = mean(x2([n,n-1,n-2]));
end
figure(5)
plot(y2)

% OR

for n=3:length(x2)
    y2(n) = (x2(n) + x2(n-1) + x2(n-2))/3;
end

hold on;
plot(y2,'r')

%% 
% 6. Create y3[n] = x2[n] ? y3[n ? 1].
% Plot y3.

y3(1) = 0;
for n = 2:length(x2)
    y3(n) = x2(n) - y3(n-1);
end
figure(6)
plot(y3)

%% 
% 7. Make a random function that returns an array of m pseudorandom 
% integers between a low and high value specified. Use the rand function in 
% your function, and call your solution myrand. Notice that myrand should 
% have three parameters: the size of the array to return (m), the lowest 
% integer value it should return, and the highest integer value it should 
% return.

lo = 0; hi = 100; n = 10;
test = round(rand(1,n)*(hi-lo)+lo);

result = myrand(lo,hi,n);

%%
% 8. With an array of 8 values between 0 and 99 (use myrand), write a 
% function to return the values in sorted order, from lowest to highest.

result = myrand(0,99,8);
sorted = sort(result)

%% Other Exercises
% 5. Using MATLAB, plot the sinusoid 2cos(2?1500t + ?/2), starting at 
% t = 0 seconds. Be sure to use enough points to make a smooth graph, and 
% only showa fewrepetitions.

fs = 1500;
t = 0:0.00001:(2*1/fs*2);
x = 2*cos(2*pi*1500*t + pi/2);
figure(8)
plot(t,x)

%% Other Exercises
% 8. The following MATLAB code should return the magnitudes and phases 
% for an array of complex numbers. The line marked problem actually has 
% two problems with it.


% a. The first problem is the ambiguity discussed in Chapter 1. Explain 
% what this problem is, and demonstrate it via examples.

%  % This function should return the magnitude and phase.
%     %
%     function [mag, phi] = my_function(y)
%     % input Y is an array of complex numbers
% 
%     yr = real(y);
%     yi = imag(y);
%     for i=1:length(yr)
%         phi(i) = atan(yi(i)/yr(i)); % problem
%         mag(i) = sqrt(yr(i)*yr(i)+yi(i)*yi(i));
%     end

% dividing by zero (the real part of the number = 0)?

% b. For the second problem, consider what happens when it processes the 
% following code, which makes it crash! What is the problem?

% it didn't crash?
