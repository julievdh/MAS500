% S5 Exercises
% 20 March 2014

% Example 7.1

% IIR filter with feedforward coefficients {4,5,6} and feedback
% coefficients {2,3}


B = [4, 5, 6];      % feedforward coefficients

% find zeros using quadratic formula
zero1 = (-B(2) + sqrt(B(2)^2 - 4*B(1)*B(3)))/(2*B(1));
zero2 = (-B(2) - sqrt(B(2)^2 - 4*B(1)*B(3)))/(2*B(1));

% find poles
A = [1,-2,-3];      % feedback coefficients
pole1 = (-A(2) + sqrt(A(2)^2 - 4*A(1)*A(3)))/(2*A(1));
pole2 = (-A(2) - sqrt(A(2)^2 - 4*A(1)*A(3)))/(2*A(1));

% visualize
zplane(A,B)


%% Exercises

% 7.1
% 	What is the z-transform of x = {3, 2, 4}?
% X(Z) = 3z^0 + 2z^-1 + 4z^-4

% 7.3 An FIR filter has the following filter coefficients: {1, ?1}.
% What is the frequency response? Where are the zeros located?
% Use MATLAB to calculate the output with an input of cos(2?100t),
% sampled at fs = 200 samples/second.


h = [1, -1];        % filter coefficients

% build impulse function
imp = zeros(1,1000);
imp(10) = 1;
% get frequency response
H = fft(conv(h, imp));
half = 1:ceil(length(H)/2);
plot(abs(H(half)))


% find zeros
% H(Z) = 1z^0 -1z^-1;
% Zero at z = 1;


% Calculate output with an input of cos(2pi100t) with fs = 200
n = 1:1000;
fs = 200;
x = cos(2*pi*100*n/200);
% Now put through filter
y = conv(x, h);
plot(n,y(1:1000))
