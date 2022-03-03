clear, clc, close all
% Part 1: Snore Signal Analysis
 
% Get a section of the sound file
data = load ('ucddb028_recm');    % Load file
inp = data.val;
fs = 128;                           % Sampling Frequency                    
N = length(inp);                    % Signal length
t = [0:1:N-1]/fs;                   % Time vector
 
% Plot the signal waveform
figure()
plot(t, inp)
xlim([0 max(t)])
ylim([-1.1*max(abs(inp)) 1.1*max(abs(inp))])
grid on
xlabel('Time, s')
ylabel('Amplitude')
title('Snore signal in the time domain')
 
% Plot the signal spectrogram
figure()
spectrogram(inp, 1024, 3/4*1024, [], fs, 'yaxis')
xlabel('Time, s')
ylabel('Frequency, Hz')
title('Spectrogram of the signal')
h = colorbar;
set(h, 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel(h, 'Magnitude, dB')
 
% Spectral analysis
w = hanning(N, 'periodic');
[X, f] = periodogram(inp, w, N, fs, 'power');
X = 20*log10(sqrt(X)*sqrt(2));
 
% Plot the signal spectrum
figure()
semilogx(f, X, 'r')
xlim([0 max(f)])
grid on
title('Amplitude spectrum of the signal')
xlabel('Frequency, Hz')
ylabel('Magnitude, dB')
 
% Plot the signal histogram
figure()
histogram(inp)
xlim([-1.1*max(abs(inp)) 1.1*max(abs(inp))])
grid on
xlabel('Signal amplitude')
ylabel('Number of samples')
title('Probability distribution of the signal')
 
%Part 2: Normalized Least Mean Square Filter
 
%Channel system order 
sysorder = 5 ; 
 
% Number of system points  
inp1 = inp'; 
n = randn(N,1);             %Create random noise
[b,a] = butter(5,0.25); 
Gz = tf(b,a,-1);            %Construct transfer function or convert to transfer function.
 
%The first sysorder weight value 
h = ldiv(b,a,sysorder)';  %This function is submitted to make inverse Z-transform (Matlab central file exchange) 
                        %The extra function ldiv is used
y = lsim(Gz,inp);       %Simulate time response of dynamic systems to arbitrary inputs.
 
%Add some noise 
d = y + n; 
totallength=size(d,1); 
 
%Take 60 points for training adaptive algorithm
N=60 ;   
 
%Begin of algorithm 
w = zeros (sysorder , 1 ) ; 
for n = sysorder : N  
    u = inp1(n:-1:n-sysorder+1) ; 
    y(n)= w' * u; 
    e(n) = d(n) - y(n) ; 
% Start with big mu for speeding the convergence then slow down to reach the correct weights 
      if n < 20 
          mu=2; 
      else 
          mu=0.5; 
      end 
% Use adaptive step to reach the solution faster mu = 0.95 * 2/M*r(0) 
  mu=0.95*2/(5*(0.001+var(u))); 
% alpha = 0.01      
    w = w + mu/(0.001+u'*u ) * u * e(n) ; 
end  
% Check of results 
for n =  N+1 : totallength 
    u = inp1(n:-1:n-sysorder+1) ; 
    y(n) = w' * u ; 
    e(n) = d(n) - y(n) ; 
end  
hold on 
subplot(2,1,1)
plot(t,d) 
title('System output') ; 
xlabel('Samples') 
ylabel('Noise output') 
subplot(2,1,2)
plot(t,-y,'r'); 
title('System output') ; 
xlabel('Samples') 
ylabel('True and estimated output') 
 
figure 
semilogy((abs(e))) ; 
title('Error curve') ; 
xlabel('Samples') 
ylabel('Error value') 
 
figure 
plot(h, 'k+') 
hold on 
plot(w, 'r*') 
legend('Actual weights','Estimated weights') 
title('Comparison of the actual weights and the estimated weights') ; 
axis([0 6 -0.5 0.5]);
 
%Part 3: Extract Formants
%Window the speech segment using a Hamming window.
x1 = inp.*hamming(length(inp))';
 
%Obtain the linear prediction coefficients. 
%To specify the model order, use the general rule that the order is two
% times the expected number of formants plus 2. 
% In the frequency range, [0,|Fs|/2], you expect three formants. 
% Therefore, set the model order equal to 8. 
% Find the roots of the prediction polynomial returned by lpc.
A = lpc(x1,14);         %Find 14th order linear predictor coefficients
rts = roots(A);         %Take square root
 
%Because the LPC coefficients are real-valued, the roots occur in complex conjugate pairs. 
%Retain only the roots with one sign for the imaginary part 
% and determine the angles corresponding to the roots.
%Because the LPC coefficients are real-valued, the roots occur in complex conjugate pairs.
% Retain only the roots with one sign for the imaginary part and determine 
% the angles corresponding to the roots.
rts = rts(imag(rts)>=0);
angz = atan2(imag(rts),real(rts));
 
%Convert the angular frequencies in rad/sample represented by the angles
% to hertz and calculate the bandwidths of the formants. 
% The bandwidths of the formants are represented by the distance of the 
% prediction polynomial zeros from the unit circle.
[frqs,indices] = sort(angz.*(fs/(2*pi)));
bw = -1/2*(fs/(2*pi))*log(abs(rts(indices)));
 
%Use the criterion that formant frequencies should be greater than 0 Hz 
% with bandwidths less than 400 Hz to determine the formants.
nn = 1;
for kk = 1:length(frqs)
    if (frqs(kk) > 0 && bw(kk) <400)
        formants(nn) = frqs(kk);
        nn = nn+1;
    end
end
formants
 


