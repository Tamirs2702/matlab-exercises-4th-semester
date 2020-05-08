%Ofek Aharony                          % Amit Shmueli
%% Default values - by homework file
% All values are changeable  
% inverse_coefficients_vector
clear ;  close all; clc;
n = -1000:1:1000;   % time domain - by defult
k = -1000:1:1000;  %Frequency domain - by defult
N = length(k);
w=2*pi./N;
input_padding = 5;
step_start = -100;
step_end = 100;
amplitude = 0.2;
set(0,'DefaultFigureWindowStyle','docked')


%% Q1 Create the Signal 'a(n)'
figure(1)
an = AllFunctions.calculate_window_function(step_start,step_end, n) ;  % Using convolution
figure(1)
title('Window Function');
xlabel('time');
ylabel('amplitude');

%% Q2 -Calculate the coefficient vector of the function an
figure(2)
ak = AllFunctions.coefficients_vector(an, n , k); 
title('Fourier Coefficients - real and image');
xlabel('k');
ylabel('amplitude');
legend('a[k]');

%% Calculate the an from the coefficevt vector  --> optional
% figure(2)
% an_chack = AllFunctions.inverse_coefficients_vector(ak, n , k);
% title('get a[n] from ak');
         
 %% Q3 - Because 'an' is real and symmetrical - 'ak' is real and symmetrical too
%Checks whether the signal an is symmetrical
symetric_bool = AllFunctions.does_the_vector_symetric(an) ;
real_bool = AllFunctions.does_the_vector_real(an);
if (real_bool == 1 && symetric_bool == 1) 
    disp('a[n] is real and symetric' )
end

%Checks whether the signal ak is symmetrical
symetric_bool = AllFunctions.does_the_vector_symetric(ak) ;
real_bool = AllFunctions.does_the_vector_real(ak);
if (real_bool == 1 && symetric_bool == 1) 
    disp('a[k] is real and symetric')
end

%% plotting of  Q3 - optional
%We can see that the imaginary part is negligible and smaller than Epsilon
figure(3)
hold on
title('Fourier Coefficients ak');
grid on
grid minor
plot(k,real(ak)); %plotting of real part of ak.
plot(k,imag(ak)); %plotting of image part of ak.
legend('Real(ak)','Imag(ak)');
hold off 

%% Q4 
figure(4)
ak_numeric = AllFunctions.coefficients_vector(an, n , k);
ak_analytic=(sin((k*pi/N)*199))./(N*sin(k*pi/N));% analytic calculate of ak using the formula. 

% plotting
plot(n, real(ak_analytic));
title('Fourier Coefficients - numeric and analitic');
hold on
xlabel('k');
ylabel('amplitude');
grid on;grid minor;
plot(k, real(ak_numeric));
legend('Analytic ak','Numeric ak');
hold off
legend('a[k] - numeric' , ' a[k] - analitic');


%% Q5
figure(5)
bk = AllFunctions.exp_multiplication(ak,150,k);
bn = AllFunctions.inverse_coefficients_vector(bk, n , k);
hold on;
grid on;grid minor;
xlabel('time');
ylabel('amplitude');
title(' b[n] compare to a[n] ');
plot(k,real(an)); %plotting of an..
legend('b[n] - shift of 150 unit time of a[n]' , ' a[n]');
hold off
title('Time Shift');




%% Q6 
figure(6)
k_to_mul = (1-exp(-1i*k*w));
ck = AllFunctions.k_multiplication(ak,k_to_mul) ;
cn = AllFunctions.inverse_coefficients_vector(ck, n , k);
grid on;grid minor;
title('Window Function Derivative');
xlabel('time');
ylabel('amplitude');
legend('C[n]');
%% Q7
figure(7)
dk = AllFunctions.speical_multiplicationl(ak,k);
hold on 
title('Convolution in the space time');
xlabel('time');
ylabel('amplitude');
dn = AllFunctions.inverse_coefficients_vector(dk, n ,k);
legend('D_n');
hold off
%% Q8

parseval_theorem = AllFunctions.parseval(dk ,dn,n);
if(parseval_theorem == 1)
        disp('parseval  was right')
end

%% Q9
figure(9)
en   = (bn.*an);
k_fix = 1:1:2001;  %% The convolution has no relative time. In order to display the functions on each other we needed to be fixed
ek = AllFunctions.coefficients_vector(en, n , k_fix) ;
fk = cconv(ak,bk,N);   %Cyclic convolution
hold on
plot(n,real(fk)); 
title('Convolution in the frequency space');
xlabel('k');
ylabel('amplitude');
legend('f[k] - by cyclic convolution' , ' e[k] = by multiply in time');

%% Q10
figure(10)
gn = an.*cos(2*pi*500*(n/N));
gk = AllFunctions.coefficients_vector(gn,n , k);

hold on
grid on;grid minor;
plot(k,real(ak));
title('The effect of multiply by cos');
xlabel('k');
ylabel('amplitude');
legend('g[k] - a[k] multiply by cos' , ' a[k]');
hold off



%% Q11   
figure(11)
m = (input_padding * n(1)):1: (input_padding * n(end)); %our new time vector.
fk = AllFunctions.zero_padding(ak, input_padding); % same as "fk=upsample(ak,input_padding);"
fk = amplitude * fk(1:end-4);
fn = AllFunctions.inverse_coefficients_vector(fk, m , k);%creating f[n] from the coefficients
hold on
grid on
plot(k,an);
title('Collapse in time');
xlabel('time');
ylabel('amplitude');
legend('f[k] - The shrinking signal' , 'a[k] - the original signal');
hold off


%% Q12
%calculate the function with different number of coefficients
figure(12)
for time_index=700:50:1000 % creating the approximations of the signal and plotting them.
    internal_time_vector = (-time_index:time_index);
    aM  = AllFunctions.inverse_coefficients_vector(ak, internal_time_vector , k);
    figure(12)
    plot(k,real(aM)); %plotting of am.
hold on
 end
grid on;grid minor;
title('Fourier approximation');
xlabel('time');
ylabel('amplitude');
legend('M = 700','M = 750','M = 800','M = 850','M = 900','M = 950','M = 1000')
hold off
%% Q13
figure(13);
    hnSin = (an.*sin(2*pi*500*n/2001));
    hkSin = imag(AllFunctions.coefficients_vector(hnSin,n,k-3)');
    
    hk_hilbert = AllFunctions.hilbert(gk);
    hn_hilbert = AllFunctions.inverse_coefficients_vector(hk_hilbert ,n, k);


plot (k, real(hkSin) , k, imag(hk_hilbert), k, real(ak));
title('The effect of multiply by sin');
grid on;
xlabel('K');
ylabel('amplitude');
legend('hk[k] - a[k] multiply by sin' ,'hk[k] - g[k] multiply by H[k]' , 'a[k]');

% The original signal is multiplied by the sine
figure(14);
plot (n, real(hnSin), n, real(hn_hilbert) );
title('Reverse Transformation of h[k]');
grid on;
xlabel('time');
ylabel('amplitude');
legend('hSin[n]' , 'hHilbert[n]');
