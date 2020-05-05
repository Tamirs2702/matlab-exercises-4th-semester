%% initial defult

clear;  
n = -1000:1:1000;   % time domain - by defult
k = -1000:1:1000;  %Frequency domain - by defult
N = length(k);
w=2*pi./N;
input_padding = 5;
step_start = -100;
step_end = 100;

%% Q1 Create the Signal 'a(n)'

an = AllFunctions.calculate_window_function(step_start,step_end, n) ;  % Using convolution

hold on 
title('a[n]')
hold off 


%% Q2 -Calculate the coefficient vector of the function an

ak = AllFunctions.calculate_coefficient_vector(an, n , k); 

hold on 
title('Fourier Coefficients ak');
hold off 

%% Calculate the an from the coefficevt vector

an_chack = AllFunctions.calculate_an_vector(ak, n , k);
hold on;
title('get a[n] from ak');
hold off;
         
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

hold on
title('Fourier Coefficients ak');
grid on
grid minor
plot(k,real(ak)); %plotting of real part of ak.
plot(k,imag(ak)); %plotting of image part of ak.
legend('Real(ak)','Imag(ak)');
hold off 

%% Q4 

ak_numeric = AllFunctions.calculate_coefficient_vector(an, n , k);
ak_analytic=(sin((k*pi/N)*199))./(N*sin(k*pi/N));% analytic calculate of ak using the formula. 

% plotting
grid on;grid minor;
plot(n, real(ak_analytic));
hold on
plot(n, real(ak_numeric));
legend('Analytic ak','Numeric ak');
hold off

%% Q5

bk = AllFunctions.exp_multiplication(ak,150,k);
bn = AllFunctions.calculate_an_vector(bk, n , k);

hold on;
title(' b[n] compare to a[n] ');
plot(k,real(an)); %plotting of an..
legend('a[n]','b[n]');
hold off


%% Q6 

k_to_mul = (1-exp(-1i*k*w));
ck = AllFunctions.k_multiplication(ak,k_to_mul) ;
cn = AllFunctions.calculate_an_vector(ck, n , k);

%% Q7

dk = AllFunctions.speical_mul(ak,k);
hold on 
dn = AllFunctions.calculate_an_vector(dk, n ,k);

%% Q8

parseval_theorem = AllFunctions.parseval(dk ,dn,n);
if(parseval_theorem == 1)
        disp('parseval  was right')
end

%% Q9
en   = (bn.*an);
ek = AllFunctions.calculate_coefficient_vector(en, n , k) ;
fk = cconv(ak,bk,length(ak));   %Cyclic convolution

hold on
plot(k,real(fk)); 
legend('ek','fk');
%% Q10
gn = an.*cos(2*pi*500*(n/N));
gk = AllFunctions.calculate_coefficient_vector(gn,n , k);

hold on
plot(k,real(ak));
legend('gk','ak');
hold off


%% Q11   
m = (input_padding * n(1)):1: (input_padding * n(end)); %time vector.
fk = AllFunctions.zero_padding(ak, input_padding); % same as fk=upsample(ak,5);
fk=fk(1:end-(input_padding-1)); % there is only one zero 
fn = AllFunctions.calculate_an_vector(fk, m , k);%creating f[n] from the coefficients
hold on
plot(k,ak);

%% Q12
%calculate the function with different number of coefficients
for time_index=700:50:1000 % creating the approximations of the signal and plotting them.
    internal_time_vector = (-time_index:time_index);
    aM  = AllFunctions.calculate_an_vector(ak, internal_time_vector , k);
    figure(12);
    plot(k,aM); %plotting of am.
hold on
 end
grid on
grid minor
xlabel('n')
ylabel('a_m [n]')
legend('M = 700','M = 750','M = 800','M = 850','M = 900','M = 950','M = 1000')
hold off
%%
