%% initial defult

clear;  
n = -1000:1:1000;   
k = 1:1:2001; 
start1 = -100;
end1 = 100;
N = length(n);

%% Q1 Create the Signal 'a(n)'

an = AllFunctions.calculate_window_function(start1,end1, n) ;

%% Q2 -Calculate the coefficient vector of the function an

ak = AllFunctions.calculate_coefficient_vector(an, n , k);

%% Calculate the an from the coefficevt vector

 an_chack = AllFunctions.calculate_an_vector(ak, n , k);

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
%% plotting  Q3
%We can see that the imaginary part is negligible and smaller than Epsilon

plot(k,real(ak)); %Some plots for fun 
 hold on
plot(k,imag(ak));
hold off

%% Q4 

ak_numeric = AllFunctions.calculate_coefficient_vector(an, n , k);
exp_by_calculation = exp(-1i*(k.')*2*pi*(n+1-1000)/2001);  % The full calculation is in the attachment
ak_analytic = (1/2001)*exp_by_calculation*an.';
hold on % plotting
plot(k, real(ak_analytic),'b');
plot(k, real(ak_numeric),'r');
hold off

%% Q5

bk = AllFunctions.exp_multiplication(ak,150,k);
hold on
bn = AllFunctions.calculate_an_vector(bk, n , k);
plot(n,real(an)) % compare to the original signal
hold off

%% Q6 

k_to_mul = (1-exp(-1j*k.'*(2*pi)/N));
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
plot(k,real(ek),k,real(fk)); 

%% Q10

gn = an.*cos(2*pi*500*(n/N));
gk = AllFunctions.calculate_coefficient_vector(gn,n , k);
plot(k,real(gk),k,real(ak));

%% Q11   

k_after_padding = 1:1:10005;
n_after_padding = -5002:1:5002;
num_of_zeros = 5;
fk = (1/5)*AllFunctions.zero_padding(ak, num_of_zeros).';
% k_after_padding = 1:1:length(fk);
% n_after_padding = circshift(n,[0,-round((length(n)/2)*num_of_zeros)]); %
% doesnt works yet 
hold on
fn = AllFunctions.calculate_an_vector(fk, n_after_padding, k_after_padding);
plot(n,real(an));
hold off

%% Q12

iterations  = 700;
number_of_examples = 5;
delta_sample = 50;
examples = [iterations:delta_sample :iterations+delta_sample * number_of_examples];
index = 0;
while index < number_of_examples
    hold on
     index = index +1; 
    AllFunctions.fourier_approx(ak, n , examples(index));
end
hold off
%%
