n= -1000:1000;
k= -1000:1000; k2 = 1:2001; newk=-5002:5002;

an = WindowFunc(1000,100);
akNumeric = (FourierCoefficients(an,k,n,2001))';
akAnalitic = (1/2000)*(exp(1i*k*66*2*pi/667)-exp(-1i*k*200*pi/2001))./(1- exp(-1i*k*2*pi/2001));

%Check the imaginary part of the Fourier coefficients we got
    akimag = imag(akNumeric);
    
    bk = TimeShift(akNumeric, 150,k,2001);
    bn = ReverseTransformation (bk ,k, n , 2001);
    
    %Derived in time by multiplying by K
    ck = DerivedByMultiplication(akNumeric , k , 2001);
    cn = ReverseTransformation (ck ,k, n , 2001);
    
    %Time convolution by frequency multiplication
    dk = 2001*akNumeric.*akNumeric;
    dn = ReverseTransformation (dk ,k, n , 2001);
    Parseval(dk ,dn);
    
    %Frequency convolution by multiplication in time
    en = an.*bn;
    
    ekFourier = (FourierCoefficients(en,k2,n,2001))';
    ekConvolution = CyclicConvolution(akNumeric , bk, 2001)';
    
    fk = (Stretch(akNumeric,5,2001));
    fn = (ReverseTransformation (fk ,newk-3, n , 2001));
  
    gn = an.*cos(2*pi*500*n/2001);  
    gk = (FourierCoefficients(gn,k,n,2001))';
    
    %6 different approximations to the original signal
    am1 = FourierApprox(akNumeric ,700, n, 2001);
    am2 = FourierApprox(akNumeric ,750, n, 2001);
    am3 = FourierApprox(akNumeric ,800, n, 2001);
    am4 = FourierApprox(akNumeric ,850, n, 2001);
    am5 = FourierApprox(akNumeric ,900, n, 2001);
    am6 = FourierApprox(akNumeric ,950, n, 2001);
    
    hnSin = (an.*sin(2*pi*500*n/2001));
    hkSin = imag(FourierCoefficients(hnSin,k-3,n,2001)');
    
    hkHilbert = Hilbert(gk, 2001);
    hnHilbert = (ReverseTransformation (hkHilbert ,k, n , 2001));

% Window Function plot
figure(1);
plot (n, an);
title('Window Function');
xlabel('time');
ylabel('amplitude');

%Real Fourier coefficients and the imaginary Fourier coefficients
figure(2);
plot (k , akNumeric , k , akimag);
title('Fourier Coefficients - real and image');
xlabel('k');
ylabel('amplitude');
legend('a[k] - real' , ' a[k] - image');

%Fourier coefficients in numerical calculation and Fourier coefficients in analytical calculation plot
figure(3);
plot (k, akNumeric , k , akAnalitic );
title('Fourier Coefficients - numeric and analitic');
xlabel('k');
ylabel('amplitude');
legend('a[k] - numeric' , ' a[k] - analitic');

%Original signal and signal after shift in time plot
figure(4);
plot (n , bn , n , an);
title('Time Shift');
xlabel('time');
ylabel('amplitude');
legend('b[n] - shift of 150 unit time of a[n]' , ' a[n]');

%The derivative of the signal in time plot
figure(5);
plot (n , cn );
title('Window Function Derivative');
xlabel('time');
ylabel('amplitude');
legend('C[n]');

%The Convolution of time signal plot
figure(6);
plot (n , dn );
title('Convolution in the space time');
xlabel('time');
ylabel('amplitude');
legend('D[n]');

%Signal received by frequency convolution versus signal received by time
% multiplication plot
figure(7);
plot (k , ekConvolution , k , ekFourier);
title('Convolution in the frequency space');
xlabel('k');
ylabel('amplitude');
legend('e[k] - by cyclic convolution' , ' e[k] = by multiply in time');

%Original signal and signal multiplied by cosine plot
figure(8);
plot (k , gk , k ,akNumeric);
title('The effect of multiply by cos');
xlabel('k');
ylabel('amplitude');
legend('g[k] - a[k] multiply by cos' , ' a[k]');

%Original and Extension Signal plot
figure(9);
plot (n, fn , n , an);
title('Collapse in time');
xlabel('time');
ylabel('amplitude');
legend('f[k] - The shrinking signal' , 'a[k] - the original signal');

%The signals received by different approximations
figure(10);
plot (n , am1 , n , am2, n , am3, n , am4, n , am5, n , am6);
title('Fourier approximation');
xlabel('time');
ylabel('amplitude');
legend('a[m]1 - 700 iterations' , 'a[m]2 - 750 iterations' , 'a[m]3 - 800 iterations', 'a[m]4 - 850 iterations' , 'a[m]5 - 900 iterations' ,'a[m]6 - 950 iterations') ;

% Fourier coefficients of the original signal are multiplied by the sine
% plot
figure(11);

plot (k, hkSin , k, imag(hkHilbert), k, akNumeric);
title('The effect of multiply by sin');
xlabel('K');
ylabel('amplitude');
legend('hk[k] - a[k] multiply by sin' ,'hk[k] - g[k] multiply by H[k]' , 'a[k]');

% The original signal is multiplied by the sine
figure(12);
plot (n, hnSin, n, hnHilbert );
title('Reverse Transformation of h[k]');
xlabel('time');
ylabel('amplitude');
legend('hSin[n]' , 'hHilbert[n]');