classdef AllFunctions
     methods(Static)
            function an = calculate_window_function(start_point, end_point, time_vector) % with convelotion
                N = length(time_vector);
                middle = round((start_point+end_point)/2,1);
                a = ones(1,(abs(start_point)+abs(end_point)));
                b=zeros(1,N);
                % now we need to chack all the options.
                if (middle==0)   %Chack if winodw is symetric
                    the_value_of_b = round(N/2, 0 );
                    b(the_value_of_b) = 1;
                else
                    if (middle>0)
                        the_value_of_b = round((N/2),0 ) - middle;
                        b(the_value_of_b) = 1;
                    else  %the middle is not zero and not greater then 0 - so must be less then zero
                        the_value_of_b = round((N/2),0 ) + middle;
                        b(the_value_of_b) = 1;
                    end
                end
                an=conv(b,a,'same');
                plot(time_vector,real(an))
            end
            
            function ak = calculate_coefficient_vector(an, time_vector , freq_vector) %Calculate the coefficient vector of the function an
                N = length(time_vector);
                w = 2*pi/N;
                exponent = exp(-1j*w*freq_vector.'*time_vector);
                ak = (exponent*an.').*(1/N);
                plot(freq_vector,real(ak))
            end
            
            function an = calculate_an_vector(ak, time_vector , freq_vector) %Calculate the coefficient vector of the function an
                N = length(time_vector);
                w = 2*pi/N;
                exponent = exp(1j*w*freq_vector.'*time_vector);
                an = (exponent.'*ak).';
                plot(time_vector,real(an));
            end
            
            
            function  symetric_bool = does_the_vector_symetric(vector)
                symetric_bool = sum(vector - flip(vector))<eps;
            end
            
            function  real_bool = does_the_vector_real(vector)
                eps_vector = eps(1:1:length(vector));
                bool_vector = sum(imag(vector) - eps_vector)<0 ;
                real_bool = all(bool_vector(:)  == 1);
            end
            
            function bk = exp_multiplication(ak,value_to_mul,freq_vector)
                bk = ak.*exp(1i*(-value_to_mul)*2*pi*(freq_vector.')/2001);
            end
            % bn = calculate_an_vector(bk, n , freq_vector);
            
            function ck = k_multiplication(ak,k_to_mul)
                ck  =(ak).*k_to_mul;
            end
            
            function dk = speical_mul(ak,k)
                N = length(k);
                dk=N*(ak).^2;
            end
            
            function was_parseval_right = parseval(frequency_domain ,time_domain, time_vector)
                N1 = length(time_vector);
                N = 2001;
                parseval_f_domain =sum((1/N)*(abs(time_domain).^2));
                parseval_t_domain =sum((abs(frequency_domain)).^2);
                max_possible_error = (length(time_domain))*(length(frequency_domain)*eps);
                the_error_itself =abs( sum((parseval_f_domain) - ( parseval_t_domain)));
                was_parseval_right = the_error_itself < max_possible_error;
            end
            
            function fk = zero_padding(original_array, padding_length)
                fk = zeros(1,padding_length*length(original_array));
                fk([1:padding_length:end])=(original_array.');
            end
            
            function aM = fourier_approx(ak, time_vector , iterations)
                N = length(time_vector);
                w = 2*pi/N;
                for index = 1:N
                    k_internal = (1:1:iterations*2+1);
                    aM(index) =(1/N) * sum(ak(k_internal).'.*exp(1j*w*(k_internal)*(time_vector(index))));
                end
                plot(time_vector,real(aM))
            end
            
            
     end
end
