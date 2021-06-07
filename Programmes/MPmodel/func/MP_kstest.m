function [rejH0,D] = MP_kstest(cdf1, size1, cdf2, size2, alpha)
%MP_KSTEST Kolmogorov-Smirnov test with CDFs as inputs
% 
% INPUTS :
% cdf1 : double array, first CDF to test
% size1 : int, size of the first sample
% cdf2 : double array, second CDF to test
% size2 : int, size of the second sample
% alpha : double, test's significance level
%
% OUTPUTS :
% rejH0 : logical, 0 if H0 can't be rejected, 1 if we can reject H0 with a
% significance level alpha
% D : double, test statistic
%

deltaCDF = abs(cdf1-cdf2);
D = max(deltaCDF);

c = sqrt(-log(alpha/2)/2);
test = sqrt((size1+size2)/(size1*size2));

rejH0 = D>c*test;
end

