function p = testMaRegress(m,sd,refVal,N,tail)
%TESTMAREGRESS test the results of the Major axis regression (t-test)
%   m  -  mean of the regression parameter (slope or intercept)
%   sd - standard deviation of the regression paramter
%   N  - number of points that went into the regression
%   ref 


if nargin <= 4
    tail = 'both';
end

v = N-1;
tval = (m-refVal) / sqrt(sd^2/N);       % Calculate T-Statistic
tdist2 = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution



switch tail
    case 'left'
        tdist1 = @(t,v) 1-(1-tdist2(t,v))/2; 
        p = abs((1-tval<0) - tdist1(tval,v));
    case 'right'
        tdist1 = @(t,v) 1-(1-tdist2(t,v))/2;              % 1-tailed t-distribution
        p = abs((1-tval>0) - tdist1(tval,v));
    case 'both'
        
        p = 1-tdist2(tval,v);
end

end

