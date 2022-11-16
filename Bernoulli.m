function y = bernoulli( P , N )
% BERNOULLI.M
% This function generates n independent draws of a Bernoulli
% random variable with probability of success P.
% first, draw N uniform random variables
x = rand( N , 1 );
% set y=1 if x is less than p. This gives probability P of success
y = ( x  <= P );
% end function definition