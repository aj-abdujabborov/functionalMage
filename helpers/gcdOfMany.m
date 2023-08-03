function [y] = gcdOfMany(x)
%GCDOFMANY Find greatest common divisor of a vector of integers
% GCD = GCDOFMANY(x), where x is a vector of integers.

assert(isequal(x, floor(x)), 'Values must all be integers.');
if length(x) == 1, y = x; return; end
assert(~isempty(x), 'Input is empty.');
y = gcdCaller(x);

function [y] = gcdCaller(x)
    xL = length(x);
    if xL > 2 % recursion
        y = gcd(x(1), gcdCaller(x(2:end)));
    elseif xL == 2 % base case
        y = gcd(x(1), x(2));
    end

