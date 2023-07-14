function [bRound] = isRound(mat)
%ISROUND Returns true if input is round
%   B = isRound(A) returns true if all values of A are whole numbers.
%%
bRound = isequal(round(mat), mat);
