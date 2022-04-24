function S = fdec2base(D,B,N)
%FDEC2BASE Convert decimal integer to base B string.
%   FDEC2BASE(D,B) returns the representation of D as a string in
%   base B.  D must be a non-negative integer array smaller than 2^52
%   and B must be an integer between 2 and 36.
%
%   FDEC2BASE(D,B,N) produces a representation with at least N digits.
%
%   FDEC2BASE is designed to be a fast replacement for dec2base. It is
%   optimized for speed and conservation of memory. It has no error checks
%   so you have to make sure your input is valid.
% 
%   Examples
%       fdec2base(23,3) returns '212'
%       fdec2base(23,3,5) returns '00212'
%
D = D(:);
D = double(D);
B = double(B);
l = length(D);
maxd = max(D);
no = max(1,round(log2(maxd+1)/log2(B)));
while B^no <= maxd, no=no+1; end
if nargin == 3 && N>no, nmax = N;else nmax = no; end
S = zeros(l,no,'uint8');
% for large arrays a lot faster than matrix multiplication
while no > 1
    no   = no - 1;
    base = B^no;
    col  = nmax-no;
    S(:,col) = uint8(D/base-0.5);
    D  = D - double(S(:,col))*base;
end
S(:,nmax) = uint8(D);
if B>10
    % conserve memory + faster than s = symbols(s+1);
    symbols = uint8('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ');
    S=S+1;
    for a = 1:10000:size(S,1)
        ind = a:min(a+9999,size(S,1));
        S(ind,:) = symbols(S(ind,:));
    end
    S = char(S);
else
    S = S+uint8('0');
    S = char(S);
end
end

