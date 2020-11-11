function res = hamming(len)
% Return a vector of length len giving the hamming window profile. Multiply element-by-element with your data vector. 
res = 0.54 + 0.46 * cos(2 * pi() * ((-len/2):1:(len/2 - 1)) / len);