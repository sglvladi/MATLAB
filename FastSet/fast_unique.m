function [c, m, n] = fast_unique(a)
%	Returns a list of ordered values appearing in both a. a does not need to be sorted. 
%	The method is analogous to MatLab's 'unique', however, it handles int64 and uint64 types correctly. MatLab seams to convert them (internally) to doubles that results in some distinct but close values being indistinguishable.
%	Performance: varies. c = fast_unique(a) is parallelized on Windows machines and achives 2-3 times the performance of the MatLab's code. 
%   Unparallelized version (Linux):  [c, ia, ib] = fast_unique(a) is slover than MatLab's by ~50%. 
%   Parallelied (Windows): [c, ia, ib] = fast_unique(a) is faster than MatLab's by ~50% for ints. Slower by 300% for int64 and doubles. 
% Performance varies by size. Small arrays - this implementation wins. Very-very large - MatLab wins. 
%   [c] = fast_unique(a): this implementation is typically 30-50% faster. 
%  
%
%	Receives:
%		a	- numeric vector of any type. Matrices ('rows' switch in union are not supported). Strings are not supported. 
%
%	Returns:
%		c	-	numeric vector of the same type as a. List of it's unique values. 
%		m  -	(optional)integer vector (int64 on x64 systems, int32 - on 32 bit systems). a(m) = c. Unlike MatLab's unique, the returned indeces m point to the first (not last!) occurence of the value in a.
%		n  -	(optional)integer vector (int64 on x64 systems, int32 - on 32 bit systems). c(n) = a
%	Notice: 
%		1. Function supports 1 or 3 returned arguments (not 2)
%
%	Compatibility:
%			Unlike MatLab's unique, the returned indeces m point to the first (not last!) occurence of the value in a.
%			As stated above, fast_unique handles 64-bit inputs properly 
%		Except for the following cases, FastSet commands are identical to the standard MatLab set commands.
%			-	input arrays must be sorted (in ascending order) and unique. Use MatLab's unique function to ensure that.
%			-	returned lists of indeces are integer formats (either 32 or 64 bits, natural for the OS)
%			-	'rows' instruction is not supported
%
%
% Example:
%	
% 	Ar = uint32(round(2e19*rand(100000,1)));
% 	tic; [x, m,n] = fast_unique(Ar); toc
%     tic; [x1,m1,n1] = unique(Ar); toc
% 
%     	Ar = uint64(round(300*rand(10000000,1)));
% 	tic; [x] = fast_unique(Ar); toc
%     tic; [x1] = unique(Ar); toc
% 
%   A = uint64(12345678901234567890);
%   B = uint64(12345678901234567891);
% 	A~=B % - A & B are different
% 	unique([A B B]) % yet, unique([A B]) returns single value (that differes from both A & B. 
% 	fast_unique([A B B]) % returns correct value
%
% Created By:
%	Dr. Lev Muchnik, LevMuchnik@gmail.com
%	Stern School Of Business, NYU
%	Jan 9, 2011