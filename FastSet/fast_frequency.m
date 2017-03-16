function [y,x] = fast_frequency(a)
%	Computes the number of occurencies of each element in a. 
%	Analogoues to [y,x] = hist(a,unique(a)); but supports any type of a (hist supports only doubles. also, consider the MatLab's bug in unique applied to int64 that may merge few close values. 
%
%	Receives:
%		a	- numeric vector of any type. Matrices ('rows' switch in union are not supported). Strings are not supported. 
%
%	Returns:
%		y  - vector  - of type uint64 on x64 platforms, uint32 - on 32-bit platforms. Each elements contains the number of the corresponding elements in a. 
%		x - vector  - (optional) same number of elements as in y, but the size is either same as y or rot90(y) depending on allignment of a. This is the same convention as in hist, unique and fast_unique
%
%	Notice: 
%		1. Function supports 1 or 2 returned arguments 
%		2. The function uses fast_unique and it's performance will be improved once fast_unique is optimized. 
%
%	Compatibility:
%			The method is faster than the built-in hist. The performance benefit grows with size. 
%			As stated above, fast_frequency handles 64-bit inputs correctly. 
%		Except for the following cases, FastSet commands are identical to the standard MatLab set commands.
%			-	input arrays must be sorted (in ascending order) and unique. Use MatLab's unique function to ensure that.
%			-	returned lists of indeces are integer formats (either 32 or 64 bits, natural for the OS)
%			-	'rows' instruction is not supported
%	To do:
%		Optimize fast_unique
%		Parallelize. (OpenMP?)
%
% Example:
%	
%	Ar = uint32(round(10000000*rand(1000000,1)));
%	tic; [y, x] = fast_frequency(Ar); toc
%
%   A = uint64(12345678901234567890);
%   B = uint64(12345678901234567891);
%	A~=B % - A & B are different
%	
%	[y x] = hist([A B B],[A B]) % fails
%	[y x] = hist(double([A B B]),[A B]) % fails
%	[y x] = hist(double([A B B]),double([A B])) % returns double that is different from A & B
%	[y x] = fast_frequency([A B B]) % returns correct values
%	
% Created By:
%	Dr. Lev Muchnik, LevMuchnik@gmail.com
%	Stern School Of Business, NYU
%	Jan 9, 2011