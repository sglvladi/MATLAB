function [c, ia, ib] = fast_setxor_sorted(a,b)
%	Find set exclusive OR of two vectors. Returns a list of values that are not in the intersection of a and b.
%	Very efficient (usually, 3-10 times faster) replacement of the standard MatLab setxor function. 
%
%	Receives:
%		a	- numeric vector of any type. Matrices ('rows' switch in union are not supported).
%		b	- numeric vector of any type. Matrices ('rows' switch in union are not supported).
%	Notice, a and b must be of the same numeric type (both doubles, int16 , ...)
%
%	Returns:
%		c	-	numeric vector of the same type as a and b. intersection of the two sets. c is sorted and unique.
%		ai  -	(optional)integer vector (int64 on x64 systems, int32 - on 32 bit systems). Lists indeces of elements in c which came from a. 
%		bi  -	(optional)integer vector (int64 on x64 systems, int32 - on 32 bit systems). Lists indeces of elements in c which came from b. 
%	Notice: 
%		1. Function supports 1 or 3 returned arguments (not 2)
%
%	Compatibility:
%		Except for the following cases, FastSet commands are identical to the standard MatLab set commands.
%			-	input arrays must be sorted (in ascending order) and unique. Use MatLab's unique function to ensure that.
%			-	returned lists of indeces are integer formats (either 32 or 64 bits, natural for the OS)
%			-	'rows' instruction is not supported
%
% Performance:
%	Due to the following reasons, the function is considerably faster than the standard union method :
%	  1. It assums that a and b are sorted sets of unique values
%	  2. It's implemented in C++
%     3. ai and bi are integres (64bit on x64 platform and 32 bit on 32 bit platform), which saves RAM and avoids conversion of integers to doubles.
%
% Example:
%
%	b = [3 : 5 : 2000000]; 
%	a = [1: 3 : 2000000]; 
%	tic, 
%	[c, ai, bi]= setxor(a,b); 
%	t1 = toc; tic, 
%	[c1,ai1,bi1] = fast_setxor_sorted(a,b); 
%	t2 = toc; 
%	disp(sprintf('union: %2.2f sec. fast_union_sorted: %2.2f (%2.2f faster)',t1,t2,t1/t2));
%	if any([nnz(c~=c1) nnz(ai~=ai1) nnz(bi~=bi1)]) 
%		error('Results differ!!!!');
%	else 
%		disp('Verified');
%	end	
%
% Created By:
%	Dr. Lev Muchnik, LevMuchnik@gmail.com
%	Stern School Of Business, NYU
%	Oct. 26, 2008