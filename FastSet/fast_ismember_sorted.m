function [tf, loc] = fast_ismember_sorted(a,b)
%	Finds wherether a elements are members of the set b. Bothe a & must be sorted.
%	Very efficient (usually, x-x times faster) replacement of the standard MatLab ismember function. 
% 
%
%	Receives:
%		a	- numeric vector of any type. Matrices ('rows' switch in union are not supported).
%		b	- numeric vector of any type. Matrices ('rows' switch in union are not supported).
%	Notice, a and b must be of the same numeric type (both doubles, int16 , ...)
%
%	Returns:
%		tf	-	logical array. same size as a. tf(i) is false if a(i) is not a member of b. tf(i) is true if a(i) is a member of b.
%		loc -  32 or 64 bit (depending on whether x64 or not) int array of the same size as a, containing location of each element of a in array b. lic(i) is 0 for all elements in a which are not part of b.
%	Notice: 
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
% 	b = [3 : 5 : 20000000]; 
% 	a = [1: 3 : 20000000]; 
% 	tic, 
% 	[c, loc]= ismember(a,b); 
% 	t1 = toc; tic, 
% 	[c1, loc1] = fast_ismember_sorted(a,b); 
% 	t2 = toc; 
% 	fprintf(1,'ismember: %2.2f sec. fast_ismember_sorted: %2.2f (%2.2f faster)\n',t1,t2,t1/t2);
% 	if any([nnz(c~=c1) nnz(loc~=loc1)]) 
% 		error('Results differ!!!!');
% 	else 
% 		disp('Verified');
% 	end	
%
% Created By:
%	Dr. Lev Muchnik, LevMuchnik@gmail.com
%	Stern School Of Business, NYU
%	Oct. 26, 2008