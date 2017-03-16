function [Indexes, Values] = fast_find_range_sorted(Data,Ranges)
%	Returns indeces and (optionally) values in the sorted Data array that fit the range in Ranges. 
%	Very efficient (usually, 3-10 times faster) replacement of the standard MatLab intersect function. 
% 
%
%	Receives:
%		Data	- numeric vector of any type. Assumed to be sorted. 
%		Ranges	- a matrix nx2 where each row represents a range. All values in column 1 must be smaller than column 2. 
%	Notice, Data and Ranges must be of the same numeric type (both doubles, int16 , ...)
%
%	Returns:
%          Indexes -  nx2 - list of indexes into Data where each raw corresponds to the appropriate range in Ranges. Each row lists the starting and ending index of the Data so that the values in between fit into the range. 
%           Values  - cell array, nx1 - (optional) List of values in Data that fit Ranges. 
%
% Performance:
%	Due to the following reasons, the function is considerably faster than the standard union method :
%	  1. It assums that a and b are sorted sets of unique values
%	  2. It's implemented in C++
%     3. ai and bi are integres (64bit on x64 platform and 32 bit on 32 bit platform), which saves RAM and avoids conversion of integers to doubles.
%
% Example:
%
%
% Created By:
%	Dr. Lev Muchnik, LevMuchnik@gmail.com
%	Stern School Of Business, NYU
%	Jul 17, 2012