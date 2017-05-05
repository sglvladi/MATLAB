GM_PHD_Filter
Version 1.05, 21 November 2013

Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au with:
- some Kalman filter update code by Tim Bailey, taken from his website http://www-personal.acfr.usyd.edu.au/tbailey/software/
- error_ellipse by AJ Johnson, taken from Matlab Central http://www.mathworks.com.au/matlabcentral/fileexchange/4705-errorellipse
Algorithm by Ba-Ngu Vo & Wing-Kin Ma in:
B.-N. Vo, W.-K. Ma, "The Gaussian Mixture Probability Hypothesis Density Filter", IEEE Transactions on Signal Processing, Vol 54, No. 11, November 2006, pp4091-4104
Optimal Subpattern Assignment (OSPA) metric by Schuhmacher, D.; Ba-Tuong Vo; Ba-Ngu Vo in:
Schuhmacher, D.; Ba-Tuong Vo; Ba-Ngu Vo, "A Consistent Metric for Performance Evaluation of Multi-Object Filters," Signal Processing, IEEE Transactions on , vol.56, no.8, pp.3447,3457, Aug. 2008
This isn't necessary for the GM-PHD filter to work but provides a nice way of visualising performance. 

This is an implementation of a gaussian mixture probability hypothesis density filter (GM-PHD) for a simulated tracking problem. The problem specification is given in the above paper - in summary, two targets move through the environment, there is a lot of clutter on the measurement, and about halfway through a third target spawns off one of the two targets.

I made a few changes, either because I couldn't understand how Vo&Ma did it, or because I wanted to make it closer to my target problem. I extended the measurement vector to include both target position AND velocity (the filter they describe tracks position and velocity but only observes position). Velocity is observed as just being dx/dt, change in position over time, between this new observation and the previous position of this target. Targets are either birthed or spawned depending on which initialisation weight function would give a higher weight; they are given the appropriate initialisation covariance, but apart from this there is no other difference in instantiation. Vo & Ma's state estimation involves repeatedly outputting targets that have a rounded weight greater than 1; I have this off by default but it can be turned on by setting OUTPUT_MULTIPLE_HIGH_WEIGHT_TARGETS and VERBOSE to 1 in GM_PHD_Initialisation.

The files are mostly implemented as scripts. The main file is GM_PHD_Filter. The structure is pretty straightforward, as it mostly just calls the other scripts in a big loop.

There are two files that initialise everything - GM_PHD_Initialise and GM_PHD_Simulate_Initialise. Changing values in these will change how the filter runs.

The other files are pretty easy to follow if you have a look in GM_PHD_Filter.m

There isn't much output on the console by default, most of it is in the graphs. Have a read of GM_PHD_Simulate_Plot to understand and edit what the different coloured dots mean. In summary, black X's are measurements, black X's with a circle around them are measurements corresponding to true targets, red/green/blue dots are true target positions, magenta circles are tracked targets, cyan triangles are the targets at the current filter iteration.

To see more filter output, set VERBOSE to 1 in GM_PHD_Initialisation. I can't remember exactly what this prints out, but it's more information about the inner working of the filter, and you can use it to add your own output statements wherever you want. I cut most of mine out because they were cluttering up the code.

To see the performance metric, set CALCULATE_OSPA_METRIC to 1 in GM_PHD_Initialisation. It appears as an extra line graph.

See license.txt for licensing information.

See ReleaseNotes.txt for changes between versions.

-- Bryan Clarke, b.clarke@acfr.usyd.edu.au