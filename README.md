# guide_bot
guide_bot is a MATLAB script that generate McStas / iFit guide optimization code from a simple input script

The use of guide_bot mainly consists of filling out a input file that defines the guide optimization task at hand.
The input files are called guide_bot_user_*.m, and includes the required input fields.
The required input consists of:
Demands, which is a description of the requested beam in terms of size, divergence, distance from source and so forth.
Requirements, which are facility related, such as closest and furthest allowed guide start and what source component should be used.
Options, which affect both assumptions made when doing some calculations for example relating to the Minimalist principle, but also contains computing options.
Defaults, which mainly describe the coating to be used in the optimized guides.

Next input strings are given that in one line describe the overall geometry of a guide to be optimized, for example:
  "S C P"
This means a Straight guide followed by a Curved guide followed by a Parabolic guide section. This corresponds to a large parameter space as no information about the geometry is given.
Constraints can be given in parenthesis behind each module, for example:
  "S(minlength=5) C(StartWidth=0.03,maxEndHeight=) P(maxlength=3)"
Which now sets a minimum length of 5 m for the straight section, a  3 cm width at the start of the curved guide and a 5 cm height at the end, and a maximum length of 3 m for the parabolic module.
Many of such input strings can be given.

When executing the guide_bot_user_*.m file, a project folder is generated containing:
  compile_all.sh : script for compiling all generated code
  launch_all.sh  : script for launching all jobs on a supported cluster
  A folder for each input string containing:
    McStas instrument files describing the desired geometry and with appropriate input.
    iFit MATLAB optimization script(s) including all parameters to be optimized and their limits
    All necessary McStas components and data files for sources
  A folder containing results including:
    Scripts for plotting characterization of each optimized guide
