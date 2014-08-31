          Streaming clustering coefficients experimental rig
          ==================================================

Author: Jason Riedy <jason.riedy@cc.gatech.edu>
Date: 2011-11-11 13:38:52 EDT


Table of Contents
=================
1 Setup 
2 Generating data 
3 Running the experimental rig 


1 Setup 
~~~~~~~~

On the Cray XMT, both `gen-streams' and `main' use the snapshot I/O
interface.  You need to run fsworker from the data file directory
and set `SWORKER\_EP' to the output id string.

2 Generating data 
~~~~~~~~~~~~~~~~~~

The primary parameters are `--scale', `--edgefact', and `--nact'.
To generate a graph of scale 24, edge factor 11 (so 2^24 vertices
and 11 * 2^24 edges), and 60 actions:
./gen-streams --scale=24 --edgefact=11 --nact=60 graph.bin actions.bin

Other parameters are given by running `./gen-streams --help'.  Note
that the data generator is only minimally parallelized.

3 Running the experimental rig 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

./main graph.bin actions.bin

For more details on the STINGER API, see ROADMAP

