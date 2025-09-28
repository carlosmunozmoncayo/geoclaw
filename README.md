#### A minimal setup to reproduce instabilities arising when solving the hyperbolic approximation of SGN by Escalante et al. (2019)

This code is behind the master branch of GeoClaw by a large number of commits.

To avoid compatibility issues, in a clawpack root, one would have to replace geoclaw by a copy of the folder in this repository
(before installation if it's not an editable one). It works for me with Clawpack 5.10.0 and 5.13.1.
Also, one would have to source the local_CLAW_env.sh file in this repository to set some environment variables.