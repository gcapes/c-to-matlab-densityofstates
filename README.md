# Request to convert C code from a thesis to MATLAB
The C code received didn't work, so I had to fix that first.

Once working, I tested the C code with the following inputs and plotted
against the figures in the thesis to confirm the code worked.

Inputs given to the C code have been 
- 5,5,0,0,1000,1000
- 9,0,0,0,1000,1000
- 10,5,0,0,1000,1000

I then rewrote the code in MATLAB so it can be modified by the user (who
doesn't know C).

It has been tested against the original C code with the following inputs:
- 5,5,0,0,1000,1000
- 5,5,0,0,10000,10000
