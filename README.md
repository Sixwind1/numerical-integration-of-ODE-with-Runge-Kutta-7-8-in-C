# numerical-integration-of-ODE-with-Runge-Kutta-7-8-in-C
This project is part of an undergraduate assignment from Universitat autonoma de barcelona (UAB). The code has 5 different parts: qr matrix decomposition, multidimensional newton method root finder, runge kutta 7-8 for a single step, numerical flow function using rk78 and an example problem solved with the previous flow function.

The original code uses catalan as the main language in variable naming and error handling.
This version is AI-translated and has no guarantee that the comments and output text are 100% correct.
However, the output file has been checked, and it gives the same result as the original code.

Changing the parameter at execution is disabled. You have to modify the files to adapt needs.

To run the program, you should compile the following c-files:
main.c flow.c rk78.c QRsolve.c basicOperation.c


After running the executable file, you should obtain a txt file named 'ophalo.txt'.

The first column is the time and the remaining columns are the positions, being the first 3 columns the 3D position of the periodic trajectory.


You can visualize the periodic trajectory through the opHalo.py file or through your favorite visualization tool.
