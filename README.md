# blocked-cannons-algorithm
These two programs perform parallel matrix mutliplication using Cannon's algorithm. The program mulBlockColCannon is a blocked implementation where matrix A is stored in row major format and matrix B in column major. mulBlockRowCannon is a naive implementation where both A and B are stored in row major format. The blocked implementation takes advantage of the CPU caches by improving temporal and spatial locality. A comparison of the two methods can be found in the report I have included with the program.

# Compilation
A makefile is included for compilation. Simply type 'make' or 'make all' in the command line to compile all source code. The executables will be stored in the newly created folder bin.

# Run
mpirun -np p ./exec N, where p is the number of processes to spawn and N the size of the matrices. 

# Example runs
$ mpirun -np 4 ./mulRowCannon 2000 <br />
serial time 23.6615 sec <br />
parallel time 3.1151 sec

$ mpirun -np 4 ./mulBlockColCannon 2000 <br />
serial time 5.0911 sec <br />
parallel time 1.3303 sec


# Additional notes
If you find the code helpful and want to reuse it for your homework or your project, you're free to do so. However, I leave no guarantee regarding correctness and functionality, and strongly discourage thoughtless copying a random piece of code. You're ultimately the one responsible for testing and making sure your program is correct. If you find a bug in the code, don't hesitate to contact me.
