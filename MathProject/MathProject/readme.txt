
-----Table of Contents-----
	1. System Requirements
	2. How to run the program
	3. Where is ___ in the program?
	4. How to input data
	5. Where is code?





1. ---------------System Requirements----------------

	-Ensure that you have the latest java jre (java runtime environment) installed on your computer
	and have the environmental variables set up, such that java files can be called in the command line.
		
	
	
2. ---------------How to Run the Program-------------------------

	-In the Program Files folder, there are three class files which are the entirety of the program
		(Driver.class, Matrix.class, MatrixMath.class).
		
	-To run, open a command line interface, with the path set to this folder. Then enter "java Driver" in order
		to begin the program.
		
	-The program will then prompt the user to enter a number corresponding to the options shown, in order to
		access certain functions of the program.
		
	-**IMPORTANT**
		-Certain functions of the program may ask that the user input a file name containing the
			Matrix which will be operated on. The file MUST be in the SAME folder as the Driver.class
			file, or else it may not be found.
		-LU and QR decomposition and Power Method will require square matrices, while Ax = b functions will
			require the square matrix A augmented with the b vector.
		-Only files that use spaces between Matrix entries can be read.


3. ---------------Where is __ in the Program?-----------------------

	-LU Decomposition is the first option available to the user.

	-QR Decomposition can be found just below that as option #2, which includes both Givens and
		Householder QR factorization.
		
	-Four Ax=b functions can be found in the 3rd option, those being:
		solving through LU decomp, solving through QR decomp,
		solving iteratively with jacobi's method, and solving iteratively with the Gauss-Seidel method.
		
	-The Power method is the 4th option.

	-Hilbert Matrix calculations, including error and solutions, can be found in the 5th option,
		which will output data to a txt file in the folder of Driver.class, rather than outputting
		to the screen.

	-Encoding and Decoding of binary streams can be found in the 6th option, which will first let the
		user determine which method they would like to decode with, then will randomly generate
		a binary stream of a user defined length, then encode and decode the stream.
		
	-In options, the user may limit the output of Matrices on the screen, such that they will be more
		readable, as they will show fewer significant digits.
		
	-The 8th option exits the program gracefully.


4. ---------------How to Input Data------------------

	-When navigating the menus of the program, the user need only enter one of the numbers shown to
		select that option.

	-When prompted for a file, the user must place the file containing the Matrix into the same folder
		as the Driver.class file, so that it may be found by the program. Then its full name must be
		input, i.e. "a.dat" "MatrixA.txt" "file.extension"

	-If the user is prompted to input an initial vector, he/she must do so by inputting values separated
		by spaces. If too many values are given, the program will use the first ones it reads. If too few,
		the user will be prompted for the rest.
		
	-If the user is prompted for a tolerance or length of a vector, he/she need only input the number.


5. ---------------Where is code?---------------------

	-The Driver class mainly contains parsers, input readers, and general functions that allow the
		user to interact with the program and see the program's output.
		**Error is always calculated with the outputs of the MatrixMath class inside of the Driver,
			however, so the functions within Driver have been named according to their function,
			so LU decomp/Power Method/etc. may easily be found.**
	
	-The Matrix class only contains functions and data for a general Matrix object, such as the addition,
		subtraction, multiplication, and transposing of Matrices. It does contain the code for generating
		a Hilbert Matrix, however.
		
	-The MatrixMath class contains the bulk of the "thinking" of the program. It performs all of the
		decompositions, iterative methods, solving methods, and encoding methods needed for the project.
		**The functions are named according to the project description.**

	-The code for the calculations of part 3 of the project description can be found in the test file
		"LeslieMatrixCalc," if needed.




