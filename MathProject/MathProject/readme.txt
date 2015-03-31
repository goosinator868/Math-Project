
-----Table of Contents-----
	1. System Requirements
	2. How to Run
	3. Where is...?
	4. How to input data





1. ---------------System Requirements----------------

	-Ensure that you have the latest java jre (java runtime environment) installed on your computer
	and have the environmental variables set up, such that java files can be called in the command line.
		
	
	
2. ---------------How to Run-------------------------

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


3. ---------------Where is...?-----------------------

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

	-When prompted for a file, the user must place the file containing the Matrix into the same folder
		as the Driver.class file, so that it may be found by the program. Then its full name must be
		input, i.e. "a.dat" "MatrixA.txt" "file.extension"

	-If the user is prompted to input an initial vector, he/she must do so by inputting values separated
		by spaces. If too many values are given, the program will use the first ones it reads. If too few,
		the user will be prompted for the rest.
		
	-If the user is prompted for a tolerance or length of a vector, he/she need only input the number.



