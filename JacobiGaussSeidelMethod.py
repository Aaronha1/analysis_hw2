"""
 * Authors: Aaron Hajaj (ID: 311338198) 
 *  
 * Jacobi Method and GaussSeidel Method
 *  
 * git: https://github.com/Aaronha1/analysis_hw2.git
"""


def PrintMatrix(matrix):
    """
    Matrix Printing Function
    :param matrix: Matrix nxn
    """
    for line in matrix:
        print('  '.join(map(str, line)))


def MaxNorm(matrix):
    """
    Function for calculating the max-norm of a matrix
    :param matrix: Matrix nxn
    :return:max-norm of a matrix
    """
    max_norm = 0
    for i in range(len(matrix)):
        norm = 0
        for j in range(len(matrix)):
            # Sum of organs per line with absolute value
            norm += abs(matrix[i][j])
        # Maximum row amount
        if norm > max_norm:
            max_norm = norm

    return max_norm


def Determinant(matrix, mul=1):
    """
    Recursive function for determinant calculation
    :param matrix: Matrix nxn
    :param mul: The double number
    :return: determinant of matrix
    """
    width = len(matrix)
    # Stop Conditions
    if width == 1:
        return mul * matrix[0][0]
    else:
        sign = -1
        det = 0
        for i in range(width):
            m = []
            for j in range(1, width):
                buff = []
                for k in range(width):
                    if k != i:
                        buff.append(matrix[j][k])
                m.append(buff)
            # Change the sign of the multiply number
            sign *= -1
            #  Recursive call for determinant calculation
            det = det + mul * Determinant(m, sign * matrix[0][i])
    return det


def MakeIMatrix(cols, rows):
    # Initialize a identity matrix
    return [[1 if x == y else 0 for y in range(cols)] for x in range(rows)]


def MultiplyMatrix(matrixA, matrixB):
    """
    Function for multiplying 2 matrices
    :param matrixA: Matrix nxn
    :param matrixB: Matrix nxn
    :return: Multiplication between 2 matrices
    """
    # result matrix initialized as singularity matrix
    result = [[0 for y in range(len(matrixB[0]))] for x in range(len(matrixA))]
    for i in range(len(matrixA)):
        # iterate through columns of Y
        for j in range(len(matrixB[0])):
            # iterate through rows of Y
            for k in range(len(matrixB)):
                result[i][j] += matrixA[i][k] * matrixB[k][j]
    return result


def InverseMatrix(matrix):
    """
    Function for calculating an inverse matrix
    :param matrix:  Matrix nxn
    :return: Inverse matrix
    """
    n = len(matrix)
    vector=[0]*n
    # Unveri reversible matrix
    if Determinant(matrix) == 0:
        raise ValueError(f"Error,Singular Matrix\n{matrix[:]}")
    if not all(len(row) == n for row in matrix):
        raise ValueError("Error, no inversion to a non-square matrix\n",matrix[:])
    # result matrix initialized as singularity matrix
    result = MakeIMatrix(n, n)
    # loop for each row
    for i in range(n):
        # turn the pivot into 1 (make elementary matrix and multiply with the result matrix )
        # pivoting process
        matrix, vector = RowXchange(matrix, vector)
        elementary = MakeIMatrix(n, n)
        elementary[i][i] = 1 / matrix[i][i]
        result = MultiplyMatrix(elementary, result)
        matrix = MultiplyMatrix(elementary, matrix)
        # make elementary loop to iterate for each row and subtracrt the number below (specific) pivot to zero  (make
        # elementary matrix and multiply with the result matrix )
        for j in range(i + 1, len(matrix)):
            elementary = MakeIMatrix(n, n)
            elementary[j][i] = -(matrix[j][i])
            matrix = MultiplyMatrix(elementary, matrix)
            result = MultiplyMatrix(elementary, result)

    # after finishing with the lower part of the matrix subtract the numbers above the pivot with elementary for loop
    # (make elementary matrix and multiply with the result matrix )
    for i in range(n - 1, 0, -1):
        for j in range(i - 1, -1, -1):
            elementary = MakeIMatrix(n, n)
            elementary[j][i] = -(matrix[j][i])
            matrix = MultiplyMatrix(elementary, matrix)
            result = MultiplyMatrix(elementary, result)

    return result


def RowXchange(matrix, vector):
    """
    Function for replacing rows with both a matrix and a vector
    :param matrix: Matrix nxn
    :param vector: Vector n
    :return: Replace rows after a pivoting process
    """

    for i in range(len(matrix)):
        max = abs(matrix[i][i])
        for j in range(i, len(matrix)):
            # The pivot member is the maximum in each column
            if abs(matrix[j][i]) > max:
                matrix[j],matrix[i]=matrix[i],matrix[j]
                vector[j],vector[i]=vector[i],vector[j]
                max = abs(matrix[i][i])

    return [matrix, vector]


def CheckDominantDiagonal(matrix):
    """
    Function for testing a dominant diagonal
    :param matrix: Matrix nxn
    :return: True or false if there is a domene diagonal
    """
    for i in range(len(matrix)):
        sum = 0
        for j in range(len(matrix)):
            if i != j:
                # Sum of line member without pivot
                sum += abs(matrix[i][j])
        # If the pivot is less than the sum of the rest of the line
        if abs(matrix[i][i]) < sum:
            return False
    return True


def DominantDiagonalFix(matrix):
    """
    Function to change a matrix to create a dominant diagonal
    :param matrix: Matrix nxn
    :return: Change the matrix to a dominant diagonal
    """
    #Check if we have a dominant for each column
    dom = [0]*len(matrix)
    result = list()
   # Find the largest organ in a row
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if (matrix[i][j] > sum(map(abs,map(float ,matrix[i])))-matrix[i][j]) :
                dom[i]=j
    for i in range(len(matrix)):
        result.append([])
        # Cannot dominant diagonal
        if i not in dom:
            print("Couldn't find dominant diagonal.")
            return matrix
    # Change the matrix to a dominant diagonal
    for i,j in enumerate(dom):
        result[j]=(matrix[i])
    return result


def minusMatrix(matrix):
    """
    :param matrix: Matrix nxn
    :return: Subtract between matrices
    """
    return [[-i for i in j] for j in matrix]


def matrixAddition(matrixA, matrixB):
    """
    :param matrix: Matrix nxn, Matrix nxn
    :return: Addition between matrices
    """
    return [[a + b for (a, b) in zip(i, j)] for (i, j) in zip(matrixA, matrixB)]



def matrixDLUdissasembly(matrix):
    """
    :param matrix: matrix nXn
    :return: Breaking down matrices L U D
    """
    D, L, U = list(), list(), list()
    for x, row in enumerate(matrix):
        D.append(list()), L.append(list()), U.append(list())
        for y, value in enumerate(row):
            # Diagonal with zeros
            if x == y:
                D[x].append(value), L[x].append(0), U[x].append(0)
            # Zeros below diagonal
            elif x < y:
                D[x].append(0), L[x].append(0), U[x].append(value)
            # Zeros above diagonal
            elif x > y:
                D[x].append(0), L[x].append(value), U[x].append(0)
    return D, L, U


def JacobiG(matrix):
    """
    :param matrix: Matrix nxn
    :return: G matrix
    """
    D, L, U = matrixDLUdissasembly(matrix)
    return MultiplyMatrix(minusMatrix(InverseMatrix(D)), matrixAddition(L, U))


def JacobiH(matrix):
    """
    :param matrix: Matrix nxn
    :return: H matrix
    """
    D, L, U = matrixDLUdissasembly(matrix)
    return InverseMatrix(D)


def GaussSeidelG(matrix):
    """
    :param matrix: Matrix nxn
    :return: G matrix
    """
    D, L, U = matrixDLUdissasembly(matrix)
    return MultiplyMatrix(minusMatrix(InverseMatrix(matrixAddition(L, D))), U)


def GaussSeidelH(matrix):
    """
    :param matrix: Matrix nxn
    :return: H matrix
    """
    D, L, U = matrixDLUdissasembly(matrix)
    return InverseMatrix(matrixAddition(L, D))



def CheckJacobiGnorm(matrix):
    return 1 > MaxNorm(JacobiG(matrix))



def CheckGaussSeidelGnorm(matrix):
    return 1 > MaxNorm(GaussSeidelG(matrix))


def InitVector(size):
    return [0 for index in range(size)]

def TestEpsilon(nxt, prv, e):
    return all(abs(n-p) < e for n, p in zip(nxt, prv))
    

def JacobiMethod(matrix, vector, epsilon, previous, counter):
    """
    Function for solving a set of equations according to the Jacobi method
    :param matrix: Matrix nxn
    :param vector: Vector n
    :param epsilon: Stop Conditions
    :param previous: Result vector
    :param counter: Number of iterations
    """

    if counter==1:
        print("\n ~~~ JacobiMethod ~~~\n")
    n=len(matrix)
    NextGuess = []
    for i in range(n):
        ins = 0
        for j in range(n):
            if i != j:
                # Insulating variables
                ins = ins + matrix[i][j]*previous[j]
        # Calculate the next iteration
        newGuess = 1/matrix[i][i]*(vector[i]-ins)
        # Result vector insertion
        NextGuess.append(newGuess)

    # Result vector insertion
    print(f"Iteration no. {counter:02}: ", ', '.join(f'{float(g):.5f}' for g in NextGuess))

    # Check stop conditions
    if TestEpsilon(NextGuess, previous, epsilon):
        return

    # Recursive call
    JacobiMethod(matrix, vector, epsilon,NextGuess,counter+1)


def GaussSeidelMethod(matrix, vector, epsilon, previous, counter):
    """
     Function for solving a set of equations according to the GaussSeidel method
     :param matrix: Matrix nxn
     :param vector: Vector n
     :param epsilon: Stop Conditions
     :param previous: Result vector
     :param counter: Number of iterations
     """

    if counter==1:
        print("\n ~~~ GaussSeidelMethod ~~~\n")
    n=len(matrix)
    NextGuess = []
    ImprovedGuess = previous.copy()
    for i in range(n):
        ins = 0
        for j in range(n):
            if i != j:
                # Insulating variables
                ins = ins + matrix[i][j]*ImprovedGuess[j]
        newGuess = 1/matrix[i][i]*(vector[i]-ins)
        # Using calculated results
        ImprovedGuess[i] = newGuess
        NextGuess.append(newGuess)

    # Result vector insertion
    print(f"Iteration no. {counter:02}: ", ', '.join(f'{float(g):.5f}' for g in NextGuess))

    # Check stop conditions
    if TestEpsilon(NextGuess, previous, epsilon):
        return

    # Recursive call
    GaussSeidelMethod(matrix, vector, epsilon,NextGuess,counter+1)


matrixA = [[4, 0, 2], [3, 12, 3], [0, 4, 8]]
b = [2,6,4]
epsilon=0.001
arg=(matrixA,b,epsilon,InitVector(len(b)),1)


try:
    if CheckDominantDiagonal(matrixA):
        print("\nThere is a dominant diagonal.")
        JacobiMethod(*arg)
        GaussSeidelMethod(*arg)
    else:
        print("There isn't a dominant diagonal.")
        print("We will try to find dominant diagonal.")
        dominantFix = DominantDiagonalFix(matrixA)
        PrintMatrix(dominantFix)
        if dominantFix != matrixA:
            arg = (dominantFix, *arg[1:])
            print("Found a dominant diagonal.")
            JacobiMethod(*arg)
            GaussSeidelMethod(*arg)
        else:
            ch = CheckGaussSeidelGnorm(matrixA)*2+CheckJacobiGnorm(matrixA)
            print("\nThe matrix", "is" if ch else "isn't","convergent.")
            
            if ch%2==1:
                JacobiMethod(*arg)
            if ch>1:
                GaussSeidelMethod(*arg)
            
                
except Exception as e:
    print("Error:", e)

