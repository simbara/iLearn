package prj.anyapp;

public class Strassen1MatrixMultiplication implements doubleMatrixMultiplication {

    public String getName () {
        return "Strassen(dynamic)";
    }

    public double[][] mult (double[][] c, double[][] a, double[][] b) {
    	return strassenMatrixMultiplication(a, b);
    }

    public double[][] strassenMatrixMultiplication(double [][] A, double [][] B) {
        int n = A.length;

        double [][] result = new double[n][n];

        if(n == 1) {
            result[0][0] = A[0][0] * B[0][0];
        } else {
            double [][] A11 = new double[n/2][n/2];
            double [][] A12 = new double[n/2][n/2];
            double [][] A21 = new double[n/2][n/2];
            double [][] A22 = new double[n/2][n/2];

            double [][] B11 = new double[n/2][n/2];
            double [][] B12 = new double[n/2][n/2];
            double [][] B21 = new double[n/2][n/2];
            double [][] B22 = new double[n/2][n/2];

            divideArray(A, A11, 0 , 0);
            divideArray(A, A12, 0 , n/2);
            divideArray(A, A21, n/2, 0);
            divideArray(A, A22, n/2, n/2);

            divideArray(B, B11, 0 , 0);
            divideArray(B, B12, 0 , n/2);
            divideArray(B, B21, n/2, 0);
            divideArray(B, B22, n/2, n/2);

            double [][] P1 = strassenMatrixMultiplication(addMatrices(A11, A22), addMatrices(B11, B22));
            double [][] P2 = strassenMatrixMultiplication(addMatrices(A21, A22), B11);
            double [][] P3 = strassenMatrixMultiplication(A11, subtractMatrices(B12, B22));
            double [][] P4 = strassenMatrixMultiplication(A22, subtractMatrices(B21, B11));
            double [][] P5 = strassenMatrixMultiplication(addMatrices(A11, A12), B22);
            double [][] P6 = strassenMatrixMultiplication(subtractMatrices(A21, A11), addMatrices(B11, B12));
            double [][] P7 = strassenMatrixMultiplication(subtractMatrices(A12, A22), addMatrices(B21, B22));

            double [][] C11 = addMatrices(subtractMatrices(addMatrices(P1, P4), P5), P7);
            double [][] C12 = addMatrices(P3, P5);
            double [][] C21 = addMatrices(P2, P4);
            double [][] C22 = addMatrices(subtractMatrices(addMatrices(P1, P3), P2), P6);

            copySubArray(C11, result, 0 , 0);
            copySubArray(C12, result, 0 , n/2);
            copySubArray(C21, result, n/2, 0);
            copySubArray(C22, result, n/2, n/2);
        }

        return result;
    }

    public double [][] addMatrices(double [][] A, double [][] B) {
        int n = A.length;

        double [][] result = new double[n][n];

        for(int i=0; i<n; i++)
        	for(int j=0; j<n; j++)
        		result[i][j] = A[i][j] + B[i][j];

        return result;
    }

    public double [][] subtractMatrices(double [][] A, double [][] B) {
        int n = A.length;

        double [][] result = new double[n][n];

        for(int i=0; i<n; i++)
            for(int j=0; j<n; j++)
                result[i][j] = A[i][j] - B[i][j];

        return result;
    }

    public void divideArray(double[][] parent, double[][] child, int iB, int jB) {
        for(int i1 = 0, i2=iB; i1<child.length; i1++, i2++)
            for(int j1 = 0, j2=jB; j1<child.length; j1++, j2++)
                child[i1][j1] = parent[i2][j2];
    }

    public void copySubArray(double[][] child, double[][] parent, int iB, int jB) {
        for(int i1 = 0, i2=iB; i1<child.length; i1++, i2++)
            for(int j1 = 0, j2=jB; j1<child.length; j1++, j2++)
                parent[i2][j2] = child[i1][j1];
    }
}
