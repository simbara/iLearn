package prj.anyapp;

public class Strassen2MatrixMultiplication implements doubleMatrixMultiplication {

    public String getName () {
        return "Strassen(cached)";
    }

    double [][] p1;
    double [][] p2;
    double [][] p3;
    double [][] p4;
    double [][] p5;
    double [][] p6;
    double [][] p7;
    double [][] t0;
    double [][] t1;

    public double[][] mult (double[][] c, double[][] a, double[][] b) {
        final int n = c.length;

        if (p1 == null || p1.length < n) {
            p1 = new double[n/2][n-1];
            p2 = new double[n/2][n-1];
            p3 = new double[n/2][n-1];
            p4 = new double[n/2][n-1];
            p5 = new double[n/2][n-1];
            p6 = new double[n/2][n-1];
            p7 = new double[n/2][n-1];
            t0 = new double[n/2][n-1];
            t1 = new double[n/2][n-1];
        }

        mult(c, a, b, 0, 0, n, 0);

        return c;
    }

    public void mult (double[][] c, double[][] a, double[][] b, int i0, int j0, int n, int offs) {
        if(n == 1) {
            c[i0][j0] = a[i0][j0] * b[i0][j0];
        } else {
            final int nBy2 = n/2;

            final int i1 = i0 + nBy2;
            final int j1 = j0 + nBy2;

            // offset applied to 'p' j index so recursive calls don't overwrite data
            final int jp0 = offs;
            final int jp1 = nBy2 + offs;

            // P1 <- (A11 + A22)(B11 + B22)
            //  T0 <- (A11 + A22), T1 <- (B11 + B22), P1 <- T0*T1
            for (int i = 0; i < nBy2; ++i) {
                for (int j = 0; j < nBy2; ++j) {
                    t0[i + i0][j + jp0] = a[i + i0][j + j0] + a[i + i1][j + j1];
                    t1[i + i0][j + jp0] = b[i + i0][j + j0] + b[i + i1][j + j1];
                }
            }

            mult(p1, t0, t1, i0, jp0, nBy2, offs + nBy2);

            // P2 <- (A21 + A22)B11
            //  T0 <- (A21 + A22), T1 <- B11, P2 <- T0*T1
            for (int i = 0; i < nBy2; ++i) {
                for (int j = 0; j < nBy2; ++j) {
                    t0[i + i0][j + jp0] = a[i + i1][j + j0] + a[i + i1][j + j1];
                    t1[i + i0][j + jp0] = b[i + i0][j + j0];
                }
            }

            mult(p2, t0, t1, i0, jp0, nBy2, offs + nBy2);

            // P3 <- A11(B12 - B22)
            //  T0 <- A11, T1 <- (B12 - B22), P3 <- T0*T1
            for (int i = 0; i < nBy2; ++i) {
                for (int j = 0; j < nBy2; ++j) {
                    t0[i + i0][j + jp0] = a[i + i0][j + j0];
                    t1[i + i0][j + jp0] = b[i + i0][j + j1] - b[i + i1][j + j1];
                }
            }

            mult(p3, t0, t1, i0, jp0, nBy2, offs + nBy2);

            // P4 <- A22(B21 - B11)
            //  T0 <- A22, T1 <- (B21 - B11), P4 <- T0*T1
            for (int i = 0; i < nBy2; ++i) {
                for (int j = 0; j < nBy2; ++j) {
                    t0[i + i0][j + jp0] = a[i + i1][j + j1];
                    t1[i + i0][j + jp0] = b[i + i1][j + j0] - b[i + i0][j + j0];
                }
            }

            mult(p4, t0, t1, i0, jp0, nBy2, offs + nBy2);

            // P5 <- (A11 + A12) B22
            //  T0 <- (A11 + A12), T1 <- B22, P5 <- T0*T1
            for (int i = 0; i < nBy2; ++i) {
                for (int j = 0; j < nBy2; ++j) {
                    t0[i + i0][j + jp0] = a[i + i0][j + j0] + a[i + i0][j + j1];
                    t1[i + i0][j + jp0] = b[i + i1][j + j1];
                }
            }

            mult(p5, t0, t1, i0, jp0, nBy2, offs + nBy2);

            // P6 <- (A21 - A11)(B11 - B12)
            //  T0 <- (A21 - A11), T1 <- (B11 - B12), P6 <- T0 * T1
            for (int i = 0; i < nBy2; ++i) {
                for (int j = 0; j < nBy2; ++j) {
                    t0[i + i0][j + jp0] = a[i + i1][j + j0] - a[i + i0][j + j0];
                    t1[i + i0][j + jp0] = b[i + i0][j + j0] - b[i + i0][j + j1];
                }
            }

            mult(p6, t0, t1, i0, jp0, nBy2, offs + nBy2);

            // P7 <- (A12 - A22)(B21 + B22)
            //  T0 <- (A12 - A22), T1 <- (B21 + B22), P7 <- T0 * T1
            for (int i = 0; i < nBy2; ++i) {
                for (int j = 0; j < nBy2; ++j) {
                    t0[i + i0][j + jp0] = a[i + i0][j + j1] - a[i + i1][j + j1];
                    t1[i + i0][j + jp0] = b[i + i1][j + j0] + b[i + i1][j + j1];
                }
            }

            mult(p7, t0, t1, i0, jp0, nBy2, offs + nBy2);

            // combine
            for (int i = 0; i < nBy2; ++i) {
                for (int j = 0; j < nBy2; ++j) {
                    // C11 = P1 + P4 - P5 + P7;
                    c[i + i0][j + j0] = p1[i + i0][j + jp0] + p4[i + i0][j + jp0] - p5[i + i0][j + jp0] + p7[i + i0][j + jp0];
                    // C12 = P3 + P5;
                    c[i + i0][j + j1] = p3[i + i0][j + jp0] + p5[i + i0][j + jp0];
                    // C21 = P2 + P4;
                    c[i + i1][j + j0] = p2[i + i0][j + jp0] + p4[i + i0][j + jp0];
                    // C22 = P1 + P3 - P2 + P6;
                    c[i + i1][j + j1] = p1[i + i0][j + jp0] + p3[i + i0][j + jp0] - p2[i + i0][j + jp0] + p6[i + i0][j + jp0];
                }
            }
        }
    }

//    void dumpInternal () {
//        System.out.println("P1");
//        TestIntMatrixMultiplication.dumpMatrix(p1);
//        System.out.println("P2");
//        TestIntMatrixMultiplication.dumpMatrix(p2);
//        System.out.println("P3");
//        TestIntMatrixMultiplication.dumpMatrix(p3);
//        System.out.println("P4");
//        TestIntMatrixMultiplication.dumpMatrix(p4);
//        System.out.println("P5");
//        TestIntMatrixMultiplication.dumpMatrix(p5);
//        System.out.println("P6");
//        TestIntMatrixMultiplication.dumpMatrix(p6);
//        System.out.println("P7");
//        TestIntMatrixMultiplication.dumpMatrix(p7);
//        System.out.println("T0");
//        TestIntMatrixMultiplication.dumpMatrix(t0);
//        System.out.println("T1");
//        TestIntMatrixMultiplication.dumpMatrix(t1);
//    }
}
