package JAMAJni;

public class LUDecomposition implements java.io.Serializable {
 private LUDecomposition() {}
 static {
    /* load library (which will contain wrapper for cblas function.)*/
    System.loadLibrary("lapack_lite_LUDecomposition");
 }
   
    /* ------------------------
     * Class variables
     * ------------------------ */
    /** Array for internal storage of decomposition.*/
    private double[] LU;

    /** Row and column dimensions, and pivot sign.*/
    private int m, n, pivsign;

    /** Internal storage of pivot vector.*/
    private int[] ipiv;

    /* ------------------------
     * Constructor
     * ------------------------ */
    public  LUDecomposition (Matrix Arg) {
	 LU = Arg.getRowPackedCopy();
	 m = Arg.getRowDimension();
	 n = Arg.getColumnDimension();
	 int matrix_layout = LUDecomposition.LAYOUT.RowMajor;
	 int [] info = new int [1];
	 int [] ipiv = new int[(m < n) ? m : n];
	 dgetrf(matrix_layout, m, n, LU, m, ipiv, info);
	}

    /* ------------------------
     * Public Methods
     * ------------------------ */

    /** Is the matrix nonsingular?*/
    public boolean isNonsingular () {
        for (int j = 0; j < n; j++) {
            if (LU[j*n + j] == 0)
                return false;
        }
        return true;
    }

   /** Return lower triangular factor */
    public Matrix getL () {
       Matrix X = new Matrix(m,n);
       double[][] L = X.getArray();
       for (int i = 0; i < m; i++) {
	       for (int j = 0; j < n; j++) {
		       if (i > j) {
			       L[i][j] = LU[i*m+j];
 	       } else if (i == j) {
			       L[i][j] = 1.0;
		       } else {
                   L[i][j] = 0.0;
		       }
	       }
       }
       return X;
    }

    /** Return upper triangular factor*/
    public Matrix getU () {
        Matrix X = new Matrix(n,n);
        double[][] U = X.getArray();
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++) {
                U[i][j] = ((i <= j) ? LU[i*m+j] : 0.0);
            }
        }
        return X;
    }
        
    /** Return pivot permutation vector */
    public int[] getPivot () {
        int[] p = new int[m];
        for (int i = 0; i < m; i++) {
            p[i] = ipiv[i];
        }
        return p;
    }
    
    /** Return pivot permutation vector as a one-dimensional double array
     @return     (double) piv
     */
    
    public double[] getDoublePivot () {
        double[] vals = new double[m];
        for (int i = 0; i < m; i++) {
            vals[i] = (double) ipiv[i];
        }
        return vals;
    }

    /** Determinant
     @return     det(A)
     @exception  IllegalArgumentException  Matrix must be square
     */
    
    public double det () {
        if (m != n) {
            throw new IllegalArgumentException("Matrix must be square.");
        }
        double d = (double) pivsign;
        for (int j = 0; j < n; j++) {
            d *= LU[j*n + j];
        }
        return d;
    }
    
    
    /** Solve A*X = B*/
    public Matrix solve (Matrix B) {
        if (B.getRowDimension() != m) {
            throw new IllegalArgumentException("Matrix row dimensions must agree.");
        }
        if (m != n) {
            throw new IllegalArgumentException("Matrix must be square.");
        }
        if (!this.isNonsingular()) {
            throw new RuntimeException("Matrix is singular.");
        }
        
        int matrix_layout = LUDecomposition.LAYOUT.RowMajor;
        char trans = LUDecomposition.TRANSPOSE.NoTrans;
        int nrhs = B.getColumnDimension();
        int ldb = nrhs;
        int lda = m;
        double[] b = B.getRowPackedCopy();
        int [] info = new int [1];
        dgetrs(matrix_layout, trans, m, nrhs, LU, lda, ipiv, b, ldb, info);
        Matrix C = new Matrix(b, B.getRowDimension());
        return C;
    }

   

    public final static class LAYOUT {
        private LAYOUT() {}
        public final static int RowMajor = 101;
        public final static int ColMajor = 102;
    }
    
    public final static class TRANSPOSE {
        private TRANSPOSE() {}
        public final static char NoTrans = 'N';         /** trans='N' */
        public final static char Trans= 'T';            /** trans='T' */
        public final static char ConjTrans= 'C';        /** trans='C' */
    }
    
    public final static class UPLO {
        private UPLO() {}
        public final static char Upper = 'U';           /** Upper triangular matrix */
        public final static char Lower= 'L';            /** Lower triangular matrix*/
    }
    
    public final static class JOBV {
        private JOBV() {}
        public final static char NoCompute = 'N';       /** eigenvectors are not computed */
        public final static char Compute= 'V';          /** eigenvectors are computed*/
    }
    
    public final static class JOB {
        private JOB() {}
        public final static char All = 'A';             /** all M columns of U are returned in array U */
        public final static char firstInU = 'S';        /** the first min(m,n) columns of U (the left singular
                                                         vectors) are returned in the array U;*/
        public final static char Overwritten = 'O';     /** the first min(m,n) columns of U (the left singular
                                                         vectors) are overwritten on the array A; */
        public final static char NoCompute = 'N';       /** no columns of U (no left singular vectors) are
                                                         computed.*/
    }
    
    public final static class ITYPE {
	private ITYPE(){}
	public final static int first = 1;
	public final static int second = 2;
	public final static int third = 3;
    }
    
    /* LU */
    public static native int dgetrf(int matrix_layout, int m, int n, double[] a,
                                    int lda, int[] ipiv, int[] info);
    
    public static native void dgetrs(int matrix_layout, char trans, int n,
                                     int nrhs, double[] a, int lda, int[] ipiv,
                                     double[] b, int ldb, int[] info);
    
}
    
    








