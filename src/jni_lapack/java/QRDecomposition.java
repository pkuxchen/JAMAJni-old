package JAMAJni;

/** QR Decomposition.
 <P>
    For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
    orthogonal matrix Q and an n-by-n upper triangular matrix R so that
    A = Q*R.
 <P>
    The QR decompostion always exists, even if the matrix does not have
    full rank, so the constructor will never fail.  The primary use of the
    QR decomposition is in the least squares solution of nonsquare systems
    of simultaneous linear equations.  This will fail if isFullRank()
    returns false.
 */

public class QRDecomposition implements java.io.Serializable {
 private QRDecomposition() {}
 static {
     
    /* load library (which will contain wrapper for cblas function.)*/
    System.loadLibrary("lapack_lite_QRDecomposition");
 }
     
    /* ------------------------
     *  Class variables
     *  ------------------------ */
    /** Array for internal storage of decomposition.*/
    private double[] QR;
    private double[] tau;
    
    /** Row and column dimensions.*/
    private int m, n;

    /* ------------------------
     Constructor
     * ------------------------ */
    public QRDecomposition (Matrix Arg) {
	    m = Arg.getRowDimension();
	    n = Arg.getColumnDimension();
/*        if (m < n) {
            throw new
            IllegalArgumentException(
            "Matrix must have number of rows larger than number of columns.");
        }
*/
        QR = Arg.getRowPackedCopy();
	    int lda = n;
/*	    int[] jpvt = new int[] {1, 2, 3}; */
        tau = new double[(m < n) ? m : n];
	    int lwork = 3*n;
	    double[] work = new double[lwork];
        int k = m;
	    int [] info = new int [1];
	    int matrix_layout = QRDecomposition.LAYOUT.RowMajor;
	    dgeqrf(matrix_layout, m, n, QR, lda, tau, work, lwork, info);
/*	    QR2 = QR;
        dorgqr(matrix_layout, m, n, k, QR2, lda, tau, work, lwork, info);
	    QR3 = Arg.getRowPackedCopy();
	    dgeqp3(matrix_layout, m, n, QR3, lda, jpvt, tau, work, lwork, info);*/
    }

     /* ------------------------
      * Public Methods
      * ------------------------ */

    /** Is the matrix full column rank?
     @return     true if R, and hence A, has full rank.
     */
    public boolean isFullRank () {
        for (int j = 0; j < n; j++) {
            if (QR[j*n + j] == 0)
                return false;
        }
        return true;
    }
       
    /** Return the upper triangular factor*/
    public Matrix getR () {
        Matrix X = new Matrix(n,n);
	    double[][] R = X.getArray();
	    for (int i = 0; i < n; i++) {
		    for (int j = 0; j < n; j++) {
                if(i <= j){
                    R[i][j] = QR[i*n+j];
                } else {
                    R[i][j] = 0.0;
                }
		    }
	    }
	    return X;
    }


   /** Generate and return the (economy-sized) orthogonal factor*/
    public Matrix getQ () {
        int lda = n;
        int lwork = 3*n;
        int k = n;
        int [] info = new int [1];
        int matrix_layout = QRDecomposition.LAYOUT.RowMajor;
        double[] work = new double[lwork];
        double[] Q = QR;
        dorgqr(matrix_layout, m, n, k, Q, lda, tau, work, lwork, info);
        return new Matrix(Q, m);
    }

    /** Least squares solution of A*X = B
     @param B    A Matrix with as many rows as A and any number of columns.
     @return     X that minimizes the two norm of Q*R*X-B.
     @exception  IllegalArgumentException  Matrix row dimensions must agree.
     @exception  RuntimeException  Matrix is rank deficient.
     */
/*
    public Matrix solve (Matrix B) {
        if (B.getRowDimension() != m) {
            throw new IllegalArgumentException("Matrix row dimensions must agree.");
        }
        if (!this.isFullRank()) {
            throw new RuntimeException("Matrix is rank deficient.");
        }
        
        // Copy right hand side
        int nx = B.getColumnDimension();
        double[][] X = B.getArrayCopy();
        
        // Compute Y = transpose(Q)*B
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < nx; j++) {
                double s = 0.0;
                for (int i = k; i < m; i++) {
                    s += QR[i][k]*X[i][j];
                }
                s = -s/QR[k][k];
                for (int i = k; i < m; i++) {
                    X[i][j] += s*QR[i][k];
                }
            }
        }
        // Solve R*X = Y;
        for (int k = n-1; k >= 0; k--) {
            for (int j = 0; j < nx; j++) {
                X[k][j] /= Rdiag[k];
            }
            for (int i = 0; i < k; i++) {
                for (int j = 0; j < nx; j++) {
                    X[i][j] -= X[k][j]*QR[i][k];
                }
            }
        }
        return (new Matrix(X,n,nx).getMatrix(0,n-1,0,nx-1));
    }
*/
 
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
    
    
    /* QR */
    public static native int dgeqrf(int matrix_layout, int m, int n, double[] a,
                                    int lda, double[] tau, double[] work, int lwork,
                                    int[] info);
    
    public static native int dorgqr(int matrix_layout, int m, int n, int k,
                                    double[] a, int lda, double[] tau, double[] work,
                                    int lwork, int[] info);
    
    public static native int dgeqp3(int matrix_layout, int m, int n, double[] a,
                                    int lda, int[] jpvt, double[] tau, double[] work,
                                    int lwork, int[] info);
    
    
    private static final long serialVersionUID = 1;
    /**inform java virtual machine that function is defined externally*/
    
}







