package JAMAJni;

public class SingularValueDecomposition implements java.io.Serializable {
 private SingularValueDecomposition() {}
 static {
     
    /* load library (which will contain wrapper for cblas function.)*/
    System.loadLibrary("lapack_lite_SingularValueDecomposition");
 
 }

 /* ------------------------
  * Class variables
  * ------------------------ */

    /** Arrays for internal storage of U and V.*/
    private double[] a,u,v;

    /** Array for internal storage of singular values.*/
    private double[] s;

    /** Row and column dimensions.*/
    private int m, n;

 /* ------------------------
  * Constructor
  * ------------------------ */
    public SingularValueDecomposition (Matrix Arg) {
        a = Arg.getRowPackedCopy();
        m = Arg.getRowDimension();
        n = Arg.getColumnDimension();
        int matrix_layout = SingularValueDecomposition.LAYOUT.RowMajor;
        char jobu = SingularValueDecomposition.JOB.All;
        char jobvt = SingularValueDecomposition.JOB.All;
        int lda = m;
        int lwork = 4*n+1;
        int ldu = n;
        int [] info = new int[1];
        int ldvt = n;
        double[] work = new double[lwork];
        double[] S = new double[n];
        double[] U = new double[m*m];
        double[] V = new double[n*n];
        double[] s = S;
        double[] u = U;
        double[] v = V;
        dgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, v, ldvt, work, lwork, info);
    
    }

    /* ------------------------
     * Public Methods
     * ------------------------ */
    /** Return the left singular vectors U*/
    public Matrix getU () {
	   Matrix X = new Matrix(m,n);
           double[][] UU = X.getArray();
	   for (int i = 0; i < m; i++) {
		   for (int j = 0; j < n; j++) {
                        UU[i][j] = this.u[i*n+j];
		   }
		  }
	   return X;
    }


    /** Return the right singular vectors V */
     public Matrix getV () {
	   Matrix X = new Matrix(n,n);
	   double[][] VV = X.getArray();
	   for (int i = 0; i < n; i++) {
		   for (int j = 0; j < n; j++) {
			   VV[i][j] = this.v[i*n+j];
		   }
	   }
		   return X;
     }

    /** Return the one-dimensional array of singular values
     @return     diagonal of S.
     */
    
    public double[] getSingularValues () {
        return s;
    }

    /** Return the diagonal matrix of singular values
     @return     S
     */
    public Matrix getS () {
        Matrix X = new Matrix(n,n);
        double[][] SS = X.getArray();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                SS[i][j] = 0.0;
            }
            SS[i][i] = this.s[i];
        }
        return X;
    }

    /** Two norm
     @return     max(S)
     */
    
    public double norm2 () {
        return s[0];
    }
    
    /** Two norm condition number
     @return     max(S)/min(S)
     */
    
    public double cond () {
        return s[0]/s[Math.min(m,n)-1];
    }
    
    public int rank () {
        double eps = Math.pow(2.0,-52.0);
        double tol = Math.max(m,n)*s[0]*eps;
        int r = 0;
        for (int i = 0; i < s.length; i++) {
            if (s[i] > tol) {
                r++;
            }
        }
        return r;
    }
    private static final long serialVersionUID = 1;

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
    
    
    /* Eigenvector and SVD */
    
    public static native int dgesvd(int matrix_layout, char jobu, char jobvt,
                                    int m, int n, double[] a, int lda, double[] s,
                                    double[] u, int ldu, double[] vt, int ldvt,
                                    double[] work, int lwork, int[] info);
    
    /**inform java virtual machine that function is defined externally*/
    
 }




