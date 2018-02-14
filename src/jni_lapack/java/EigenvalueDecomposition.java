package JAMAJni;

/*import JAMAJni.jni_blas.*;

/** Eigenvalues and eigenvectors of a real matrix.
 <P>
    If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
    diagonal and the eigenvector matrix V is orthogonal.
    I.e. A = V.times(D.times(V.transpose())) and
    V.times(V.transpose()) equals the identity matrix.
 <P>
    If A is not symmetric, then the eigenvalue matrix D is block diagonal
    with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
    lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda].  The
    columns of V represent the eigenvectors in the sense that A*V = V*D,
    i.e. A.times(V) equals V.times(D).  The matrix V may be badly
    conditioned, or even singular, so the validity of the equation
    A = V*D*inverse(V) depends upon V.cond().
 **/

public class EigenvalueDecomposition  implements java.io.Serializable {
 private EigenvalueDecomposition() {}
 static {
     
    /* load library (which will contain wrapper for cblas function.)*/
    System.loadLibrary("lapack_lite_EigenvalueDecomposition");
 
 }

   /* ------------------------
    * Class variables
    * ------------------------ */
    /** Row and column dimension (square matrix).*/
    private int n;

    /** Symmetry flag.*/
    private boolean issymmetric;

    /**Array for elements storage of A */
    private double[] a;

    /**Array for eigenvalue of A */
    private double[] w;

  /* ------------------------
   * Constructor
   * ------------------------ */
    
    /** Use upper triangular part of a symmetric matrix,
     to construct the eigenvalue decomposition
     Structure to access D and V.
     @param A    Square matrix
     */
    
    public EigenvalueDecomposition (Matrix A) {
      a = A.getRowPackedCopy();
      n = A.getColumnDimension();
      char jobvl = EigenvalueDecomposition.JOBV.Compute;
      char uplo = EigenvalueDecomposition.UPLO.Upper;
      int matrix_layout =  EigenvalueDecomposition.LAYOUT.RowMajor;
      int lda = n;
      double[] w = new double[n];
      int lwork = 4*n+1;
      double[] work = new double[lwork];
      int[] info = new int[]{0};
      dsyev(matrix_layout, jobvl, uplo, n, a, lda, w, work, lwork, info);
    }

    /* ------------------------
     * Public Methods
     * ------------------------ */
      /** Return the eigenvector matrix */
      public Matrix getV () {
	     Matrix X = new Matrix(n,n);
         double[][] V = X.getArray();
	     for (int i = 0; i < n; i++) {
		     for (int j = 0; j < n; j++) {
			     V[i][j] = a[i*n+j];
             }
	     }
	    return X;
      }

      /**Return the singular values of A */
      public double[] getD () {
             return w;
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
    
    
    
    public static native void dsyev(int matrix_layout, char jobz, char uplo, int n,
                                    double[] a, int lda, double[] w, double[] work,
                                    int lwork, int[] info);
 
    /**inform java virtual machine that function is defined externally*/
    
}







