package JAMAJni;

import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.text.FieldPosition;
import java.io.PrintWriter;
import java.io.BufferedReader;
import java.io.StreamTokenizer;
import JAMAJni.util.*;


/**CBLAS.java*/

public class Matrix implements Cloneable, java.io.Serializable {
 private Matrix() {}
    static {
     
    /* load library (which will contain wrapper for cblas function.)*/
    System.loadLibrary("cblas_Matrix");
 
 }
    
    /**inform java virtual machine that function is defined externally*/
    
    
    /* -----------------------
     * Class variables
     * ---------------------- */
    
    private double[][] A;
    
    private int m, n;
    
    /* -----------------------
     * Constructors
     * ----------------------- */
    
    /** Construct an m-by-n matrix of zeros. */
    public Matrix (int m, int n) {
        this.m = m;
        this.n = n;
        A = new double[m][n];
    }
    
    /** Construct an m-by-n constant matrix. */
    public Matrix (int m, int n, double s) {
        this.m = m;
        this.n = n;
        A = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = s;
            }
        }
    }
    
    /** Construct a matrix from a 2-D array. */
    public Matrix (double[][] A) {
        m = A.length;
        n = A[0].length;
        for (int i = 0; i < m; i++) {
            if (A[i].length != n) {
                throw new IllegalArgumentException("All rows must have the same length.");
            }
        }
        this.A = A;
    }
    
    /** Construct a matrix quickly without checking arguments. */
    public Matrix (double[][] A, int m, int n) {
        this.A = A;
        this.m = m;
        this.n = n;
    }
    
    /** Construct a matrix from a one-dimensional packed array
     @param vals One-dimensional array of doubles, packed by columns (ala Fortran).
     @param m    Number of rows.
     @exception  IllegalArgumentException Array length must be a multiple of m.
     */
    
    public Matrix (double vals[], int m) {
        this.m = m;
        n = (m != 0 ? vals.length/m : 0);
        if (m*n != vals.length) {
            throw new IllegalArgumentException("Array length must be a multiple of m.");
        }
        A = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = vals[i+j*m];
            }
        }
    }
    
    /* ------------------------
     *    Public Methods
     * ------------------------ */
    
    /** Construct a matrix from a copy of a 2-D array.*/
    
    public static Matrix constructWithCopy(double[][] A) {
        int m = A.length;
        int n = A[0].length;
        Matrix X = new Matrix(m,n);
        double[][] C = X.getArray();
        for (int i = 0; i < m; i++) {
            if (A[i].length != n) {
                throw new IllegalArgumentException
                ("All rows must have the same length.");
            }
            for (int j = 0; j < n; j++) {
                C[i][j] = A[i][j];
            }
        }
        return X;
    }
    
    
    /** Make a deep copy of a matrix
     */
    
    public Matrix copy () {
        Matrix X = new Matrix(m,n);
        double[][] C = X.getArray();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                C[i][j] = A[i][j];
            }
        }
        return X;
    }
    
    /** Clone the Matrix object.
     */
    
    public Object clone () {
        return this.copy();
    }
    
    /** Access the internal two-dimensional array. */
    public double[][] getArray () {
        return A;
    }

     /** Copy the internal two-dimensional array.*/
     public double[][] getArrayCopy () {
	 double[][] C = new double[m][n];
         for (int i = 0; i < m; i++) {	 
		 for (int j = 0; j < n; j++) {
			 C[i][j] = A[i][j];
		 }
	 }
	 return C;
     }

    /** Make a one-dimensional column packed copy of the internal array.*/
    public double[] getColumnPackedCopy () {
        double[] vals = new double[m*n];
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                vals[i + j*m] = A[i][j];
            }
        }
        return vals;
    }
    
    /** Make a one-dimensional row packed copy of the internal array. */
    public double[] getRowPackedCopy() {
        double[] vals = new double[m*n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                vals[i*n + j] = A[i][j];
            }
        }
        return vals;
    }

    /** Get row dimension.*/
    public int getRowDimension () {
        return m;
    }
    /** Get column dimension.*/
    public int getColumnDimension () {
        return n;
    }
    /** Get a single element.*/
    public double get (int i, int j) {
        return A[i][j];
    }
    
    public Matrix getMatrix (int i0, int i1, int j0, int j1) {
        Matrix X = new Matrix(i1-i0+1,j1-j0+1);
        double[][] B = X.getArray();
        try {
            for (int i = i0; i <= i1; i++) {
                for (int j = j0; j <= j1; j++) {
                    B[i-i0][j-j0] = A[i][j];
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }

    
    /** Get a submatrix.
     @param r    Array of row indices.
     @param c    Array of column indices.
     @return     A(r(:),c(:))
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */
    
    public Matrix getMatrix (int[] r, int[] c) {
        Matrix X = new Matrix(r.length,c.length);
        double[][] B = X.getArray();
        try {
            for (int i = 0; i < r.length; i++) {
                for (int j = 0; j < c.length; j++) {
                    B[i][j] = A[r[i]][c[j]];
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }
    
   
    /** Get a submatrix.
     @param i0   Initial row index
     @param i1   Final row index
     @param c    Array of column indices.
     @return     A(i0:i1,c(:))
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */
    
    public Matrix getMatrix (int i0, int i1, int[] c) {
        Matrix X = new Matrix(i1-i0+1,c.length);
        double[][] B = X.getArray();
        try {
            for (int i = i0; i <= i1; i++) {
                for (int j = 0; j < c.length; j++) {
                    B[i-i0][j] = A[i][c[j]];
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }
    
    /** Get a submatrix.
     @param r    Array of row indices.
     @param j0   Initial column index
     @param j1   Final column index
     @return     A(r(:),j0:j1)
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */
    
    public Matrix getMatrix (int[] r, int j0, int j1) {
        Matrix X = new Matrix(r.length,j1-j0+1);
        double[][] B = X.getArray();
        try {
            for (int i = 0; i < r.length; i++) {
                for (int j = j0; j <= j1; j++) {
                    B[i][j-j0] = A[r[i]][j];
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }
    
    /** Set a single element.
     @param i    Row index.
     @param j    Column index.
     @param s    A(i,j).
     @exception  ArrayIndexOutOfBoundsException
     */
    
    public void set (int i, int j, double s) {
        A[i][j] = s;
    }
    
    /** Set a submatrix.
     @param i0   Initial row index
     @param i1   Final row index
     @param j0   Initial column index
     @param j1   Final column index
     @param X    A(i0:i1,j0:j1)
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */
    
    public void setMatrix (int i0, int i1, int j0, int j1, Matrix X) {
        try {
            for (int i = i0; i <= i1; i++) {
                for (int j = j0; j <= j1; j++) {
                    A[i][j] = X.get(i-i0,j-j0);
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
    }
    
    /** Set a submatrix.
     @param r    Array of row indices.
     @param c    Array of column indices.
     @param X    A(r(:),c(:))
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */
    
    public void setMatrix (int[] r, int[] c, Matrix X) {
        try {
            for (int i = 0; i < r.length; i++) {
                for (int j = 0; j < c.length; j++) {
                    A[r[i]][c[j]] = X.get(i,j);
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
    }
    
    /** Set a submatrix.
     @param r    Array of row indices.
     @param j0   Initial column index
     @param j1   Final column index
     @param X    A(r(:),j0:j1)
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */
    
    public void setMatrix (int[] r, int j0, int j1, Matrix X) {
        try {
            for (int i = 0; i < r.length; i++) {
                for (int j = j0; j <= j1; j++) {
                    A[r[i]][j] = X.get(i,j-j0);
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
    }
    
    /** Set a submatrix.
     @param i0   Initial row index
     @param i1   Final row index
     @param c    Array of column indices.
     @param X    A(i0:i1,c(:))
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */
    
    public void setMatrix (int i0, int i1, int[] c, Matrix X) {
        try {
            for (int i = i0; i <= i1; i++) {
                for (int j = 0; j < c.length; j++) {
                    A[i][c[j]] = X.get(i-i0,j);
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
    }
    
    /** Matrix transpose.
     @return    A'
     */
    
    public Matrix transpose () {
        Matrix X = new Matrix(n,m);
        double[][] C = X.getArray();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                C[j][i] = A[i][j];
            }
        }
        return X;
    }
    
    /** One norm
     @return    maximum column sum.
     */
    
    public double norm1 () {
        double f = 0;
        for (int j = 0; j < n; j++) {
            double s = 0;
            for (int i = 0; i < m; i++) {
                s += Math.abs(A[i][j]);
            }
            f = Math.max(f,s);
        }
        return f;
    }
    
    /** Two norm
     @return    maximum singular value.
     */
    
    public double norm2 () {
        return (new SingularValueDecomposition(this).norm2());
    }

    /** Infinity norm
     @return    maximum row sum.
     */
    
    public double normInf () {
        double f = 0;
        for (int i = 0; i < m; i++) {
            double s = 0;
            for (int j = 0; j < n; j++) {
                s += Math.abs(A[i][j]);
            }
            f = Math.max(f,s);
        }
        return f;
    }
    
    /** Frobenius norm */
    
    public double normF () {
        double f = 0;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                f = Maths.hypot(f,A[i][j]);
            }
        }
        return f;
    }
    
    /**  Unary minus
     @return    -A
     */
    
    public Matrix uminus ( ) {
        double[] a = this.getColumnPackedCopy();
        dscal(m * n, -1, a, 1);
        Matrix C = new Matrix(a, m);
        return C;
    }
    
    /* No function in Level 1-3 involves operation of addition, operation of matrix
     addition is written without using functions in Level 1-3        */
    /* C=A+B     */
    public Matrix plus (Matrix B){
        checkMatrixDimensions(B);
        double[] a = this.getColumnPackedCopy();
        double[] b = B.getColumnPackedCopy();
        daxpy(m * n, 1, b, 1, a, 1);
        Matrix C = new Matrix (a , m);
        return C;
    }
    
    
    /** A = A + B
     @param B    another matrix
     @return     A + B
     */
    
    public Matrix plusEquals (Matrix B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = A[i][j] + B.A[i][j];
            }
        }
        return this;
    }

    /* C = A - B     */
    public Matrix minus (Matrix B){
        checkMatrixDimensions(B);
        double[] a = this.getColumnPackedCopy();
        double[] b = B.getColumnPackedCopy();
        daxpy(m * n, -1.0, b, 1, a, 1);
        Matrix X = new Matrix (a , m);
        return X;
    }
    
    /** A = A - B
     @param B    another matrix
     @return     A - B
     */
    
    public Matrix minusEquals (Matrix B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = A[i][j] - B.A[i][j];
            }
        }
        return this;
    }

    
    /** Element-by-element multiplication, C = A.*B
     @param B    another matrix
     @return     A.*B
     */
    
    public Matrix arrayTimes (Matrix B) {
        checkMatrixDimensions(B);
        Matrix X = new Matrix(m,n);
        double[][] C = X.getArray();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                C[i][j] = A[i][j] * B.A[i][j];
            }
        }
        return X;
    }
    
    
    
    /** Element-by-element multiplication in place, A = A.*B
     @param B    another matrix
     @return     A.*B
     */
    
    public Matrix arrayTimesEquals (Matrix B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = A[i][j] * B.A[i][j];
            }
        }
        return this;
    }
    
    /** Element-by-element right division, C = A./B
     @param B    another matrix
     @return     A./B
     */
    
    public Matrix arrayRightDivide (Matrix B) {
        checkMatrixDimensions(B);
        Matrix X = new Matrix(m,n);
        double[][] C = X.getArray();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                C[i][j] = A[i][j] / B.A[i][j];
            }
        }
        return X;
    }
    
    /** Element-by-element right division in place, A = A./B
     @param B    another matrix
     @return     A./B
     */
    
    public Matrix arrayRightDivideEquals (Matrix B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = A[i][j] / B.A[i][j];
            }
        }
        return this;
    }
    
    /** Element-by-element left division, C = A.\B
     @param B    another matrix
     @return     A.\B
     */
    
    public Matrix arrayLeftDivide (Matrix B) {
        checkMatrixDimensions(B);
        Matrix X = new Matrix(m,n);
        double[][] C = X.getArray();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                C[i][j] = B.A[i][j] / A[i][j];
            }
        }
        return X;
    }
    
    /** Element-by-element left division in place, A = A.\B
     @param B    another matrix
     @return     A.\B
     */
    
    public Matrix arrayLeftDivideEquals (Matrix B) {
        checkMatrixDimensions(B);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = B.A[i][j] / A[i][j];
            }
        }
        return this;
    }
    
    
    /*Declaration
     B = alpha * A
     */
    
    public Matrix times (double alpha){
        double[] a = this.getColumnPackedCopy();
        dscal(m * n, alpha, a, 1);
        Matrix X = new Matrix(a, m);
        return X;
    }
    
    /** Multiply a matrix by a scalar in place, A = s*A
     @param s    scalar
     @return     replace A by s*A
     */
    
    public Matrix timesEquals (double s) {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = s*A[i][j];
            }
        }
        return this;
    }
    
    /*Declaration
     X = A * B
     */
    
    public Matrix times (Matrix B) {
        double[] a = this.getColumnPackedCopy();
        double[] b = B.getColumnPackedCopy();
        double[] c = new double[m * B.getColumnDimension()];
        dgemm(Matrix.LAYOUT.ColMajor, Matrix.TRANSPOSE.NoTrans, Matrix.TRANSPOSE.NoTrans,
              m, B.getColumnDimension(), n, 1, a, m, b, n, 0, c, m);
        Matrix X = new Matrix(c, m);
        return X;
    }
    
    
    
    /* Tentatively made by Diyang */
    
    
    /* Using daxpy, constants times a matrix plus a matrix
     C = A + alpha * B
     */
    public Matrix plus (Matrix B, double alpha){
        checkMatrixDimensions(B);
        double[] a = this.getColumnPackedCopy();
        double[] b = B.getColumnPackedCopy();
        daxpy(m * n, alpha, b, 1, a, 1);
        Matrix C = new Matrix (a , m);
        return C;
    }
    
    
    /* Using dtrmm, C=alpha*op(A)*B or  C=alpha*B*op(A)        */
    public  Matrix tritimes (Matrix B, double alpha){
        double[] a = this.getColumnPackedCopy();
        double[] b = B.getColumnPackedCopy();
        /* need to check if A is square matrix*/
        dtrmm(Matrix.LAYOUT.ColMajor, Matrix.SIDE.Left, Matrix.UPLO.Upper, Matrix.TRANSPOSE.NoTrans, Matrix.DIAG.NonUnit, B.getRowDimension(), B.getColumnDimension(), alpha, a, B.getRowDimension(), b, B.getColumnDimension());
        Matrix C = new Matrix (b, B.getRowDimension());
        return C;
    }
    
    
    /* Instead of using ddot, which is to compute the dot product of two vectors, computation
     of the Frobenius norm of a matrix     */
    public double getnorm (){
        double f = 0;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                f = Math.hypot(f,A[i][j]);
            }
        }
        return f;
    }
    
 
    /**   Below is matrix division. Don't know if it should be here because division is somewhat not useful compared to multiplication
    */
    

    
    /*****************************************************************/
    /****************************LAPACK*******************************/
    /*****************************************************************/
    
    /** LU Decomposition
     @return     LUDecomposition
     @see LUDecomposition
     */
    
    public LUDecomposition lu () {
        return new LUDecomposition(this);
    }
    
    /** QR Decomposition
     @return     QRDecomposition
     @see QRDecomposition
     */
    
    public QRDecomposition qr () {
        return new QRDecomposition(this);
    }
    
    /** Cholesky Decomposition
     @return     CholeskyDecomposition
     @see CholeskyDecomposition
     */
    
    public CholeskyDecomposition chol () {
        return new CholeskyDecomposition(this);
    }
    
    /** Singular Value Decomposition
     @return     SingularValueDecomposition
     @see SingularValueDecomposition
     */
    
    public SingularValueDecomposition svd () {
        return new SingularValueDecomposition(this);
    }
    
    /** Eigenvalue Decomposition
     @return     EigenvalueDecomposition
     @see EigenvalueDecomposition
     */
    
    public EigenvalueDecomposition eig () {
        return new EigenvalueDecomposition(this);
    }
    
    /** Solve A*X = B
     @param B    right hand side
     @return     solution if A is square, least squares solution otherwise
     */
    
    public Matrix solve (Matrix B) {
        return (m == n ? (new LUDecomposition(this)).solve(B) :
                (new QRDecomposition(this)).solve(B));
    }
    
    /** Solve X*A = B, which is also A'*X' = B'
     @param B    right hand side
     @return     solution if A is square, least squares solution otherwise.
     */
    
    public Matrix solveTranspose (Matrix B) {
        return transpose().solve(B.transpose());
    }
    
    /** Matrix inverse or pseudoinverse
     @return     inverse(A) if A is square, pseudoinverse otherwise.
     */
    
    public Matrix inverse () {
        return solve(identity(m,m));
    }
    
    /** Matrix determinant
     @return     determinant
     */
    
    public double det () {
        return new LUDecomposition(this).det();
    }
    
    /** Matrix rank
     @return     effective numerical rank, obtained from SVD.
     */
    
    public int rank () {
        return new SingularValueDecomposition(this).rank();
    }
    
    /** Matrix condition (2 norm)
     @return     ratio of largest to smallest singular value.
     */
    
    public double cond () {
        return new SingularValueDecomposition(this).cond();
    }
    
    /** Matrix trace.
     @return     sum of the diagonal elements.
     */
    
    public double trace () {
        double t = 0;
        for (int i = 0; i < Math.min(m,n); i++) {
            t += A[i][i];
        }
        return t;
    }
    
    /*****************************************************************/
    /****************************  end  *******************************/
    /*****************************************************************/
    
    
    
    /** Generate matrix with random elements
     @param m    Number of rows.
     @param n    Number of colums.
     @return     An m-by-n matrix with uniformly distributed random elements.
     */
    
    public static Matrix random (int m, int n) {
        Matrix A = new Matrix(m,n);
        double[][] X = A.getArray();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                X[i][j] = Math.random();
            }
        }
        return A;
    }

    
    /** Generate identity matrix
     @param m    Number of rows.
     @param n    Number of colums.
     @return     An m-by-n matrix with ones on the diagonal and zeros elsewhere.
     */
    
    public static Matrix identity (int m, int n) {
        Matrix A = new Matrix(m,n);
        double[][] X = A.getArray();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                X[i][j] = (i == j ? 1.0 : 0.0);
            }
        }
        return A;
    }
    
    
    /** Print the matrix to stdout.   Line the elements up in columns
     * with a Fortran-like 'Fw.d' style format.
     @param w    Column width.
     @param d    Number of digits after the decimal.
     */
    
    public void print (int w, int d) {
        print(new PrintWriter(System.out,true),w,d); }
    
    /** Print the matrix to the output stream.   Line the elements up in
     * columns with a Fortran-like 'Fw.d' style format.
     @param output Output stream.
     @param w      Column width.
     @param d      Number of digits after the decimal.
     */
    
    public void print (PrintWriter output, int w, int d) {
        DecimalFormat format = new DecimalFormat();
        format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
        format.setMinimumIntegerDigits(1);
        format.setMaximumFractionDigits(d);
        format.setMinimumFractionDigits(d);
        format.setGroupingUsed(false);
        print(output,format,w+2);
    }
    
    /** Print the matrix to stdout.  Line the elements up in columns.
     * Use the format object, and right justify within columns of width
     * characters.
     * Note that is the matrix is to be read back in, you probably will want
     * to use a NumberFormat that is set to US Locale.
     @param format A  Formatting object for individual elements.
     @param width     Field width for each column.
     @see java.text.DecimalFormat#setDecimalFormatSymbols
     */
    
    public void print (NumberFormat format, int width) {
        print(new PrintWriter(System.out,true),format,width); }
    
    // DecimalFormat is a little disappointing coming from Fortran or C's printf.
    // Since it doesn't pad on the left, the elements will come out different
    // widths.  Consequently, we'll pass the desired column width in as an
    // argument and do the extra padding ourselves.
    
    /** Print the matrix to the output stream.  Line the elements up in columns.
     * Use the format object, and right justify within columns of width
     * characters.
     * Note that is the matrix is to be read back in, you probably will want
     * to use a NumberFormat that is set to US Locale.
     @param output the output stream.
     @param format A formatting object to format the matrix elements
     @param width  Column width.
     @see java.text.DecimalFormat#setDecimalFormatSymbols
     */
    
    public void print (PrintWriter output, NumberFormat format, int width) {
        output.println();  // start on new line.
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                String s = format.format(A[i][j]); // format the number
                int padding = Math.max(1,width-s.length()); // At _least_ 1 space
                for (int k = 0; k < padding; k++)
                    output.print(' ');
                output.print(s);
            }
            output.println();
        }
        output.println();   // end with blank line.
    }
    
    /** Read a matrix from a stream.  The format is the same the print method,
     * so printed matrices can be read back in (provided they were printed using
     * US Locale).  Elements are separated by
     * whitespace, all the elements for each row appear on a single line,
     * the last row is followed by a blank line.
     @param input the input stream.
     */
    
    public static Matrix read (BufferedReader input) throws java.io.IOException {
        StreamTokenizer tokenizer= new StreamTokenizer(input);
        
        // Although StreamTokenizer will parse numbers, it doesn't recognize
        // scientific notation (E or D); however, Double.valueOf does.
        // The strategy here is to disable StreamTokenizer's number parsing.
        // We'll only get whitespace delimited words, EOL's and EOF's.
        // These words should all be numbers, for Double.valueOf to parse.
        
        tokenizer.resetSyntax();
        tokenizer.wordChars(0,255);
        tokenizer.whitespaceChars(0, ' ');
        tokenizer.eolIsSignificant(true);
        java.util.Vector<Double> vD = new java.util.Vector<Double>();
        
        // Ignore initial empty lines
        while (tokenizer.nextToken() == StreamTokenizer.TT_EOL);
        if (tokenizer.ttype == StreamTokenizer.TT_EOF)
            throw new java.io.IOException("Unexpected EOF on matrix read.");
        do {
            vD.addElement(Double.valueOf(tokenizer.sval)); // Read & store 1st row.
        } while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);
        
        int n = vD.size();  // Now we've got the number of columns!
        double row[] = new double[n];
        for (int j=0; j<n; j++)  // extract the elements of the 1st row.
            row[j]=vD.elementAt(j).doubleValue();
        java.util.Vector<double[]> v = new java.util.Vector<double[]>();
        v.addElement(row);  // Start storing rows instead of columns.
        while (tokenizer.nextToken() == StreamTokenizer.TT_WORD) {
            // While non-empty lines
            v.addElement(row = new double[n]);
            int j = 0;
            do {
                if (j >= n) throw new java.io.IOException
                    ("Row " + v.size() + " is too long.");
                row[j++] = Double.valueOf(tokenizer.sval).doubleValue();
            } while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);
            if (j < n) throw new java.io.IOException
                ("Row " + v.size() + " is too short.");
        }
        int m = v.size();  // Now we've got the number of rows.
        double[][] A = new double[m][];
        v.copyInto(A);  // copy the rows out of the vector
        return new Matrix(A);
    }

    
    /* ------------------------
     Private Methods
     * ------------------------ */
    
    /** Check if size(A) == size(B) **/
    
    private void checkMatrixDimensions (Matrix B) {
        if (B.m != m || B.n != n) {
            throw new IllegalArgumentException("Matrix dimensions must agree.");
        }
    }
    
    private static final long serialVersionUID = 1;
    
    
    public static class LAYOUT {
        private LAYOUT() {}
        /** row-major arrays */
        public final static int RowMajor= 101;
        /** column-major arrays */
        public final static int ColMajor= 102;
    }
    
    public static class TRANSPOSE {
        private TRANSPOSE() {}
        /** trans = 'N' */
        public final static int NoTrans = 111;
        /** trans = 'T' */
        public final static int Trans= 112;
        /** trans = 'C'*/
        public final static int ConjTrans= 113;
    }

    public static class UPLO {
        private UPLO() {}
        /** Upper triangular matrix */
        public final static int Upper = 121;
        /** Lower triangular matrix*/
        public final static int Lower = 122;
    }
    
    public static class DIAG {
        private DIAG() {}
        /** not assumed to be unit  */
        public final static int NonUnit = 131;
        /** assumed to be unit */
        public final static int Unit = 132;
    }
    
    public static class SIDE {
        private SIDE() {}
        /** B := alpha*op( A )*B */
        public final static int Left = 141;
        /** B := alpha*B*op( A ) */
        public final static int Right = 142;
    }
    
/*
 * ===========================================================================
 * Prototypes for level 1 BLAS functions (complex are recast as routines)
 * ===========================================================================
 */
public static native float  sdsdot(int N, float alpha, float[] X, int incX, float[] Y, int incY);

public static native double dsdot(int N, float[] X, int incX, float[] Y, int incY);

public static native float  sdot(int N, float[] X, int incX, float[] Y, int incY);

public static native double ddot(int N, double[] X, int incX, double[] Y, int incY);

/*
 * Functions having prefixes Z and C only
 */

public static native void   cdotu_sub(int N, float[] x, int incX, float[] y, int incY, float[] dotu);

public static native void   cdotc_sub(int N, float[] x, int incX, float[] y, int incY, float[] dotc);

public static native void   zdotu_sub(int N, double[] x, int incX, double[] y, int incY, double[] dotu);

public static native void   zdotc_sub(int N, double[] x, int incX, double[] y, int incY, double[] dotc);


/*
 * Functions having prefixes S D SC DZ
 */

public static native float  snrm2(int N, float[] X, int incX);

public static native float  sasum(int N, float[] X, int incX);

public static native double dnrm2(int N, double[] X, int incX);

public static native double dasum(int N, double[] X, int incX);

public static native float  scnrm2(int N, float[] X, int incX);

public static native float  scasum(int N, float[] X, int incX);

public static native double dznrm2(int N, double[] X, int incX);

public static native double dzasum(int N, double[] X, int incX);



/*
 * ===========================================================================
 * Prototypes for level 1 BLAS routines
 * ===========================================================================
 */

/* 
 * Routines with standard 4 prefixes (s, d, c, z)
 */


public static native void sswap(int N, float[] X, int incX, float[] Y, int incY);

public static native void scopy(int N, float[] X, int incX, float[] Y, int incY);

public static native void saxpy(int N, float alpha, float[] X, int incX, float[] Y, int incY);

public static native void dswap(int N, double[] X, int incX, double[] Y, int incY);

public static native void dcopy(int N, double[] X, int incX, double[] Y, int incY);

public static native void daxpy(int N, double alpha, double[] X, int incX, double[] Y, int incY);

public static native void cswap(int N, float[] x, int incX, float[] y, int incY);

public static native void ccopy(int N, float[] x, int incX, float[] y, int incY);

public static native void caxpy(int N, float[] alpha, float[] x, int incX, float[] y, int incY);

public static native void zswap(int N, double[] x, int incX, double[] y, int incY);

public static native void zcopy(int N, double[] x, int incX, double[] y, int incY);

public static native void zaxpy(int N, double[] alpha, double[] x, int incX, double[] y, int incY);


/* 
 * Routines with S and D prefix only
 */

public static native void srotg(float[] a, float[] b, float[] c, float[] s);

public static native void srotmg(float[] d1, float[] d2, float[] b1, float b2, float[] p);

public static native void srot(int N, float[] x, int incX, float[] y, int incY, float c, float s);

public static native void srotm(int N, float[] x, int incX, float[] y, int incY, float[] p);

public static native void drotg(double[] a, double[] b, double[] c, double[] s);

public static native void drotmg(double[] d1, double[] d2, double[] b1, double b2, double[] p);

public static native void drot(int N, double[] x, int incX, double[] y, int incY, double c, double  s);

public static native void drotm(int N, double[] x, int incX, double[] y, int incY, double[] p);


/* 
 * Routines with S D C Z CS and ZD prefixes
 */

public static native void sscal(int N, float alpha, float[] x, int incX);

public static native void dscal(int N, double alpha, double[] x, int incX);

public static native void cscal(int N, float[] alpha, float[] x, int incX);

public static native void zscal(int N, double[] alpha, double[] x, int incX);

public static native void csscal(int N, float alpha, float[] x, int incX);

public static native void zdscal(int N, double alpha, double[] x, int incX);

/*
 * ===========================================================================
 * Prototypes for level 2 BLAS
 * ===========================================================================
 */

/* 
 * Routines with standard 4 prefixes (S, D, C, Z)
 */

public static native void sgemv(int order, int TransA, int M, int N, float alpha, float[] a, int lda, float[] x, int incX, float beta, float[] y, int incY);

public static native void sgbmv(int order, int TransA, int M, int N, int KL, int KU, float alpha, float[] a, int lda, float[] x, int incX, float beta, float[] y, int incY);

public static native void strmv(int order, int Uplo, int TransA, int Diag, int N, float[] a, int lda, float[] x, int incX);

public static native void stbmv(int order, int Uplo, int TransA, int Diag, int N, int K, float[] a, int lda, float[] x, int incX);

public static native void stpmv(int order, int Uplo, int TransA, int Diag, int N, float[] Ap, float[] x, int incX);

public static native void strsv(int order, int Uplo, int TransA, int Diag, int N, float[] A, int lda, float[] X, int incX);

public static native void stbsv(int order, int Uplo, int TransA, int Diag, int N, int K, float[] A, int lda, float[] X, int incX);

public static native void stpsv(int order, int Uplo, int TransA, int Diag, int N, float[] Ap, float[] X, int incX);

public static native void dgemv(int order, int TransA, int M, int N, double alpha, double[] A, int lda, double[] X, int incX, double beta, double[] Y, int incY);

public static native void dgbmv(int order, int TransA, int M, int N, int KL, int KU, double alpha, double[] A, int lda, double[] X, int incX, double beta, double[] Y, int incY);

public static native void dtrmv(int order, int Uplo, int TransA, int Diag, int N, double[] A, int lda, double[] X, int incX);

public static native void dtbmv(int order, int Uplo, int TransA, int Diag, int N, int K, double[] A, int lda, double[] X, int incX);

public static native void dtpmv(int order, int Uplo, int TransA, int Diag, int N, double[] Ap, double[] X, int incX);

public static native void dtrsv(int order, int Uplo, int TransA, int Diag, int N, double[] A, int lda, double[] X, int incX);

public static native void dtbsv(int order, int Uplo, int TransA, int Diag, int N, int K, double[] A, int lda, double[] X, int incX);

public static native void dtpsv(int order, int Uplo, int TransA, int Diag, int N, double[] Ap, double[] X, int incX);

public static native void cgemv(int order, int TransA, int M, int N, float[] alpha, float[] A, int lda, float[] X, int incX, float[] beta, float[] Y, int incY);

public static native void cgbmv(int order, int TransA, int M, int N, int KL, int KU, float[] alpha, float[] A, int lda, float[] X, int incX, float[] beta, float[] Y, int incY);

public static native void ctrmv(int order, int Uplo, int TransA, int Diag, int N, float[] A, int lda, float[] X, int incX);

public static native void ctbmv(int order, int Uplo, int TransA, int Diag, int N, int K, float[] A, int lda, float[] X, int incX);

public static native void ctpmv(int order, int Uplo, int TransA, int Diag, int N, float[] Ap, float[] X, int incX);

public static native void ctrsv(int order, int Uplo, int TransA, int Diag, int N, float[] A, int lda, float[] X, int incX);

public static native void ctbsv(int order, int Uplo, int TransA, int Diag, int N, int K, float[] A, int lda, float[] X, int incX);

public static native void ctpsv(int order, int Uplo, int TransA, int Diag, int N, float[] Ap, float[] X, int incX);

public static native void zgemv(int order, int TransA, int M, int N, double[] alpha, double[] A, int lda,double[] X, int incX, double[] beta, double[] Y, int incY);

public static native void zgbmv(int order, int TransA, int M, int N, int KL, int KU, double[] alpha, double[] A, int lda, double[] X, int incX, double[] beta, double[] Y, int incY);

public static native void ztrmv(int order, int Uplo, int TransA, int Diag, int N, double[] A, int lda, double[] X, int incX);

public static native void ztbmv(int order, int Uplo, int TransA, int Diag, int N, int K, double[] A, int lda, double[] X, int incX);

public static native void ztpmv(int order, int Uplo, int TransA, int Diag, int N, double[] Ap, double[] X, int incX);

public static native void ztrsv(int order, int Uplo, int TransA, int Diag, int N, double[] A, int lda, double[] X, int incX);

public static native void ztbsv(int order, int Uplo, int TransA, int Diag, int N, int K, double[] A, int lda, double[] X, int incX);

public static native void ztpsv(int order, int Uplo, int TransA, int Diag, int N, double[] Ap, double[] X, int incX);


/* 
 * Routines with S and D prefixes only
 */

public static native void ssymv(int order, int Uplo, int N, float alpha, float[] A, int lda, float[] X, int incX, float beta, float[] Y, int incY);

public static native void ssbmv(int order, int Uplo, int N, int K, float alpha, float[] A, int lda, float[] X, int incX, float beta, float[] Y, int incY);

public static native void sspmv(int order, int Uplo, int N, float alpha, float[] Ap, float[] X, int incX, float beta, float[] Y, int incY);

public static native void sger(int order, int M, int N, float alpha, float[] X, int incX, float[] Y, int incY, float[] A, int lda);

public static native void ssyr(int order, int Uplo, int N, float alpha, float[] X, int incX, float[] A, int lda);

public static native void sspr(int order, int Uplo, int N, float alpha, float[] X, int incX, float[] Ap);

public static native void ssyr2(int order, int Uplo, int N, float alpha, float[] X, int incX, float[] Y, int incY, float[] A, int lda);

public static native void sspr2(int order, int Uplo, int N, float alpha, float[] X, int incX, float[] Y, int incY, float[] A);

public static native void dsymv(int order, int Uplo, int N, double alpha, double[] A, int lda, double[] X, int incX, double beta, double[] Y, int incY);

public static native void dsbmv(int order, int Uplo, int N, int K, double alpha, double[] A, int lda, double[] X, int incX, double beta, double[] Y, int incY);

public static native void dspmv(int order, int Uplo, int N, double alpha, double[] Ap, double[] X, int incX, double beta, double[] Y, int incY);

public static native void dger(int order, int M, int N, double alpha, double[] X, int incX, double[] Y, int incY, double[] A, int lda);

public static native void dsyr(int order, int Uplo, int N, double alpha, double[] X, int incX, double[] A, int lda);

public static native void dspr(int order, int Uplo, int N, double alpha, double[] X, int incX, double[] Ap);

public static native void dsyr2(int order, int Uplo, int N, double alpha, double[] X, int incX, double[] Y, int incY, double[] A, int lda);

public static native void dspr2(int order, int Uplo, int N, double alpha, double[] X, int incX, double[] Y, int incY, double[] A);


/* 
 * Routines with C and Z prefixes only
 */

public static native void chemv(int order, int Uplo, int N, float[] alpha, float[] A, int lda, float[] X, int incX,float[] beta, float[] Y, int incY);

public static native void chbmv(int order, int Uplo, int N, int K, float[] alpha, float[] A, int lda, float[] X, int incX,float[] beta, float[] Y, int incY);

public static native void chpmv(int order, int Uplo, int N, float[] alpha, float[] Ap, float[] X, int incX,float[] beta, float[] Y, int incY);

public static native void cgeru(int order, int M, int N, float[] alpha, float[] X, int incX, float[] Y, int incY, float[] A, int lda);

public static native void cgerc(int order, int M, int N, float[] alpha, float[] X, int incX, float[] Y, int incY, float[] A, int lda);

public static native void cher(int order, int Uplo, int N, float alpha, float[] X, int incX, float[] A, int lda);

public static native void chpr(int order, int Uplo, int N, float alpha, float[] X, int incX, float[] A);

public static native void cher2(int order, int Uplo, int N, float[] alpha, float[] X, int incX, float[] Y, int incY, float[] A, int lda);

public static native void chpr2(int order, int Uplo, int N, float[] alpha, float[] X, int incX, float[] Y, int incY, float[] Ap);

public static native void zhemv(int order, int Uplo, int N, double[] alpha, double[] A, int lda, double[] X, int incX, double[] beta, double[] Y, int incY);

public static native void zhbmv(int order, int Uplo, int N, int K, double[] alpha, double[] A, int lda, double[] X, int incX, double[] beta, double[] Y, int incY);

public static native void zhpmv(int order, int Uplo, int N, double[] alpha, double[] Ap, double[] X, int incX, double[] beta, double[] Y, int incY);

public static native void zgeru(int order, int M, int N, double[] alpha, double[] X, int incX, double[] Y, int incY, double[] A, int lda);

public static native void zgerc(int order, int M, int N, double[] alpha, double[] X, int incX, double[] Y, int incY, double[] A, int lda);

public static native void zher(int order, int Uplo, int N, double alpha, double[] X, int incX, double[] A, int lda);

public static native void zhpr(int order, int Uplo, int N, double alpha, double[] X, int incX, double[] A);

public static native void zher2(int order, int Uplo, int N, double[] alpha, double[] X, int incX, double[] Y, int incY, double[] A, int lda);

public static native void zhpr2(int order, int Uplo, int N, double[] alpha, double[] X, int incX, double[] Y, int incY, double[] Ap);

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/* 
 * Routines with standard 4 prefixes (S, D, C, Z)
 */

public static native void sgemm(int Order, int TransA, int TransB, int M, int N, int K, float alpha, float[] A, int lda, float[] B, int ldb, float beta, float[] C, int ldc);

public static native void ssymm(int Order, int Side, int Uplo, int M, int N, float alpha, float[] A, int lda, float[] B, int ldb, float beta, float[] C, int ldc);

public static native void ssyrk(int Order, int Uplo, int Trans, int N, int K, float alpha, float[] A, int lda, float beta, float[] C, int ldc);

public static native void ssyr2k(int Order, int Uplo, int Trans, int N, int K, float alpha, float[] A, int lda, float[] B, int ldb,  float beta, float[] C, int ldc);

public static native void strmm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, float alpha, float[] A, int lda, float[] B, int ldb);

public static native void strsm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, float alpha, float[] A, int lda, float[] B, int ldb);

public static native void dgemm(int Order, int TransA, int TransB, int M, int N, int K, double alpha, double[] A, int lda, double[] B, int ldb, double beta, double[] C, int ldc);

public static native void dsymm(int Order, int Side, int Uplo, int M, int N, double alpha, double[] A, int lda, double[] B, int ldb, double beta, double[] C, int ldc);

public static native void dsyrk(int Order, int Uplo, int Trans, int N, int K, double alpha, double[] A, int lda, double beta, double[] C, int ldc);

public static native void dsyr2k(int Order, int Uplo, int Trans, int N, int K, double alpha, double[] A, int lda, double[] B, int ldb, double beta, double[] C, int ldc);

public static native void dtrmm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, double alpha, double[] A, int lda, double[] B, int ldb);

public static native void dtrsm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, double alpha, double[] A, int lda, double[] B, int ldb);

public static native void cgemm(int Order, int TransA, int TransB, int M, int N, int K, float[] alpha, float[] A, int lda, float[] B, int ldb, float[] beta, float[] C, int ldc);

public static native void csymm(int Order, int Side, int Uplo, int M, int N, float[] alpha, float[] A, int lda, float[] B, int ldb, float[] beta, float[] C, int ldc);

public static native void csyrk(int Order, int Uplo, int Trans, int N, int K, float[] alpha, float[] A, int lda, float[] beta, float[] C, int ldc);

public static native void csyr2k(int Order, int Uplo, int Trans, int N, int K, float[] alpha, float[] A, int lda, float[] B, int ldb, float[] beta, float[] C, int ldc);

public static native void ctrmm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, float[] alpha, float[] A, int lda, float[] B, int ldb);

public static native void ctrsm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, float[] alpha, float[] A, int lda, float[] B, int ldb);

public static native void zgemm(int Order, int TransA, int TransB, int M, int N, int K, double[] alpha, double[] A, int lda, double[] B, int ldb, double[] beta, double[] C, int ldc);

public static native void zsymm(int Order, int Side, int Uplo, int M, int N, double[] alpha, double[] A, int lda, double[] B, int ldb, double[] beta, double[] C, int ldc);

public static native void zsyrk(int Order, int Uplo, int Trans, int N, int K, double[] alpha, double[] A, int lda, double[] beta, double[] C, int ldc);

public static native void zsyr2k(int Order, int Uplo, int Trans, int N, int K, double[] alpha, double[] A, int lda, double[] B, int ldb, double[] beta, double[] C, int ldc);

public static native void ztrmm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, double[] alpha, double[] A, int lda, double[] B, int ldb);

public static native void ztrsm(int Order, int Side, int Uplo, int TransA, int Diag, int M, int N, double[] alpha, double[] A, int lda, double[] B, int ldb);


/* 
 * Routines with prefixes C and Z only
 */

public static native void chemm(int Order, int Side, int Uplo, int M, int N, float[] alpha, float[] A, int lda, float[] B, int ldb, float[] beta, float[] C, int ldc);

public static native void cherk(int Order, int Uplo, int Trans, int N, int K, float alpha, float[] A, int lda, float beta, float[] C, int ldc);

public static native void cher2k(int Order, int Uplo, int Trans, int N, int K, float[] alpha, float[] A, int lda, float[] B, int ldb, float beta, float[] C, int ldc);

public static native void zhemm(int Order, int Side, int Uplo, int M, int N, double[] alpha, double[] A, int lda, double[] B, int ldb, double[] beta, double[] C, int ldc);

public static native void zherk(int Order, int Uplo, int Trans, int N, int K, double alpha, double[] A, int lda, double beta, double[] C, int ldc);

public static native void zher2k(int Order, int Uplo, int Trans, int N, int K, double[] alpha, double[] A, int lda, double[] B, int ldb, double beta, double[] C, int ldc);


       
    
}







