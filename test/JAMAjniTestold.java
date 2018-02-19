import JAMAJni.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

public class JAMAjniTest {
    public static void main(String[] args) {
        Matrix A,B,C,Z,O,I,R,S,X,SUB,M,T,SQ,DEF,SOL, UU;
        
        //
        // Prepare the matrices and other parameters
        //
        System.out.println("###   Testfile for JAMAJni   ### ");
        int errorCount=0;
        int warningCount=0;
        double tmp, s;
        double[] columnwise = {1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.};
        double[] rowwise = {1.,4.,7.,10.,2.,5.,8.,11.,3.,6.,9.,12.};
        double[][] avals = {{1.,4.,7.,10.},{2.,5.,8.,11.},{3.,6.,9.,12.}};
        double[][] rankdef = avals;
        double[][] tvals =  {{1.,2.,3.},{4.,5.,6.},{7.,8.,9.},{10.,11.,12.}};
        double[][] subavals = {{5.,8.,11.},{6.,9.,12.}};
        double[][] rvals = {{1.,4.,7.},{2.,5.,8.,11.},{3.,6.,9.,12.}};
        double[][] pvals = {{4.,1.,1.},{1.,2.,3.},{1.,3.,6.}};
        double[][] ivals = {{1.,0.,0.,0.},{0.,1.,0.,0.},{0.,0.,1.,0.}};
        double[][] evals =
        {{0.,1.,0.,0.},{1.,0.,2.e-7,0.},{0.,-2.e-7,0.,1.},{0.,0.,1.,0.}};
        double[][] square = {{166.,188.,210.},{188.,214.,240.},{210.,240.,270.}};
        double[][] sqSolution = {{13.},{15.}};
        double[][] condmat = {{1.,3.},{7.,9.}};
        double[][] badeigs = {{0,0,0,0,0}, {0,0,0,0,1},{0,0,0,1,0},
            {1,1,0,0,1},{1,0,1,0,1}};
        int rows=3,cols=4;
        int invalidld=5;/* should trigger bad shape for construction with val */
        int raggedr=0; /* (raggedr,raggedc) should be out of bounds in ragged array */
        int raggedc=4;
        int validld=3; /* leading dimension of intended test Matrices */
        int nonconformld=4; /* leading dimension which is valid, but nonconforming */
        int ib=1,ie=2,jb=1,je=3; /* index ranges for sub Matrix */
        int[] rowindexset = {1,2};
        int[] badrowindexset = {1,3};
        int[] columnindexset = {1,2,3};
        int[] badcolumnindexset = {1,2,4};
        double columnsummax = 33.;
        double rowsummax = 30.;
        double sumofdiagonals = 15;
        double sumofsquares = 650;

        //
        // Set Option   (better to change into matrix.... not LUDecomposition)
        //
        int matrix_layout = LUDecomposition.LAYOUT.RowMajor;
        char Trans = LUDecomposition.TRANSPOSE.NoTrans;
        char uplo = LUDecomposition.UPLO.Upper;
        char Jobvl = LUDecomposition.JOBV.Compute;
        char Jobvr = LUDecomposition.JOBV.Compute;
        char Jobu = LUDecomposition.JOB.All;
        char Jobvt = LUDecomposition.JOB.All;
        char Jobz = LUDecomposition.JOB.Overwritten;
        int itype = LUDecomposition.ITYPE.first;
        int result;
        
        SUB = new Matrix(subavals);
        Z = new Matrix(rows, cols);
        /**
         LA methods:
         transpose
         times
         cond
         rank
         det
         trace
         norm1
         norm2
         normF
         normInf
         solve
         solveTranspose
         inverse
         chol
         eig
         lu
         qr
         svd
         **/
        
        print("\nTesting linear algebra methods...\n");
        A = new Matrix(columnwise,3);
        T = new Matrix(tvals);
        T = A.transpose();
        try {
            check(A.transpose(),T);
            try_success("transpose...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"transpose()...","transpose unsuccessful");
        }
        A.transpose();
        try {
            check(A.norm1(),columnsummax);
            try_success("norm1...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"norm1()...","incorrect norm calculation");
        }
        try {
            check(A.normInf(),rowsummax);
            try_success("normInf()...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"normInf()...","incorrect norm calculation");
        }
 /*       try {
            check(A.normF(),Math.sqrt(sumofsquares));
            try_success("normF...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"normF()...","incorrect norm calculation");
        }*/
        try {
            check(A.trace(),sumofdiagonals);
            try_success("trace()...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"trace()...","incorrect trace calculation");
        }
        try {
            check(A.getMatrix(0,A.getRowDimension()-1,0,A.getRowDimension()-1).det(),0.);
            try_success("det()...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"det()...","incorrect determinant calculation");
        }
        SQ = new Matrix(square);
        try {
            check(A.times(A.transpose()),SQ);
            try_success("times(Matrix)...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"times(Matrix)...","incorrect Matrix-Matrix product calculation");
        }
        try {
            check(A.times(0.),Z);
            try_success("times(double)...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"times(double)...","incorrect Matrix-scalar product calculation");
        }
        
        A = new Matrix(columnwise,4);
        QRDecomposition QR = A.qr();
        R = QR.getR();
        try {
            check(A,QR.getQ().times(R));
            try_success("QRDecomposition...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"QRDecomposition...","incorrect QR decomposition calculation");
        }
        
        SingularValueDecomposition SVD = A.svd();
        try {
            check(A,SVD.getU().times(SVD.getS().times(SVD.getV().transpose())));
            try_success("SingularValueDecomposition...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"SingularValueDecomposition...","incorrect singular value decomposition calculation");
        }
        DEF = new Matrix(rankdef);
        try {
            check(DEF.rank(), Math.min(DEF.getRowDimension(),DEF.getColumnDimension())-1);
            try_success("rank()...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"rank()...","incorrect rank calculation");
        }
        B = new Matrix(condmat);
        SVD = B.svd();
        double [] singularvalues = SVD.getSingularValues();
        try {
            check(B.cond(),singularvalues[0]/singularvalues[Math.min(B.getRowDimension(),B.getColumnDimension())-1]);
            try_success("cond()...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"cond()...", "incorrect condition number calculation");
        }
        int n = A.getColumnDimension();
        A = A.getMatrix(0, n-1, 0, n-1);
        A.set(0,0,0.);
        LUDecomposition LU = A.lu();
        int[] pv = {2, 0, 1}; /**pv has problem....*/
        try {
            check(A.getMatrix(pv, 0, n-1), LU.getL().times(LU.getU()));
            try_success("LUDecomposition...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"LUDecomposition...","incorrect LU decomposition calculation");
        }
        X = A.inverse();
        try {
            check(A.times(X),Matrix.identity(3,3));
            try_success("inverse()...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"inverse()...","incorrect inverse calculation");
        }
        O = new Matrix(SUB.getRowDimension(),1,1.0);
        SOL = new Matrix(sqSolution);
        SQ = SUB.getMatrix(0,SUB.getRowDimension()-1,0,SUB.getRowDimension()-1);
        try {
            check(SQ.solve(SOL),O);
            try_success("solve()...","");
        } catch ( java.lang.IllegalArgumentException e1 ) {
            errorCount = try_failure(errorCount,"solve()...",e1.getMessage());
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"solve()...",e.getMessage());
        }
        A = new Matrix(pvals);
        CholeskyDecomposition Chol = A.chol();
        Matrix L = Chol.getL();
        try {
            check(A, L.times(L.transpose()));
            try_success("CholeskyDecomposition...", "");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount, "CholeskyDecomposition...","incorrect Cholesky decomposition calculation");
        }
        X = Chol.solve(Matrix.identity(3,3));
        try {
            check(A.times(X), Matrix.identity(3,3));
            try_success("CholeskyDecomposition solve()...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"CholeskyDecomposition solve()...","incorrect Choleskydecomposition solve calculation");
        }
        EigenvalueDecomposition Eig = A.eig();
        Matrix D = Eig.getD();
        Matrix V = Eig.getV();
        try {
            check(A.times(V),V.times(D));
            try_success("EigenvalueDecomposition (symmetric)...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"EigenvalueDecomposition (symmetric)...","incorrect symmetric Eigenvalue decomposition calculation");
        }
        A = new Matrix(evals);
        Eig = A.eig();
        D = Eig.getD();
        V = Eig.getV();
        try {
            check(A.times(V),V.times(D));
            try_success("EigenvalueDecomposition (nonsymmetric)...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"EigenvalueDecomposition (nonsymmetric)...","incorrect nonsymmetric Eigenvalue decomposition calculation");
        }
        
        try {
            print("\nTesting Eigenvalue; If this hangs, we've failed\n");
            Matrix bA = new Matrix(badeigs);
            EigenvalueDecomposition bEig = bA.eig();
            try_success("EigenvalueDecomposition (hang)...","");
        } catch ( java.lang.RuntimeException e ) {
            errorCount = try_failure(errorCount,"EigenvalueDecomposition (hang)...",
                                     "incorrect termination");
        }
        
        
        print("\nTestMatrix completed.\n");
        print("Total errors reported: " + Integer.toString(errorCount) + "\n");
        print("Total warnings reported: " + Integer.toString(warningCount) + "\n");
        
        
        
        
        
    }
   
     //print function//
     /** Print the matrix X. */
      private static void printMatrix2(String prompt, Matrix A) {
	    System.out.println(prompt);
	    for (int i=0; i<A.getRowDimension(); i++) {
		    for (int j=0; j<A.getColumnDimension(); j++)
			    System.out.print("\t"+string(A.get(i,j)));
		    System.out.println();
	    }
    }

   
    /** Print the matrix X. */
     private static void printMatrix(String prompt, int layout, double[] X, int I, int J) {
        System.out.println(prompt);
        if (layout == LUDecomposition.LAYOUT.ColMajor) {
            for (int i=0; i<I; i++) {
                for (int j=0; j<J; j++)
                    System.out.print("\t" + string(X[j*I+i]));
                System.out.println();
            }
        }	
        else if (layout == LUDecomposition.LAYOUT.RowMajor){
            for (int i=0; i<I; i++) {
                for (int j=0; j<J; j++)
                    System.out.print("\t" + string(X[i*J+j]));
                System.out.println();
            }
        }
        else{System.out.println("** Illegal layout setting");}
    }
    
    private static void printIntArray(String prompt, int[] X, int L) {
        System.out.println(prompt);
        for (int i=0; i<L; i++) {
                System.out.print("\t" + string(X[i]));
        }
        System.out.println();
    }
    
    /** Shorter string for real number. */
    private static String string(double re) {
        String s="";
        if (re == (long)re)
            s += (long)re;
        else
            s += re;
        return s;
    }
    
    /*******************************************************/
    /*******************************************************/
    /*******************************************************/
    /** private utility routines **/
    
    /** Check magnitude of difference of scalars. **/
    
    private static void check(double x, double y) {
        double eps = Math.pow(2.0,-52.0);
        if (x == 0 & Math.abs(y) < 10*eps) return;
        if (y == 0 & Math.abs(x) < 10*eps) return;
        if (Math.abs(x-y) > 10*eps*Math.max(Math.abs(x),Math.abs(y))) {
            throw new RuntimeException("The difference x-y is too large: x = " + Double.toString(x) + "  y = " + Double.toString(y));
        }
    }
    
    /** Check norm of difference of "vectors". **/
    
    private static void check(double[] x, double[] y) {
        if (x.length == y.length ) {
            for (int i=0;i<x.length;i++) {
                check(x[i],y[i]);
            }
        } else {
            throw new RuntimeException("Attempt to compare vectors of different lengths");
        }
    }
    
    /** Check norm of difference of arrays. **/
    
    private static void check(double[][] x, double[][] y) {
        Matrix A = new Matrix(x);
        Matrix B = new Matrix(y);
        check(A,B);
    }
    
    /** Check norm of difference of Matrices. **/
    
    private static void check(Matrix X, Matrix Y) {
        double eps = Math.pow(2.0,-52.0);
        if (X.norm1() == 0. & Y.norm1() < 10*eps) return;
        if (Y.norm1() == 0. & X.norm1() < 10*eps) return;
        if (X.minus(Y).norm1() > 1000*eps*Math.max(X.norm1(),Y.norm1())) {
            throw new RuntimeException("The norm of (X-Y) is too large: " +  Double.toString(X.minus(Y).norm1()));
        }
    }
    
    /** Shorten spelling of print. **/
    
    private static void print (String s) {
        System.out.print(s);
    }
    
    /** Print appropriate messages for successful outcome try **/
    
    private static void try_success (String s,String e) {
        print(">    " + s + "success\n");
        if ( e != "" ) {
            print(">      Message: " + e + "\n");
        }
    }
    /** Print appropriate messages for unsuccessful outcome try **/
    
    private static int try_failure (int count, String s,String e) {
        print(">    " + s + "*** failure ***\n>      Message: " + e + "\n");
        return ++count;
    }
    
    /** Print appropriate messages for unsuccessful outcome try **/
    
    private static int try_warning (int count, String s,String e) {
        print(">    " + s + "*** warning ***\n>      Message: " + e + "\n");
        return ++count;
    }
    
    /** Print a row vector. **/
    
    private static void print(double[] x, int w, int d) {
        // Use format Fw.d for all elements.
        System.out.print("\n");
        new Matrix(x,1).print(w,d);
        print("\n");
    }
    
    /*******************************************************/
    /*******************************************************/
    /*******************************************************/
    
}



