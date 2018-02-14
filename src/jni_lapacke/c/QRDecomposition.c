#include <jni.h>
#include <assert.h>
#include <lapacke.h>

extern void dtrmm_(char *side, char *uplo, char *transa, char *diag, int *m,
                   int *n, double *alpha, double *A, int *lda, double *B,
                   int *ldb);

JNIEXPORT jint Java_JAMAJni_QRDecomposition_dgeqrf (JNIEnv *env, jclass klass, jint matrix_layout,
                                                    jint m, jint n, jdoubleArray a, jint lda,
                                                    jdoubleArray tau){
    
    double *aElems, *tauElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    tauElems = (*env)-> GetDoubleArrayElements (env, tau, NULL);
    
    assert(aElems && tauElems);
    
    info = LAPACKE_dgeqrf ((int) matrix_layout, (lapack_int) m, (lapack_int) n, aElems, (lapack_int) lda, tauElems);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, tau, tauElems, 0);
    
    return info;
    
}

JNIEXPORT jint Java_JAMAJni_QRDecomposition_dorgqr (JNIEnv *env, jclass klass, jint matrix_layout,
                                                    jint m, jint n, jint k, jdoubleArray a, jint lda,
                                                    jdoubleArray tau){
    
    double *aElems, *tauElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    tauElems = (*env)-> GetDoubleArrayElements (env, tau, NULL);
    
    assert(aElems && tauElems);
    
    info = LAPACKE_dorgqr((int) matrix_layout, (lapack_int) m, (lapack_int) n, (lapack_int) k, aElems, (lapack_int) lda, tauElems);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, tau, tauElems, 0);
    
    return info;
    
}


JNIEXPORT jint Java_JAMAJni_QRDecomposition_dgeqp3  (JNIEnv *env, jclass klass, jint matrix_layout,
                                                     jint m, jint n, jdoubleArray a, jint lda,
                                                     jintArray jpvt, jdoubleArray tau){
    
    double *aElems, *tauElems;
    lapack_int *jpvtElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    tauElems = (*env)-> GetDoubleArrayElements (env, tau, NULL);
    jpvtElems = (*env)-> GetIntArrayElements (env, jpvt, NULL);
    
    assert(aElems && tauElems && jpvtElems);
    
    info = LAPACKE_dgeqp3 ((int) matrix_layout, (lapack_int) m, (lapack_int) n, aElems, (lapack_int) lda, jpvtElems, tauElems);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, tau, tauElems, 0);
    (*env)-> ReleaseIntArrayElements (env, jpvt, jpvtElems, 0);
    
}

JNIEXPORT jint Java_JAMAJni_QRDecomposition_dormqr
(JNIEnv *env, jclass klass, jint matrix_layout, jchar side, jchar trans,
 jint m, jint n, jint k, jdoubleArray a, jint lda, jdoubleArray tau,
 jdoubleArray c, jint ldc){
    
    double *aElems, *tauElems, *cElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    tauElems = (*env)-> GetDoubleArrayElements (env, tau, NULL);
    cElems = (*env)-> GetDoubleArrayElements (env, c, NULL);
    
    assert(aElems && tauElems && cElems);
    
    info = LAPACKE_dormqr ((int) matrix_layout, (char) side, (char) trans,
                           (lapack_int) m, (lapack_int) n, (lapack_int) k,
                           aElems, (lapack_int) lda, tauElems, cElems,
                           (lapack_int) ldc);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, tau, tauElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, c, cElems, JNI_ABORT);
    
    return info;
    
}

JNIEXPORT void Java_JAMAJni_QRDecomposition_dtrsm
(JNIEnv *env, jclass klass, jchar side, jchar uplo, jchar transa,
 jchar diag, jint m, jint n, jdouble alpha, jdoubleArray A, jint lda,
 jdoubleArray B, jint ldb){
    
    /*  DTRSM  solves one of the matrix equations
     
     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
     
     where alpha is a scalar, X and B are m by n matrices, A is a unit, or
     non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     
     op( A ) = A   or   op( A ) = A**T.
     
     The matrix X is overwritten on B.*/
    
    double *aElems, *bElems;
    
    aElems = (*env)-> GetDoubleArrayElements (env,A, NULL);
    bElems = (*env)-> GetDoubleArrayElements (env,B, NULL);
    
    assert(aElems && bElems);
    
    dtrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, aElems, &lda,
           bElems, &ldb);
    
    (*env)-> ReleaseDoubleArrayElements (env, B, bElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, A, aElems, JNI_ABORT);
    
}

