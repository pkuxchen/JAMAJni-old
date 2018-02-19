#include <jni.h>
#include <assert.h>
#include <lapacke.h>

extern void dtrmv_(char *uplo, char *trans, char *diag, int *n, double *A,
                   int *lda, double *x, int *incx);

#define jniRowMajor 101
#define jniColMajor 102

#define jniNoTrans 111
#define jniTrans   112
#define jniConjTrans    113

#define jniUpper   121
#define jniLower   122

#define jniNonUnit 131
#define jniUnit    132

#define jniLeft    141
#define jniRight   142


JNIEXPORT jint Java_JAMAJni_EigenvalueDecomposition_dgeev
(JNIEnv *env, jclass klass, jint matrix_layout, jchar jobvl, jchar jobvr,
 jint n, jdoubleArray a, jint lda, jdoubleArray wr, jdoubleArray wi,
 jdoubleArray vl, jint ldvl, jdoubleArray vr, jint ldvr){
    
    double *aElems, *wrElems, *wiElems, *vlElems, *vrElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    wrElems = (*env)-> GetDoubleArrayElements (env, wr, NULL);
    wiElems = (*env)-> GetDoubleArrayElements (env, wi, NULL);
    vlElems = (*env)-> GetDoubleArrayElements (env, vl, NULL);
    vrElems = (*env)-> GetDoubleArrayElements (env, vr, NULL);
    
    assert(aElems && wrElems && wiElems && vlElems && vrElems);
    
    info = LAPACKE_dgeev((int) matrix_layout, (char) jobvl, (char) jobvr, (lapack_int) n, aElems, (lapack_int) lda, wrElems, wiElems, vlElems, (lapack_int) ldvl, vrElems, (lapack_int) ldvr);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, vl, vlElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, vr, vrElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, wr, wrElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, wi, wiElems, 0);
    
    return info;
    
}


JNIEXPORT void JNICALL Java_JAMAJni_EigenvalueDecomposition_dsyev
(JNIEnv *env, jclass obj, jint layout, jchar jobz, jchar uplo, jint n,
 jdoubleArray ja, jint lda, jdoubleArray jw)
{
    jdouble *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    jdouble *w = (*env)-> GetDoubleArrayElements (env, jw, NULL);
    
    assert(a && w);
    
    LAPACKE_dsyev((int) layout, (char) jobz, (char) uplo, (lapack_int) n, a, (lapack_int) lda, w);
    
    (*env)-> ReleaseDoubleArrayElements (env, ja, a, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jw, w, 0);
}

JNIEXPORT void JNICALL Java_JAMAJni_EigenvalueDecomposition_dtrsm
(JNIEnv *env, jclass obj, jint layout, jchar jobz, jchar uplo, jint n,
 jdoubleArray ja, jint lda, jdoubleArray jw)
{
    jdouble *a = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    jdouble *w = (*env)-> GetDoubleArrayElements (env, jw, NULL);
    
    assert(a && w);
    
    LAPACKE_dsyev((int) layout, (char) jobz, (char) uplo, (lapack_int) n, a, (lapack_int) lda, w);
    
    (*env)-> ReleaseDoubleArrayElements (env, ja, a, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jw, w, 0);
}


