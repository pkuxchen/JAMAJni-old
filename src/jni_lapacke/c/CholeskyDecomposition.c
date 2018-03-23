#include <jni.h>
#include <assert.h>
#include <lapacke.h>

/* Cholesky */

JNIEXPORT jint Java_JAMAJni_CholeskyDecomposition_dpotrf (JNIEnv *env, jclass klass, jint matrix_layout, jchar uplo, jint n, jdoubleArray a, jint lda){
    
    double *aElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    
    assert(aElems);
    
    info = LAPACKE_dpotrf((int) matrix_layout, (char) uplo, (lapack_int) n, aElems, (lapack_int) lda);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    
    return info;
    
}


JNIEXPORT jint Java_JAMAJni_CholeskyDecomposition_dpotri(JNIEnv *env, jclass klass, jint matrix_layout, jchar uplo, jint n, jdoubleArray a, jint lda){
    
    double *aElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    
    assert(aElems);
    
    info = LAPACKE_dpotri((int) matrix_layout, (char) uplo, (lapack_int) n, aElems, (lapack_int) lda);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    
    return info;
    
}

JNIEXPORT jint Java_JAMAJni_CholeskyDecomposition_dpotrs (JNIEnv *env, jclass klass, jint matrix_layout, jchar uplo, jint n, jint nrhs, jdoubleArray a, jint lda, jdoubleArray b, jint ldb){
    
    double *aElems, *bElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    bElems = (*env)-> GetDoubleArrayElements (env, b, NULL);
    
    assert(aElems && bElems);
    
    info = LAPACKE_dpotrs ((int) matrix_layout, (char) uplo, (lapack_int) n, (lapack_int) nrhs, aElems, (lapack_int) lda, bElems, (lapack_int) ldb);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, b, bElems, 0);
    
    return info;
    
}




