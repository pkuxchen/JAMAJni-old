#include <jni.h>
#include <assert.h>
#include <lapacke.h>

/* LU */

JNIEXPORT jint Java_JAMAJni_LUDecomposition_dgetrf (JNIEnv *env, jclass klass, jint matrix_layout, jint m, jint n, jdoubleArray a, jint lda, jintArray ipiv){
    
    double *aElems;
    lapack_int *ipivElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    ipivElems = (*env)-> GetIntArrayElements (env, ipiv, NULL);
    
    assert(aElems && ipivElems);
    
    info = LAPACKE_dgetrf ((int) matrix_layout, (lapack_int) m, (lapack_int) n, aElems, (lapack_int) lda, ipivElems);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseIntArrayElements (env, ipiv, ipivElems, 0);
    
    return info;
}

JNIEXPORT jint Java_JAMAJni_LUDecomposition_dgetrs (JNIEnv *env, jclass klass, jint matrix_layout, jchar trans, jint n, jint nrhs, jdoubleArray a, jint lda, jintArray ipiv, jdoubleArray b, jint ldb){
    
    double *bElems;
    double *aElems;
    lapack_int *ipivElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    bElems = (*env)-> GetDoubleArrayElements (env, b, NULL);
    ipivElems = (*env)-> GetIntArrayElements (env, ipiv, NULL);
    
    assert(aElems && bElems && ipivElems);

    info = LAPACKE_dgetrs ((int) matrix_layout, (char) trans, (lapack_int) n, (lapack_int) nrhs, aElems, (lapack_int) lda, ipivElems, bElems, (lapack_int) ldb);

    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, b, bElems, 0);
    (*env)-> ReleaseIntArrayElements (env, ipiv, ipivElems, JNI_ABORT);
    
    return info;
}






