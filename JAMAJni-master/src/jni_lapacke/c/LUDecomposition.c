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
    
    printf(" before dgetrf ipiv: %d %d %d \n", ipivElems[0], ipivElems[1], ipivElems[2]);
    info = LAPACKE_dgetrf ((int) matrix_layout, (lapack_int) m, (lapack_int) n, aElems, (lapack_int) lda, ipivElems);
    printf(" after dgetrf ipiv: %d %d %d \n", ipivElems[0], ipivElems[1], ipivElems[2]);
    
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
    printf(" before dgetrs ipiv: %d %d %d \n", ipivElems[0], ipivElems[1], ipivElems[2]);
    info = LAPACKE_dgetrs ((int) matrix_layout, (char) trans, (lapack_int) n, (lapack_int) nrhs, aElems, (lapack_int) lda, ipivElems, bElems, (lapack_int) ldb);
    printf(" after dgetrs ipiv: %d %d %d \n", ipivElems[0], ipivElems[1], ipivElems[2]);
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, b, bElems, 0);
    (*env)-> ReleaseIntArrayElements (env, ipiv, ipivElems, JNI_ABORT);
    
    return info;
}






