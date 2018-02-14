#include <jni.h>
#include <assert.h>
#include <lapacke.h>

JNIEXPORT jint Java_JAMAJni_SingularValueDecomposition_dgesvd (JNIEnv *env, jclass klass, jint matrix_layout, jchar jobu, jchar jobvt, jint m, jint n, jdoubleArray a, jint lda, jdoubleArray s, jdoubleArray u, jint ldu, jdoubleArray vt, jint ldvt, jdoubleArray superb){
    
    //superb contains the unconverged superdiagonal elements of an upper bidiagonal matrix B whose diagonal is in S (not necessarily sorted).
    
    double *aElems, *sElems, *uElems, *vtElems, *superbElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    sElems = (*env)-> GetDoubleArrayElements (env, s, NULL);
    uElems = (*env)-> GetDoubleArrayElements (env, u, NULL);
    vtElems = (*env)-> GetDoubleArrayElements (env, vt, NULL);
    superbElems = (*env)-> GetDoubleArrayElements (env, superb, NULL);
    
    assert(aElems && sElems && uElems && vtElems && superbElems);
    
    info = LAPACKE_dgesvd((int) matrix_layout, (char) jobu, (char) jobvt, (lapack_int) m, (lapack_int) n, aElems, (lapack_int) lda, sElems, uElems, (lapack_int) ldu, vtElems, (lapack_int) ldvt, superbElems);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, s, sElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, u, uElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, vt, vtElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, superb, superbElems, 0);
    
    return info;
}

JNIEXPORT jint Java_JAMAJni_SingularValueDecomposition_dgesdd
(JNIEnv *env, jclass klass, jint matrix_layout, jchar jobz, jint m, jint n,
 jdoubleArray a, jint lda, jdoubleArray s, jdoubleArray u, jint ldu,
 jdoubleArray vt, jint ldvt){
    
    double *aElems, *sElems, *uElems, *vtElems;
    int info;
    
    aElems = (*env)-> GetDoubleArrayElements (env, a, NULL);
    sElems = (*env)-> GetDoubleArrayElements (env, s, NULL);
    uElems = (*env)-> GetDoubleArrayElements (env, u, NULL);
    vtElems = (*env)-> GetDoubleArrayElements (env, vt, NULL);
    
    assert(aElems && sElems && uElems && vtElems);
    
    info = LAPACKE_dgesdd((int) matrix_layout, (char) jobz, (lapack_int) m, (lapack_int) n, aElems, (lapack_int) lda, sElems, uElems, (lapack_int) ldu, vtElems, (lapack_int) ldvt);
    
    (*env)-> ReleaseDoubleArrayElements (env, a, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, s, sElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, u, uElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, vt, vtElems, 0);
    
    return info;
}

JNIEXPORT jint Java_JAMAJni_SingularValueDecomposition_dgeev (JNIEnv *env, jclass klass, jint matrix_layout, jchar jobvl, jchar jobvr, jint n, jdoubleArray a, jint lda, jdoubleArray wr, jdoubleArray wi, jdoubleArray vl, jint ldvl, jdoubleArray vr, jint ldvr){
    
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


