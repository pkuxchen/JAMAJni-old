#include <assert.h>
#include <jni.h>
#include "cblas.h"


// sdsdot: sdsdot = alpha+ inner product of vector x and y
JNIEXPORT jfloat JNICALL Java_JAMAJni_Matrix_sdsdot
  (JNIEnv *env, jclass obj, jint n, jfloat alpha, jfloatArray jx, jint incx, jfloatArray jy, jint incy)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && y);

  float result = cblas_sdsdot(n, alpha, x, incx, y, incy);
  
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,JNI_ABORT);

  return result;
}


// dsdot: dsdot = inner product of vector x and y
// x and y are single precision, the output is double precision
JNIEXPORT jdouble JNICALL Java_JAMAJni_Matrix_dsdot
  (JNIEnv *env, jclass obj, jint n, jfloatArray jx, jint incx, jfloatArray jy, jint incy)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && y);

  double result = cblas_dsdot(n, x, incx, y, incy);
  
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,JNI_ABORT);

  return result;
}

// sdot: sdot = inner product of vector x and y
JNIEXPORT jfloat JNICALL Java_JAMAJni_Matrix_sdot
  (JNIEnv *env, jclass obj, jint n, jfloatArray jx, jint incx, jfloatArray jy, jint incy)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && y);

  float result = cblas_sdot(n, x, incx, y, incy);
  
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,JNI_ABORT);

  return result;
}

//note: x, y, and the output are all double precision
JNIEXPORT jdouble JNICALL Java_JAMAJni_Matrix_ddot
  (JNIEnv *env, jclass obj, jint n, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(x && y);

  double result = cblas_ddot(n, x, incx, y, incy);
  
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,JNI_ABORT);

  return result;
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_cdotu_1sub
  (JNIEnv *env, jclass obj, jint n, jfloatArray jx, jint incx, jfloatArray jy, jint incy, jfloatArray jdotu)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  jfloat *dotu = (*env)->GetFloatArrayElements(env,jdotu,NULL);
  assert(x && y && dotu);
  
  cblas_cdotu_sub(n,x,incx,y,incy,dotu);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jdotu,dotu,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_cdotc_1sub
  (JNIEnv *env, jclass obj, jint n, jfloatArray jx, jint incx, jfloatArray jy, jint incy, jfloatArray jdotc)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  jfloat *dotc = (*env)->GetFloatArrayElements(env,jdotc,NULL);
  assert(x && y && dotc);

  cblas_cdotc_sub(n,x,incx,y,incy,dotc);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jdotc,dotc,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zdotu_1sub
  (JNIEnv *env, jclass obj, jint n, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy, jdoubleArray jdotu)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  jdouble *dotu = (*env)->GetDoubleArrayElements(env,jdotu,NULL);
  assert(x && y && dotu);

  cblas_zdotu_sub(n,x,incx,y,incy,dotu);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jdotu,dotu,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zdotc_1sub
  (JNIEnv *env, jclass obj, jint n, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy, jdoubleArray jdotc)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  jdouble *dotc = (*env)->GetDoubleArrayElements(env,jdotc,NULL);
  assert(x && y && dotc);
  
  cblas_zdotc_sub(n,x,incx,y,incy,dotc);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jdotc,dotc,0);
}


JNIEXPORT jfloat JNICALL Java_JAMAJni_Matrix_snrm2
  (JNIEnv *env, jclass obj, jint n, jfloatArray jx, jint incx)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  assert(x);
  
  float result = cblas_snrm2(n, x, incx);
  
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  
  return result;
}


JNIEXPORT jfloat JNICALL Java_JAMAJni_Matrix_sasum
  (JNIEnv *env, jclass obj, jint n, jfloatArray jx, jint incx)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  assert(x);

  float result = cblas_sasum(n, x, incx);
  
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  
  return result;
}


JNIEXPORT jdouble JNICALL Java_JAMAJni_Matrix_dnrm2
  (JNIEnv *env, jclass obj, jint n, jdoubleArray jx, jint incx)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(x);
  
  double result = cblas_dnrm2(n, x, incx);
  
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  
  return result;
}


JNIEXPORT jdouble JNICALL Java_JAMAJni_Matrix_dasum
  (JNIEnv *env, jclass obj, jint n, jdoubleArray jx, jint incx)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(x);
  
  double result = cblas_dasum(n, x, incx);
  
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  
  return result;
}


JNIEXPORT jfloat JNICALL Java_JAMAJni_Matrix_scnrm2
  (JNIEnv *env, jclass obj, jint n, jfloatArray jx, jint incx)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  assert(x);
  
  float result = cblas_scnrm2(n,x,incx);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  
  return result;
}


JNIEXPORT jfloat JNICALL Java_JAMAJni_Matrix_scasum
  (JNIEnv *env, jclass obj, jint n, jfloatArray jx, jint incx)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  assert(x);
  
  float result = cblas_scasum(n,x,incx);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  
  return result;
}


JNIEXPORT jdouble JNICALL Java_JAMAJni_Matrix_dznrm2
  (JNIEnv *env, jclass obj, jint n, jdoubleArray jx, jint incx)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(x);
  
  double result = cblas_dznrm2(n,x,incx);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  
  return result;
}


JNIEXPORT jdouble JNICALL Java_JAMAJni_Matrix_dzasum
  (JNIEnv *env, jclass obj, jint n, jdoubleArray jx, jint incx)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(x);
  
  double result = cblas_dzasum(n,x,incx);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  
  return result;
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_sswap
  (JNIEnv *env, jclass obj, jint n, jfloatArray jx, jint incx, jfloatArray jy, jint incy)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && y);

  cblas_sswap(n,x,incx,y,incy);

  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_scopy
  (JNIEnv *env, jclass obj, jint n, jfloatArray jx, jint incx, jfloatArray jy, jint incy)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && y);

  cblas_scopy(n,x,incx,y,incy);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_saxpy
  (JNIEnv *env, jclass obj, jint n, jfloat a, jfloatArray jx, jint incx, jfloatArray jy, jint incy)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);

  cblas_saxpy(n,a,x,incx,y,incy);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dswap
  (JNIEnv *env, jclass obj, jint n, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(x && y);

  cblas_dswap(n,x,incx,y,incy);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dcopy
  (JNIEnv *env, jclass obj, jint n, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(x && y);

  cblas_dcopy(n,x,incx,y,incy);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_daxpy
  (JNIEnv *env, jclass obj, jint n, jdouble a, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(x && y);

  cblas_daxpy(n,a,x,incx,y,incy);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_cswap
  (JNIEnv *env, jclass obj, jint n, jfloatArray jx, jint incx, jfloatArray jy, jint incy)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && y);
  
  cblas_cswap(n,x,incx,y,incy);

  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ccopy
  (JNIEnv *env, jclass obj, jint n, jfloatArray jx, jint incx, jfloatArray jy, jint incy)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && y);
  
  cblas_ccopy(n,x,incx,y,incy);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
}

JNIEXPORT void JNICALL Java_JAMAJni_Matrix_caxpy
  (JNIEnv *env, jclass obj, jint n, jfloatArray ja, jfloatArray jx, jint incx, jfloatArray jy, jint incy)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  assert(x && y && a);
  
  cblas_caxpy(n,a,x,incx,y,incy);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zswap
  (JNIEnv *env, jclass obj, jint n, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(x && y);

  cblas_zswap(n,x,incx,y,incy);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zcopy
  (JNIEnv *env, jclass obj, jint n, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(x && y);

  cblas_zcopy(n,x,incx,y,incy);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zaxpy
  (JNIEnv *env, jclass obj, jint n, jdoubleArray ja, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  assert(x && y && a);

  cblas_zaxpy(n,a,x,incx,y,incy);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_srotg
  (JNIEnv *env, jclass obj, jfloatArray ja, jfloatArray jb, jfloatArray jc, jfloatArray js)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *b = (*env)->GetFloatArrayElements(env,jb,NULL);
  jfloat *c = (*env)->GetFloatArrayElements(env,jc,NULL);
  jfloat *s = (*env)->GetFloatArrayElements(env,js,NULL);
  assert(a && b && c && s);
  
  cblas_srotg(a,b,c,s);

  (*env)->ReleaseFloatArrayElements(env,ja,a,0);
  (*env)->ReleaseFloatArrayElements(env,jb,b,0);
  (*env)->ReleaseFloatArrayElements(env,jc,c,0);
  (*env)->ReleaseFloatArrayElements(env,js,s,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_srotmg
  (JNIEnv *env, jclass obj, jfloatArray jsd1, jfloatArray jsd2, jfloatArray jsx1, jfloat sy1, jfloatArray jsparam)
{
  jfloat *sd1 = (*env)->GetFloatArrayElements(env,jsd1,NULL);
  jfloat *sd2 = (*env)->GetFloatArrayElements(env,jsd2,NULL);
  jfloat *sx1 = (*env)->GetFloatArrayElements(env,jsx1,NULL);
  jfloat *sparam = (*env)->GetFloatArrayElements(env,jsparam,NULL);
  assert(sd1 && sd2 && sx1 && sparam);
  
  cblas_srotmg(sd1,sd2,sx1,sy1,sparam);

  (*env)->ReleaseFloatArrayElements(env,jsd1,sd1,0);
  (*env)->ReleaseFloatArrayElements(env,jsd2,sd2,0);
  (*env)->ReleaseFloatArrayElements(env,jsx1,sx1,0);
  (*env)->ReleaseFloatArrayElements(env,jsparam,sparam,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_srot
  (JNIEnv *env, jclass obj, jint n, jfloatArray jx, jint incx, jfloatArray jy, jint incy, jfloat c, jfloat s)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && y);

  cblas_srot(n,x,incx,y,incy,c,s); 

  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_srotm
  (JNIEnv *env, jclass obj, jint n, jfloatArray jx, jint incx, jfloatArray jy, jint incy, jfloatArray jsparam)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  jfloat *sparam = (*env)->GetFloatArrayElements(env,jsparam,NULL);
  assert(x && y && sparam);

  cblas_srotm(n,x,incx,y,incy,sparam); 

  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
  (*env)->ReleaseFloatArrayElements(env,jsparam,sparam,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_drotg
  (JNIEnv *env, jclass obj, jdoubleArray ja, jdoubleArray jb, jdoubleArray jc, jdoubleArray js)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *b = (*env)->GetDoubleArrayElements(env,jb,NULL);
  jdouble *c = (*env)->GetDoubleArrayElements(env,jc,NULL);
  jdouble *s = (*env)->GetDoubleArrayElements(env,js,NULL);
  assert(a && b && c && s);
  
  cblas_drotg(a,b,c,s);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,0);
  (*env)->ReleaseDoubleArrayElements(env,jb,b,0);
  (*env)->ReleaseDoubleArrayElements(env,jc,c,0);
  (*env)->ReleaseDoubleArrayElements(env,js,s,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_drotmg
  (JNIEnv *env, jclass obj, jdoubleArray jdd1, jdoubleArray jdd2, jdoubleArray jdx1, jdouble dy1, jdoubleArray jdparam)
{
  jdouble *dd1 = (*env)->GetDoubleArrayElements(env,jdd1,NULL);
  jdouble *dd2 = (*env)->GetDoubleArrayElements(env,jdd2,NULL);
  jdouble *dx1 = (*env)->GetDoubleArrayElements(env,jdx1,NULL);
  jdouble *dparam = (*env)->GetDoubleArrayElements(env,jdparam,NULL);
  assert(dd1 && dd2 && dx1 && dparam);
  
  cblas_drotmg(dd1,dd2,dx1,dy1,dparam);

  (*env)->ReleaseDoubleArrayElements(env,jdd1,dd1,0);
  (*env)->ReleaseDoubleArrayElements(env,jdd2,dd2,0);
  (*env)->ReleaseDoubleArrayElements(env,jdx1,dx1,0);
  (*env)->ReleaseDoubleArrayElements(env,jdparam,dparam,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_drot
  (JNIEnv *env, jclass obj, jint n, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy, jdouble c, jdouble s)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(x && y);

  cblas_drot(n,x,incx,y,incy,c,s); 

  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_drotm
  (JNIEnv *env, jclass obj, jint n, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy, jdoubleArray jdparam)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  jdouble *dparam = (*env)->GetDoubleArrayElements(env,jdparam,NULL);
  assert(x && y && dparam);
  
  cblas_drotm(n,x,incx,y,incy,dparam); 

  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
  (*env)->ReleaseDoubleArrayElements(env,jdparam,dparam,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_sscal
  (JNIEnv *env, jclass obj, jint n, jfloat a, jfloatArray jx, jint incx)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  assert(x);
  
  cblas_sscal(n,a,x,incx);

  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dscal
  (JNIEnv *env, jclass obj, jint n, jdouble a, jdoubleArray jx, jint incx)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(x);
  
  cblas_dscal(n,a,x,incx);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_cscal
  (JNIEnv *env, jclass obj, jint n, jfloatArray ja, jfloatArray jx, jint incx)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  assert(x && a);

  cblas_cscal(n,a,x,incx);

  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zscal
  (JNIEnv *env, jclass obj, jint n, jdoubleArray ja, jdoubleArray jx, jint incx)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  assert(x && a);
  
  cblas_zscal(n,a,x,incx);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_csscal
  (JNIEnv *env, jclass obj, jint n, jfloat a, jfloatArray jx, jint incx)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  assert(x);

  cblas_csscal(n,a,x,incx);

  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zdscal
  (JNIEnv *env, jclass obj, jint n, jdouble a, jdoubleArray jx, jint incx)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(x);
  
  cblas_zdscal(n,a,x,incx);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_sgemv
  (JNIEnv *env, jclass obj, jint order, jint transA, jint m, jint n, jfloat alpha, jfloatArray ja, jint lda, jfloatArray jx, jint incx, jfloat beta, jfloatArray jy, jint incy)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && a && y);
  
  cblas_sgemv((CBLAS_ORDER) order,(CBLAS_TRANSPOSE) transA, m, n, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_sgbmv
  (JNIEnv *env, jclass obj, jint order, jint transA, jint m, jint n, jint kl, jint ku, jfloat alpha, jfloatArray ja, jint lda, jfloatArray jx, jint incx, jfloat beta, jfloatArray jy, jint incy)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && a && y);
  
  cblas_sgbmv((CBLAS_ORDER) order,(CBLAS_TRANSPOSE) transA, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_strmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jfloatArray ja, jint lda, jfloatArray jx, jint incx)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  assert(x && a);

  cblas_strmv((CBLAS_ORDER) order,(CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, lda, x, incx);

  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_stbmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jint k, jfloatArray ja, jint lda, jfloatArray jx, jint incx)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  assert(x && a);

  cblas_stbmv((CBLAS_ORDER) order,(CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, k, a, lda, x, incx);

  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_stpmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jfloatArray ja, jfloatArray jx, jint incx)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  assert(x && a);

  cblas_stpmv((CBLAS_ORDER) order,(CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, x, incx);

  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_strsv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jfloatArray ja, jint lda, jfloatArray jx, jint incx)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  assert(x && a);

  cblas_strsv((CBLAS_ORDER) order,(CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, lda, x, incx);

  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_stbsv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jint k, jfloatArray ja, jint lda, jfloatArray jx, jint incx)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  assert(x && a);

  cblas_stbsv((CBLAS_ORDER) order,(CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, k, a, lda, x, incx);

  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_stpsv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jfloatArray ja, jfloatArray jx, jint incx)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  assert(x && a);

  cblas_stpsv((CBLAS_ORDER) order,(CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, x, incx);

  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dgemv
  (JNIEnv *env, jclass obj, jint order, jint transA, jint m, jint n, jdouble alpha, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx, jdouble beta, jdoubleArray jy, jint incy)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(a && x && y);
  
  cblas_dgemv((CBLAS_ORDER) order, (CBLAS_TRANSPOSE) transA, m, n, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dgbmv
  (JNIEnv *env, jclass obj, jint order, jint transA, jint m, jint n, jint kl, jint ku, jdouble alpha, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx, jdouble beta, jdoubleArray jy, jint incy)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(a && x && y);
  
  cblas_dgbmv((CBLAS_ORDER) order, (CBLAS_TRANSPOSE) transA, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dtrmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_dtrmv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, lda, x, incx);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dtbmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jint k, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_dtbmv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, k, a, lda, x, incx);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dtpmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jdoubleArray ja, jdoubleArray jx, jint incx)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_dtpmv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, x, incx);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dtrsv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_dtrsv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, lda, x, incx);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dtbsv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jint k, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_dtbsv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, k, a, lda, x, incx);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dtpsv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jdoubleArray ja, jdoubleArray jx, jint incx)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_dtpsv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, x, incx);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_cgemv
  (JNIEnv *env, jclass obj, jint order, jint transA, jint m, jint n, jfloatArray jalpha, jfloatArray ja, jint lda, jfloatArray jx, jint incx, jfloatArray jbeta, jfloatArray jy, jint incy)
{
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *beta = (*env)->GetFloatArrayElements(env,jbeta,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(alpha && a && x && beta && y);
  
  cblas_cgemv((CBLAS_ORDER) order, (CBLAS_TRANSPOSE) transA, m, n, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);  
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_cgbmv
  (JNIEnv *env, jclass obj, jint order, jint transA, jint m, jint n, jint kl, jint ku, jfloatArray jalpha, jfloatArray ja, jint lda, jfloatArray jx, jint incx, jfloatArray jbeta, jfloatArray jy, jint incy)
{
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *beta = (*env)->GetFloatArrayElements(env,jbeta,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(alpha && a && x && beta && y);
  
  cblas_cgbmv((CBLAS_ORDER) order, (CBLAS_TRANSPOSE) transA, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);  
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ctrmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jfloatArray ja, jint lda, jfloatArray jx, jint incx)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_ctrmv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, lda, x, incx);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ctbmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jint k, jfloatArray ja, jint lda, jfloatArray jx, jint incx)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_ctbmv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, k, a, lda, x, incx);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ctpmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jfloatArray ja, jfloatArray jx, jint incx)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_ctpmv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, x, incx);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ctrsv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jfloatArray ja, jint lda, jfloatArray jx, jint incx)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_ctrsv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, lda, x, incx);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ctbsv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jint k, jfloatArray ja, jint lda, jfloatArray jx, jint incx)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_ctbsv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, k, a, lda, x, incx);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ctpsv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jfloatArray ja, jfloatArray jx, jint incx)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_ctpsv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, x, incx);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zgemv
  (JNIEnv *env, jclass obj, jint order, jint transA, jint m, jint n, jdoubleArray jalpha, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx, jdoubleArray jbeta, jdoubleArray jy, jint incy)
{
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *beta = (*env)->GetDoubleArrayElements(env,jbeta,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(alpha && a && x && beta && y);
  
  cblas_zgemv((CBLAS_ORDER) order, (CBLAS_TRANSPOSE) transA, m, n, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);  
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zgbmv
  (JNIEnv *env, jclass obj, jint order, jint transA, jint m, jint n, jint kl, jint ku, jdoubleArray jalpha, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx, jdoubleArray jbeta, jdoubleArray jy, jint incy)
{
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *beta = (*env)->GetDoubleArrayElements(env,jbeta,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(alpha && a && x && beta && y);
  
  cblas_zgbmv((CBLAS_ORDER) order, (CBLAS_TRANSPOSE) transA, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);  
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ztrmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_ztrmv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, lda, x, incx);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ztbmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jint k, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_ztbmv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, k, a, lda, x, incx);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ztpmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jdoubleArray ja, jdoubleArray jx, jint incx)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_ztpmv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, x, incx);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ztrsv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_ztrmv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, lda, x, incx);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ztbsv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jint k, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_ztbsv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, k, a, lda, x, incx);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ztpsv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint transA, jint diag, jint n, jdoubleArray ja, jdoubleArray jx, jint incx)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_ztpsv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, n, a, x, incx);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ssymv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jfloat alpha, jfloatArray ja, jint lda, jfloatArray jx, jint incx, jfloat beta, jfloatArray jy, jint incy)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && a && y);
  
  cblas_ssymv((CBLAS_ORDER) order,(CBLAS_UPLO) uplo, n, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ssbmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jint k, jfloat alpha, jfloatArray ja, jint lda, jfloatArray jx, jint incx, jfloat beta, jfloatArray jy, jint incy)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && a && y);
  
  cblas_ssbmv((CBLAS_ORDER) order,(CBLAS_UPLO) uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_sspmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jfloat alpha, jfloatArray ja, jfloatArray jx, jint incx, jfloat beta, jfloatArray jy, jint incy)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && a && y);

  cblas_sspmv((CBLAS_ORDER) order,(CBLAS_UPLO) uplo, n, alpha, a, x, incx, beta, y, incy);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_sger
  (JNIEnv *env, jclass obj, jint order, jint m, jint n, jfloat alpha, jfloatArray jx, jint incx, jfloatArray jy, jint incy, jfloatArray ja, jint lda)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && a && y);
  
  cblas_sger((CBLAS_ORDER) order, m, n, alpha, x, incx, y, incy, a, lda);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,ja,a,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ssyr
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jfloat alpha, jfloatArray jx, jint incx, jfloatArray ja, jint lda)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  assert(x && a);

  cblas_ssyr((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, a, lda);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,ja,a,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_sspr
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jfloat alpha, jfloatArray jx, jint incx, jfloatArray ja)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  assert(x && a);
  
  cblas_sspr((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, a);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,ja,a,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ssyr2
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jfloat alpha, jfloatArray jx, jint incx, jfloatArray jy, jint incy, jfloatArray ja, jint lda)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && a && y);
  
  cblas_ssyr2((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, y, incy, a, lda);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,ja,a,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_sspr2
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jfloat alpha, jfloatArray jx, jint incx, jfloatArray jy, jint incy, jfloatArray ja)
{
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(x && a && y);
  
  cblas_sspr2((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, y, incy, a);

  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,ja,a,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dsymv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jdouble alpha, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx, jdouble beta, jdoubleArray jy, jint incy)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(x && a && y);
  
  cblas_dsymv((CBLAS_ORDER) order,(CBLAS_UPLO) uplo, n, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dsbmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jint k, jdouble alpha, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx, jdouble beta, jdoubleArray jy, jint incy)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(x && a && y);
  
  cblas_dsbmv((CBLAS_ORDER) order,(CBLAS_UPLO) uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dspmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jdouble alpha, jdoubleArray ja, jdoubleArray jx, jint incx, jdouble beta, jdoubleArray jy, jint incy)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(x && a && y);
  
  cblas_dspmv((CBLAS_ORDER) order,(CBLAS_UPLO) uplo, n, alpha, a, x, incx, beta, y, incy);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dger
  (JNIEnv *env, jclass obj, jint order, jint m, jint n, jdouble alpha, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy, jdoubleArray ja, jint lda)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(x && a && y);
  
  cblas_dger((CBLAS_ORDER) order, m, n, alpha, x, incx, y, incy, a, lda);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dsyr
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jdouble alpha, jdoubleArray jx, jint incx, jdoubleArray ja, jint lda)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  assert(x && a);
  
  cblas_dsyr((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, a, lda);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,0);
}

JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dspr
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jdouble alpha, jdoubleArray jx, jint incx, jdoubleArray ja)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  assert(x && a);
  
  cblas_dspr((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, a);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dsyr2
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jdouble alpha, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy, jdoubleArray ja, jint lda)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(x && a && y);
  
  cblas_dsyr2((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, y, incy, a, lda);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dspr2
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jdouble alpha, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy, jdoubleArray ja)
{
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(x && a && y);
  
  cblas_dspr2((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, y, incy, a);

  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,0);
}



JNIEXPORT void JNICALL Java_JAMAJni_Matrix_chemv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jfloatArray jalpha, jfloatArray ja, jint lda, jfloatArray jx, jint incx, jfloatArray jbeta, jfloatArray jy, jint incy)
{
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *beta = (*env)->GetFloatArrayElements(env,jbeta,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(alpha && a && x && beta && y);
  
  cblas_chemv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_chbmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jint k, jfloatArray jalpha, jfloatArray ja, jint lda, jfloatArray jx, jint incx, jfloatArray jbeta, jfloatArray jy, jint incy)
{
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *beta = (*env)->GetFloatArrayElements(env,jbeta,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(alpha && a && x && beta && y);

  cblas_chbmv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_chpmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jfloatArray jalpha, jfloatArray ja, jfloatArray jx, jint incx, jfloatArray jbeta, jfloatArray jy, jint incy)
{
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *beta = (*env)->GetFloatArrayElements(env,jbeta,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(alpha && a && x && beta && y);
  
  cblas_chpmv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, a, x, incx, beta, y, incy);

  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_cgeru
  (JNIEnv *env, jclass obj, jint order, jint m, jint n, jfloatArray jalpha, jfloatArray jx, jint incx, jfloatArray jy, jint incy, jfloatArray ja, jint lda)
{
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(alpha && a && x && y);
  
  cblas_cgeru((CBLAS_ORDER) order, m, n, alpha, x, incx, y, incy, a, lda);

  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,ja,a,0);
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_cgerc
  (JNIEnv *env, jclass obj, jint order, jint m, jint n, jfloatArray jalpha, jfloatArray jx, jint incx, jfloatArray jy, jint incy, jfloatArray ja, jint lda)
{
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  assert(alpha && a && x && y);

  cblas_cgerc((CBLAS_ORDER) order, m, n, alpha, x, incx, y, incy, a, lda);

  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,ja,a,0);
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_cher
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jfloat alpha, jfloatArray jx, jint incx, jfloatArray ja, jint lda)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_cher((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, a, lda);

  (*env)->ReleaseFloatArrayElements(env,ja,a,0);
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_chpr
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jfloat alpha, jfloatArray jx, jint incx, jfloatArray ja)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  assert(a && x);
  
  cblas_chpr((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, a);

  (*env)->ReleaseFloatArrayElements(env,ja,a,0);
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_cher2
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jfloatArray jalpha, jfloatArray jx, jint incx, jfloatArray jy, jint incy, jfloatArray ja, jint lda)
{
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  assert(alpha && x && y && a);
  
  cblas_cher2((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, y, incy, a, lda);

  (*env)->ReleaseFloatArrayElements(env,ja,a,0);
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_chpr2
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jfloatArray jalpha, jfloatArray jx, jint incx, jfloatArray jy, jint incy, jfloatArray ja)
{
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  jfloat *x = (*env)->GetFloatArrayElements(env,jx,NULL);
  jfloat *y = (*env)->GetFloatArrayElements(env,jy,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  assert(alpha && x && y && a);
  
  cblas_chpr2((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, y, incy, a);

  (*env)->ReleaseFloatArrayElements(env,ja,a,0);
  (*env)->ReleaseFloatArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jy,y,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zhemv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jdoubleArray jalpha, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx, jdoubleArray jbeta, jdoubleArray jy, jint incy)
{
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *beta = (*env)->GetDoubleArrayElements(env,jbeta,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(alpha && a && x && beta && y);

  cblas_zhemv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zhbmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jint k, jdoubleArray jalpha, jdoubleArray ja, jint lda, jdoubleArray jx, jint incx, jdoubleArray jbeta, jdoubleArray jy, jint incy)
{
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *beta = (*env)->GetDoubleArrayElements(env,jbeta,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(alpha && a && x && beta && y);
  
  cblas_zhbmv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);

  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zhpmv
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jdoubleArray jalpha, jdoubleArray ja, jdoubleArray jx, jint incx, jdoubleArray jbeta, jdoubleArray jy, jint incy)
{
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *beta = (*env)->GetDoubleArrayElements(env,jbeta,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(alpha && a && x && beta && y);
  
  cblas_zhpmv((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, a, x, incx, beta, y, incy);

  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zgeru
  (JNIEnv *env, jclass obj, jint order, jint m, jint n, jdoubleArray jalpha, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy, jdoubleArray ja, jint lda)
{
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(alpha && a && x && y);

  cblas_zgeru((CBLAS_ORDER) order, m, n, alpha, x, incx, y, incy, a, lda);

  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,0);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zgerc
  (JNIEnv *env, jclass obj, jint order, jint m, jint n, jdoubleArray jalpha, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy, jdoubleArray ja, jint lda)
{
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  assert(alpha && a && x && y);
  
  cblas_zgerc((CBLAS_ORDER) order, m, n, alpha, x, incx, y, incy, a, lda);

  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,0);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zher
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jdouble alpha, jdoubleArray jx, jint incx, jdoubleArray ja, jint lda)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(a && x);

  cblas_zher((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, a, lda);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,0);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zhpr
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jdouble alpha, jdoubleArray jx, jint incx, jdoubleArray ja)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  assert(a && x);

  cblas_zhpr((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, a);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,0);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zher2
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jdoubleArray jalpha, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy, jdoubleArray ja, jint lda)
{
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  assert(alpha && x && y && a);
  
  cblas_zher2((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, y, incy, a, lda);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,0);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zhpr2
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint n, jdoubleArray jalpha, jdoubleArray jx, jint incx, jdoubleArray jy, jint incy, jdoubleArray ja)
{
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
  jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  assert(alpha && x && y && a);
  
  cblas_zhpr2((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, n, alpha, x, incx, y, incy, a);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,0);
  (*env)->ReleaseDoubleArrayElements(env,jx,x,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jy,y,JNI_ABORT);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_sgemm
  (JNIEnv *env, jclass obj, jint order, jint transA, jint transB, jint m, jint n, jint k, jfloat alpha, jfloatArray ja, jint lda, jfloatArray jb, jint ldb, jfloat beta, jfloatArray jc, jint ldc)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *b = (*env)->GetFloatArrayElements(env,jb,NULL);
  jfloat *c = (*env)->GetFloatArrayElements(env,jc,NULL);
  assert(a && b && c);
  
  cblas_sgemm((CBLAS_ORDER) order, (CBLAS_TRANSPOSE) transA, (CBLAS_TRANSPOSE) transB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ssymm
  (JNIEnv *env, jclass obj, jint order, jint side, jint uplo, jint m, jint n, jfloat alpha, jfloatArray ja, jint lda, jfloatArray jb, jint ldb, jfloat beta, jfloatArray jc, jint ldc)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *b = (*env)->GetFloatArrayElements(env,jb,NULL);
  jfloat *c = (*env)->GetFloatArrayElements(env,jc,NULL);
  assert(a && b && c);
  
  cblas_ssymm((CBLAS_ORDER) order, (CBLAS_SIDE) side, (CBLAS_UPLO) uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ssyrk
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint trans, jint n, jint k, jfloat alpha, jfloatArray ja, jint lda, jfloat beta, jfloatArray jc, jint ldc)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *c = (*env)->GetFloatArrayElements(env,jc,NULL);
  assert(a && c);
  
  cblas_ssyrk((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) trans, n, k, alpha, a, lda, beta, c, ldc);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ssyr2k
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint trans, jint n, jint k, jfloat alpha, jfloatArray ja, jint lda, jfloatArray jb, jint ldb, jfloat beta, jfloatArray jc, jint ldc)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *b = (*env)->GetFloatArrayElements(env,jb,NULL);
  jfloat *c = (*env)->GetFloatArrayElements(env,jc,NULL);
  assert(a && b && c);
  
  cblas_ssyr2k((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_strmm
  (JNIEnv *env, jclass obj, jint order, jint side, jint uplo, jint transA, jint diag, jint m, jint n, jfloat alpha, jfloatArray ja, jint lda, jfloatArray jb, jint ldb)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *b = (*env)->GetFloatArrayElements(env,jb,NULL);
  assert(a && b);
  
  cblas_strmm((CBLAS_ORDER) order, (CBLAS_SIDE) side, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, m, n, alpha, a, lda, b, ldb);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jb,b,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_strsm
  (JNIEnv *env, jclass obj, jint order, jint side, jint uplo, jint transA, jint diag, jint m, jint n, jfloat alpha, jfloatArray ja, jint lda, jfloatArray jb, jint ldb)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *b = (*env)->GetFloatArrayElements(env,jb,NULL);
  assert(a && b);
  
  cblas_strsm((CBLAS_ORDER) order, (CBLAS_SIDE) side, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, m, n, alpha, a, lda, b, ldb);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jb,b,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dgemm
  (JNIEnv *env, jclass obj, jint order, jint transA, jint transB, jint m, jint n, jint k, jdouble alpha, jdoubleArray ja, jint lda, jdoubleArray jb, jint ldb, jdouble beta, jdoubleArray jc, jint ldc)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *b = (*env)->GetDoubleArrayElements(env,jb,NULL);
  jdouble *c = (*env)->GetDoubleArrayElements(env,jc,NULL);
  assert(a && b && c);
  
  cblas_dgemm((CBLAS_ORDER) order, (CBLAS_TRANSPOSE) transA, (CBLAS_TRANSPOSE) transB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dsymm
  (JNIEnv *env, jclass obj, jint order, jint side, jint uplo, jint m, jint n, jdouble alpha, jdoubleArray ja, jint lda, jdoubleArray jb, jint ldb, jdouble beta, jdoubleArray jc, jint ldc)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *b = (*env)->GetDoubleArrayElements(env,jb,NULL);
  jdouble *c = (*env)->GetDoubleArrayElements(env,jc,NULL);
  assert(a && b && c);
  
  cblas_dsymm((CBLAS_ORDER) order, (CBLAS_SIDE) side, (CBLAS_UPLO) uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dsyrk
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint trans, jint n, jint k, jdouble alpha, jdoubleArray ja, jint lda, jdouble beta, jdoubleArray jc, jint ldc)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *c = (*env)->GetDoubleArrayElements(env,jc,NULL);
  assert(a && c);
  
  cblas_dsyrk((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) trans, n, k, alpha, a, lda, beta, c, ldc);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dsyr2k
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint trans, jint n, jint k, jdouble alpha, jdoubleArray ja, jint lda, jdoubleArray jb, jint ldb, jdouble beta, jdoubleArray jc, jint ldc)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *b = (*env)->GetDoubleArrayElements(env,jb,NULL);
  jdouble *c = (*env)->GetDoubleArrayElements(env,jc,NULL);
  assert(a && b && c);
  
  cblas_dsyr2k((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dtrmm
  (JNIEnv *env, jclass obj, jint order, jint side, jint uplo, jint transA, jint diag, jint m, jint n, jdouble alpha, jdoubleArray ja, jint lda, jdoubleArray jb, jint ldb)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *b = (*env)->GetDoubleArrayElements(env,jb,NULL);
  assert(a && b);
  
  cblas_dtrmm((CBLAS_ORDER) order, (CBLAS_SIDE) side, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, m, n, alpha, a, lda, b, ldb);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jb,b,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_dtrsm
  (JNIEnv *env, jclass obj, jint order, jint side, jint uplo, jint transA, jint diag, jint m, jint n, jdouble alpha, jdoubleArray ja, jint lda, jdoubleArray jb, jint ldb)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *b = (*env)->GetDoubleArrayElements(env,jb,NULL);
  assert(a && b);
  
  cblas_dtrsm((CBLAS_ORDER) order, (CBLAS_SIDE) side, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, m, n, alpha, a, lda, b, ldb);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jb,b,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_cgemm
  (JNIEnv *env, jclass obj, jint order, jint transA, jint transB, jint m, jint n, jint k, jfloatArray jalpha, jfloatArray ja, jint lda, jfloatArray jb, jint ldb, jfloatArray jbeta, jfloatArray jc, jint ldc)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *b = (*env)->GetFloatArrayElements(env,jb,NULL);
  jfloat *c = (*env)->GetFloatArrayElements(env,jc,NULL);
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  jfloat *beta = (*env)->GetFloatArrayElements(env,jbeta,NULL);
  assert(a && b && c && alpha && beta);
  
  cblas_cgemm((CBLAS_ORDER) order, (CBLAS_TRANSPOSE) transA, (CBLAS_TRANSPOSE) transB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_csymm
  (JNIEnv *env, jclass obj, jint order, jint side, jint uplo, jint m, jint n, jfloatArray jalpha, jfloatArray ja, jint lda, jfloatArray jb, jint ldb, jfloatArray jbeta, jfloatArray jc, jint ldc)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *b = (*env)->GetFloatArrayElements(env,jb,NULL);
  jfloat *c = (*env)->GetFloatArrayElements(env,jc,NULL);
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  jfloat *beta = (*env)->GetFloatArrayElements(env,jbeta,NULL);
  assert(a && b && c && alpha && beta);
  
  cblas_csymm((CBLAS_ORDER) order, (CBLAS_SIDE) side, (CBLAS_UPLO) uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_csyrk
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint trans, jint n, jint k, jfloatArray jalpha, jfloatArray ja, jint lda, jfloatArray jbeta, jfloatArray jc, jint ldc)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *c = (*env)->GetFloatArrayElements(env,jc,NULL);
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  jfloat *beta = (*env)->GetFloatArrayElements(env,jbeta,NULL);
  assert(a && c && alpha && beta);
  
  cblas_csyrk((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) trans, n, k, alpha, a, lda, beta, c, ldc);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_csyr2k
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint trans, jint n, jint k, jfloatArray jalpha, jfloatArray ja, jint lda, jfloatArray jb, jint ldb, jfloatArray jbeta, jfloatArray jc, jint ldc)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *b = (*env)->GetFloatArrayElements(env,jb,NULL);
  jfloat *c = (*env)->GetFloatArrayElements(env,jc,NULL);
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  jfloat *beta = (*env)->GetFloatArrayElements(env,jbeta,NULL);
  assert(a && b && c && alpha && beta);
  
  cblas_csyr2k((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ctrmm
  (JNIEnv *env, jclass obj, jint order, jint side, jint uplo, jint transA, jint diag, jint m, jint n, jfloatArray jalpha, jfloatArray ja, jint lda, jfloatArray jb, jint ldb)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *b = (*env)->GetFloatArrayElements(env,jb,NULL);
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  assert(a && b && alpha);
  
  cblas_ctrmm((CBLAS_ORDER) order, (CBLAS_SIDE) side, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, m, n, alpha, a, lda, b, ldb);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jb,b,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ctrsm
  (JNIEnv *env, jclass obj, jint order, jint side, jint uplo, jint transA, jint diag, jint m, jint n, jfloatArray jalpha, jfloatArray ja, jint lda, jfloatArray jb, jint ldb)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *b = (*env)->GetFloatArrayElements(env,jb,NULL);
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  assert(a && b && alpha);
  
  cblas_ctrsm((CBLAS_ORDER) order, (CBLAS_SIDE) side, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, m, n, alpha, a, lda, b, ldb);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jb,b,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zgemm
  (JNIEnv *env, jclass obj, jint order, jint transA, jint transB, jint m, jint n, jint k, jdoubleArray jalpha, jdoubleArray ja, jint lda, jdoubleArray jb, jint ldb, jdoubleArray jbeta, jdoubleArray jc, jint ldc)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *b = (*env)->GetDoubleArrayElements(env,jb,NULL);
  jdouble *c = (*env)->GetDoubleArrayElements(env,jc,NULL);
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  jdouble *beta = (*env)->GetDoubleArrayElements(env,jbeta,NULL);
  assert(a && b && c && alpha && beta);
  
  cblas_zgemm((CBLAS_ORDER) order, (CBLAS_TRANSPOSE) transA, (CBLAS_TRANSPOSE) transB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zsymm
  (JNIEnv *env, jclass obj, jint order, jint side, jint uplo, jint m, jint n, jdoubleArray jalpha, jdoubleArray ja, jint lda, jdoubleArray jb, jint ldb, jdoubleArray jbeta, jdoubleArray jc, jint ldc)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *b = (*env)->GetDoubleArrayElements(env,jb,NULL);
  jdouble *c = (*env)->GetDoubleArrayElements(env,jc,NULL);
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  jdouble *beta = (*env)->GetDoubleArrayElements(env,jbeta,NULL);
  assert(a && b && c && alpha && beta);
  
  cblas_zsymm((CBLAS_ORDER) order, (CBLAS_SIDE) side, (CBLAS_UPLO) uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zsyrk
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint trans, jint n, jint k, jdoubleArray jalpha, jdoubleArray ja, jint lda, jdoubleArray jbeta, jdoubleArray jc, jint ldc)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *c = (*env)->GetDoubleArrayElements(env,jc,NULL);
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  jdouble *beta = (*env)->GetDoubleArrayElements(env,jbeta,NULL);
  assert(a && c && alpha && beta);
  
  cblas_zsyrk((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) trans, n, k, alpha, a, lda, beta, c, ldc);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zsyr2k
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint trans, jint n, jint k, jdoubleArray jalpha, jdoubleArray ja, jint lda, jdoubleArray jb, jint ldb, jdoubleArray jbeta, jdoubleArray jc, jint ldc)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *b = (*env)->GetDoubleArrayElements(env,jb,NULL);
  jdouble *c = (*env)->GetDoubleArrayElements(env,jc,NULL);
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  jdouble *beta = (*env)->GetDoubleArrayElements(env,jbeta,NULL);
  assert(a && b && c && alpha && beta);
  
  cblas_zsyr2k((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ztrmm
  (JNIEnv *env, jclass obj, jint order, jint side, jint uplo, jint transA, jint diag, jint m, jint n, jdoubleArray jalpha, jdoubleArray ja, jint lda, jdoubleArray jb, jint ldb)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *b = (*env)->GetDoubleArrayElements(env,jb,NULL);
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  assert(a && b && alpha);
  
  cblas_ztrmm((CBLAS_ORDER) order, (CBLAS_SIDE) side, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, m, n, alpha, a, lda, b, ldb);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jb,b,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_ztrsm
  (JNIEnv *env, jclass obj, jint order, jint side, jint uplo, jint transA, jint diag, jint m, jint n, jdoubleArray jalpha, jdoubleArray ja, jint lda, jdoubleArray jb, jint ldb)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *b = (*env)->GetDoubleArrayElements(env,jb,NULL);
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  assert(a && b && alpha);
  
  cblas_ztrsm((CBLAS_ORDER) order, (CBLAS_SIDE) side, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA, (CBLAS_DIAG) diag, m, n, alpha, a, lda, b, ldb);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jb,b,0);
}



JNIEXPORT void JNICALL Java_JAMAJni_Matrix_chemm
  (JNIEnv *env, jclass obj, jint order, jint side, jint uplo, jint m, jint n, jfloatArray jalpha, jfloatArray ja, jint lda, jfloatArray jb, jint ldb, jfloatArray jbeta, jfloatArray jc, jint ldc)
{
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *b = (*env)->GetFloatArrayElements(env,jb,NULL);
  jfloat *beta = (*env)->GetFloatArrayElements(env,jbeta,NULL);
  jfloat *c = (*env)->GetFloatArrayElements(env,jc,NULL);
  assert(alpha && a && b && beta && c);
  
  cblas_chemm((CBLAS_ORDER) order, (CBLAS_SIDE)side, (CBLAS_UPLO) uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_cherk
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint trans, jint n, jint k, jfloat alpha, jfloatArray ja, jint lda, jfloat beta, jfloatArray jc, jint ldc)
{
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *c = (*env)->GetFloatArrayElements(env,jc,NULL);
  assert(a && c);

  cblas_cherk((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) trans, n, k, alpha, a, lda, beta, c, ldc);

  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_cher2k
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint trans, jint n, jint k, jfloatArray jalpha, jfloatArray ja, jint lda, jfloatArray jb, jint ldb, jfloat beta, jfloatArray jc, jint ldc)
{
  jfloat *alpha = (*env)->GetFloatArrayElements(env,jalpha,NULL);
  jfloat *a = (*env)->GetFloatArrayElements(env,ja,NULL);
  jfloat *b = (*env)->GetFloatArrayElements(env,jb,NULL);
  jfloat *c = (*env)->GetFloatArrayElements(env,jc,NULL);
  assert(alpha && a && b && c);
  
  cblas_cher2k((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseFloatArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseFloatArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zhemm
  (JNIEnv *env, jclass obj, jint order, jint side, jint uplo, jint m, jint n, jdoubleArray jalpha, jdoubleArray ja, jint lda, jdoubleArray jb, jint ldb, jdoubleArray jbeta, jdoubleArray jc, jint ldc)
{
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *b = (*env)->GetDoubleArrayElements(env,jb,NULL);
  jdouble *beta = (*env)->GetDoubleArrayElements(env,jbeta,NULL);
  jdouble *c = (*env)->GetDoubleArrayElements(env,jc,NULL);
  assert(alpha && a && b && beta && c);

  cblas_zhemm((CBLAS_ORDER) order, (CBLAS_SIDE)side, (CBLAS_UPLO) uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jbeta,beta,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zherk
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint trans, jint n, jint k, jdouble alpha, jdoubleArray ja, jint lda, jdouble beta, jdoubleArray jc, jint ldc)
{
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *c = (*env)->GetDoubleArrayElements(env,jc,NULL);
  assert(a && c);
  
  cblas_zherk((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) trans, n, k, alpha, a, lda, beta, c, ldc);

  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jc,c,0);
}


JNIEXPORT void JNICALL Java_JAMAJni_Matrix_zher2k
  (JNIEnv *env, jclass obj, jint order, jint uplo, jint trans, jint n, jint k, jdoubleArray jalpha, jdoubleArray ja, jint lda, jdoubleArray jb, jint ldb, jdouble beta, jdoubleArray jc, jint ldc)
{
  jdouble *alpha = (*env)->GetDoubleArrayElements(env,jalpha,NULL);
  jdouble *a = (*env)->GetDoubleArrayElements(env,ja,NULL);
  jdouble *b = (*env)->GetDoubleArrayElements(env,jb,NULL);
  jdouble *c = (*env)->GetDoubleArrayElements(env,jc,NULL);
  assert(alpha && a && b && c);
  
  cblas_zher2k((CBLAS_ORDER) order, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

  (*env)->ReleaseDoubleArrayElements(env,jalpha,alpha,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,ja,a,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jb,b,JNI_ABORT);
  (*env)->ReleaseDoubleArrayElements(env,jc,c,0);
}
