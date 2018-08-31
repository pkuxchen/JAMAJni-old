DEST = .
CC = gcc
CFLAGS = $(INCLUDE) $(LIB) 

JAVA_HOME = ..

JNI_LAPACK_SRC = $(JAVA_HOME)/src/jni_lapack
JNI_CBLAS_SRC = $(JAVA_HOME)/src/jni_cblas
JNI_LAPACKE_SRC = $(JAVA_HOME)/src/jni_lapacke
JNI_UTIL = $(JAVA_HOME)/src/util

INCLUDE = -I$(JAVA_HOME)/include

LIB = -L/usr/local/lib #-L/usr/local/lib64 -L/usr/lib 

C_SRC_LAPACK1 = $(JNI_LAPACKE_SRC)/c/LUDecomposition.c 
C_SRC_LAPACK2 = $(JNI_LAPACKE_SRC)/c/CholeskyDecomposition.c
C_SRC_LAPACK3 = $(JNI_LAPACKE_SRC)/c/QRDecomposition.c
C_SRC_LAPACK4 = $(JNI_LAPACKE_SRC)/c/EigenvalueDecomposition.c
C_SRC_LAPACK5 = $(JNI_LAPACKE_SRC)/c/SingularValueDecomposition.c
C_SRC_CBLAS = $(JNI_CBLAS_SRC)/c/Matrix.c
DYLIB_LAPACK1 = $(JAVA_HOME)/lib/liblapacke_LUDecomposition.so
DYLIB_LAPACK2 = $(JAVA_HOME)/lib/liblapacke_CholeskyDecomposition.so
DYLIB_LAPACK3 = $(JAVA_HOME)/lib/liblapacke_QRDecomposition.so
DYLIB_LAPACK4 = $(JAVA_HOME)/lib/liblapacke_EigenvalueDecomposition.so
DYLIB_LAPACK5 = $(JAVA_HOME)/lib/liblapacke_SingularValueDecomposition.so
DYLIB_CBLAS = $(JAVA_HOME)/lib/libcblas_Matrix.so
	

# For package JAMAJni
JC = javac
JFLAGS = -d
JAVA_SRC = $(JNI_UTIL)/Maths.java $(JNI_LAPACKE_SRC)/java/LUDecomposition.java $(JNI_LAPACKE_SRC)/java/CholeskyDecomposition.java $(JNI_LAPACKE_SRC)/java/QRDecomposition.java $(JNI_LAPACKE_SRC)/java/EigenvalueDecomposition.java $(JNI_LAPACKE_SRC)/java/SingularValueDecomposition.java $(JNI_CBLAS_SRC)/java/Matrix.java 
PACKAGE = $(JAVA_HOME)/class/*


all: $(PACKAGE) $(DYLIB_LAPACK1) $(DYLIB_LAPACK2) $(DYLIB_LAPACK3) $(DYLIB_LAPACK4) $(DYLIB_LAPACK5) $(DYLIB_CBLAS)

$(PACKAGE): $(JAVA_SRC)
	@mkdir -p $(@D)
	$(JC) $(JFLAGS) $(JAVA_HOME)/class $(JAVA_SRC)

$(DYLIB_LAPACK1): $(C_SRC_LAPACK1)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_LAPACK1) -llapacke -lcblas -o $(DYLIB_LAPACK1)

$(DYLIB_LAPACK2): $(C_SRC_LAPACK2)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_LAPACK2) -llapacke -lcblas -o $(DYLIB_LAPACK2)
$(DYLIB_LAPACK3): $(C_SRC_LAPACK3)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_LAPACK3) -llapacke -lcblas -o $(DYLIB_LAPACK3)

$(DYLIB_LAPACK4): $(C_SRC_LAPACK4)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_LAPACK4) -llapacke -lcblas -o $(DYLIB_LAPACK4)

$(DYLIB_LAPACK5): $(C_SRC_LAPACK5)
	@mkdir -p $(@D)
	 $(CC) $(CFLAGS) -shared -fPIC $(C_SRC_LAPACK5) -llapacke -lcblas -o $(DYLIB_LAPACK5)

$(DYLIB_CBLAS): $(C_SRC_CBLAS) 
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC_CBLAS) -lcblas -o $(DYLIB_CBLAS)

clean: 
	rm -fr $(JAVA_HOME)/class
	rm -fr $(JAVA_HOME)/lib
	rm -f $(JAVA_HOME)/lib/*.so
	rm -f $(JAVA_HOME)/test/*.class
	rm -f $(JAVA_HOME)/test/*.txt
	rm -f $(JAVA_HOME)/test/*.out
	rm -f $(JAVA_HOME)/test/*.serial
