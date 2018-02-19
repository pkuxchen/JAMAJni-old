JAVA_HOME = ..

LAPACK_CLASS = lapacktest.class
LAPACK_RUN = lapack_test_result.txt
BLAS_CLASS = blastest.class
BLAS_RUN = blas_test_result.txt
JAMAJNI_CLASS = JAMAjniTest.class
JAMAJNI_RUN = JAMAjni_test_result.txt

all:  $(JAMAJNI_CLASS) $(JAMAJNI_RUN)
# $(BLAS_CLASS) $(BLAS_RUN) $(LAPACK_CLASS) $(LAPACK_RUN)

#$(LAPACK_CLASS):
#	javac -classpath $(JAVA_HOME)/class:. lapacktest.java 

#$(LAPACK_RUN):
#	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. lapacktest > lapack_test_result.txt

#$(BLAS_CLASS):
#	javac -classpath $(JAVA_HOME)/class:. blastest.java 

#$(BLAS_RUN):
#	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. blastest #> blas_test_result.txt


$(JAMAJNI_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JAMAjniTest.java 

$(JAMAJNI_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JAMAjniTest > JAMAjni_test_result.txt

clean: 
	rm -f *.txt *.class
