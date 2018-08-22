JAVA_HOME = ..

JAMAJNI_CLASS = JAMAjniTest.class
JAMAJNI_RUN = JAMAjni_test_result.txt
CBLAS_CLASS = JAMAJniExamplesCBLAS.class
CBLAS_RUN = CBLAS_examples_result.txt
LAPACKE_CLASS = JAMAJniExamplesLAPACKE.class
LAPACKE_RUN = LAPACKE_examples_result.txt
JAMA_CLASS = JAMAExamples.class
JAMA_RUN = JAMA_examples_result.txt

all:  $(JAMAJNI_CLASS) $(JAMAJNI_RUN) $(CBLAS_CLASS) $(CBLAS_RUN) $(LAPACKE_CLASS) $(LAPACKE_RUN) $(JAMA_CLASS) $(JAMA_RUN)

$(JAMAJNI_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JAMAjniTest.java 

$(JAMAJNI_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JAMAjniTest > JAMAjni_test_result.txt

$(CBLAS_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JAMAJniExamplesCBLAS.java 

$(CBLAS_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JAMAJniExamplesCBLAS > CBLAS_examples_result.txt

$(LAPACKE_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JAMAJniExamplesLAPACKE.java 

$(LAPACKE_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JAMAJniExamplesLAPACKE > LAPACKE_examples_result.txt

$(JAMA_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JAMAExamples.java 

$(JAMA_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JAMAExamples > JAMA_examples_result.txt

clean: 
	rm -f *.txt *.class *.out *.serial
