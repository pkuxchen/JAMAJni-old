JAVA_HOME = ..

JAMAJNI_CLASS = JAMAjniTest.class
JAMAJNI_RUN = JAMAjni_test_result.txt

all:  $(JAMAJNI_CLASS) $(JAMAJNI_RUN)

$(JAMAJNI_CLASS):
	javac -classpath $(JAVA_HOME)/class:. JAMAjniTest.java 

$(JAMAJNI_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. JAMAjniTest > JAMAjni_test_result.txt

clean: 
	rm -f *.txt *.class *.out *.serial
