# JAMAJni
* JAMAJni is a JAVA package providing a java interface for lapack and blas libraries and using the classes defined by JAMA Package. It's a tool for numerical linear algebra computations. JAMAJni calls an intermediate c interface - cblas and - lapacke which requires previous installation. (Detailed installation instructions are listed below). 

* JAMAJni has a sibling package called JAMAJniLite, which calls lapack and blas libararies. You can find it on GitHub here.


Build Instructions
------------------

* JAMAJni requires the installation of lapacke and cblas. (The command in Ubuntu is "sudo apt-get install libatlas-base-dev liblapacke-dev"). We recommend you to use cblas library in "libatlas-base-dev" for some other cblas libraries may fail to compile.

* To compile the package, enter src directory and execute "make". Notice that you may have to change the extension of generated libraries in the Makefile based on your operating system. On OS X you have to change all the extensions of dynamic library to .dylib while on Linux the corresponding extensions are .so. 

* To clean generated file, type “make clean” on the command line.  


Running the tests
-----------------
* For testing, enter test directory and execute “make” . If you want to clean testing results and all class files, type "make clean".  

* There are four test files. The "JAMAjniTest.java" will test all the methods in JAMAJni package and report the errors. The "JAMAExamples.java" will provide specific examples for basic linear algebra operations. It can clearly show you how to use methods defined in JAMAJni. If you are interested in how to use functions in blas and lapack library to do matrix operations, the "JAMAJniExamplesCBLAS.java" and "JAMAJniExamplesLAPACKE.java" will give you specific examples. However, it is not necessary to go into blas and lapack if you just want to be a user of JAMAJni.

Notes
---------
This package is intended for some basic problems we encounter when solving linear algebra problems. So we only include several most basic and widely used routines of the blas and lapack library.


Source Repository
-----------------
JAMAJni's source-code repository is hosted here on GitHub.


Authors
---------

| Name   | Email       |              |
|:------ |:----------- | :----------- |
| Xiang Chen (maintainer)| pkuchenxiang@pku.edu.cn   | Visiting student, Department of Biostatistics  UCLA|
| Lu Zhang | lu.zhang@ucla.edu    | PhD student, Department of Biostatistics UCLA  |                             
| Sudipto Banerjee | sudipto@ucla.edu   | Professor, Department of Biostatistics  UCLA |
| Diyang Wu | wudiyangabc@hotmail.com    | MS student, Department of Biostatistics  UCLA|
| Zijie Zhao | zijzhao@ucla.edu    | MS student, Department of Biostatistics Harvard University|
<!--- --->
                             


Licensing
---------
JAMAJni is licensed under the Creative Commons Attribution License. 



