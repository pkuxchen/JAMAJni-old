<b>JAMAJni</b> is a JAVA package providing a java interface for lapack and blas library and using the classes defined by JAMA Package. It's a tool for numerical linear algebra computations.


Build Instructions
------------------

* JAMAJni requires the installation of lapack, blas, lapacke and cblas.

* To compile the package, enter src directory and execute "make". Notice that you may have to change the extension of generated libraries in the Makefile base on your operating system. On OS X you have to change all the extensions of dynamic library to .dylib while on Linux the corresponding extensions are .so or .a. 

* To clean generated file, type “make clean” on the command line.  


Running the tests
-----------------
For testing, enter test directory and execute “make” . If you want to clean testing results and all class files, type "make clean".  
There are four test files. The "JAMAjniTest.java" will test all the methods in JAMAJni package and report the errors. The "JAMAExamples.java" will provide specific examples for basic linear algebra operations. It can clearly show you how to use methods defined in JAMAJni. If you are interested in how to use functions in blas and lapack library to do matrix operations, the "JAMAJniExamplesCBLAS.java" and "JAMAJniExamplesLAPACKE.java" will give you specific examples. However, it is not necessary to go into blas and lapack if you just want to use JAMAJni to do some operations.

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
| Xiang Chen | pkuchenxiang@pku.edu.cn   | Visiting student, Department of Biostatistics  UCLA|
| Lu Zhang (maintainer)| lu.zhang@ucla.edu    | PhD student, Department of Biostatistics UCLA  |                             
| Sudipto Banerjee | sudipto@ucla.edu   | Professor, Department of Biostatistics  UCLA |
| Diyang Wu | wudiyangabc@hotmail.com    | MS student, Department of Biostatistics  UCLA|
| Zijie Zhao | zijzhao@ucla.edu    | MS student, Department of Biostatistics Havard University|
<!--- --->
                             


Licensing
---------
JAMAJni is licensed under the Creative Commons Attribution License. 



