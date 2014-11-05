ECHO "Running on 64-bit windows version"
@ECHO OFF
set CLASSPATH=.
set CLASSPATH=%CLASSPATH%;src;lib/jhdf.jar;lib/jhdf5.jar;lib/jhdf5obj.jar;lib/jhdfobj.jar;examples

java -Djava.library.path=lib/native/windows/64 -cp %CLASSPATH% ReadGctxExample
