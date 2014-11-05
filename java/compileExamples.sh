#!/bin/bash
export CLASSPATH="src:lib/jhdf.jar:lib/jhdf5.jar:lib/jhdf5obj.jar:lib/jhdfobj.jar"

javac -cp $CLASSPATH examples/ReadGctxExample.java

# -Djava.library.path=lib/macosx/64 
