#!/bin/bash
echo "Running on 64 bit Linux"
export CLASSPATH="src:lib/jhdf.jar:lib/jhdf5.jar:lib/jhdf5obj.jar:lib/jhdfobj.jar:examples"
java -Djava.library.path=lib/native/linux/64 -cp $CLASSPATH ReadGctxExample
