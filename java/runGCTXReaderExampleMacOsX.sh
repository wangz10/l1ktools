#!/bin/bash
echo "Running on 64 bit MacOsX"
java -Djava.library.path=lib/native/macosx/64 -cp $CLASSPATH ReadGctxExample
