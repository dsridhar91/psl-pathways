#!/bin/bash

mvn compile > /dev/null
mvn package > /dev/null
JAR=$(ls target/edu.ucsc.cs-*.jar | head -n1)
CLASSPATH=$(mvn dependency:build-classpath | grep -v "\[INFO\]")	
java -cp $CLASSPATH:$JAR $*
