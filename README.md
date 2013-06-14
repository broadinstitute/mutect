Overview
--------------------------------------------------------
For complete details see

http://www.broadinstitute.org/cancer/cga/mutect

And for support please visit

http://gatkforums.broadinstitute.org/categories/mutect


How to build MuTect
--------------------------------------------------------

PREREQUISITES

To compile MuTect you must be using Java 1.6 and Ant 1.8.  In addition, In addition, it is necessary to download the Apache BCEL Library and put it in either your system classpath, or your Ant installation lib directory.  After building successfully, the MuTect JAR file will be in mutect-src/gatk-protected/dist/packages/muTect-*/muTect.jar 

BUILD STEPS

    # make a new source directory (e.g. mutect-src)
    mkdir mutect-src
    cd mutect-src

    # get MuTect source
    git clone git@github.com:broadinstitute/mutect.git
    cd ..

    # get the GATK source and set to the latest tested version
    git clone git@github.com:broadgsa/gatk-protected.git
    cd gatk-protected
    git reset --hard 2.5
    
    # build
    ant -Dexternal.dir=`pwd`/../mutect-src -Dexecutable=mutect package
