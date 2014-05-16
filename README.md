Overview
--------------------------------------------------------
For complete details see

http://www.broadinstitute.org/cancer/cga/mutect

And for support please visit

http://gatkforums.broadinstitute.org/categories/mutect


How to build MuTect
--------------------------------------------------------

PREREQUISITES

To compile MuTect you must be using Java 1.7 and Maven 3.0+.  After building successfully, the MuTect JAR file will be in mutect/target/mutect-*.jar

BUILD STEPS

    # make a new source directory (e.g. mutect-src)
    mkdir mutect-src
    cd mutect-src

    # get MuTect source
    git clone git@github.com:broadinstitute/mutect.git

    # get the GATK source and set to the latest tested version
    git clone git@github.com:broadgsa/gatk-protected.git
    cd gatk-protected
    git reset --hard 3.1 
    
    # build the GATK first and install it to the local mvn repo.  Once GATK publishes to a public repo this will be much simpler
    mvn -Ddisable.queue install

	# build MuTect and run unit tests (the target jar will be in target/mutect-*.jar)
	cd ../mutect
	mvn verify
	
	# run integration tests, if you like
	./run_regression.pl
	./tieout_regression.pl

