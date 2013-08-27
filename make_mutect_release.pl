#! /usr/bin/perl -w
use strict;

if (scalar(@ARGV) < 3) {
    die("usage: make_mutect_release.pl <tmp-dir> <mutect-tag> <gatk-tag>\n  for example: make_walker_release.pl /tmp 1.1.5 2.5\n");
}
my ($TMP_DIR, $mutect_tag, $gatk_tag) = @ARGV;

my $BASE_DIR = "$TMP_DIR/mutect-dist";
my $GATK_DIR = "$BASE_DIR/gatk-protected";
my $MUTECT_DIR = "$BASE_DIR/mutect-src";
my $cwd = `pwd`;
chomp($cwd);

my $TMP_DIST = "/tmp/mutect-dist-zip";
if (-e $TMP_DIST) { `rm -rf $TMP_DIST`; }
`mkdir -p $TMP_DIST`;

if (-e $BASE_DIR) { `rm -rf $BASE_DIR`; }
`mkdir -p $MUTECT_DIR`;


# update CGA and get revision info
system("cd $MUTECT_DIR && git clone git\@github.com:broadinstitute/mutect.git") == 0 or die();

# check to see if this tag exists
my $cnt = `cd $MUTECT_DIR/mutect && git ls-remote --tags -q | grep refs/tags/$mutect_tag | wc -l`;
chomp($cnt);
if ($cnt == 0) { die("ERROR: release tag $mutect_tag does not exist!\n"); }
system("cd $MUTECT_DIR/mutect && git reset --hard $mutect_tag") == 0 or die();
`cd $MUTECT_DIR/mutect && git describe --tags | awk '{ print "MuTect Revision: " \$0 }' > $TMP_DIST/version.txt`;

my $outputZip = "$cwd/muTect-$mutect_tag-bin.zip";
if (-e $outputZip) { die("release $outputZip already exists!!\n"); }

# update GATK and get revision info
chdir($BASE_DIR);
system("git clone git\@github.com:broadgsa/gatk-protected.git") == 0 or die();
chdir($GATK_DIR);
system("git reset --hard $gatk_tag") == 0 or die();
`git describe --tags | awk '{ print "GATK Revision: " \$0 }' >> $TMP_DIST/version.txt`;

# do a clean build
system("cd $GATK_DIR && ant clean") == 0 or die();
system("cd $GATK_DIR && ant -Dexternal.dir=$MUTECT_DIR -Dexecutable=mutect package") == 0 or die();

# move the executable over to the release directory
system("cp $GATK_DIR/dist/packages/muTect-*/muTect.jar $TMP_DIST/muTect-$mutect_tag.jar") == 0 or die();

# move the license over to the release directory
system("cp $MUTECT_DIR/mutect/mutect.LICENSE.TXT $TMP_DIST/LICENSE.TXT") == 0 or die();

# zip it up
chdir($cwd);
system("zip -j $outputZip $TMP_DIST/*") == 0 or die();




