#! /usr/bin/perl -w
use strict;

if (scalar(@ARGV) < 3) {
    die("usage: make_walker_release.pl <mutect-dir> <gatk-dir> <mutect-tag> <gatk-tag>\n  for example: make_walker_release.pl  /xchip/cga1/kcibul/cga/analysis_pipeline/gatk /xchip/cga1/kcibul/sting 1.1.5 2.5\n");
}
my ($CGA_DIR, $GATK_DIR, $tag, $gatk_tag) = @ARGV;
my $tool = "muTect";

my $cwd = `pwd`;
chomp($cwd);

my $tmp = "/tmp/$tool-dist-tmp";
if (-e $tmp) { `rm -rf $tmp`; }
`mkdir -p $tmp`;

# check to see if this tag exists
my $cnt = `git ls-remote --tags -q | grep refs/tags/$tag | wc -l`;
chomp($cnt);
if ($cnt == 0) { die("ERROR: release tag $tag does not exist!\n"); }

# update CGA and get revision info
chdir($CGA_DIR);
system("git pull") == 0 or die();
system("git reset --hard $tag") == 0 or die();
`git describe | awk '{ print "$tool Revision: " \$0 }' > $tmp/version.txt`;

my $outputZip = "$cwd/muTect-$tag-bin.zip";
if (-e $outputZip) { die("release $outputZip already exists!!\n"); }

# update GATK and get revision info
chdir($GATK_DIR);
system("git pull") == 0 or die();
system("git reset --hard $gatk_tag") == 0 or die();
`git describe | awk '{ print "GATK Revision: " \$0 }' >> $tmp/version.txt`;

# do a clean build
chdir($GATK_DIR);
system("ant clean") == 0 or die();
system("ant -Dexternal.dir=$CGA_DIR/.. -Dexecutable=" . lc($tool) . " package") == 0 or die();

# move the executable over to the release directory
system("cp $GATK_DIR/dist/packages/$tool-*/$tool.jar $tmp/$tag.jar") == 0 or die();

# move the license over to the release directory
system("cp $CGA_DIR/" . lc($tool) . ".LICENSE.TXT $tmp/LICENSE.TXT") == 0 or die();

# zip it up
chdir($cwd);
system("zip -j $outputZip $tmp/*") == 0 or die();

# tag subversion with the version
#system("svn copy https://svnrepos/CancerGenomeAnalysis/trunk/analysis_pipeline/gatk https://svnrepos/CancerGenomeAnalysis/tags/$tag -m \"Tagging $tag\"") == 0 or die();




