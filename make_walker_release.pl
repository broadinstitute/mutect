#! /usr/bin/perl -w
use strict;

if (scalar(@ARGV) < 4) {
    die("usage: make_walker_release.pl <package-name> <cga-gatk-dir> <sting-dir> <tag>\n  for example: make_walker_release.pl muTect /xchip/cga1/kcibul/cga/analysis_pipeline/gatk /xchip/cga1/kcibul/sting muTect-0.8.0\n");
}
my ($tool, $CGA_DIR, $GATK_DIR, $tag) = @ARGV;

my $cwd = `pwd`;
chomp($cwd);

my $tmp = "/tmp/$tool-dist-tmp";
if (-e $tmp) { `rm -rf $tmp`; }
`mkdir -p $tmp`;

# update CGA and get revision info
chdir($CGA_DIR);
system("svn up > /dev/null") == 0 or die();
`svn info | grep "Revision:" | awk '{ print "$tool Revision: " \$2 }' > $tmp/version.txt`;
my $version = `cat $tmp/version.txt | grep -i $tool | awk '{print \$3 }'`;
chomp($version);

# check to see if this tag exists
system("svn info https://svnrepos/CancerGenomeAnalysis/tags/$tag &> /dev/null") != 0 or die("Tag for $tag already exists in subversion!\n"); 

my $outputZip = "$cwd/$tag-bin.zip";
if (-e $outputZip) { die("release $outputZip already exists!!\n"); }

# update GATK and get revision info
chdir($GATK_DIR);
system("git pull > /dev/null") == 0 or die();
`git describe | awk '{ print "GATK Revision: " \$0 }' >> $tmp/version.txt`;

# do a clean build
#chdir($GATK_DIR);
#system("ant clean") == 0 or die();
#system("ant -Dexternal.dir=$CGA_DIR");

#system("ant dist") == 0 or die();
system("ant clean") == 0 or die();
system("ant -Dexternal.dir=$CGA_DIR -Dexecutable=" . lc($tool) . " package") == 0 or die();

# move the executable over to the release directory
system("cp $GATK_DIR/dist/packages/$tool-*/$tool.jar $tmp/$tag.jar") == 0 or die();

# move the license over to the release directory
system("cp $CGA_DIR/walkers/" . lc($tool) . ".LICENSE.TXT $tmp/LICENSE.TXT") == 0 or die();

# zip it up
chdir($cwd);
system("zip -j $outputZip $tmp/*") == 0 or die();

# tag subversion with the version
system("svn copy https://svnrepos/CancerGenomeAnalysis/trunk/analysis_pipeline/gatk https://svnrepos/CancerGenomeAnalysis/tags/$tag -m \"Tagging $tag\"") == 0 or die();




