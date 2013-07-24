#!/usr/bin/perl -w
use strict;

my $gold = "testdata/gold";

my $size = shift() || "micro";
my $new = shift() || "/tmp/regression_trial";

# if you're adding new columns, or expect values to change, you can ignore
# them temporarily by adding them to @ignore.  However, once tied out the gold
# results should be updated and @ignore set back to empty for the commit
my @ignore = ();
#my @ignore = ("failure_reasons", "total_reads", "total_pairs","improper_pairs", "tumor_ref_qscores", "tumor_alt_qscores", "normal_ref_qscores", "normal_alt_qscores");

# similar to the above, if keywords are changing a global SED can be applied to the
# the new data, but upon commit this should be reset
my $NEW_SED = "";

foreach my $prefix ("ccle.$size") {
    # if the gold  file does not exist, check if a gzip version does and uncompress it to tmp
    my $goldCs = "$gold/$prefix.call_stats.txt";
    if (-e "$goldCs.gz") {
        system("mkdir -p /tmp/regression_gold && gunzip -c $goldCs.gz > /tmp/regression_gold/$prefix.call_stats.txt");
        $goldCs = "/tmp/regression_gold/$prefix.call_stats.txt";
    }

    my $newCs = "$new/$prefix.call_stats.txt";

    my $goldCut = "/tmp/gold.cut";
    my $newCut = "/tmp/new.cut";

    # figure out what columns to ignore
    my $CUT_CMD = getCutCommand($goldCs, @ignore);
    my $NEW_CUT_CMD = getCutCommand($newCs, @ignore);
    
    print "Comparing call stats (ignoring variable columns with " . join(",", @ignore) . " )\n";

    system("cat $goldCs | grep -v \"#\" $CUT_CMD > $goldCut") == 0 or die();
    system("cat $newCs | grep -v \"#\" $NEW_CUT_CMD $NEW_SED > $newCut") == 0 or die();

    system("diff $goldCut $newCut") == 0 or print STDERR "FAILURE: new call stats differs from gold file\n";

    my $lineCount = `wc -l $goldCut`;
    chomp($lineCount);

    print "SUCCESS! compared $lineCount lines\n";

    system("gunzip -c $gold/$prefix.vcf.gz | grep -v MuTect > /tmp/gold.vcf");
    system("cat $new/$prefix.vcf | grep -v MuTect > /tmp/new.vcf");

    system("diff /tmp/gold.vcf /tmp/new.vcf") == 0 or print STDERR "FAILURE: new VCF differs from gold file\n";

    $lineCount = `wc -l /tmp/gold.vcf`;
    chomp($lineCount);

    print "SUCCESS! compared VCFs with $lineCount lines\n";


    foreach my $suf ("tumor.depth.wig.txt","normal.depth.wig.txt") {
        print "Comparing $suf files...";
        
		my $fullgold = "$gold/$prefix.$suf";
        if (-e "$gold/$prefix.$suf.gz") {
	        $fullgold = "/tmp/gold.$prefix.$suf";
			system("gunzip -c $gold/$prefix.$suf > $fullgold") == 0 or die();
		}
		
        if (system("diff $fullgold $new/$prefix.$suf") == 0) {
            print "SUCCESS!\n";
        } else {
            print STDERR "FAILURE: new $suf wiggle differs from gold file\n";
        }
    }


}

sub getCutCommand {
    my $cs = shift();
    my @ignore = @_;

    # figure out what columns to ignore
    my %goldMap = ();
    my $goldHeader = `tail -n+2 $cs | head -1`;
    chomp($goldHeader);
    my @parts = split("\t", $goldHeader);
    for(my $i=0; $i<scalar(@parts); $i++) { $goldMap{$parts[$i]} = 1 + $i;}
    
    my $CUT_CMD = "";
    if (scalar(@ignore) > 0) {
        foreach my $name (@ignore) {
            if (exists $goldMap{$name}) { 
                delete $goldMap{$name}; 
            }
        }
    }


    my @idx = ();
    while ( my ($key, $value) = each(%goldMap) ) {
        push(@idx, $value);
    }
    my @sIdx = sort {$a <=> $b} @idx;
    $CUT_CMD = "| cut -f" . join(",",@sIdx); 
    
    return ($CUT_CMD);
}
