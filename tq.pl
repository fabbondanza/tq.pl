#!/usr/bin/perl -w

use strict;
use Getopt::Std;

##########################
#Provide an array to hold all the lines in from association file.

my @assfile;

my %opts;

getopt('qtr', \%opts);

##########################
#Open input filehandle for the qtdt association file.

if (defined ($opts{q})) {

	print "\nYour file containing the association p-values is called: $opts{q}\n";

	open (ASSFILE, "<$opts{q}");
}

else {

	print "\nPlease enter the filename for the file containing your association p-values: ";

	my $assfile = <STDIN>;

	chomp($assfile);

	open (ASSFILE, "<$assfile");

}

##########################
#Open output filehandle for the table to be drawn to.

if (defined ($opts{t})) {

	print "\nYour file containing the new table will be called: $opts{t}\n";

	open (TABLEFILE, ">$opts{t}");
}

else {

	print "\nPlease enter the filename for me to write the new table to: ";

	my $tablefile = <STDIN>;

	chomp($tablefile);

	open (TABLEFILE, ">$tablefile");
}

##########################
#Open input filehandle for regress.tbl file.

if (defined ($opts{r})) {

        print "\nYour file containing the regress.tbl is called: $opts{r}\n";

        open (REGRESSFILE, "<$opts{r}");
}

#else {

#        print "\nPlease enter the filename for the file containing your regress.tbl: ";

#        my $regressfile = <STDIN>;

#        chomp($regressfile);
 
#        open (REGRESSFILE, "<$regressfile");
#}

##########################
#Read the assfile data into the assfile array.

my $line;
my $count = 0;

while(defined($line = <ASSFILE>)) {

		$assfile[$count] = $line;
		++$count;
}

##########################
#Read the regressfile data into the regress array if required.

$count = 0;

my @regressFile;

if (defined ($opts{r})) {

	while(defined($line = <REGRESSFILE>)) {

		$regressFile[$count] = $line;
		++$count;
	}

	close (REGRESSFILE);
}

##########################
#Close the ASSFILE filehandle

close (ASSFILE);


##########################
#Read through each line and identify the trait name(s), which will become keys to the 'results' hash, and marker name(s), that will become keys to hashes of the 'results' hash.

my %results;
my %table;

my @markerOrder;
my $tempTrait;
my $tempMarker;
my $tempAllele;
my $tempAll = 0;

foreach $line (@assfile){

	if($line =~ /^Testing trait:\s+(\S+)\n/) {

		$results{$1} = {};
		$tempTrait = $1;

	}

        if($line =~ /^Testing marker:\s+(\S+)\n/) {
        
                $results{$tempTrait}{$1} = [];
		$tempMarker = $1;
		$tempAllele = 0;

		if (keys %results == 1) {

			push @markerOrder, $1;
		}

        }

	if($line =~ /\s+(\S+)\s+\*\*\* not tested \*\*\*\s+\(.+probands\)/) {

		if(($tempAllele + 1) eq $1) {

			$results{$tempTrait}{$tempMarker}[$tempAllele] = "?";
			$table{$tempMarker}[$tempAllele]{$tempTrait} = "?";

		$tempAllele++;
		}
	}

	else {

 	       if($line =~ /\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)?\s+\(.+probands\)/) {

			if(($tempAllele + 1) eq $1) {

			        if(defined($2)) {

					$results{$tempTrait}{$tempMarker}[$tempAllele] = $2;
					$table{$tempMarker}[$tempAllele]{$tempTrait} = $2;
				}

				else {

					$results{$tempTrait}{$tempMarker}[$tempAllele] = ".";
                        	        $table{$tempMarker}[$tempAllele]{$tempTrait} = ".";
				}

				$tempAllele++;

			}

			else {

				if($1 eq "All") {

					$tempAll = 1;

		                        if(defined($2)) {
        
                		                $results{$tempTrait}{$tempMarker}[$tempAllele] = $2;
        		                        $table{$tempMarker}[$tempAllele]{$tempTrait} = $2;
                        		}
                                
		                        else {
                                
        	        	                $results{$tempTrait}{$tempMarker}[$tempAllele] = "-";
                		                $table{$tempMarker}[$tempAllele]{$tempTrait} = "-";
                        		}


				}

				else {
					print "Something's gone wrong here - go call Tom for help.";

				}

			}
         
        	}

	}

}

if ($tempAll == 0) {

	print TABLEFILE "Marker(s)\tAlleles(s)\t";
	}

else {
	
        print TABLEFILE "Marker(s)\tGlobal\t";
        }


foreach	$tempTrait (sort (keys %results)) {

	print TABLEFILE "$tempTrait\t";
}

print TABLEFILE "\n";

my $Allele;


##########################
#First figure out the marker order

my $tempPlusTrait;
my $tempPlusMarker;

foreach $tempMarker (@markerOrder) {

	$Allele = 0;

	print TABLEFILE "$tempMarker";

	if (defined @{$table{$tempMarker}}) {

		foreach $tempAllele (@{$table{$tempMarker}}) {

			if ($tempAll == 0) {

				print TABLEFILE "\t" . (++$Allele);
				}

			else {

        	                print TABLEFILE "\tAll";
                	        }

			foreach $tempTrait (sort (keys %{$tempAllele})) {

				my $direction = "";

				if (defined ($opts{r})) {

				my $regress;

				$count = 0;

				$tempPlusTrait = $tempTrait;
                                $tempPlusMarker = $tempMarker;

				$tempPlusTrait =~ s/\+/\\\+/g;
                                $tempPlusMarker =~ s/\+/\\\+/g;

					while (defined ($regress = $regressFile[$count])) {

						++$count;

#						print "$tempPlusTrait\n";
#						print "$tempPlusMarker\n";

						if ($regress =~ /Trait:\s$tempPlusTrait\s*Marker:\s$tempPlusMarker\s*Allele:\s$Allele/) {

							while (defined ($regress = $regressFile[$count])) {

							++$count;

								if ($regress =~ /FULL HYPOTHESIS/) {

									while (defined ($regress = $regressFile[$count])) {

										++$count;

										if ($regress =~ /\s+means\s:\s+\S+\s+(-)?\S/) {

											if (defined ($1)) {
												$direction = $1;
											}

											else {

											$direction = "+";
											}

										last;
										}
									}
								last;
								}
							}
						last;
						}
					}
				}

				print TABLEFILE "\t" . $direction . ${$tempAllele}{$tempTrait};

			}

			print TABLEFILE "\n";
		}
	}

	else {

        print TABLEFILE "\t?" x (keys %results) . "\t\?\n";

	}
}

close (TABLEFILE);
