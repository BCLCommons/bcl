#/usr/bin/perl
use Getopt::Long;

my $usage =
"
\t from a loop file generated e.g. by make_all_loops.pl create all locator and site definitions required as input by
   BCL::FusionProtein

\t-l or -loop_file     \tthe loop file
\t-m or -min_length    \tthe minimal length of the loops to consider
\t-o or -output_prefix \tthe output prefix
\t-v or -verbose       \tverbose mode

author \tNils Woetzel
date   \t11/06/2013\n
";

## initialize the commandline variables
$loop_file     = "";
$min_length    = 0;
$output_prefix = "";

## print help if no argument is given
die $usage if( $#ARGV == -1);

## set the options form the commandline
GetOptions
(
  "l|loop_file=s"     => \$loop_file,
  "m|min_length=i"    => \$min_length,
  "o|output_prefix=s" => \$output_prefix,
  "v|verbose"         => \$verbose
) or die "invalid command line!\n";

## display the variables
if( defined $verbose)
{
  print "l|loop_file     = $loop_file\n";
  print "m|min_length    = $min_length\n";
  print "o|output_prefix = $output_prefix\n";
}

my $locators_path = $output_prefix."locators";
if( !-d $locators_path)
{
  mkdir $locators_path or die "cannot create directory for locators: $locators_path\n";
}

my $sites_path = $output_prefix."sites";
if( !-d $sites_path)
{
  mkdir $sites_path or die "cannot create directory for sites: $sites_path\n";
}

# open the loop file
open( IN, "<", $loop_file) or die "cannot open file $loop_file for reading!\n";

# iterate over all lines in the file
while( $line = <IN>)
{
  ( my $pdb, my $chain_id, my $loop1, my $sse1, my $res_name1, my $res_id1, my $loop2, my $sse2, my $res_name2, my $res_id2) = split( /\s+/, $line);
  if( $res_id2 - $res_id1 < $min_length)
  {
    next;
  }
  my $loop = $res_id1."_".$res_id2;
  my $filename_base = $pdb."_".$loop1.$sse1.$loop2.$sse2.'_'.$loop;
  
  my $filename_site = $sites_path."/".$filename_base.".site";
  open( OUT, ">", $filename_site) or die "cannot open file for writing: $filename_site\n";

  print OUT "bcl::pdb::Site\n";
  print OUT "\"AC1\"\n";
  print OUT "SOFTWARE\n";
  print OUT "\"BINDING $loop\"\n";
  print OUT "bcl::storage::List<bcl::pdb::ResidueSimple>\n";
  print OUT "  2\n";
  print OUT "  bcl::pdb::ResidueSimple\n";
  print OUT "    \"$res_name1\"     $chain_id       $res_id1      ' '\n";
  print OUT "  bcl::pdb::ResidueSimple\n";
  print OUT "    \"$res_name2\"     $chain_id       $res_id2      ' '\n";
  print OUT "bcl::storage::List<bcl::pdb::ResidueSimple>\n";
  print OUT "  0\n";
  print OUT "bcl::util::ShPtr<bcl::pdb::Ligand>\n";
  print OUT "  0\n";
  print OUT "  NULL\n";

  close OUT;
  
  my $filename_locator = $locators_path."/". $filename_base.".locators";
  open( OUT, ">", $filename_locator) or die "cannot open file for writing: $filename_locator\n";

  print OUT "bcl::assemble::LocatorAA\n";
  print OUT "  LocatorAA(locator_chain=$chain_id,seq_id=".($res_id1-1).",use_pdb_id=1)\n";
  print OUT "bcl::assemble::LocatorAA\n";
  print OUT "  LocatorAA(locator_chain=$chain_id,seq_id=".($res_id2+1).",use_pdb_id=1)\n";

  close OUT;
}
close( IN);
