#!/usr/bin/perl
use Getopt::Long;

my $usage =
"
\t try to identify all possible loops from a given pdb file for a given chain

\t-p or -pdb     \tthe pdb file
\t-c or -chain   \tthe one letter chain identifier
\t-v or -verbose \tverbose mode

author \tNils Woetzel
date   \t11/06/2013\n
";

## initialize the commandline variables
$pdb_file  = "";
$chain_id  = " ";

## print help if no argument is given
die $usage if( $#ARGV == -1);

## set the options form the commandline
GetOptions
(
  "p|pdb=s"      => \$pdb_file,
  "c|chain=s"    => \$chain_id,
  "v|verbose"    => \$verbose
) or die "invalid command line!\n";

## display the variables
if( defined $verbose)
{
  print "p|pdb   = $pdb_file\n";
  print "c|chain = $chain_id\n";
}

open( IN, "<", $pdb_file) or die "cannot open file $pdb_file for reading!\n";
my $pdb_name = substr( $pdb_file, 0, length( $pdb_file)-4);

# store the residue as keys and the values are for the SSE type
my %end_residues   = {};
my %start_residues = {};

# collect all possible start and end residues of loops by parsing HELIX and SHEET entries and substracting/adding one to
# start and end residue seqids respectively
while( $line = <IN>)
{
  if( $line =~ /^SHEET/)
  {
    my $current_chain_id = substr( $line, 21, 1);
    # test the chain id
    ( $current_chain_id eq $chain_id) or next;
    my $res_end = substr( $line, 33, 4) + 1;
    my $res_begin   = substr( $line, 22, 4) - 1;
    $start_residues{ $res_end} = 'S';
    $end_residues{ $res_begin} = 'S';
#    print "$res_end\t$res_begin\n";
    next;
  }
  elsif( $line =~ /^HELIX/)
  {
    my $current_chain_id = substr( $line, 19, 1);
    # test the chain id
    ( $current_chain_id eq $chain_id) or next;

    my $res_end = substr( $line, 33, 4) + 1;
    my $res_begin   = substr( $line, 21, 4) - 1;
    $start_residues{ $res_end} = 'H';
    $end_residues{ $res_begin} = 'H';
#    print "$res_end\t$res_begin\n";
    next;
  }
  elsif( $line =~ /^ATOM/)
  {
    last;
  }
}

# parse all atom lines and get the full residue information required by the res locators
my @loops;
my $count = 0;
while( $line = <IN>)
{
  if( $line =~ /^ATOM/)
  {
    my $current_chain_id = substr( $line, 21, 1);
    # test the chain id
    if( $current_chain_id ne $chain_id)
    {
      next;
    }

    my $res = substr( $line, 22, 4) + 0;
    if( exists $start_residues{$res})
    {
      ++$count;
      my $res_name = substr( $line, 17, 3);

      my $loop_def = "$pdb_name $chain_id ".sprintf( "%02d", $count)." $start_residues{$res} $res_name ".sprintf( "%4d", $res);
      push( @loops, $loop_def);

      delete $start_residues{$res};
    }
    elsif( exists $end_residues{$res})
    {
      my $res_name = substr( $line, 17, 3);
      my $loop_def_end = " $end_residues{$res} $res_name ".sprintf( "%4d", $res);
      foreach my $loop_def_start (@loops)
      {
        print $loop_def_start." ".sprintf( "%02d", $count).$loop_def_end."\n";
      }
      delete $end_residues{$res};
    }
  }
}

close( IN);
