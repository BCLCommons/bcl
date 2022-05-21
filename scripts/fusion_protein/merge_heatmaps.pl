#!/usr/bin/perl
use Getopt::Long;
use List::Util qw[min max];

my $usage =
"
\t from multiple heatmap files that only have a small tile for a given problem, generate the complete heatmap
   heatmap files are located through: {prefix}{scaffold}_*_*_{donor}_*_*score_heatmap_{scheme}.gnuplot

\t-p  or -prefix   \tthe path and prefix for all gnuplot file that are to be merged
\t-s  or -scaffold \tthe scaffold name
\t-d  or -donor    \tthe donor name
\t-sc or -scheme   \tthe scoring scheme to merge
\t-pl pr -plot     \tplot the gnuplot at the end
\t-v  or -verbose  \tverbose mode

author \tNils Woetzel
date   \t11/04/2013\n
";

## initialize the commandline variables
$prefix   = "";
$scaffold = "";
$donor    = "";
$scheme   = "";

## print help if no argument is given
die $usage if( $#ARGV == -1);

## set the options form the commandline
GetOptions
(
  "p|prefix=s"   => \$prefix,
  "sc|scheme=s"  => \$scheme,
  "s|scaffold=s" => \$scaffold,
  "d|donor=s"    => \$donor,
  "pl|plot"      => \$plot,
  "v|verbose"    => \$verbose
) or die "invalid command line!\n";

## display the variables
if( defined $verbose)
{
  print "p|prefix   = $prefix\n";
  print "sc|scheme  = $scheme\n";
  print "s|scaffold = $scaffold\n";
  print "d|donor    = $donor\n";
}

# list all heatmap files
@files = `ls ${prefix}${scaffold}_*_*_${donor}_*_*score_heatmap_${scheme}.gnuplot`;

my $header;
my $set_object_rest;
my @heatmap;
my $min_x =  10000;
my $max_x =      0;
my $min_y =  10000;
my $max_y =      0;
my $min_z =  10000;
my $max_z = -10000;

my @recs;
$pallet_low  = "set palette defined (0 1 1 1, 0.00001 0 0 1, 1 0 1 0, 2 1 0 0)";
$pallet_high = "set palette defined (0 0 0 1, 1 0 1 0, 2 1 0 0, 2.00001 1 1 1)";

# store the coordinates of the lowest/highest energy
$low_best = 0;
$best_x = 0;
$best_y = 0;

foreach my $file (@files)
{
  chomp( $file);
  # open current heatmap file for reading
  open( IN, "<", $file) or die "cannot open gnuplot file: $file\n";

  my $x_offset;
  my $y_offset;
  my $rec;
  if( not defined( $header))
  {
    while( my $line = <IN>)
    {
      if( $line =~ /range/){ next;}
      if( $line =~ /set terminal/)
      {
        $header .= "set terminal png enhanced transparent font \"Arial,12\" size 4000,4000\n";
        next;
      }
      if( $line =~ /^set output/)
      {
        $header .= "set output \"${prefix}${scaffold}_${donor}_merge_score_heatmap_${scheme}.png\"\n";
        next;
      }
      if( $line =~ /^set xtics.+\(\"\s*(\d+)\s*\"/)
      {
        $x_offset = $1;
        next;
      }
      if( $line =~ /^set ytics \(\"\s*(\d+)\s*\"/)
      {
        $y_offset = $1;
        next;
      }
      if( $line =~ /^set palette defined \(0 (\d+.\d+.\d+)/)
      {
        $low_best = ( $1 eq "1 1 1");
        next;
      }
      if( $line =~ /^set object rect from (\d+.\d+,\d+.\d+) to \d+.\d+,\d+.\d+\s*(.+)$/)
      {
        $rec = $1;
        $set_object_rest = $2;
        last;
      }
      $header .= $line;
    }
  }

  while( my $line = <IN>)
  {
    if( $line =~ /^set xtics.+\(\"\s*(\d+)\s*\"/)
    {
      $x_offset = $1;
      next;
    }
    if( $line =~ /^set ytics \(\"\s*(\d+)\s*\"/)
    {
      $y_offset = $1;
      next;
    }
    if( $line =~ /^set object rect from (\d+.\d+,\d+.\d+) to \d+.\d+,\d+.\d+\s*(.+)$/)
    {
      $rec = $1;
      $set_object_rest = $2;
      next;
    }

    # store the actual value
    if( $line =~ /^(\d+)\s+(\d+)\s+(.+)$/)
    {
      my $x = $1;
      my $y = $2;
      my $z = $3;
      if( $z == 0){ next;}
      $x += $x_offset;
      $y += $y_offset;
      $min_x = min( $x, $min_x);
      $max_x = max( $x, $max_x);
      $min_y = min( $y, $min_y);
      $max_y = max( $y, $max_y);
      $min_z = min( $z, $min_z);
      $max_z = max( $z, $max_z);
      if( ($low_best && $min_z == $z) || (!$low_best && $max_z == $z))
      {
        $best_x = $x;
        $best_y = $y;
      }
      $heatmap[$x][$y] = $z;
    }
  }
  
  # correct the lower and upper rect
  ( my $rec_x, $rec_y) = split( /,/, $rec);
  $rec_x += $x_offset;
  $rec_y += $y_offset;
  push( @recs, "$rec_x,$rec_y");
  close( IN);
}

$min_x -= 1;
$max_x += 1;
$min_y -= 1;
$max_y += 1;

foreach my $rec (@recs)
{
  ( my $rec_x, my $rec_y) = split( /,/, $rec);
  $rec_x -= $min_x;
  $rec_y -= $min_y;
  $header .= "set object rect from $rec_x,$rec_y to ".($rec_x+1).",".($rec_y+1)." $set_object_rest\n";
}
if( 0)
{
  $best_x -= $min_x;
  $best_y -= $min_y;
  $best_x += 0.5;
  $best_y += 0.5;
  $header .= "set object rect from $best_x,$best_y to ".($best_x+1).",".($best_y+1)." front fillcolor rgb \"red\" fs empty border 1 linewidth 2\n";
}

if( $low_best)
{
  $header .= "$pallet_low\n";
}
else
{
  $header .= "$pallet_high\n";
}

my $output_filename = "${prefix}${scaffold}_${donor}_merge_score_heatmap_${scheme}.gnuplot";
$verbose && print "printing heatmap of size: ".($max_x-$min_x+1)." x ".($max_y-$min_y+1)." to $output_filename\n";
$verbose && print "x: [$min_x,$max_x]\ty: [$min_y,$max_y]\tz: [$min_z,$max_z]\n";
open( OUT, ">", $output_filename) or die "cannot open outputfile: $output_filename\n";
print "writing merged gnuplot to: $output_filename\n";
print OUT $header;
print OUT "set size ratio ".($max_y-$min_y)/($max_x-$min_x)."\n";
print OUT "set xrange [ -0.5 : ".($max_x-$min_x+0.5)." ]\n";
print OUT "set yrange [ -0.5 : ".($max_y-$min_y+0.5)." ]\n";
print OUT "set cbrange [ $min_z : $max_z ]\n";

print OUT "set xtics rotate by 90 (";
for( my $x = $min_x; $x <= $max_x; ++$x)
{
  if( $x != $min_x){ print OUT ",";}
  print OUT "\"".sprintf( "%4d", $x)."\" ".($x-$min_x); 
}
print OUT ")\n";

print OUT "set ytics (";
for( my $y = $min_y; $y <= $max_y; ++$y)
{
  if( $y != $min_y){ print OUT ",";}
  print OUT "\"".sprintf( "%4d", $y)."\" ".($y-$min_y);
}
print OUT ")\n";

$empty = 0;
$diff = ($max_z - $min_z)/1000;

if( $low_best)
{
  $empty = $min_z-$diff;
}
else
{
  $empty = $max_z+$diff;
}

print OUT "splot '-' using 1:2:3 with image\n";

for( my $x = $min_x; $x <= $max_x; ++$x)
{
  my $row = $heatmap[$x-$min_x];

  for( my $y = $min_y; $y <= $max_y; ++$y)
  {
    print OUT ($x-$min_x)."\t".($y-$min_y)."\t";
    if( defined( $heatmap[$x][$y]))
    {
      print OUT "$heatmap[$x][$y]";
      if( $heatmap[$x][$y] == $max_z)
      {
        print OUT "#max";
      }
      elsif( $heatmap[$x][$y] == $min_z)
      {
        print OUT "#min";
      }
    }
    else
    {
      print OUT "$empty";
    }
    print OUT "\n";
  }
  print OUT "\n";
}

print OUT "e";
close( OUT);

if( defined( $plot))
{
  $verbose && print "plotting the merged plot to file: ${prefix}${scaffold}_${donor}_merge_score_heatmap_${scheme}.png\n";
  `gnuplot $output_filename`;
}
