#!/usr/bin/perl
use Getopt::Long;
use File::Basename;

########################################################################################################################
## visualize_sspred.pl
##  
## This perl script allows visualizaton of secondary structure predictions as well as native or pool SSE definitions  
##  @author Mert Karakas
##  @date 04/08/2011
##  
########################################################################################################################

$usage = 
"
usage:

\t for an amino acid sequence, different secondary structure prediction methods and secondary structure element pools
\t can be plotted in a gnuplot generated png file.

\t-pdb_id       \t 5letter pdb ide
\t-pools_native \t list of native pools (no filtering will be applied); use \"empty\" is an separator line is desired 
\t-pools        \t list of pools files, they will be filtered according to minimal helix and strand lengths; use \"empty\" is an separator line is desired
\t-methods      \t list of ssmethods (file extensions) to plot
\t-methods_dir  \t directory where sspred files are stored
\t-min_helix    \t minimal helix length (default is 5)
\t-min_strand   \t minimal strand length (default is 3)
\t-out_dir      \t directory where the pngs should be output
\t-h -help      \t print help
\t-v -verbose   \t run in verbose mode

author\tMert Karakas, Nils Woetzel
date\t 04/08/2011\n
";

$pdb_id;
@cmd_pools_native;
@cmd_pools;
@methods;
$methods_dir;
$min_helix = 5;
$min_strand = 3;
$verbose;
$help;
$out_dir;

die $usage if( $#ARGV == -1);

$cl_check =
GetOptions
(
  "pdb_id=s"       => \$pdb_id,
  "pools_native=s{,}" => \@cmd_pools_native,
  "pools=s{,}"      => \@cmd_pools,
  "methods=s{,}"    => \@methods,
  "methods_dir=s"   => \$methods_dir,
  "min_helix=i"     => \$min_helix,
  "min_strand=i"    => \$min_strand,
  "out_dir=s"       => \$out_dir,
  "h|help"          => \$help,
  "v|verbose"       => \$verbose
);

## check parameters
die $usage if( $help || !$cl_check);
die "there needs to be at least one method provided!\n" if( ! defined @methods);

@pools = ();
push( @pools, @cmd_pools_native) if( defined @cmd_pools_native);
push( @pools, @cmd_pools) if( defined @cmd_pools);

$nr_pools = scalar( @pools);
$nr_pools_native = scalar( @cmd_pools_native);

$nr_methods = scalar( @methods);
$nr_res = 0;
for( $i = 0; $i < $nr_methods; ++$i)
{
  $file = $methods_dir.$pdb_id.".".$methods[ $i];
  print "reading file $file\n";
  open( IN, $file) or die "can't open $file!\n";
  $nr_lines = 0;
  while( $line = <IN>)
  {
    chomp( $line);
    next if $line =~ /^\#/ || $line !~ /\w/;
    ++$nr_lines;
  }
  # if first one
  if( $i == 0){ print "$nr_lines lines found\n"; $nr_res = $nr_lines;}
  elsif($nr_res != $nr_lines)
  {
    die "the number of lines do not match ".$nr_res." in ".$methods[$i-i]." vs ".$nr_lines." in ".$methods[$i]."!\n";
  }
}


# initialize gnuplot lines for native structures
$pool_lines = "";
$obj_ctr = 0;

# iterate over native files
for( $i = 0; $i < $nr_pools; ++$i)
{
  $this_file = $pools[ $i];
  if( $this_file eq "empty")
  {
    ++$obj_ctr;
    $y_min = 1.02 + ($nr_pools - $i - 1) * 0.1;
    $y_max = $y_min + 0.1;
   $pool_lines .=
     "set obj ".$obj_ctr." rectangle from ".(0.5).",".$y_min.
     " to ".($nr_res+0.5).",".$y_max." fc rgb \"black\" fs solid 0.25\n";
    next;
  }
  @suffixes = ( ".pdb", ".pool");
  $this_tag = fileparse( $this_file, @suffixes);
  $this_tag =~ s/${pdb_id}//;
  $this_tag =~ s/SSPred//;
  $tag4 = substr($pdb_id, 0, 4);
  $this_tag =~ s/$tag4//;
  $this_tag =~ s/^\.//g;
  if( length( $this_tag) == 0){ $this_tag = $pdb_id}
  print "file is: $this_file\n";
  print "tag is:  $this_tag\n";
  $y_min = 1.02 + ($nr_pools - $i - 1) * 0.1;
  $y_max = $y_min + 0.1;
  $is_native = $i < $nr_pools_native;

  $pool_lines .= "set label \"".$this_tag."\" at ".($nr_res+1).", ".($y_min + 0.06)." textcolor lt ".($is_native + 1)."\n";
  open( IN, $this_file) or die "can't open file #".$i." at ".$this_file."!\n";
  $border_offset = 0.25;
  while( $line = <IN>)
  {
    chomp( $line);
    if( $line =~ /^HELIX/)
    {
      $start= substr( $line, 22, 3);
      $end = substr( $line, 34, 3);
      $length = $end - $start + 1;
      $helix_type = substr( $line, 38 ,2);
      $color = "blue";
      # add the rectangle
      
      # if not a native pool and the length does not match the minimal then skip this definition
      if( !$is_native && $length < $min_helix)
      {next;}

      ++$obj_ctr;
      $pool_lines .=
        "set obj ".$obj_ctr." rectangle from ".($start - $border_offset).",".$y_min.
        " to ".($end + $border_offset).",".$y_max." fc rgb \"".$color."\" fs solid 0.8\n";
      # add label that represents the helix type
      if( $helix_type != 1)
      {
        $pool_lines .=
        "set label \"".$helix_type."\" at ".(($start + $end - 1.5)/2).", ".($y_min + 0.05)." textcolor lt 7\n";
      }
    }
    elsif( $line =~ /^SHEET/)
    {
      $start = substr( $line, 23,3);
      $end = substr( $line, 34, 3);
      $length = $end - $start + 1;

      # if not a native pool and the length does not match the minimal then skip this definition
      if( !$is_native && $length < $min_strand)
      {next;}
      ++$obj_ctr;
      $pool_lines .=
      "set obj ".$obj_ctr." rectangle from ".($start - $border_offset).",".$y_min.
      " to ".($end + $border_offset).",".$y_max." fc rgb \"red\" fs solid 0.8\n";
    }
  }
}

print "$nr_res residues found in sequence\n";
$size_x = 1500 * $nr_res / 100;
$size_y = 300;

$method_pngs_string; 

# iterate over each method
for( $x = 0; $x < $nr_methods; ++$x)
{
  #variables
  $method = $methods[ $x];
  print "plotting ".$method."\n";
  $extension = $method;
  $gnuplot_file = $out_dir.$pdb_id."_".$method.".gnuplot";
  $png_file = $out_dir.$pdb_id."_".$method.".png";
  $method_pngs_string .= $png_file." ";
  $data_file = $methods_dir.$pdb_id.".".$extension;
  $y_range = 1.02 + $nr_pools * 0.1;
  
  # open and output the gnuplot file
  open( OUT, ">".$gnuplot_file) or die "cannot open output file $gnuplot_file\n";
  print OUT "set terminal png font \"arialbd,12\" size ".$size_x.",".$size_y."\n";
  print OUT "set output \"$png_file\"\n";
  print OUT "set key outside vertical bottom right\n";
  # labels
  print OUT "set label \"".$pdb_id."\" at ".($nr_res+1.0).", 0.9 \n";
  print OUT "set label \"".$method."\" at ".($nr_res+1.0).", 0.8 \n";
  # ranges
  print OUT "set xrange [ 0.5: ".($nr_res + 0.5)." ]\n";
  print OUT "set yrange [ 0 : $y_range]\n";
  # add native lines
  print OUT $pool_lines;
  # xtics and ytics
  print OUT "set xtics ( \"0\" -0.5 ";
  for( $i = 1; $i <= int( $nr_res / 5); ++$i) {  $id = $i * 5; print OUT ", \"$id\" $id";  }
  print OUT ") mirror out\n";
  print OUT "set ytics (\"0\" 0, \"0.2\" 0.2, \"0.4\" 0.4, \"0.6\" 0.6, \"0.8\" 0.8, \"1.0\" 1.0)\n"; 
  # plot lines
  print OUT "plot \"".$data_file."\" u 1:5 title \"H\" w l lc rgb \"blue\" lw 2, ";
  print OUT     " \"".$data_file."\" u 1:6 title \"S\" w l lc rgb \"red\" lw 2, ";
  print OUT     " \"".$data_file."\" u 1:4 title \"C\" w l lc rgb \"green\" lw 2, ";
#  print OUT     " \"".$data_file."\" u 1:5:6 w filledcu notitle,";
  print OUT             " 0.5 w l lw 1.5 lt -1 notitle";
  for( $i = 0; $i < $nr_pools; ++$i)
  {
    $line_y = 1.02 + $i * 0.1;
    print OUT             ", ".$line_y." w l notitle lw 1.5 lt -1";
  }
  close OUT;
  
  # call gnuplot
  system( "gnuplot ".$gnuplot_file);
  $verbose || system( "rm -f $gnuplot_file"); # do not remove in verbose mode  
  #system( "convert -trim ".$png_file." ".$png_file);
}

system( "convert -append ".$method_pngs_string." ".$out_dir.$pdb_id.".png");
system( "rm -f ".$method_pngs_string);
