#!/usr/bin/perl

use Time::HiRes qw( gettimeofday tv_interval );
use threads;
use threads::shared;

# This script gets the list of available examples from bcl executable and runs them one by one
# @author Mert Karakas, 22 Oct 2008; Sten Heinze (Threading), Jeff Mendenhall (Threading; Improved output)

if($#ARGV != 1 && $#ARGV != 2 && $#ARGV != 3)
{
  die "usage: run_individual_examples.pl <bcl_nightly_build_results_folder> <bcl_executable_name> [max_examples_to_run_simultaneously] [Exclude apps]\nIf not given, max_examples_to_run_simultaneously defaults to 6";
}

## set maximum thread count
my $max_thread_count = 6;

## script start time
my $script_begin_time = gettimeofday();

## store the executable name
$results_folder=$ARGV[0];
$executable=$ARGV[1];
my $exclude_apps=0;
## if a maximum thread was given, use it
if( $#ARGV >= 3)
{
  $exclude_apps = int($ARGV[3]);
}

@examples_list = ();
share(@examples_list);
my( $command, $needed_virtual_screen, @exampleslist)=GetCommandAndExamples( $executable, $exclude_apps);
@examples_list = @exampleslist;

## if a maximum thread was given, use it
if( $#ARGV >= 2)
{
  $max_thread_count = int($ARGV[2]);
  if( $max_thread_count < 2)
  {
  	# always need at least two threads in the current setup; one to launch jobs, one to run them
  	$max_thread_count = 2;
  }
}

$total_run_time = 0;
share($total_run_time);
$total_run_time_with_startup = 0;
share($total_run_time_with_startup);
    
#print "EXAMPLES_LIST\n";

#number of errors encountered
$errors_found = 0;
share($errors_found);

if( scalar( @examples_list) == 0)
{
  print "Executable does not appear to have any examples!";
  $errors_found = -1;
}

%all_examples = ();
share(%all_examples);

$example_index=0;
share($example_index);

$error_string="";
share($error_string);

open( OUT, ">$results_folder/failed_examples.txt") or die "cannot open failed_examples.txt for writing!";

## iterate over examples
$num_threads = scalar( @examples_list) < $max_thread_count ? scalar( @examples_list) : $max_thread_count;
# create a number of threads 
for( $count = 0; $count < $num_threads; $count++)
{
  # create thread and test the result; stop creating more threads in case of an error
  my( $thread_process) = threads->create(\&run_example);
  if( ! $thread_process) 
  {
    printf( "Warning: Ran out of threads!\n");
    last;
  }
}

# wait until all threads are done  
foreach my $thr( threads->list()) 
{
  $thr->join();
}
if( $errors_found > 0)
{
  print OUT $error_string;
}
close OUT;

## prepare the output page
&prepare_output_page(\%all_examples, $script_begin_time);

print "Errors: ".$errors_found;

if( $needed_virtual_screen == 1) # have to kill the Xvfb instance that was started earlier
{
  $user_name = $ENV{'USER'};
  system("killall -u $user_name Xvfb");
}
          
open( OUT, ">$results_folder/example_errors.txt") or die "cannot open example_errors.txt for writing!";

print OUT $errors_found;

close OUT;

if( $errors_found > 0)
{
  &prepare_error_output_page(\%all_examples);
}
else
{
  open( OUT, ">$results_folder/example_errors.html") or die "cannot open example_errors.html for writing!";
  close OUT;
}

exit;

################################################################################
## subroutine for getting the next index out of the example
################################################################################
sub get_next_example
{
  my $next_example = "";
  
  {
    lock($example_index);
    if( $example_index < scalar(@examples_list))
    {
      $next_example = $examples_list[$example_index];
      print $next_example." ";
      $example_index += 1;
      if( !( $example_index % 6))
      {
        print "\n";
      }
    }
  }
  return $next_example;
}

################################################################################
## subroutine for executing one example command
################################################################################
sub run_example
{
  my $thread_error_string = "";
  my $thread_errors = 0;
  my $thread_runtime = 0;
  my $thread_runtime_with_startup = 0;
  
  for( $example_name = get_next_example(); length($example_name) > 0; $example_name = get_next_example())
  {
    ## record starting time
    my $begin_time = gettimeofday();
  
    my $example_fixed_name = $example_name;
    $example_fixed_name =~ s/</\\</g;
    $example_fixed_name =~ s/>/\\>/g;
    $output_filename=$results_folder."/output_".$example_fixed_name.".txt";
    $command_1 = "$command $example_fixed_name >& ".$output_filename;
  
    ## execute that example
    system($command_1);
    
    ## record ending time
    $end_time = gettimeofday();
    
    ## parse the output file, look for important lines
    open( IN, "<", $output_filename) or die "cannot open $output_filename";
    $time_in_seconds=-1;
    $reported_errors=-1;
    $max_virtual_memory=-1;
    $max_ram=-1;
    $found_end_examples=0;
    while( <IN>)
    {
      if( $found_end_examples == 0)
      {
        if( $_ =~ /BCL Example *\| END *: All Examples/) 
        {
          $found_end_examples=1;
        }
      }
      elsif( $_ =~ /^total /)
      {
        ( $reported_errors )= ( $_ =~ /^total *[0-9]* *([0-9]*)/);
      }
      elsif( $_ =~ /bcl has run for /)
      {
        ( $hours, $minutes, $seconds, $peakvm, $peakram) = ( $_ =~ /bcl has run for ([0-9]{2}):([0-9]{2}):([0-9]{2})(, peak virtual memory used: [0-9]* MB)?(, peak physical RAM used: [0-9]* MB)?/);
        $time_in_seconds = $seconds + 60 * ( $minutes + 60 * $hours);
        if( $peakvm)
        {
          ( $max_virtual_memory ) = ( $peakvm =~ /, peak virtual memory used: ([0-9]*)/);
        }
        if( $peakram)
        {
          ( $max_ram ) = ( $peakram =~ /, peak physical RAM used: ([0-9]*)/);
        }
      }
    }
    close IN;
  
    ## calculate interval in microseconds and update total run time
    $run_time = sprintf( "%.3f", ($end_time - $begin_time));
    $thread_runtime_with_startup += $run_time;
    if( $time_in_seconds == -1)
    {
      $time_in_seconds = $run_time;
    }
    $thread_runtime += $time_in_seconds;
    
    if( $reported_errors > 0)
    {
      $thread_error_string .= "$example_name failed $reported_errors example checks\n";
    }
    elsif( $reported_errors < 0)
    {
      $thread_error_string .= "$example_name crashed\n";
    }
    $thread_errors += abs( $reported_errors);
    
    @output = ( $time_in_seconds, $run_time, $reported_errors, $max_virtual_memory, $max_ram);
      
    {
      lock(%all_examples);
      $all_examples{ $example_name} = shared_clone([ @output]);
    }
  }
  
  if($thread_errors)
  {
    lock($errors_found);
    $errors_found += $thread_errors;
  }
  {
    lock($total_run_time);
    $total_run_time += $thread_runtime;
  }
  {
    lock($total_run_time_with_startup);
    $total_run_time_with_startup += $thread_runtime_with_startup;
  }
  if( length($thread_error_string))
  {
    lock($error_string);
    $error_string .= $thread_error_string;
  }
}

################################################################################
## subroutine for preparing the output html page text that includes
## the table for all examples their runtimes and outputs and error states etc.
################################################################################
sub prepare_output_page
{
  my(%all_examples) = %{$_[0]};
  my($start_time) = $_[1];
  
  # determine whether ram was determined for the architecture
  my $have_ram = 0;
  
  # determine whether vram was available
  my $have_vm = 0;
  my $min_vm_value = -1; # store the first value seen for virtual memory on the architecture
  my $max_vm_value = -1;
  my $min_ram_value = -1;
  my $max_ram_value = -1;
  foreach $example( sort keys %all_examples)
  {
    if( $all_examples{$example}[2] == 0)
    {
      if( $all_examples{$example}[4] != -1)
      {
        $have_ram = 1;
        $min_ram_value = $all_examples{$example}[4];
      }
      if( $all_examples{$example}[3] != -1)
      {
      	$have_vm = 1;
      	$min_vm_value = $all_examples{$example}[3];
      }
    }
  }
  
  # OS's often allocate large amount of base vram to the bcl and other large processes (partly a result of libraries) 
  # So take the minimum vram and subtract it from all reported figures
  if( $have_vm > 0)
  {
  	$max_vm_value = $min_vm_value;
  	foreach $example( keys %all_examples)
    {
      if( $all_examples{$example}[2] == 0)
      {
      	$vm = $all_examples{$example}[3];
        if( $vm < $min_vm_value)
        {
          $min_vm_value = $vm;
        }
        elsif( $vm > $max_vm_value)
        {
          $max_vm_value = $vm;
        }
      }
    }
  }
  
  # OS's often allocate some base ram to the bcl (to load the executable into memory)
  # So take the minimum ram and subtract it from all reported figures
  if( $have_ram > 0)
  {
  	$max_ram_value = $min_ram_value;
    foreach $example( keys %all_examples)
    {
      if( $all_examples{$example}[2] == 0)
      {
        $ram = $all_examples{$example}[4];
        if( $ram < $min_ram_value)
        {
          $min_ram_value = $ram;
        }
        elsif( $ram > $max_ram_value)
        {
          $max_ram_value = $ram;
        }
      }
    }
  }
  
  ## open a php file for writing
  open( OUT, ">$results_folder/example_log.html") or die "cannot open example_log
        ( $hours, $minutes, $seconds, $peakvm, $peakram) = ( $_ =~ /bcl has.html for writing!";
  
  ## get the formatted day and time
  ($sec,$min,$hour,$mday,$mon,$year)=localtime(time);
  
  $today= sprintf( "%4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec);

  print OUT "<HTML>\n";
  # make the table sortable if the sorttable javascript is available
  print OUT "<HEAD><script src=\"sorttable.js\"></script></HEAD>";
  print OUT "<TITLE> Example run results from ".$today."</TITLE>\n";
  print OUT
 '
<BODY>
        <DIV STYLE="width: 850px; overflow: auto;">
';


  print OUT "<b>".$today." EXAMPLE RUN RESULTS"."</b><br>\n";
 
  ## record ending time
  $script_run_time = sprintf( "%.3f", (gettimeofday() - $start_time));
 
  # print out the total run time information and the number of failed tests
  print OUT "Total run time for all examples: ".$total_run_time." seconds<BR>\n";
  print OUT "Total run time for all examples, counting executable startup: ".$total_run_time_with_startup." seconds<BR>\n";
  print OUT "Script runtime: ".$script_run_time." seconds<BR>\n";
  if( $have_vm == 1)
  {
  	print OUT "Base virtual memory: ".$min_vm_value." MB<BR>\n";
  	print OUT "Max additional virtual memory used: ".( $max_vm_value - $min_vm_value)." MB<BR>\n";
  }
  if( $have_ram == 1)
  {
  	print OUT "Base physical ram: ".$min_ram_value." MB<BR>\n";
  	print OUT "Max extra physical ram used: ".( $max_ram_value - $min_ram_value)." MB<BR>\n";
  }
  print OUT "Total number of failed tests: ".$errors_found."<BR>\n";
  print "To sort columns in the example_log.html file, you will need to copy sorttable.js into the same directory\n";
  print "sorttable.js can be obtained from http://www.kryogenix.org/code/browser/sorttable/\n";
  print "sorttable.js is also available at /blue/meilerlab/apps/Linux2/noarch/shell_scripts/bcl_nightly_build_and_check/sorttable.js\n";
  print OUT "Click a column to sort<BR>\n";

  print OUT "<TABLE CELLSPACING=3 CELLPADDING=10 STYLE=\"font-size: 8pt;\" CLASS=\"sortable\" >\n";
  print OUT "  <THEAD BGCOLOR=#D8D8D8\">\n";
  print OUT "    <TH> Example Name </TH><TH> Runtime(s)</TH><TH> Runtime with startup</TH>";
  if( $have_vm == 1)
  {
    print OUT "<TH> Virtual Memory Increase (MB)</TH>";
  }
  if( $have_ram == 1)
  {
    print OUT "<TH> Physical RAM Increase Over Base (MB)</TH>";
  }
  print OUT "<TH> Status</TH><TH> # Failed Tests</TH><TH> Output</TH>\n";
  print OUT "  </THEAD>\n";
  print OUT "  <TBODY>\n";

  foreach $example( keys %all_examples)
  {
    $ex_run_time = $all_examples{$example}[0];
    $ex_run_time_with_startup = $all_examples{$example}[1];
    $ex_error= $all_examples{$example}[2];
    $max_vm= $all_examples{$example}[3] - $min_vm_value;
    $max_ram= $all_examples{$example}[4] - $min_ram_value;
    
    ## if the example took to much time to run put a different background color
    if( $ex_error)
    {
      $max_vm = $max_ram = 0;
      print OUT "  <TR BGCOLOR=#FF0000>\n ";
    }
    elsif( $ex_run_time >= 10)
    {
      print OUT "  <TR BGCOLOR=#006666>\n ";
    }
    else
    {
      print OUT "  <TR BGCOLOR=#E8E8E8>\n ";
    }
    
    $error_description="";
    if( $ex_error == 0)
    {
      $error_description="passed";
    }
    elsif( $ex_error == -1)
    {
      $error_description="crashed";
      $ex_error=1;
    }
    else
    {
      $error_description="failed";
    }
    
    print OUT "    <TD> ".$example."</TD>\n";
    print OUT "    <TD> ".$ex_run_time."</TD>\n";
    print OUT "    <TD> ".$ex_run_time_with_startup."</TD>\n";
    if( $have_vm == 1)
    {
      print OUT "    <TD>".$max_vm."</TD>\n";
    }
    if( $have_ram == 1)
    {
      print OUT "    <TD>".$max_ram."</TD>\n";
    }
    print OUT "    <TD> ".$error_description." </TD>\n";
    print OUT "    <TD> ".$ex_error."</TD>\n";
    print OUT "    <TD><A HREF=\""."output_".$example.".txt"."\">output</A></TD>\n";
    print OUT "  </TR>\n";

  }
  print OUT "</TBODY>\n";
  print OUT "</TABLE>\n";

  print OUT
'

        </DIV>
</BODY>
</HTML>
                               
';
  close OUT;
}

################################################################################
## subroutine for preparing the output html page text that includes
## the table for failed examples and outputs and error states etc.
################################################################################
sub prepare_error_output_page
{
  my(%all_examples) = %{$_[0]};
  ## open a php file for writing
  open( OUT, ">$results_folder/example_errors.html") or die "cannot open example_errors.html for writing!";

  print OUT "<table cellspacing=3 style=\"font-size: 10pt;\" width=400px>\n";
  print OUT "  <tbody>\n";
  $number_examples=scalar keys %all_examples;
  if( $number_examples == 0)
  {
    # handle the case where no examples are present
    print OUT "<tr><td>Executable was not built or could not be run</td></tr>";
  }
  else
  {
    # handle the case where there were examples, some of which failed
  foreach $example( sort keys %all_examples)
  {
    $ex_error= $all_examples{$example}[2];
    if( $ex_error != 0) ## do not print the total or completed tests
    {
      # compose an error string to represent the problem
      $error_string="";
      if( $ex_error == -1)
      {
        # crash
        $error_string = "crashed";
      }
      else
      {
        # example check failure
        $error_string = "failed ".( $ex_error == 1 ? "1 example check " : $ex_error." example checks");
      }
      print OUT "<tr><td><a href=\"output_".$example.".txt"."\">".$example."</a>: ";
      print OUT $error_string."</td></tr>\n";
    }
    }
  }
  print OUT "</tbody>\n";
  print OUT "</table>\n";
  close OUT;
}

#####################################################################################
## subroutine for initialization; determines proper command line, opens screen buffer
## if necessary, and returns a list of examples
#####################################################################################
sub GetCommandAndExamples
{
  my $executable= $_[0];
  my $exclude_apps = $_[1];
  my $command=$executable." Examples -exec ";

  ## run the bcl example with a bogus example name to get list of examples
  my $error_number=system("$command -help >& examples.out");
  
  ## flag to indicate whether we needed to start a virtual screen because the machine running this command was headless
  ## if we started such a virtual screen, we need to kill it before closing
  my $needed_virtual_screen=0;
  
  if( $error_number != 0)
  {
    # may have been the wrong architecture.  Look in the output file for /architecture/ or /format/
    # Try prepending the command with wine to run the executable
    open( IN, "examples.out") or die "cannot open examples.out";
    $wrong_architecture=0;
    while( <IN>)
    {
      if( $_ =~ /archictecture/ || $_ =~ /format/ || $_ =~ /binary/) 
      {
        $wrong_architecture=1;
      }
    }
    close IN;
    if( $wrong_architecture == 1)
    {
      print "Running executable without wine failed, trying wine\n";
      # try to run under wine
      $command="wine ".$executable." Examples -exec ";
      $wrong_architecture=system("$command -help >& examples.out");
      if( $wrong_architecture != 0)
      {
        print "Running executable with wine but without a virtual screen buffer failed, making a virtual screen buffer\n";
        # wine cannot run on a headless machine if X11 connections have not been forwarded
        # which won't happen in the nightly build, so we need to create a virtual screen buffer
        
        # first kill any existing virtual screen buffers for this user (can't have multiple virtual screen buffers with the same id)
        $user_name=$ENV{'USER'};
        print "Killing any virtual screen buffers already running for $user_name\n";
        system("killall -u $user_name Xvfb; sleep 1"); # also sleep for a second to give the currently running instances time to close
        
        # now make the virtual screen buffer
        $had_display=system("Xvfb :8 -fp /usr/share/X11/fonts/misc -screen 0 800x640x16 &");
        if( $had_display == 0)
        {
          $needed_virtual_screen=1;
          print "Making a virtual screen buffer succeeded, exporting display\n";
          #export the screen buffer as the display
          $command="export DISPLAY=localhost:8; wine ".$executable." Examples -exec ";
          $wrong_architecture=system("$command -help >& examples.out");
        }
        else
        {
          #the command failed, most likely Xvfb is not installed
          print "Making a virtual screen buffer failed, please ensure that Xvfb is installed and in the default path\n";
        }
        if( $wrong_architecture != 0)
        {
          system( "cat examples.out");
          unlink( "examples.out");
          if( $needed_virtual_screen == 1) # have to kill the Xvfb instance that was started earlier
          {
            system("killall -u $user_name Xvfb");
          }
          die $executable." appears to be in the wrong architecture for this computer and could not be run under wine, even with a virtual screen buffer";
        }
      }
      else
      {
        print "Succeeded.  Running examples under wine\n";
      }
    }
  }
  
  ## open the examples out for parsing
  open( IN, "examples.out") or die "cannot open examples.out";
  
  ## while reading
  while( <IN>)
  {
    ## if the line with list of examples is found
    if( $_ =~ /<example>/)
    {
      ## store it and break
      $examples_line = $_;
      $examples_line =~ s/^\s+//;
      $examples_line =~ s/\s+$//;
      $examples_line .= ' ';
      if( !( $_ =~ /}/))
      {
        while( <IN>)
        {
          $trimmed = $_;
          $trimmed =~ s/^\s+//;
            $trimmed =~ s/\s+$//;
            $trimmed .= ' ';
          $examples_line .= $trimmed;
          if( $_ =~ /}/)
          {
            last;
          }
        }
      }
    }
  }
  ## close filestream
  close IN;
  unlink( "examples.out");
  
  ## remove new lines
  $examples_line =~ tr/\n//d;
  
  ## now parse the examples line
  $examples_line =~ s/.*\{(.+)\}.*/$1/;
  
  ## remove commas
  $examples_line =~ tr/,//d;
  
  ## split by whitespace
  my @exampleslist = split(/\W/,$examples_line,0);

  if( $exclude_apps > 0)
  {
    @exampleslist_no_apps = ();
    
    foreach $example (@exampleslist)
    {
      if( !( $example =~ /^App[A-Z]/))
      {
        push( @exampleslist_no_apps, $example);
      }
    }
    return ( $command, $needed_virtual_screen, @exampleslist_no_apps);
  }
  return ( $command, $needed_virtual_screen, @exampleslist);
}
