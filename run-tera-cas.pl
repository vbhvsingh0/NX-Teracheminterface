#!/usr/bin/env perl
#
#====================================================================================
#
# This file is a Perl template to build an interface between 
# Newton-X CS and a third-party program (TPP) to run surface 
# hopping dynamics.
#
#====================================================================================
#
use strict;
use warnings;
use File::Copy;
use File::Path;
use lib join('/',$ENV{"NX"},"lib");
use colib_perl;

# We are using "strict." Thus, all variables mus be declared.
my ($mld,$mdle,$BASEDIR,$DEBUG,$JAD,$JND,$ctd,$fpar);
my ($scr,$scr1,$scr2,$grad,$out);
my ($nat,$istep,$nstat,$nstatdyn,$ndamp,$kt,$dt,$t,$tmax,$nintc);
my ($nxrestart,$thres,$killstat,$timekill,$prog,$lvprt,$etot_jump,$etot_drift);
my (@symb,@zn,@x,@y,@z,@Mass,$ia);
my ($typeofinput,$type);
my (@epot,$ns);
my ($vdoth,$run_flag);
my ($ncoup,$nc);
my (@gx,@gy,@gz,@hx,@hy,@hz);
my ($kross_d,$cascade_d,$current_d,$never_state_d,$include_pair_d);
my ($e_ci_d,$ci_cons_d,$cio_options_d,$cisc_options_d,$idalton_d);
my ($cprog_d,$coptda_d,$ncore_d,$ndisc_d,$blasthread_d);
my ($blasthread,$CIO_OPTIONS,$CISC_OPTIONS,$CPROG);
my (%progconf);
#Declaration of some variables which were not done earlier:
my ($mem,$keyword1_d,$keyword2_d,$keyword1,$keyword2,$file,$retval,$mldcio);

# Define basic variables
define_variables();

# Write debug messages
write_debug1();

# Read current geometry
read_geom();

# Read TPP parameters
#read_fpar();  

# Prepare TPP input
prepare_input();

# Run TPP
run_program();

# Read and write energies
read_energy();

# Read and write gradients
read_gradients();

# Nonadiabatic couplings
treat_nacme();

#
#====================================================================================
#
# 				START SUBROUTINES
#
#====================================================================================
#

sub define_variables{
#
#====================================================================================
#
# Define basic variables to be used everywhere.
#
#------------------------------------------------------------------------------------
#

  $mld    = $ENV{"NX"};
  $mdle   = "run-tera-cas.pl:";  
  # There are 3 scr.geom files because subroutine callprog doesn't work if directory is changed and hence terachem scr directories are saved in TEMP rather than in GRAD, NAD_N1N2, NAD_N1N0! 
  $scr    = "scr.geom";      #this will contain GRAD
  $scr1    = "scr.geom.1";    #this will contain NAD_N1N0, except nstatdyn == 1, here it will contain NAD_N1N2
  $scr2    = "scr.geom.2";    #this will contain NAD_N1N2
  $grad   = "grad.xyz";
  $out    = "geom.out";
  $BASEDIR=`pwd`;
  chomp ($BASEDIR);

  $DEBUG  = "DEBUG";
  $JAD    = "JOB_AD";
  $JND    = "JOB_NAD";
  $ctd    = "control.d";

  # Load current dynamics status (read control.d)
  ($nat,$istep,$nstat,$nstatdyn,$ndamp,$kt,$dt,$t,$tmax,$nintc,
   $mem,$nxrestart,$thres,$killstat,$timekill,$prog,$lvprt,
   $etot_jump,$etot_drift)=load_status($ctd,$mdle); 

  %progconf = prog_config($prog);
  $fpar = $progconf{parfile};  
}

sub write_debug1{
#
#====================================================================================
#
# Write debug messages.
#
#------------------------------------------------------------------------------------
#
  if ($lvprt>=3) {print_STDOUT("$mdle has taken over\n",$istep,$kt);}
  if ($lvprt>=3) {print_STDOUT("$mdle running in $BASEDIR\n",$istep,$kt);}
  if ($lvprt>=3) {print_STDOUT("$mdle beginning here\n",$istep,$kt);}
  if ($lvprt>=3) {
     print_STDOUT("$mdle Dynamics control: \n",$istep,$kt);
     print_STDOUT("$nat $istep $nstat $nstatdyn $ndamp $kt $dt $t $tmax \n",$istep,$kt);
     print_STDOUT("$nintc $mem $nxrestart $thres $killstat $timekill $prog $lvprt\n",$istep,$kt);
     print_STDOUT("$etot_jump $etot_drift\n\n",$istep,$kt);
  }
}

sub read_geom{
#
#====================================================================================
#
# This subroutine reads the current geometry and populates the
# following vectors:
# @symb      - atomic symbol
# @zn        - atomic number
# @x, @y, @z - cartesian coordinates (au)
# @Mass      - atomic mass (amu)
# The vector index runs from $ia = 0 to $ia = $nat-1.
# Thus, coordinate x of forst atom in in $x[0], for instance.
#
#------------------------------------------------------------------------------------
#

  if ($lvprt>=3) {print_STDOUT("$mdle Reading current geometry.\n",$istep,$kt);}

  open(GM,"geom") or die "$mdle Cannot open geom";
  $ia = 0;
  while(<GM>){
    chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
    ($symb[$ia],$zn[$ia],$x[$ia],$y[$ia],$z[$ia],$Mass[$ia])=split(/\s+/,$_);
    $ia++;
  }
  close(GM);
}

#sub read_fpar{
#
#
# Any information about the TPP can be given through $fpar file, which is
# defined in the %progconf hash colibperl,pm.
#
# For instance, with Gaussian, $fpar="gau.par". We can tel NX the Gaussian version 
# to be used through the variable $g_vers. Thus, if gau.par contains the line
# g_vers = 16
# NX will read it through
# $g_vers = getkeyword($fpar,"g_vers",$g_vers_d);
# where $g_vers_d is the default value of $g_vers.
#
# The commands below will read keyword1 and keyword2.
# You should change to appropriate variable names.
# You can add as manys variables as you need.
#
# 1. Declare the variables in the beggining of the code.
# 2. Define the default values in the load_defaults_tpp() subroutine
#    in lib/colib_perl.pm.
# 3. Document the variables.
#
#------------------------------------------------------------------------------------
#
 
#  if ($lvprt>=3) {print_STDOUT("$mdle Reading parameters from $fpar.\n",$istep,$kt);}

#  ($keyword1_d,$keyword2_d) = load_defaults_tpp($prog,"","");

#  $keyword1 = getkeyword($fpar,"keyword1",$keyword1_d);
#  $keyword2 = getkeyword($fpar,"keyword2",$keyword2_d);

#}

sub prepare_input{
#
# This subroutine prepares the TPP input.
# Template files (given in JOB_AB/ or JOB_NAD/) were already
# copied to $BASEDIR/. Now, you should change these files to
# run TPP with the proper geometry and other options.
#
# Remember:
# 1. The current geometry is stored in the vectors
# @symb,@zn,@x,@y,@z,@Mass.
#
# 2. If you need to convert the geometry into Angstrom,
# the variable $au2ang will be useful.
# ($au2ang = 0.52917720859)
#
# 3. NX instructions to TPP are stored in the variables
# read by the subroutine read_fpar.
#
# 4. All specific info about the dynamics are in the
# variables read by load_status. Some of those that you
# may need are:
# $nstat    - number of state
# $nstatdyn - current state
#
#------------------------------------------------------------------------------------
#

  my $au2ang = units("au2ang");

  if ($lvprt>=3) {print_STDOUT("$mdle Preparing input for third-party program.\n",$istep,$kt);}

  find_typeofinput();

  if ($typeofinput eq "jad"){
	# Read the original XYZ file content
	my $input_filename = "geom.xyz";  # Replace with your input file name
	open my $input_fh, '<', $input_filename or die "Could not open $input_filename: $!";
	my @input_lines = <$input_fh>;
	close $input_fh;
	# Create the new XYZ file with updated coordinates
	my $output_filename = "geom.xyz";  # Replace with your output file name
	open my $output_fh, '>', $output_filename or die "Could not open $output_filename: $!";

	# Write the total number of atoms
	print $output_fh scalar(@symb), "\n";

	# Write a blank line
	print $output_fh "\n";

	# Write the updated atomic symbols and coordinates
	for my $i (0..$#symb) {
    	print $output_fh "$symb[$i]   ".$x[$i]*$au2ang."   ".$y[$i]*$au2ang."   ".$z[$i]*$au2ang."\n";
	}

	close $output_fh;

  }
  elsif($typeofinput eq "jnd"){
    # Prepare input for nonadiabatic run (energy + gradients + NACME like in JOB_NAD/)
    # If TPP does not provide NACME, you do not need to enter anything here.
	if ($nstatdyn == 1) {
		`cp NAD_N1N2/geom.xyz .`; # or die "Is geom file saved as geom.xyz: $!";
		my $input_filename = "geom.xyz";  # Replace with your input file name
        	open my $input_fh, '<', $input_filename or die "Could not open $input_filename: $!";
        	my @input_lines = <$input_fh>;
        	close $input_fh;
        # Create the new XYZ file with updated coordinates
        	my $output_filename = "geom.xyz";  # Replace with your output file name
        	open my $output_fh, '>', $output_filename or die "Could not open $output_filename: $!";

        # Write the total number of atoms
        	print $output_fh scalar(@symb), "\n";

        # Write a blank line
        	print $output_fh "\n";

        # Write the updated atomic symbols and coordinates
        	for my $i (0..$#symb) {
        	print $output_fh "$symb[$i]   ".$x[$i]*$au2ang."   ".$y[$i]*$au2ang."   ".$z[$i]*$au2ang."\n";
        }

        	close $output_fh; 
		`cp geom.xyz NAD_N1N2/.`;
		`cp geom.xyz GRAD/.`;
		chdir('NAD_N1N2') or die "Could not change directory:NAD_N1N2 $!";
		# Open the 'start.sp' file for reading
		open my $input, '<', 'start.sp' or die "Could not open 'start.sp' for reading: $!";
		# Read the file content into an array
		my @file_content = <$input>;
		# Close the input file
		close $input;

		# Open the 'start.sp' file for writing
		open my $output, '>', 'start.sp' or die "Could not open 'start.sp' for writing: $!";
		# Iterate through the file content, update nacstate1 and nacstate2
		foreach my $line (@file_content) {
    			if ($line =~ /^nacstate1\s+(\d+)/) {
        			$line = "nacstate1               " . ($nstatdyn - 1) . "\n";
    			} 
			elsif ($line =~ /^nacstate2\s+(\d+)/) {
        			$line = "nacstate2               $nstatdyn\n";
    			}
			# Write the updated line to the output file
    			print $output $line;
		}
		close $output;
		chdir('../GRAD') or die "Could not change directory:GRAD $!";
		# Open the 'start.sp' file for reading
                open my $input1, '<', 'start.sp' or die "Could not open 'start.sp' for reading: $!";
		# Read the file content into an array
                my @file_content1 = <$input1>;
                # Close the input file
                close $input1;

		# Open the 'start.sp' file for writing
                open my $output1, '>', 'start.sp' or die "Could not open 'start.sp' for writing: $!";
                # Iterate through the file content, update current state
		foreach my $line (@file_content1) {
			if ($line =~ /^castarget\s+(\d+)/) {
                                $line = "castarget               " . ($nstatdyn - 1) . "\n";
                        }
			# Write the updated line to the output file
                	print $output1 $line;
		}
		close $output1;
		chdir('../');
	}
	elsif ($nstatdyn > 1 and $nstatdyn != $nstat) {
                `cp NAD_N1N2/geom.xyz .` ;
		my $input_filename = "geom.xyz";  # Replace with your input file name
        	open my $input_fh, '<', $input_filename or die "Could not open $input_filename: $!";
        	my @input_lines = <$input_fh>;
        	close $input_fh;
        # Create the new XYZ file with updated coordinates
        	my $output_filename = "geom.xyz";  # Replace with your output file name
        	open my $output_fh, '>', $output_filename or die "Could not open $output_filename: $!";

        # Write the total number of atoms
        	print $output_fh scalar(@symb), "\n";

        # Write a blank line
        	print $output_fh "\n";

        # Write the updated atomic symbols and coordinates
        	for my $i (0..$#symb) {
        	print $output_fh "$symb[$i]   ".$x[$i]*$au2ang."   ".$y[$i]*$au2ang."   ".$z[$i]*$au2ang."\n";
        }

        	close $output_fh;
	        `cp geom.xyz NAD_N1N0/.`;	
		`cp geom.xyz NAD_N1N2/.`;
                `cp geom.xyz GRAD/.`;
		#`cp NAD_N1N2/start.sp NAD_N1N0/.`; #This is to make sure start.sp is present in NAD_N1N0
		chdir('NAD_N1N0') or die "Could not change directory:NAD_N1N0 $!";
		# Open the 'start.sp' file for reading
                open my $input, '<', 'start.sp' or die "Could not open 'start.sp' for reading: $!";
                # Read the file content into an array
                my @file_content = <$input>;
                # Close the input file
                close $input;

                # Open the 'start.sp' file for writing
                open my $output, '>', 'start.sp' or die "Could not open 'start.sp' for writing: $!";
                # Iterate through the file content, update nacstate1 and nacstate2
                foreach my $line (@file_content) {
                        if ($line =~ /^nacstate1\s+(\d+)/) {
                                $line = "nacstate1               " . ($nstatdyn - 2) . "\n";
                        }
                        elsif ($line =~ /^nacstate2\s+(\d+)/) {
                                $line = "nacstate2               " . ($nstatdyn - 1) . "\n";
                        }
                        # Write the updated line to the output file
                        print $output $line;
                         }
                close $output;
                chdir('../NAD_N1N2') or die "Could not change directory:NAD_N1N2 $!";
                # Open the 'start.sp' file for reading
                open my $input1, '<', 'start.sp' or die "Could not open 'start.sp' for reading: $!";
                # Read the file content into an array
                my @file_content1 = <$input1>;
                # Close the input file
                close $input1;

                # Open the 'start.sp' file for writing
                open my $output1, '>', 'start.sp' or die "Could not open 'start.sp' for writing: $!";
                # Iterate through the file content, update nacstate1 and nacstate2
                foreach my $line (@file_content1) {
                        if ($line =~ /^nacstate1\s+(\d+)/) {
                                $line = "nacstate1               " . ($nstatdyn - 1) . "\n";
                        }
                        elsif ($line =~ /^nacstate2\s+(\d+)/) {
                                $line = "nacstate2               $nstatdyn\n";
                        }
                        # Write the updated line to the output file
                        print $output1 $line;	
			 }
                close $output1;
		chdir('../GRAD') or die "Could not change directory:GRAD $!";
                # Open the 'start.sp' file for reading
                open my $input2, '<', 'start.sp' or die "Could not open 'start.sp' for reading: $!";
                # Read the file content into an array
                my @file_content2 = <$input2>;
                # Close the input file
                close $input2;

                # Open the 'start.sp' file for writing
                open my $output2, '>', 'start.sp' or die "Could not open 'start.sp' for writing: $!";
                # Iterate through the file content, update current state
                foreach my $line (@file_content2) {
                        if ($line =~ /^castarget\s+(\d+)/) {
                                $line = "castarget               " . ($nstatdyn - 1) . "\n";
                        }
                        # Write the updated line to the output file
                        print $output2 $line;
                }
                close $output2;
		chdir('../');
	}
	elsif ($nstatdyn > 1 and $nstatdyn == $nstat) {
		#`cp NAD_N1N0/geom.xyz .` or die "Is geom file saved as geom.xyz: $!";
		`cp NAD_N1N0/geom.xyz .` ;
		my $input_filename = "geom.xyz";  # Replace with your input file name
                open my $input_fh, '<', $input_filename or die "Could not open $input_filename: $!";
                my @input_lines = <$input_fh>;
                close $input_fh;
        # Create the new XYZ file with updated coordinates
                my $output_filename = "geom.xyz";  # Replace with your output file name
                open my $output_fh, '>', $output_filename or die "Could not open $output_filename: $!";

        # Write the total number of atoms
                print $output_fh scalar(@symb), "\n";

        # Write a blank line
                print $output_fh "\n";

        # Write the updated atomic symbols and coordinates
                for my $i (0..$#symb) {
                print $output_fh "$symb[$i]   ".$x[$i]*$au2ang."   ".$y[$i]*$au2ang."   ".$z[$i]*$au2ang."\n";
	        }
		close $output_fh;
		`cp geom.xyz NAD_N1N0/.`;
		`cp geom.xyz GRAD/.`;
		#`cp NAD_N1N2/start.sp NAD_N1N0/.` ;#This is to make sure start.sp is present in NAD_N1N0
		chdir('NAD_N1N0') or die "Could not change directory:NAD_N1N0 $!";
		open my $input, '<', 'start.sp' or die "Could not open 'start.sp' for reading: $!";
		my @file_content = <$input>;
		close $input;
		open my $output, '>', 'start.sp' or die "Could not open 'start.sp' for writing: $!";
		foreach my $line (@file_content) {
                        if ($line =~ /^nacstate1\s+(\d+)/) {
                                $line = "nacstate1               " . ($nstatdyn - 2) . "\n";
                        }
                        elsif ($line =~ /^nacstate2\s+(\d+)/) {
                                $line = "nacstate2               " . ($nstatdyn - 1) . "\n";
                        }
                        # Write the updated line to the output file
                        print $output $line;
                         }
                close $output;
		chdir('../GRAD') or die "Could not change directory GRAD: $!";
                # Open the 'start.sp' file for reading
                open my $input1, '<', 'start.sp' or die "Could not open 'start.sp' for reading: $!";
                # Read the file content into an array
                my @file_content1 = <$input1>;
                # Close the input file
                close $input1;

                # Open the 'start.sp' file for writing
                open my $output1, '>', 'start.sp' or die "Could not open 'start.sp' for writing: $!";
                # Iterate through the file content, update current state
                foreach my $line (@file_content1) {
                        if ($line =~ /^castarget\s+(\d+)/) {
                                $line = "castarget               " . ($nstatdyn - 1) . "\n";
                        }
                        # Write the updated line to the output file
                        print $output1 $line;
                }
                close $output1;
                chdir('../');
	}
  }

}

sub find_typeofinput{
#
#====================================================================================
#
# Check type of dynamics.
#
#------------------------------------------------------------------------------------
#
 
  if ($lvprt>=3) {print_STDOUT("$mdle Checking type of dynamics.\n",$istep,$kt);}
 
  # Read $type from type_of_dyn.out 
  my $typeout="type_of_dyn.out";
  $type = 2;
  if (-s $typeout){
    open(INP,$typeout) or die "$mdle Cannot open $typeout";
    $_=<INP>;
    chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
    my $type = $_;
    close(INP);
  }

  # Define type of input
  if ((-e $JAD) and (!-e $JND)){
    $typeofinput = "jad";
  }elsif((!-e $JAD) and (-e $JND)){
    $typeofinput = "jnd";
  }elsif((-e $JAD) and (-e $JND)){
    if ($type == 1){
      $typeofinput = "jad";
    }elsif($type != 1){
      $typeofinput = "jnd";
    }
  }

}


#the subroutine below will be called for running program
sub process_directory {
    my $tpp = $ENV{"TeraChem"};
    my ($directory) = @_;
    use Cwd;

    #my $current_dir = getcwd();  used to DEBUG=UTD
    #my $dir1 = "$current_dir/$directory"; used to DEBUG=UTD
    #print_STDOUT("$current_dir\n"); #UTD
    #chdir($dir1) or die "Could not change directory: $directory $!";#uncomment of you want to use this. was giving error!! #UTD
    #my $current_dir1 = getcwd(); #UTD
    #print_STDOUT("$current_dir1\n"); #UTD

    #opendir(my $dh, '.') or die "Could not open current directory: $!";
    #while (my $item = readdir($dh)) {
    #	next if $item =~ /^\./;  # Skip . and ..
    #	print_STDOUT("$item\n");
    #	}
    #closedir($dh);
    # Remove old results containing folder
    #if ($istep != 0){
    #	rmtree($scr, { safe => 1 }) or die "No $scr directory found: $!";
    #	rmtree($scr1, { safe => 1 }) or die "No $scr directory found: $!";
    #	rmtree($scr2, { safe => 1 }) or die "No $scr directory found: $!";}
    #system("rm -rf $scr") or die "No $scr directory found: $!";

    # Remove old geom.out file
    # system("rm geom.out"); # Uncomment this line if needed

    # Run the third-party program
    `rm $directory/$out`;
    $retval = callprog($tpp, "/bin/terachem $directory/start.sp >> $directory/$out", $mdle);

    if ($retval != 0) {
        print_STDOUT("Error in the third-party program execution!\n", $istep, $kt);
    }


    #chdir($current_dir);
}

sub run_program{
# 
#====================================================================================
#
# This subroutine executes the TPP.
# It uses the subroutine callprog("path","execution command",$mdle).
#
# Let me give an example. Columbus is executed with the command line:
# $COLUMBUS/runc -m 1600 > runls
# With callprog, we do:
# $retval = callprog(\$COLUMBUS, "runc -m 1600 > runls", $mdle);
#
# If the execution is fine, $retval = 0.
# Note the \ in \$COLUMBUS. If you don't add it, $COLUMBUS will be
# interpreted as a perl variable, not as an environment path.
# Alternatively, you can load the environment path into a perl variable.
# In the case of COLUMBUS, for example:
# my $cbus = $ENV{"COLUMBUS"};
# $retval = callprog($cbus, "runc -m 1600 > runls", $mdle);
#
#------------------------------------------------------------------------------------
#
   
  if ($lvprt>=3) {print_STDOUT("$mdle Executing third-party program.\n",$istep,$kt);}

  my $tpp = $ENV{"TeraChem"}; 
  # I am supposing that TPP has different calls for when computing 
  # NACME (jnd) or not (jad). If it is not the case, you can simplify 
  # the code below, removing the conditional.
  #
  if ($typeofinput eq "jad"){
	  #system("rm -rf scr.geom") or die "No $scr directory found: $!"
    #remove old results containing folder 
    rmtree($scr, { safe => 1 }) ;
    `rm $out`;
    #system("rm geom.out"); 
    $retval = callprog($tpp, "/bin/terachem start.sp >> $out", $mdle);
    if ($retval != 0){print_STDOUT("Error in the third-party program execution!\n",$istep,$kt);}
  }
  if($typeofinput eq "jnd"){
	if ($nstatdyn == 1) {
		if ($istep != 0){
               		`rm -rf scr.geom*`; }
		process_directory('NAD_N1N2');
		#`cp -rf scr.geom NAD_N1N2/.` or die "No $scr directory found: $!";
		#system("rm -rf scr.geom");
		process_directory('GRAD');
		#`cp -rf scr.geom GRAD/.` or die "No $scr directory found: $!";
		#system("rm -rf scr.geom");
	}
	elsif ($nstatdyn > 1) {
	   if ($nstatdyn != $nstat){
	   	 if ($istep != 0){
                        `rm -rf scr.geom*`; }
		process_directory('NAD_N1N0');
		#`cp -rf scr.geom NAD_N1N0/.` or die "No $scr directory found: $!";
		#system("rm -rf scr.geom");
		process_directory('NAD_N1N2');
		#`cp -rf scr.geom NAD_N1N2/.` or die "No $scr directory found: $!";
		#system("rm -rf scr.geom");
                process_directory('GRAD');
		#`cp -rf scr.geom GRAD/.` or die "No $scr directory found: $!";
		#system("rm -rf scr.geom");
           }
	   elsif ($nstatdyn == $nstat){
		    if ($istep != 0){
                        `rm -rf scr.geom*`; }
		process_directory('NAD_N1N0');
		#`cp -rf scr.geom NAD_N1N0/.` or die "No $scr directory found: $!";
		#system("rm -rf scr.geom");
                process_directory('GRAD');
		#`cp -rf scr.geom GRAD/.` or die "No $scr directory found: $!";
		#system("rm -rf scr.geom");
           }		   
	}
  }
  # The commandline for JAD or JNAD are the same now, but however, can be changed later on, if needed.
  # You may consider adding a subroutine here to check convergence and other possible 
  # issues in the TPP execution.
  #
  # check_for_problems();

}

sub read_energy{
#
#====================================================================================
#
# This subroutine reads TPP output to get the potential energies.
# These energies are written to the vector @epot, with indexes
# running from $ns = 0 to $nstat-1. Thus, the ground state energy
# should be in $epot[0]. Enegies must be in Hartree, with the
# maximum precision provided by TPP.
#
# Physical constants and conversion factors are available in
# lib/colib_perl.pm. For instance, to convert from eV to au,
# you may define:
# my $au2ev = units("au2ev");
# where $au2ev = 27.21138386.
#
#------------------------------------------------------------------------------------
#
if ($typeofinput eq "jad"){
  if ($lvprt>=3) {print_STDOUT("$mdle Reading Epot.\n",$istep,$kt);}
	
  open(OUT,"$out") or die "Cannot find TPP output";
  $ns = 0;
  #number of lines where energies are stored
  my $nlines = $nstat + 2 ;
  # Extract the portion containing the data
  my $data_section = `grep -A $nlines '  Root   Mult.   Total Energy (a.u.)' $out`;
  #Saving everything as a list or @array
  
  if ($data_section) {
	  # Read line by line in the above data and capture the total energies
	my @lines = split("\n",$data_section);
	my @total_energies;
	foreach my $line (@lines) {
		if ($line =~ /^\s+\d+\s+\S+\s+([\d.-]+)/) {
			push @total_energies, $1;
		}
	}

	# Save extracted Total Energy values
	for (my $i = 0; $i < scalar(@total_energies); $i++) {
		my $root = $i + 1;
		$epot[$i] = $total_energies[$i];
	}
	} else {
		print "Energies section not found in the file.\n";
	}
}
elsif ($typeofinput eq "jnd"){
  if ($lvprt>=3) {print_STDOUT("$mdle Reading Epot.\n",$istep,$kt);}
  chdir('GRAD') or die "Could not change directory GRAD: $!";
  open(OUT,"$out") or die "Cannot find TPP output";
  $ns = 0;
  #number of lines where energies are stored
  my $nlines = $nstat + 2 ;

  # Extract the portion containing the data

  my $data_section = `grep -A $nlines '  Root   Mult.   Total Energy (a.u.)' $out`;

  #Saving everything as a list or @array

  if ($data_section) {
          # Read line by line in the above data and capture the total energies
        my @lines = split("\n",$data_section);
        my @total_energies;
        foreach my $line (@lines) {
                if ($line =~ /^\s+\d+\s+\S+\s+([\d.-]+)/) {
                        push @total_energies, $1;
                }
        }

        # Print extracted Total Energy values
        for (my $i = 0; $i < scalar(@total_energies); $i++) {
                my $root = $i + 1;
                $epot[$i] = $total_energies[$i];
        }
        } else {
                print "Data section not found in the file.\n";
        }
    chdir('../');
  }
  write_energy();

}

sub write_energy{
#
#====================================================================================
#
# This subroutine writes the potential energies to epot file.
# It supposes that @epot vector contains these energies.
# See subroutine read_energies();
# oldepot and new epot files, used for interpolation are also
# written.
#
#------------------------------------------------------------------------------------
#
  
  if ($lvprt>=3) {print_STDOUT("$mdle Writing Epot.\n",$istep,$kt);}

  $file = "epot";
  if (-s $file){
    copy($file,"oldepot") or die "Copy failed: $!";
  } 
  open(OUT,">$file") or die "Cannot write to $file";
  foreach(@epot){
    print OUT "$_\n"; 
  }
  close(OUT);
  copy($file,"newepot") or die "Copy failed: $!";
}

sub read_gradients{
# 
#====================================================================================
#
# This subroutine reads the energy gradients from TPP outputs.
# The gradient components are stored in the @gx, @gy, and @gz arrays.
# $gx[$ns][$ia] with indexes $ns = 0 to $nstat-1 and $ia = 0 to $nat-1.
# Thus, gx[0][1] is the component x of the gradient of the ground state
# for atom 2.
#
# Gradients must be converted to au (Hartree/Bohr) and stored
# with the highest precision provided by the TPP.
#
# For many programs, only the gradient of the current state ($nstatdyn)
# is computed. In this cased, @gx, @gy, and @gz of the other states must
# be given as zeros.
#
#------------------------------------------------------------------------------------
#
if ($typeofinput eq "jad"){
  if ($lvprt>=3) {print_STDOUT("$mdle Reading Gradients.\n",$istep,$kt);}
  
  # Initialize array with zeros
  for ($ns = 0; $ns <= $nstat-1; $ns++){
    for ($ia = 0; $ia <= $nat-1; $ia++){
      $gx[$ns][$ia] = 0.0;
      $gy[$ns][$ia] = 0.0;
      $gz[$ns][$ia] = 0.0;
    }
  }
 
  # Copying grad file from scr.geom directory
  `cp $scr/$grad .`; 
  # Opening Gradient file 
  open my $fhh, '<', $grad or die "Couldn't open file to read gradients" ;

  #Skipping the first 2 lines
  <$fhh>;
  <$fhh>;
  #Storing gradients in arrays 
  my (@gx_new, @gy_new, @gz_new);
  while (my $line = <$fhh>) {
	  chomp $line;
	  my ($crp, $symbol, $gxx, $gyy, $gzz) = split /\s+/, $line;
	  push @gx_new, $gxx;
	  push @gy_new, $gyy;
	  push @gz_new, $gzz;
  }

  close $fhh;

  for ($ia = 0; $ia <= $nat-1; $ia++) {
	  $gx[$nstatdyn-1][$ia] = $gx_new[$ia];
	  $gy[$nstatdyn-1][$ia] = $gy_new[$ia];
	  $gz[$nstatdyn-1][$ia] = $gz_new[$ia];
  }
}

elsif ($typeofinput eq "jnd"){
  if ($lvprt>=3) {print_STDOUT("$mdle Reading Gradients.\n",$istep,$kt);}

  # Initialize array with zeros
  for ($ns = 0; $ns <= $nstat-1; $ns++){
    for ($ia = 0; $ia <= $nat-1; $ia++){
      $gx[$ns][$ia] = 0.0;
      $gy[$ns][$ia] = 0.0;
      $gz[$ns][$ia] = 0.0;
    }
  }

  # Copying grad file from scr.geom directory
  `cp $scr/$grad .`;
  # Opening Gradient file
  open my $fhh, '<', $grad or die "Couldn't open file to read gradients" ;

  #Skipping the first 2 lines
  <$fhh>;
  <$fhh>;
  #Storing gradients in arrays
  my (@gx_new, @gy_new, @gz_new);
  while (my $line = <$fhh>) {
          chomp $line;
          my ($crp, $symbol, $gxx, $gyy, $gzz) = split /\s+/, $line;
          push @gx_new, $gxx;
          push @gy_new, $gyy;
          push @gz_new, $gzz;
  }

  close $fhh;

  for ($ia = 0; $ia <= $nat-1; $ia++) {
          $gx[$nstatdyn-1][$ia] = $gx_new[$ia];
          $gy[$nstatdyn-1][$ia] = $gy_new[$ia];
          $gz[$nstatdyn-1][$ia] = $gz_new[$ia];
  }
}

  # This is just a pseudocode suggestion. It must be adapted to the actual case.
  #open(IN,"my_TPP.output") or die "Cannot find TPP output";
  #$ns = 0;
  #$ia = 0;
  #while(<IN>){
    # Search for Gradient of state $ns+1.
    # If it is found:
    # $gx[$ns][$ia] = $value;
    # $gy[$ns][$ia] = $value;
    # $gz[$ns][$ia] = $value;
    # $ns++;
    # $ia++;
    #
    #}
    #close(IN);
  
  write_grad();

  }

sub write_grad{
#
#====================================================================================
#
# This subroutine writes the energy gradients to grad and grad.all files.
#
#------------------------------------------------------------------------------------
#

  if ($lvprt>=3) {print_STDOUT("$mdle Writing Gradients.\n",$istep,$kt);}
 
  open(OUT1,">grad.all") or die "Cannot write to grad.all";
  open(OUT2,">grad")     or die "Cannot write to grad";

  for ($ns = 0; $ns <= $nstat-1; $ns++){
    for ($ia = 0; $ia <= $nat-1; $ia++){
      print OUT1 "$gx[$ns][$ia]  $gy[$ns][$ia]  $gz[$ns][$ia]\n";
      if ($ns == $nstatdyn-1){
        print OUT2 "$gx[$ns][$ia]  $gy[$ns][$ia]  $gz[$ns][$ia]\n";
      }
    }
  }

  close(OUT1);
  close(OUT2);
  }

sub treat_nacme{
#
#====================================================================================
#
# This subroutine processes the NACME. 
# If it was computed by TPP, it reads and writes them.
# If not, it may call cioverlap to get the couplings.  
#
#------------------------------------------------------------------------------------
#
 
  if ($lvprt>=3) {print_STDOUT("$mdle Treating couplings.\n",$istep,$kt);}

  $run_flag=read_single_value("which_run_is_that","");

  # $run_flag = "second run" is a new program call after hopping.
 
  if ($run_flag ne "second run"){

    # ncoup is half of the number of coupling vectors
    $ncoup = $nstat*($nstat-1)/2;

    # Intialize arrays
    for ($nc = 0; $nc <= $ncoup-1; $nc++){
      for ($ia = 0; $ia <= $nat-1; $ia++){
         $hx[$nc][$ia] = 0.0;
         $hy[$nc][$ia] = 0.0;
         $hz[$nc][$ia] = 0.0;
      }
    }

    # Create or update oldh (coupling in the previous step) for interpolation.
    if (-s "nad_vectors"){
      copy("nad_vectors","oldh") or die "Copy of nad_vectors to oldh failed: $!";
    }elsif(!-s "nad_vectors"){
      write_nacme("oldh");
    }

    $vdoth = getkeyword("sh.inp","vdoth",0);
    if ($vdoth == 0){
       read_nacme();
    }elsif(($vdoth == 1) or ($vdoth < 0)){
       run_cioverlap();
    }  # Even if you don't define any coupling here, you can still use TD-BA

    # Write nonadiabatic couplings
    write_nacme("nad_vectors");

    # Correct nonadiabatic coupling phases
    adjust_phase();
 
    # Create or update newh (coupling in the current step) for interpolation.
    copy("nad_vectors","newh") or die "Copy failed: $!";

  }
 }

sub read_nacme{
#
#====================================================================================
#
# This subroutine reads nonadiabatic coupling vectors from TPP output.
# The couplings components are stored in the @hx, @hy, and @hz arrays
# written to the nad_vectors file.
#
# $hx[$nc][$ia] has indexes $nc = 0 to $ncoup-1 and $ia = 0 to $nat-1.
# ncoup = nstat*(nstat-1)/2 is half of the number of coupling vectors.
# For instance, for nstat = 4, we have ncoup = 6 coupling vectors.
#
# The order in which couplings are stored is that of the elements of a
# lower triangular matrix. For example:
#     1    2    3    4
# 1
# 2  h21
# 3  h31  h32
# 4  h41  h42  h43
# leads to the order:
# nc       =  0    1    2    3    4    5
# coupling = h21  h31  h32  h41  h42  h43
#
# Thus, hx[0][1] is the component x of the h21 coupling for atom 2.
#
# Nonadiabatic couplings must be converted to au (1/Bohr) and stored
# with the highest precision provided by the TPP.
#
# If one or more couplings are not computed, they must be given as zeros.
#
# Note the writing convention has (higher-state,lower-state)-index order.
# If TPP computes the opposite, you should invert the indexes (remember
# that hkj = -hjk). For istance, if TPP computes h12 instead of h21, just
# write h21 = -h12.
#
#------------------------------------------------------------------------------------
  my (@nac01x, @nac01y, @nac01z);
  my (@nac12x, @nac12y, @nac12z);
  my $grad01 = ("grad01.xyz");
  my $grad12 = ("grad12.xyz");
  my ($ncc1, $ncc2);

  if ($lvprt>=3) {print_STDOUT("$mdle Reading coupling vectors.\n",$istep,$kt);}
  
  if ($nstatdyn > 1) {
    if ($nstatdyn != $nstat){	  
  	`cp $scr1/$grad $grad01`;
  	`cp $scr2/$grad $grad12`; #grad01 and grad12 are NACS for the current state (1) and the state below (0) and the state above (2)
	open my $nac10, '<', $grad01 or die "Couldn't open file to read NACs" ;
	<$nac10>;
  	<$nac10>;
	#Storing NAC_state1_state0 in arrays
	while (my $line = <$nac10>) {
		chomp $line;
		my ($crp, $symbol, $nadx, $nady, $nadz) = split /\s+/, $line;
		push @nac01x, $nadx;
		push @nac01y, $nady;
		push @nac01z, $nadz;
	}
	#storing NAC_state1_state2 in arrays
  	open my $nac21, '<', $grad12 or die "Couldn't open file to read NACs" ;  
	<$nac21>;
  	<$nac21>;

  	while (my $line = <$nac21>) {
        	chomp $line;
        	my ($crp1, $symbol1, $nadx1, $nady1, $nadz1) = split /\s+/, $line;
        	push @nac12x, $nadx1;
        	push @nac12y, $nady1;
        	push @nac12z, $nadz1;
        }
     }
     elsif ($nstatdyn == $nstat){
     	`cp $scr1/$grad $grad01`; #only state below will be counted
	open my $nac10, '<', $grad01 or die "Couldn't open file to read NACs";
        <$nac10>;
        <$nac10>;
        #Storing NAC_state1_state0 in arrays
        while (my $line = <$nac10>) {
                chomp $line;
                my ($crp, $symbol, $nadx, $nady, $nadz) = split /\s+/, $line;
                push @nac01x, $nadx;
                push @nac01y, $nady;
                push @nac01z, $nadz;
        }
     }
  }
	
  elsif ($nstatdyn == 1) {
        `cp $scr1/$grad $grad12`;
	 #storing NAC_state1_state2 in arrays
  	open my $nac21, '<', $grad12 or die "Couldn't open file to read NACs" ;
  	<$nac21>;
  	<$nac21>;

  	while (my $line = <$nac21>) {
        	chomp $line;
        	my ($crp, $symbol, $nadx, $nady, $nadz) = split /\s+/, $line;
        	push @nac12x, $nadx;
        	push @nac12y, $nady;
        	push @nac12z, $nadz;
        }
    }

  $ncc1 = ($nstatdyn * ($nstatdyn - 1)/2)-1;
  $ncc2 = (($nstatdyn * ($nstatdyn - 1)/2)-1)+$nstatdyn;
    for ($ia = 0; $ia <= $nat-1; $ia++){
        if ($nstatdyn == 1) {
                $hx[$ncc2][$ia] = -$nac12x[$ia];
                $hy[$ncc2][$ia] = -$nac12y[$ia];
                $hz[$ncc2][$ia] = -$nac12z[$ia];
        }
        elsif ($nstatdyn > 1){
                if ($nstatdyn != $nstat){
                        $hx[$ncc1][$ia] = -$nac01x[$ia];
                        $hy[$ncc1][$ia] = -$nac01y[$ia];
                        $hz[$ncc1][$ia] = -$nac01z[$ia];
                        $hx[$ncc2][$ia] = -$nac12x[$ia];
                        $hy[$ncc2][$ia] = -$nac12y[$ia];
                        $hz[$ncc2][$ia] = -$nac12z[$ia];
                }
                elsif ($nstatdyn == $nstat){
                        $hx[$ncc1][$ia] = -$nac01x[$ia];
                        $hy[$ncc1][$ia] = -$nac01y[$ia];
                        $hz[$ncc1][$ia] = -$nac01z[$ia];
                }
        }
  }
		        

}


sub write_nacme{
#
#====================================================================================
#
# This subroutine writes @h to a file. The definition and order  of @h are explained
# in subroutine @read_nacme. The output file is passed as an argument of the
# subroutine.
#
#------------------------------------------------------------------------------------
#
 
  if ($lvprt>=3) {print_STDOUT("$mdle Writing coupling vectors.\n",$istep,$kt);}

  ($file) = @_;  
  open(OUT,">$file") or die "Cannot write to $file";

  for ($nc = 0; $nc <= $ncoup-1; $nc++){
    for ($ia = 0; $ia <= $nat-1; $ia++){
      print OUT "$hx[$nc][$ia]  $hy[$nc][$ia]  $hz[$nc][$ia]\n";
    }
  } 

  close(OUT);  

}

sub adjust_phase{
#
#====================================================================================
#
# This subroutine computes cos(q) where q is the angle between couplings at t and t-Dt.
# If cos(q) >= 0, phase =1
#           < 0, phase = -1
# The coupling h(t) is rewritten as h(t)*phase.
#
#------------------------------------------------------------------------------------
#

  if ($lvprt>=3) {print_STDOUT("$mdle Fixing the coupling phase.\n",$istep,$kt);}

  write_warning1();

  # Create or update newh (coupling in the current step) for interpolation.
  copy("nad_vectors","newh") or die "Copy failed: $!";
  
  $retval = callprogsp( $mld, "escalar", $mdle );
  if ($retval != 0){
      die "$mdle is dying now\n";
  }
  copy("nadv","nad_vectors") or die "Copy failed: $!";

  write_debug2();

}

sub write_warning1{
#
#====================================================================================
#
# Wrting warning messages concerning phases.
#
#------------------------------------------------------------------------------------
#

  my  $getphase = getkeyword("sh.inp","getphase","1");
  $vdoth    = getkeyword("sh.inp","vdoth","0");

  if ($lvprt >= 3){
    if ($vdoth != 0){
      print_STDOUT("\n");
      print_STDOUT("    *********************************************\n");
      print_STDOUT("    WARNING: phase adjustment is not implemented \n");
      print_STDOUT("    for vdoth=$vdoth. However, cioverlap takes   \n");
      print_STDOUT("    care of the phase.                           \n");
      print_STDOUT("    *********************************************\n");
      print_STDOUT("\n");
    }
  }

  if ($getphase != 1){
    print_STDOUT("\n");
    print_STDOUT("    *********************************************\n");
    print_STDOUT("    WARNING: phase adjustment is only implemented \n");
    print_STDOUT("    for getphase = 1 \n");
    print_STDOUT("    *********************************************\n");
    print_STDOUT("\n");
  }

}

sub write_debug2{
#
#====================================================================================
#
# This subroutine write debug messages concerning the phase.
#
#------------------------------------------------------------------------------------
#

  if ($lvprt >= 3){copy("escalar.log","../$DEBUG/.") or die "Copy failed: $!";}

  if ($lvprt >= 2){
    open( EP, "escph" ) or die "Cannot open escph to read.";
    print_STDOUT("The global phase(s) is(are):",$istep,$kt);
    while(<EP>){
      chomp;$_ =~ s/^\s*//;$_ =~ s/\s*$//;
      print_STDOUT(" $_  ",$istep,$kt);
    }
    print_STDOUT("\n\n",$istep,$kt);
    close(EP);
  }

}

sub run_cioverlap{
#
#====================================================================================
#
# This subroutine is going to control the wavefunction overlap calculation.
# This template is just a draft given the general outline of the
# execution.
#
# You can use the overlap programs integrated into Newton-X or using your own.
#
# cioverlap is a set of programs to build the wavefunction and compute their overlaps.
# See NX docs for info about these programs.
#  
# First step: prepare the wavefunction for each state.
# For linear-response, 
# i) cis_casida creates a CIS-like wavefunction from LR coeff.
# ii) cis_slatergen generates a CIS-like determinant basis.
# See NX docs for info about these programs.
#
# Second step: after having the wavefuntion, you need overlap integrals
# between atomic orbitals computed at the geometry of time t (current) and
# t-Dt (previous step). The eay we usually do that is creating a double
# molecule composed of the coordinates at t + coordinates at t+Dt.
# Then, we ask the TPP to compute the integrals for this double molecule.
# The geometry in t is writen in geom, that in t-Dt is writen in geom.old.
#  
# Third step: compute wavefunction overlaps. 
# Here you have three options (CPROG keyword):
#
# i) cioverlap DD ($NX/cioverlap-64/cioverlap program).
# It is general and accounts for determinants with multiple excitations.
#
# ii) cioverlap OD ($NX/cioverlap-64/cioverlap.od program).
# It's very fast, but it works only for CIS wavefunctions.
#
# iii) Alternatively, you may use any some external program to compute the
# overlap. In this case, instead of the previous code (if ($CPROG == ...)),
# you should call here this external program and write the output to a file
# named run_ciovelap.log, containing the overlap matrix in the format:
# Nstat Nstat
# C(1,1)     C(1,2)     ...  C(1,Nstat)
# ...
# C(Nstat,1) C(Nstat,2) ...  C(Nstat,Nstat)
# where Nstat is the number of states and 
# C(m,n) = < psi_m(t) | psi_n(t-dt) >
# For example, for three states, the overlap marix ay look like:
# 3 3
# 0.995441495621 -0.002585762807 0.000241713276
# 0.003428804668 0.981313114771 -0.157099932930
# 0.000553617088 0.157249101829 0.981336068872
#
#------------------------------------------------------------------------------------
#
 
  if ($lvprt>=3) {print_STDOUT("$mdle Computing state-overlap matrix.\n",$istep,$kt);}
 
  # The main variable determining the behavior of the cioverlap 
  # programs are in jiri.inp file.
  #
  if ($lvprt>=3) {print_STDOUT("$mdle Computing state-overlap matrix.\n",$istep,$kt);}

  ($kross_d,$cascade_d,$current_d,$never_state_d,$include_pair_d,$e_ci_d,$ci_cons_d,
   $cio_options_d,$cisc_options_d,$idalton_d,$cprog_d,$coptda_d,
   $ncore_d,$ndisc_d,$blasthread_d)=load_defaults("jiri.inp",$prog,$vdoth);

  $CPROG = getkeyword("jiri.inp","cprog",$cprog_d);

  if     ($CPROG == 1) {
    $retval = callprog($mldcio,"cioverlap $CIO_OPTIONS < cioverlap.input >cioverlap.out",$mdle,$BASEDIR);
    if ($retval != 0){die "$mdle is dying now (cioverlap)\n";}
  }elsif ($CPROG == 2) {
    $retval = callprog($mldcio,"cioverlap.od ",$mdle,$BASEDIR);
    if ($retval != 0){die "$mdle is dying now (cioverlap.od)\n";}
  }elsif ($CPROG == 3) {
    # Call your program to compute wavefunction overlaps.
  } 
 
  # Read cioverlap results and write them to nad_vectors or use them for local diabtization
  #
  if ($vdoth == 1){
    if ($lvprt>=3) {print_STDOUT("Overlap matrix will be used to compute time-derivative couplings.\n",$istep,$kt);}
    $retval = callprogsp($mld, "read_cioverlap.pl", $mdle);
    if ($retval != 0){die "$mdle is dying now (read_cioverlap.pl)\n";}
  }elsif($vdoth < 0){
    if ($lvprt>=3) {print_STDOUT("Overlap matrix will be used to do local diabatization.\n",$istep,$kt);}
    # In this case, nothing needs to be done. sh program will automatically read run_cioverlap.log.
  }

}
