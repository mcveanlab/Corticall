package DM;
use strict;

use version 0.77;
our $VERSION = qv('0.2.3');
use 5.006;
use warnings;
use File::Temp qw/tempdir tempfile/;
use File::Basename;
use Carp;

=head1 NAME

DM - Distributed Make: A perl module for running pipelines

=head1 VERSION

0.2.3

=head1 SYNOPSIS

use DM 0.002001;

# create a DM object

my $dm = DM->new(
    dryRun => 0,
    numJobs => 5
)

# add rule with target, prerequisite, and command to use to update the target

$dm->addRule( 'targetFile', 'prerequisiteFile', 'touch targetFile' );

# add more rules ...

# executed the pipeline 

$dm->execute();

=head1 DESCRIPTION

DM is a perl module for running pipelines.  DM is based on GNU make.  Currently, DM supports running on a single computer or an SGE managed cluster.

=head1 GOOD PRACTICE

=over

=item * 

Never make a directory a dependency. DM creates directories as it needs them.

=item *

Never create rules that delete files. Delete files by hand instead. Chances are, You will be sorry otherwise.

=item *

make runs in dryRun mode by default (this is for your own safety!).  Pass in 'dryRun => 0' to new() to run.

=back

=head1 OPTIONS

Any of the following options can be passed in to a call to new() in order to change the defaults on how make is run by DM. The default value is listed behind the option name.

=head2 GNU make specific options

=over

=item dryRun 1

show what will be run, but don't actually run anything. Corresponds to -n option in GNU make.

=item numJobs 1

maximum number of jobs to run, or "" for maximum concurrency permitted by dependencies. Applicable to queue and non-queue situations. Corresponds to -j option in GNU make.

=item keepGoing 0

if any job returns a non-zero exit status, the default behaviour is not to submit any further jobs and wait for the others to finish.  If this option is true, then any jobs that do not depend on the failed job(s) will still be submitted. Corresponds to -k option in GNU make.

=item alwaysMake 0

Unconditionally make all targets. Corresponds to -B option in GNU make.

=item touch 0

If true, touch all targets such that make will think that all files have been made successfully. This is only partially supported, as touch will not create any prerequisite directories. Corresponds to -t option in GNU make. 

=item ignoreErrors 0

Corresponds to -i option in GNU make.

=back

=head2 SGE specific options

These options are passed to qsub for submitting jobs on an SGE cluster

=over

=item cluster undef

Type of cluster (localhost, SGE, PBS, LSF).  Is detected automagically by DM.

=item queue undef

Corresponds to -q.

=item projectName undef

Corresponds to -P.

=item PE { name => undef, range => undef }

Anonymous hash reference. Corresponds to -pe option

=item name

Corresponds to -N option.

=back

=head2 other options

=over

=item globalTmpDir undef

Directory for storing temporary files that can be accessed by every node on the cluster (usually not /tmp)

=item tmpdir /tmp

Directory for storing temporary files.

=back

=head1 GENERAL FUNCTIONS

=head2 new()

Returns a DM object.  Options (see the Options section) can be passed to new() as key value pairs.

=over

=item Required Arguments

none

=item Returns

DM object

=back

=cut

sub new {
    my ( $class, %args ) = @_;

    my %self = (

        # Make options
        'dryRun'         => 1, # show what will be run, but don't actually run anything
        'numJobs'        => 1
        , # maximum number of jobs to run, or "" for maximum concurrency permitted by dependencies
          # Applicable to queue and non-queue situations
        'keepGoing'      => 0,
        'alwaysMake'     => 0,
        'debugging'      => 0,
        'ignoreErrors'   => 0,
        'printDirectory' => 0,
        'touch'          => 0,
        'unlink'         => 1,    # 0 = don't clean tmp file

        # Cluster engine options
        'queue'          => undef,
        'cluster'        => undef,
        'PE'             => { name => undef, range => undef },  # parallel environment
        'memRequest'     => 4,                                  # in gigabytes
        'rerunnable'     => 0,
        'name'           => undef,
        'projectName'    => undef,
        'outputFile'     => 'distributedmake.log',
        'extra'          => '',
        supportedEngines => { SGE => 1, localhost => 1, LSF => 0, PBS => 0 },

        # make options
        'tmpdir'         => '/tmp',
        'target'         => 'all',
        'targets'        => [],

        # job array related information
        # check for undef to determine if job array has been
        # started but not ended
        'currentJobArrayObject' => undef,
        'globalTmpDir'          => undef,    # necessary for running job arrays

        ## other attributes...
        %args,
    );

    $self{'makefile'} = new File::Temp(
        TEMPLATE => "$self{'tmpdir'}/DM_XXXXXX",
        SUFFIX   => ".makefile",
        UNLINK   => $self{'unlink'}
    );

    chomp( my $sge_qmaster = qx(which sge_qmaster 2>/dev/null) );
    chomp( my $pbsdsh      = qx(which pbsdsh 2>/dev/null) );
    chomp( my $bsub        = qx(which bsub 2>/dev/null) );

    unless ( defined $self{cluster} ) {

        if ( -e $sge_qmaster ) { $self{'cluster'} = 'SGE'; }

  #    elsif ( -e $pbsdsh )      { $self{'cluster'} = 'PBS'; } not supported yet
        elsif ( -e $bsub ) { $self{'cluster'} = 'LSF'; }
        else               { $self{'cluster'} = 'localhost'; }

    }

    bless \%self, $class;

    my $self = \%self;
    $self->_check_arg_consistency(%self);

    return \%self;
}

=head2 addRule()

This function creates a basic dependency between a prerequisite, target and command.  The prerequisite is a file that is required to exist in order to create the target file.  The command is used to create the target file from the prerequisite.

=over

=item Required Arguments

<string>  -  target file

<string or ref to array of strings>  -  prerequisite file(s)

<string>  -  command

=item Returns

none

=back

=cut

sub addRule {
    my ( $self, $targetsref, $dependenciesref, $cmdsref, %batchjoboverrides ) =
      @_;
    my @targets =
      ( ref($targetsref) eq 'ARRAY' ) ? @$targetsref : ($targetsref);
    my @dependencies =
      ( ref($dependenciesref) eq 'ARRAY' )
      ? @$dependenciesref
      : ($dependenciesref);
    my @cmds = ( ref($cmdsref) eq 'ARRAY' ) ? @$cmdsref : ($cmdsref);

    my %bja = (
        'cluster'     => $self->{'cluster'},
        'queue'       => $self->{'queue'},
        'PE'          => $self->{'PE'},
        'memRequest'  => $self->{'memRequest'},
        'rerunnable'  => $self->{'rerunnable'},
        'name'        => $self->{'name'},
        'projectName' => $self->{'projectName'},
        'outputFile'  => $self->{'outputFile'},
        'extra'       => $self->{'extra'},
        %batchjoboverrides,
    );

    # we should really be checking arguments in one central place,
    # instead of ad-hoc throughout the script - winni
    $self->_check_arg_consistency(%bja);

# Setup the pre-commands (things like pre-making directories that will hold log files and output files)
    my @precmds;
    my $logdir = dirname( $bja{'outputFile'} );
    if ( !-e $logdir ) {
        my $mklogdircmd = "\@test \"! -d $logdir\" && mkdir -p $logdir";
        push( @precmds, $mklogdircmd );
    }

    foreach my $target (@targets) {
        my $rootdir = dirname($target);

        my $mkdircmd = "\@test \"! -d $rootdir\" && mkdir -p $rootdir";
        push( @precmds, $mkdircmd );
    }

# Setup the user's commands, taking care of imposing memory limits and adding in cluster prefix commands
    for ( my $i = 0 ; $i <= $#cmds ; $i++ ) {
        if (   $cmds[$i] =~ /^java /
            && $cmds[$i] =~ / -jar /
            && $cmds[$i] !~ / -Xmx/ )
        {
            $cmds[$i] =~ s/^java /java -Xmx$bja{'memLimit'}g /;
        }
    }

    if ( !defined( $bja{'name'} ) ) {
        my $firstcmd = $cmds[0];
        my $name     = "unknown";
        if ( $firstcmd =~ /java/ && $firstcmd =~ /-jar/ ) {
            ($name) = $firstcmd =~ /-jar\s+(.+?)\s+/;
        }
        else {
            $firstcmd =~ m/([A-Za-z0-9\_\/]+)/;
            $name = $1;
        }

        $bja{'name'} = &basename($name);
    }

    # -l is not supported on all SGE systems
    my $memRequest = $bja{'memRequest'};

    my $cmdprefix  = "";
    my $cmdpostfix = "";

    if ( $bja{'cluster'} eq 'SGE' ) {
        $cmdprefix = "qsub -sync y -cwd -V -b yes -j y"
          . (
            defined $memRequest
            ? qq/ -l mem_free=${memRequest}G,h_vmem=${memRequest}G/
            : q//
          ) . " -o $bja{'outputFile'} -N $bja{'name'}";
        $cmdprefix .=
          ( defined( $bja{'projectName'} ) )
          ? " -P $bja{'projectName'}"
          : "";
        $cmdprefix .= ( $bja{'rerunnable'} == 1 ) ? " -r yes" : " -r no";
        $cmdprefix .=
          defined( $bja{'queue'} )
          ? " -q $bja{'queue'}"
          : "";
        $cmdprefix .=
          defined( $bja{PE}->{name} )
          ? " -pe " . $bja{PE}->{name} . q/ / . $bja{PE}->{range}
          : "";
        $cmdprefix .= $bja{'extra'};
    }
    elsif ( $bja{'cluster'} eq 'PBS' ) {

    }
    elsif ( $bja{'cluster'} eq 'LSF' ) {

#$cmdprefix = "bsub -q $bja{'queue'} -M $memCutoff -P $bja{'projectName'} -o $bja{'outputFile'} -u $bja{'mailTo'} -R \"rusage[mem=$integerMemRequest]\" $wait $rerunnable $migrationThreshold $bja{'extra'}";
    }
    elsif ( $bja{'cluster'} eq 'localhost' ) {
        #$cmdpostfix = "| tee -a $bja{'outputFile'}";
    }
    else {
        croak
"programming error. unknkown engine $bja{cluster} was not caught by argument checking";
    }

    my @modcmds;

    foreach my $cmd (@cmds) {
        my $modcmd = $cmd;

        # protect single quotes if running on SGE
        # perhaps this could be an issue with one-liners
        #using double quotes? -- winni
        if ( $bja{cluster} eq q/SGE/ ) {
            $modcmd =~ s/'/"'/g;
            $modcmd =~ s/'/'"/g;
            $modcmd =~ s/\$/\$\$/g;
        }

        # protect $ signs from make by turning them into $$
        if ( $bja{cluster} eq q/localhost/ ) {
            $modcmd =~ s/\$/\$\$/g;
        }

        push( @modcmds, "$cmdprefix   $modcmd   $cmdpostfix" );
    }

# Setup the post-commands (touching output files to make sure the timestamps don't get screwed up by clock skew between cluster nodes).
    my @postcmds;
    foreach my $target (@targets) {
        push( @postcmds, "\@touch -c $target" );
    }

    # Emit the makefile commands
    print { $self->{'makefile'} } "$targets[0]: "
      . join( " ",    @dependencies ) . "\n\t"
      . join( "\n\t", @precmds ) . "\n\t"
      . join( "\n\t", @modcmds ) . "\n\t"
      . join( "\n\t", @postcmds ) . "\n\n";

    push( @{ $self->{'targets'} }, $targets[0] );
}

=head2 execute()

This method is called after all rules have been defined in order to write the make file and execute it.  No mandatory options. Takes only overrides.

=over

=item Required Arguments

none

=item Returns

none

=back

=cut

sub execute {
    my ( $self, %overrides ) = @_;

    # checking to make sure all started job arrays were ended.
    die "Need to end all started job arrays with endJobArray()"
      if defined $self->{currentJobArrayObject};

    print { $self->{'makefile'} } "all: "
      . join( " ", @{ $self->{'targets'} } ) . "\n\n";
    print { $self->{'makefile'} } ".DELETE_ON_ERROR:\n";

    my %makeargs = (
        'dryRun'         => $self->{'dryRun'},
        'numJobs'        => $self->{'numJobs'},
        'keepGoing'      => $self->{'keepGoing'},
        'alwaysMake'     => $self->{'alwaysMake'},
        'debugging'      => $self->{'debugging'},
        'ignoreErrors'   => $self->{'ignoreErrors'},
        'printDirectory' => $self->{'printDirectory'},
        'touch'          => $self->{'touch'},
        'target'         => $self->{'target'},
        'touchFiles'     => $self->{'touchFiles'},
        %overrides,
    );

    my $numjobs = $makeargs{'numJobs'};

    my $makecmd = "make"
      . ( $makeargs{'dryRun'}         ? " -n" : "" )
      . ( $makeargs{'keepGoing'}      ? " -k" : "" )
      . ( $makeargs{'alwaysMake'}     ? " -B" : "" )
      . ( $makeargs{'ignoreErrors'}   ? " -i" : "" )
      . ( $makeargs{'printDirectory'} ? " -w" : "" )
      . ( $makeargs{'touch'}          ? " -t" : "" )
      . (
        $makeargs{'debugging'} =~ /[abvijm]+/
        ? " --debug=$makeargs{'debugging'}"
        : ""
      )
      . (    $makeargs{'debugging'} =~ /\d+/
          && $makeargs{'debugging'} == 1 ? " -d" : "" )
      . " -j $numjobs" . " -f "
      . $self->{'makefile'}->filename
      . " $makeargs{'target'}";

    print "$makecmd\n";
    system($makecmd);
    print "$makecmd\n";
}

=head1 JOB ARRAY FUNCTIONS

=head2 Workflow

First, initialize a job array with startJobArray().  Add rules to the job array with addJobArrayRule().  Last, call endJobArray() to signal that no more rules will be added to this particular job array. Multiple job arrays can be defined after each other in this manner. execute() can only be called if the most recently started job array has been completed with endJobArray.

On SGE, the job array will only be started once the prerequisites of all job array rules have been updated.  On other platforms, each job will start once its prerequisite has been updated.  However, on all platforms, the job array target will only be updated once all rules of that job array have completed successfully.  

Only the target specified in startJobArray() should be used as a prerequisite for other rules.  The targets specified through addJobArrayRule() should never be used as prerequisites for other rules.

=head2 startJobArray()

daes nothing unless 'cluster' eq 'SGE'.
Requires 'target' and 'globalTmpDir' to be specified as key value pairs:
    startJobArray(target=>$mytarget, globalTmpDir=>$mytmpdir)

=cut

sub startJobArray {
    my ( $self, %overrides ) = @_;

    die "startJobArray was called before endJobArray"
      if defined $self->{currentJobArrayObject};

    my %args = (
        commandsFile => undef,
        targetsFile  => undef,
        prereqsFile  => undef,
        target       => undef,
        globalTmpDir => undef,
        name         => undef,
        %overrides,
    );

    die "startJobArray needs a target to be specified"
      unless defined $args{target};

    # pull object tmp dir if none was passed in
    $args{globalTmpDir} =
      defined $args{globalTmpDir} ? $args{globalTmpDir} : $self->{globalTmpDir};

    # make sure globalTmpDir is defined and exists
    die
"startJobArray needs a global temporary directory to be specified with globalTmpDir and for that direcory to exist"
      unless defined $args{globalTmpDir} && -d $args{globalTmpDir};

    # definition of jobArrayObject
    my $jobArrayObject = {
        fileHandles => {},
        files       => {},

        # final target to touch when all job array rules completed successfully
        target => $args{target},

        # lists of all targets and prereqs of all rules added to job array
        arrayTargets => [],
        arrayPrereqs => [],

        # name of job array
        name => $args{name},
    };

    ## initialize files to hold targets, commands and prereqs for job array
    # open file handles
    for my $name (qw(commands targets prereqs)) {
        (
            $jobArrayObject->{fileHandles}->{$name},
            $jobArrayObject->{files}->{$name}
          )
          = tempfile(
            $name . '_XXXX',
            DIR    => $args{globalTmpDir},
            UNLINK => 1
          );
    }

    # save new object
    $self->{currentJobArrayObject} = $jobArrayObject;
    return $jobArrayObject;
}

=head2 addJobArrayRule()

This structure is designed to work with SGE's job array functionality.  Any rules added to a jobArray structure will be treated as simple add rules when running on localhost, LSF or PBS, but will be executed as a jobArray on SGE.

=head3 Required Arguments

takes three inputs: target, prereqs, command as key value pairs:

addJobArrayRule( 
    target  => $mytarget, 
    prereqs => \@myprereqs, 
    command => $mycommand 
);

prereqs may also be a scalar (string).  

The target is only for internal updating by the job array.  The target may not be used as a prerequisite for another rule.  Use the job array target instead.

=head3 Returns

none

=cut

sub addJobArrayRule {
    my $self = shift;

    # get input
    my %args = @_;

    # check to make sure startJobArray() has been run
    die "need to run startJobArray() first"
      unless defined $self->{currentJobArrayObject};

    # check required args.
    foreach my $arg (qw/target prereqs command/) {
        die "need to define $arg" unless defined $args{$arg};
    }

    # keep track of all rule targets
    my @targets =
      ref( $args{target} ) eq 'ARRAY' ? @{ $args{target} } : ( $args{target} );
    push @{ $self->{currentJobArrayObject}->{arrayTargets} }, $targets[0];

    # keep track of all rule prereqs
    my @prereqs = (
        ref( $args{prereqs} ) eq 'ARRAY'
        ? @{ $args{prereqs} }
        : $args{prereqs}
    );
    push @{ $self->{currentJobArrayObject}->{arrayPrereqs} }, @prereqs;

    # just use addRule unless we are in an SGE cluster
    unless ( $self->{cluster} eq 'SGE'
        || ( defined $args{cluster} && $args{cluster} eq 'SGE' ) )
    {
        $self->addRule( \@targets, $args{prereqs}, $args{command}, %args );
    }
    else {

        ### Add target, prereqs and command to respective files
        # TARGET
        print { $self->{currentJobArrayObject}->{fileHandles}->{targets} }
          $targets[0] . "\n";

        # COMMAND
        # need to make sure target directory exists
        my @precmds;
        foreach my $target (@targets) {
            my $rootdir  = dirname($target);
            my $mkdircmd = "test \"! -d $rootdir\" && mkdir -p $rootdir";
            push( @precmds, $mkdircmd );
        }

        print { $self->{currentJobArrayObject}->{fileHandles}->{commands} }
          join( q/ && /, ( @precmds, $args{command} ) ) . "\n";

        # PREREQS - also add prereqs to job array prereqs file
        print { $self->{currentJobArrayObject}->{fileHandles}->{prereqs} }
          join( q/ /, @prereqs ) . "\n";
    }
}

=head2 endJobArray()

Adds the rule that kicks off the job array. See Workflow for further description.

=head3 Requried Arguments

none

=head3 Returns

The target of the job array

=cut

sub endJobArray {

    my $self = shift;

    # close all file handles
    map { close( $self->{currentJobArrayObject}->{fileHandles}->{$_} ) }
      keys %{ $self->{currentJobArrayObject}->{fileHandles} };

    # determine how many tasks to kick off in job array
    my $arrayTasks = @{ $self->{currentJobArrayObject}->{arrayTargets} };

    # add job array rule
    #  makes sure target is touched when everything ran through successfully
    my $target = $self->{currentJobArrayObject}->{target};
    if ( $self->{cluster} eq 'SGE' ) {
        $self->addRule(
            $self->{currentJobArrayObject}->{target},
            $self->{currentJobArrayObject}->{arrayPrereqs},
            " -t 1-$arrayTasks:1  sge_job_array.pl  -t "
              . $self->{currentJobArrayObject}->{files}->{targets} . ' -p '
              . $self->{currentJobArrayObject}->{files}->{prereqs} . ' -c '
              . $self->{currentJobArrayObject}->{files}->{commands}
              . " && touch $target",
            name => $self->{currentJobArrayObject}->{name}
        );
    }
    else {
        $self->addRule(
            $self->{currentJobArrayObject}->{target},
            $self->{currentJobArrayObject}->{arrayTargets},
            "touch $target",
            name => $self->{currentJobArrayObject}->{name}
        );

    }

    $self->{currentJobArrayObject} = undef;
}

=head1 INTERNAL FUNCTIONS

=head2 _check_arg_consistency()

This function is used for checking the consistency of arguments passed to a DM object.

=cut

sub _check_arg_consistency {

    my ( $self, %overrides ) = @_;
    my %bja = ( %{$self}, %overrides );

    if ( defined $bja{PE}->{name} xor defined $bja{PE}->{name} ) {
        croak
          "both 'name' and 'range' need to specified when using the PE option";
    }

    # testing supported clusters
    my $isSupported = 0;
    foreach my $engine ( sort keys %{ $bja{supportedEngines} } ) {
        $isSupported = 1
          if ( $bja{cluster} eq $engine && $bja{supportedEngines}->{$engine} );
    }
    croak "$bja{cluster} is not a supported engine" unless $isSupported;

    # cluster related tests
    my $error = q//;
    if ( $bja{cluster} ne 'localhost' ) {
        $error .= "cluster is not 'localhost'\n\t";

        unless ( defined $bja{PE}->{name} || defined $bja{queue} ) {
            carp $error. "either 'queue' or 'PE' or both need to be defined";
        }
    }
}

=head1 AUTHORS

Kiran V Garimella <kiran@well.ox.ac.uk> and Warren W. Kretzschmar <warren.kretzschmar@well.ox.ac.uk>

=head1 BUGS

Please report any bugs or feature requests to C<bug-dm at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=DM>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc DM


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=DM>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/DM>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/DM>

=item * Search CPAN

L<http://search.cpan.org/dist/DM/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2012 Kiran V Garimella.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1;    # End of DM

