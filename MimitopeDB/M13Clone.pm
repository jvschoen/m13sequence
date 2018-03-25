package M13Clone;

use strict;
use warnings;

# Must determine a standard naming convention where I can pull
# the clone name from the sequence file name that our sequencing
# service provides

#long_name is the name of the file

sub new{
	my $class = shift;
	my $long_name = shift;
	my $self = {};
	bless $self, $class;
	
	$self->_initialize($long_name);
	
	return $self;
	
}

sub _initialize{
	my $self = shift;
	my $long_name = shift;
	
	$self->{id} = undef;
	$self->{type} = undef;
	$self->{name} = undef;
	$self->{short_name} = undef;
	$self->{insert} = undef;
	$self->{aa_seq} = undef;
	$self->{dna_obj} = undef;
	
	$self->{seq_files} = undef;
	$self->{date} = undef;
	
	
	$self->_setName($long_name);
	$self->_setDate();
	$self->_setType();
	
	return $self;
	
}


# Work on setting a type either plasmid or phage
sub _setType{
	my $self = shift;
	my $long_name = $self->{name};
	if($long_name =~ /.*p.*_M13/){
		$self->{type} = "plasmid";
	}elsif($long_name =~ /.*M13.*M13/){
		$self->{type} = "M13_phage";
	}
	
}

# Uses Regex to parse name from sequence file name
# long_name is the name of the file. The name is captured
# which is everything starting with lowercase "p" (denoting plasmid)
# before "_M13R"
# Macrogen Example: 171002-06_A01_pVRC01_A4_ec1_092817_M13R.ab1
# GeneWiz Example: pVRC01_A4_ec5_M13R-M13R_C08.ab1
sub _setName{
	my $self = shift;
	my $long_name = shift;
	
	$self->{name} = $long_name;
	
	#capture if plasmid
	$long_name =~ /.*(p.*)_M13/;
	my $short_name = $1;
 	#capture if M13 clone
 	if (!defined($short_name)){
 		$long_name =~ /.*(M13_.*)_96/; #primer name is 96GIII_M13Rev_Primer
 		$short_name = $1;
 	}
	$self->{short_name} = $short_name;
	
}

sub _setDate{
	my $self = shift;
	my $short_name = $self->{short_name};
	
	$short_name =~ /.*(\d{6})$/;
	my $date = $1;
	
	$self->{date} = $date;
}

sub translate{
	my $self = shift;
	my $dna_obj = shift;
	
	my $aa_seq;
	my $insert;
	my $aa_obj;
	
	if($aa_obj = _findFrame($dna_obj)){
		#$aa_seq = findFrame($node)->seq;
		$aa_seq = $aa_obj->primary_seq->seq;
		$insert = _getInsert(_findFrame($dna_obj));
		
	}else{ 
		$aa_seq = "null";
		$insert = "unable to determine"
		
		};			
	$self->{dna_obj} = $dna_obj;
	$self->{aa_seq} = $aa_seq;
	$self->{insert} = $insert;
	
}

# Finds the frame which has FYSH or ETVE found in the sequence
# If foward frames fail, a reverse complement is made and checked
sub _findFrame{
	my $dna_obj = shift;
	my $sequence_name = $dna_obj->primary_id;
	# print "Sequence: $sequence_name\n";
	for (my $frame = 0; $frame < 3; $frame++){
		# print "-Trying forward reading frame $frame-\n";
		my $aa_obj = $dna_obj->translate(-frame => $frame);
		my $sequence = $aa_obj->seq;
		
		# print "\n" . $sequence . "\n";
		
		if ($sequence =~ m/FYSH|ETVE/){
			#print "\nSequence Found: $sequence\n\n";
			return $aa_obj;
			
		}else{
			# print "-Trying reverse reading frame $frame-\n";
			my $dna_obj_rev = _reverseSeq($dna_obj);
			my $aa_obj_rev = $dna_obj_rev->translate(-frame => $frame);
			$sequence = $aa_obj_rev->seq;
			if ($sequence =~ m/FYSH|ETVE/){
				#print "Sequence Found: $sequence\n\n";
				return $aa_obj_rev;
			}
		}
		
	}
	print "Searched Sequence not found in any frame for:\n".
			"$sequence_name\n\n";
	return 0;
	
}

# Reverse Complement when using reverse primer
sub _reverseSeq{
	my $dna_obj = shift;
	#print "Before Reverse: " . $dna_obj->seq . "\n";
	
	my $revcomp = $dna_obj->revcom;
	my $revseq = $revcomp->seq;
	$dna_obj->seq($revseq);
	
	#print "After Reverse: " . $dna_obj->seq . "\n";
	
	return $dna_obj;
}

# Uses regex to extract the insert sequence from the
# Full AA sequence. 
# Requires that sequencing quality is good enough to
# have complete flanking sequences FYSHS-12-GGG or
# FYSHA-c7c-GGG
sub _getInsert{
	my $aa_obj = shift;
	my $aa_seq = $aa_obj->seq;
	my $sample_name = $aa_obj->primary_id;
	
	my $insert;
	
	if ($aa_seq =~ /FYSHS(.{12})GGG/){
		$insert = $1;
	}elsif ($aa_seq =~ /FYSHSA(.{9})GGG/){
		$insert = $1;
	}else { 
		$insert = "NoMatch";
	};
	
	return $insert;
}

1;
