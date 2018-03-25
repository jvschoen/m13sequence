#!/usr/bin/perl

# This program takes in file with 1+ sequences and attempts to
# extract the insert sequence that has been engineered between
# known flanking regions of M13 phage genome on pIII protein

# want to put this information into an SQLlite database where
# the sequence field is UNIQUE since we only want one clone per
# insert.
# The database should house the clone name, Insert seq, DNAseq,
# AAseq, sequence file name


# Design a matrix that has clone names based on Macrogen plate
# Readout. Must input another metadata file to match.

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Data::Dumper; ##must be sudo to install Simple Package
use DBI;
use XML::Simple;
use FindBin;
use M13Clone;
use Table;

use Cwd;
use File::Copy;
use File::Basename;

my $cwd = getcwd;
my $format = "fasta";

# Get files an DNA sequences
my $dir_path = $cwd . "/sequence_files";
my @seq_files = getFiles($dir_path);
my $dnas = getDNAs($dir_path, \@seq_files);


#make our hash table holding information to store in a database
my $clones = makeCloneHash($dnas);

# connect to database
my $db_file = "mimitopes_log.db";
my $db_source = "dbi:SQLite:dbname=$db_file";
my $dbh = DBI->connect($db_source);

# Add tables to database
my @table_names = ("unique_clones", "all_clones");
my @field_list = ("id", "type", "name", "short_name", "insert_aa", "aa_seq", "dna_seq","date");
my $tables = makeTables(\@table_names, \@field_list);
makeSQLtables($dbh, $tables);
insertData($dbh, $tables, $clones);

outputSQLtables($dbh, $tables, $cwd);

moveFiles($dir_path,\@seq_files); # Move processed files to another directory



# Outputs the tables to TSV files with name being the Table Name
sub outputSQLtables{
	my $dbh = shift;
	my $tables = shift;
	my $cwd = shift;
	my $out_dir = $cwd . "/output_files/";
	
	foreach my $table(@$tables){
		my $header = $table->getHeader();
		my $fields = $table->{fields};
		my $fields_holder = join(",", @$fields);
		
		my $name = $table->{name};
		
		my $sql = "SELECT rowid, * FROM $name";
		my $sth = $dbh->prepare($sql);
		$sth->execute();
			
		my $out_filename = $out_dir . "output_" . $name .".tsv";
		open my $out_file, '>', $out_filename;
		
		print $out_file "row_id\t" .$header . "\n";

		while(my @row = $sth->fetchrow_array){
			my $line = join("\t", @row);
			print $out_file $line . "\n";
		}
		close $out_file;
			
	}
}

sub insertData{
	my $dbh = shift;
	my $tables = shift;
	my $clones_hash = shift;
	
	
	foreach my $table(@$tables){
		my $table_name = $table->{name};
		my $field_list = $table->{fields};
		my $num_fields = scalar (@$field_list);
		
		my $fields_holder = join(",", @$field_list);
		 
		
		my $values_holder = ();
		for (my $i = 0; $i<$num_fields; $i++){
			push @$values_holder, "?";	
		}
		
		$values_holder = join(",", @$values_holder);
		
		my $clones = $clones_hash->{$table_name};
		
		foreach my $key(keys %$clones){
			my $clone = $clones->{$key};
			
			my $id 			= $clone->{id};
			#if (!defined($id)){$id = "null"}
			
			my $type 		= $clone->{type};
			#if (!defined($type)){$type = "null"}
			
			my $name 		= $clone->{name};
			my $short_name 	= $clone->{short_name};
			my $insert 		= $clone->{insert};
			my $aa_seq 		= $clone->{aa_seq};
			my $dna_seq		= $clone->{dna_obj}->{primary_seq}->{seq};
			my $date 		= $clone->{date};
			#if (!defined($date)){$date = "null"}
			
			my @values = ($id, $type, $name, $short_name, $insert, $aa_seq, $dna_seq, $date);
			
			my $sql = "REPLACE INTO $table_name($fields_holder)
						VALUES($values_holder)";
			my $sth = $dbh->prepare($sql);
			$sth->execute(@values);
			
		}
		
		
	}
} 

sub makeSQLtables{
	my $dbh = shift;
	my $tables = shift;
	
	foreach my $table(@$tables){
		my $fields = $table->{fields};
		my $name = $table->{name};
		
		my $sql = "
				CREATE TABLE IF NOT EXISTS $name(
				id			VARCHAR,
				type			VARCHAR,
				name			VARCHAR PRIMARY KEY,
				short_name	VARCHAR,
				insert_aa	VARCHAR,
				aa_seq		VARCHAR,
				dna_seq		VARCHAR,
				date			VARCHAR
				)";
		my $sth = $dbh->prepare($sql);
		$sth->execute();
	}
	return $dbh;
	
}

sub makeTables{
	my $table_names = shift;
	my $field_list = shift;
	my $tables = ();
	
	foreach my $name(@$table_names){
		my $table = Table->new($name);
		$table->{fields} = $field_list;
		push @$tables, $table; 
		
	}
	return $tables;
	
}

# Makes a hash full of M13Clone Objects
# Works through translating the DNA sequence, and filling
# elements of the M13Clone Object.
# unique clones housed in a hash,
# all clones in an array.
sub makeCloneHash{
	
	my $dnas = shift;
	my $clone_str = {};
	my $unique_clones_node = {};
	my $all_clones_node = {};
	
	# pass in first dna_object from main @dnas
	foreach my $dna_obj(@$dnas){
		
		my $mimName = $dna_obj->primary_id;
		
		my $clone = M13Clone->new($mimName);
		$clone->translate($dna_obj);
		my $insert = $clone->{insert};
		$unique_clones_node->{$insert} = $clone;
		$all_clones_node->{$mimName} = $clone;
	}
	
	$clone_str->{unique_clones} = $unique_clones_node;
	$clone_str->{all_clones} = $all_clones_node;
	
	return $clone_str;
}

# Moves file to subdirectory called 'processed'
sub moveFiles{
	my $dir_path = shift;
	my $seq_files = shift;
	my $processed_dir = $dir_path . "/processed";
	
	foreach my $file(@$seq_files){
		my $old_file_path = $dir_path . "/$file";
		my $new_file_path = $processed_dir . "/$file";
		move($old_file_path, $new_file_path);
	}
	
}

# Pulls all the DNA sequences from a file and pushes 
# into array @dnas
sub getDNAobjects{
	my $file = shift;
	my $format = shift;
	my $dnas = ();
	
	# Make a Sequence file reader
	my $seqin_obj = Bio::SeqIO->new(
		-file => "$file",
		-format => $format,
		-alphabet => 'dna'
		);
		
	while(my $dna_obj = $seqin_obj->next_seq()){	
		push(@$dnas, $dna_obj);
	}
	
	
	return $dnas;
	
}

# Gets all the DNA objects from all the files in the directory
sub getDNAs{
	my $dir_path = shift;
	my $seq_files = shift;
	my $dnas = ();

	foreach my $file (@$seq_files){
		# Fill dnas with dna objects from file;
		my $file_path = $dir_path . "/$file";
		my $dna_objects = getDNAobjects($file_path, $format);
		push @$dnas, @$dna_objects;
	}
	return $dnas;

}

#for getting all files from directory
sub getFiles{
	my $dir_path = shift;
	opendir my $dir, $dir_path or die "Cannot open directory: $!";
		my @seq_files = grep {/\.txt$/ || /\.fasta$/} readdir $dir;
	closedir $dir;
	return @seq_files;
}
