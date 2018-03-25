package Table;

use strict;
use warnings;

# Tables in an SQL database

sub new{
	my $class = shift;
	my $name = shift;
	
	my $self = {};
	bless $self, $class;
	
	$self->_initialize($name);
	
	return $self;
}

sub _initialize{
	my $self = shift;
	my $name = shift;
	
	$self->{name} = $name;
	$self->{fields} = ();
}

sub getHeader{
	my $self = shift;
	
	my $fields = $self->{fields};
	my $column_header = join("\t",@$fields);
	
	return $column_header;
}
+
1;