use Bio::TreeIO;
use strict;
use warnings;
use LWP::Simple;

if (scalar @ARGV != 2){
  #print STDOUT "\nError : Nombre d'arguments invalide ! \nUsage : perl $0 langues.new hwt.txt\n\n";
  #exit;
}

my $fichier_langue 	= "gray-atkinson_num3.new";
my $input           = new Bio::TreeIO(-file   => "$fichier_langue", -format => "newick");
my $tree = $input->next_tree;

my @ids_dest = ();
my @ids_source = ();
my $flag_source = 0;
my @ids = ();
my @ida = ();
my @idd = ();

foreach my $id (@ARGV){
  chomp($id);
  
  if ( $id =~ /^[0-9][0-9]/){
    $id =~ s/-[0-9]//g;
    my @nodes = $tree->find_node(-id => "$id");
    if( $flag_source == 0){
      push(@ids_source,$id); 
		  push(@ids,$nodes[0]);
      push(@ida,$nodes[0]);
    }
    elsif ($flag_source == 1){
      push(@ids_dest,$id);     
		  push(@idd,$nodes[0]);
      push(@ida,$nodes[0]);
    }
  }
  elsif ($id eq "-"){
    $flag_source = 1;
  }

}

my $d_node = $idd[0];
$d_node = $tree->get_lca(@idd) if ( scalar (@idd) > 1);
my $s_node = $ids[0];
$s_node = $tree->get_lca(@ids) if ( scalar (@ids) > 1);

my $a_node = $tree->get_lca(@ida);

#print STDOUT "\n" . join(",",@ids_source) . ":" . $s_node->height;
#print STDOUT "\n" . join(",",@ids_dest)   . ":" . $d_node->height;
#print STDOUT "\nall:" . $a_node->height;

print STDOUT ($a_node->height- $s_node->height + $a_node->height-$d_node->height);
