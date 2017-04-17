use Bio::TreeIO;
use strict;
use warnings;
use LWP::Simple;

my $fichier_langue 	= $ARGV[0];
my $input           = new Bio::TreeIO(-file   => "$fichier_langue", -format => "newick");
my $tree            = $input->next_tree;

$tree->reroot(trouverRacine("Root",$tree));

my @ids_dest = ();
my @ids_source = ();
my $flag = "";
my @ids = ();
my @ida = ();
my @idd = ();

foreach my $id (@ARGV){
  chomp($id);
  if ( $id eq "_src_" or $id eq "_dest_" ){
    $flag = $id;
  }
  elsif ( $flag eq "_src_" or $flag eq "_dest_" ){
    my @nodes = $tree->find_node(-id => "$id");
    if( $flag eq "_src_"){
      #print "\ntrouve(src)=$id";
      push(@ids_source,$id);
		  push(@ids,$nodes[0]);
      push(@ida,$nodes[0]);
    }
    elsif ($flag eq "_dest_"){
      #print "\ntrouve(dest)=$id";
      push(@ids_dest,$id);
		  push(@idd,$nodes[0]);
      push(@ida,$nodes[0]);
    }
  }
}

my $d_node = $idd[0];
$d_node = $tree->get_lca(@idd) if ( scalar (@idd) > 1);
my $s_node = $ids[0];
$s_node = $tree->get_lca(@ids) if ( scalar (@ids) > 1);

my $a_node = $tree->get_lca(@ida);

print STDOUT ($a_node->height- $s_node->height + $a_node->height-$d_node->height);


sub trouverRacine{

  my ($racine,$tree) = @_;
  my @nodes = $tree->find_node(-id => "$racine");
  # print "\nTrouve : " . $nodes[0]->id_output;
  return $nodes[0];
}
