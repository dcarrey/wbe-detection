#!/usr/bin/perl

use strict;
use File::Copy qw/ copy /;
use Bio::TreeIO;
use PostTraitement;
use Data::Dumper;


my $header = <<'END_HEADER';
===============================================
= WBE-Detection
= Authors : Alix Boc and Vladimir Makarenkov
= Version : 1.0
= Date : January 2017
===============================================
END_HEADER

#======================================================
#= VERIFICATION DES ARGUMENTS DE LA LIGNE DE COMMANDE
#======================================================
if( scalar @ARGV < 1){
    print "\nErreur\noptions :";
    print STDOUT "\n-inputfile\t  : see README for details";
    print STDOUT "\n-c1\t  : see README for details";
    print STDOUT "\n-c2\t  : see README for details";
    print STDOUT "\n-blk\t  : see README for details";
    print STDOUT "\n-mininternalnodes : minimum internal nodes (see README for details)";
    print STDOUT "\n-minexternalnodes : minimum external nodes (see README for details)";
    print STDOUT "\n-verbose\t  : see README for details";
    print STDOUT "\n-blk\t  : see README for details";
    print STDOUT "\n\n";
    exit 0;
}

my %parameters = (
  outputfile => {value => "output.txt", description =>"" },
  mininternalnodes => {value => 3, description =>"" },
  minexternalnodes => {value => 2, description =>"" },
  blk => {value => 0.20, description =>"" },
  c1 => {value => 1505, description =>"" },
  c2 => {value => 1.5, description =>"" },
  clean => {value => "yes", description =>"" },
  verbose => {value => "0", description =>"" },
  path => {value => "./", description =>"" },
  inputfile => {value => "", description =>""}
);

readParameters(\%parameters,@ARGV);

my $results    = $parameters{"path"}{"value"} . "results.txt";
my $hgtplus    = $parameters{"path"}{"value"} . "hgtplus.txt";
my $tmp_input  = $parameters{"path"}{"value"} . "_wbe_input.txt";
my $returnFile = $parameters{"path"}{"value"} . "return.txt";
my $logFile    = $parameters{"path"}{"value"} . "wbe.log";
my $filteredTree = $parameters{"path"}{"value"} . "_filteredLangueTree.new";
my $inputFile = $parameters{"path"}{"value"}.$parameters{"inputfile"}{"value"};
my $outputFile = $parameters{"path"}{"value"}.$parameters{"outputfile"}{"value"};
my $translationsFile = $parameters{"path"}{"value"}."_translations.txt";
my $c1 = $parameters{"c1"}{"value"};
my $c2 = $parameters{"c2"}{"value"};
my $blk = $parameters{"blk"}{"value"};
my $cmd = "./hgt ";

unlink($outputFile);
unlink($hgtplus);

#== read input file
my $languageTree = getTree($inputFile,"language_tree");
my $wordTree = getTree($inputFile,"word_tree");
my %groups = getGroups($inputFile,"group_content");
my $translations = getTranslations($inputFile,"translations",newickToBioTree($wordTree));

save_to_file($languageTree . "\n" . $wordTree, $tmp_input);
save_to_file($translations,$translationsFile);

#=== ADD ROOT TO TREES ===
$cmd = "./hgt  -inputfile=$tmp_input -addroot=yes -speciesroot=midpoint -generoot=midpoint > $logFile";
execute_hgt($cmd);
$wordTree = getTree("_wbe_word.new","");

save_to_file($languageTree . "\n" . $wordTree, $tmp_input);

#===========================================================================
#======================== EXECUTION DU PROGRAMME ===========================
#===========================================================================
$cmd = "./hgt  -inputfile=$tmp_input -outputfile=$outputFile -translationsfile=$translationsFile -constraints=3 -blk=$blk -c1=$c1 -c2=$c2 -speciesroot=midpoint -generoot=midpoint >> $logFile";

#print STDERR "\nPERL : $cmd";
execute_hgt($cmd);
my $filteredlanguageTree = getTree("_wbe_languages.new","");
my $postTraitement = new PostTraitement(newickToBioTree($languageTree),
                                        newickToBioTree($filteredlanguageTree),
                                        \%groups,
                                        $parameters{"mininternalnodes"}{"value"},
                                        $parameters{"minexternalnodes"}{"value"},
                                        $results,
                                        $hgtplus,
                                        );
$postTraitement->findAdditionnalsWBE();
saveResultats($hgtplus, $outputFile);

if ($parameters{"clean"} eq "yes"){
  deleteTempFiles(qw/speciesRoot.txt input_.txt geneRoot.txt _wbe_*/);
}
exit_program(0,$returnFile,"");

#===============================================================================
#=============================== FUNCTIONS =====================================
#===============================================================================
sub par_num {return $a <=> $b}

sub deleteTempFiles{
  foreach my $file (@_){
    unlink($file);
  }
}

sub read_line{
  my ($line) = @_;
  chomp($line);
  return $line;
}

sub exit_program{
  my($val,$file,$message) = @_;
  open(RET,">$file") || die "Cannot open $file";
  print RET $val;
  close(RET);
}

sub execute_hgt{
    my ($cmd) = @_;
    my $retour = 0;
    print STDOUT "\n$cmd";
    system($cmd);
}

#
# Apres une exécution de hgt-detection, on filtre les résultats
# en considérant les tranferts ajoutés au niveau des noeuds internes
#
sub saveResultats {
  my ($hgtplus,$outputFile) = @_;
  return if(! -f $hgtplus);
  printToFile($header, $outputFile, ">");
  printToFile("\nList of word borrowing events found :\n\n", $outputFile, ">>");
  open (IN,$hgtplus) or die($! . "($hgtplus)");
  while( my $ligne =<IN>) {
    chomp($ligne);
    printToFile("$ligne\n", $outputFile, ">>") if ($ligne !~ /^$/);
  }
  close (IN);
}

sub printToFile{
  my ($content, $file, $mode) = @_;
  open(OUT, "$mode$file") or die($! . "($results)");
  print OUT $content;
  close(OUT);
}

sub getTrees{
  my $file = $_[0];
  open(my $io,$file) or die("Cannot open $file");
  my $treeio = Bio::TreeIO->new(-format => 'newick', -fh => $io);
  my $tree1 = $treeio->next_tree;
  my $tree2 = $treeio->next_tree;
  close($io);
  return ($tree1,$tree2);
}

sub newickToBioTree{
  my $tree = $_[0];
  open(my $io,'<',\$tree) or die("Cannot read $tree");
  my $treeio = Bio::TreeIO->new(-format => 'newick', -fh => $io);
  my $tree1 = $treeio->next_tree;
  close($io);
  return $tree1;
}

sub getTree{
  my ($file,$key) = @_;
  my $tree="";
  open(my $io,$file) or die("Cannot open $file");
  if ($key eq ""){
    $tree = <$io>;
    chomp($tree);
  }
  else{
    while(my $ligne =<$io>){
      chomp($ligne);
      if( $ligne eq "$key:"){
        $tree = <$io>;
        chomp($tree);
      }
    }
  }
  close($io);
  return $tree;
}

sub filterTrees{
  my ($t1,$t2) = @_;

  my $tree1 = $t1->clone();
  my $tree2 = $t2->clone();

  my @nodes = ();
  foreach my $node ($tree1->get_nodes()){
    if($node->is_Leaf){
      if($tree2->findnode_by_id($node->id()) eq ""){
        push @nodes, $node;
      }
    }
  }
  foreach my $node (@nodes){
    $tree1->remove_Node($node);
    $tree1->contract_linear_paths();
  }
  return ($tree1,$tree2);
}

sub getTranslations{
  my ($file,$key,$tree) = @_;
  my $content  = "";
  my $nbLeaves = scalar ( $tree->get_leaf_nodes());
  my $cpt=0;
  my @leaves = ();
  my $start_reading = 0;
  open(IN,$file) or die("Cannot open $file");
  while(my $line =<IN>){
    chomp($line);
    if($start_reading==1 and $line !~ /^$/){
      my @tmp = split("=",$line);
      $content .= $tmp[0] . " " . $tmp[1] . "\n";
      push @leaves, $tmp[0];
    }
    $start_reading = 1 if($line eq "$key:");
    $start_reading = 0 if($line =~ /^$/);
  }
  for my $node ($tree->get_leaf_nodes()){
    my $trouve = 0;
    for my $leaf (@leaves){
      if($leaf eq $node->id_output){
        $trouve=1;
      }
    }
    if ( $trouve == 0 and $node->id_output ne "Root"){
      die("Error : No translation for word leaf " . $node->id_output . "\n");
    }
  }
  return $content;
}

sub getGroups{
  my ($file,$key) = @_;
  my %groups = ();
  my $start_reading = 0;
  open(IN,$file) or die("Cannot open $file");
  while(my $line =<IN>){
    chomp($line);
    if($start_reading==1 and $line !~ /^$/){
      my @tmp = split("=",$line);
      my $groupId = $tmp[0];
      chomp($groupId);
      for my $elt (split(",",$tmp[1])){
        $groups{$elt} = $groupId;
      }
    }
    $start_reading = 1 if($line eq "$key:");
    $start_reading = 0 if($line =~ /^$/);
  }
  return %groups;
}

sub save_to_file{
  my ($content,$file) = @_;
  open(OUT,">$file") or  die ("Cannot open $file");
  print OUT $content;
  close(OUT);
}

sub readParameters{
  my ($parameters,@line) = @_;
  foreach my $elt (@line){
    my ($key,$value) = ($elt =~ /-([^=]+)=([^\s]+)/);
    die("\nunknown option : $key\n") if (!exists $$parameters{$key});
    $$parameters{$key}{"value"} = $value;
  }
}
