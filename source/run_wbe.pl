#!/usr/bin/perl

use strict;
use File::Copy qw/ copy /;
use Bio::TreeIO;
use PostTraitement;

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
    print STDOUT "\n-mininternalnodes : minimum internal nodes (see README for details)";
    print STDOUT "\n-minexternalnodes : minimum external nodes (see README for details)";
    print STDOUT "\n\n";
    exit 0;
}

my $cmd = "./hgt ";
my $inputFile = "";
my $bootstrap = "no";
my $path = "./";
my $viewtree="no";
my @tmp_tab;
my @tmp_tab_init;
my %hgt;
my $ligne;
my $nbLines = 5;
my %hgt_number_tab;
my %hgt_description_tab;
my %hgt_compteur_tab;
my %hgt_criterion_tab;
my %hgt_nbHGT_tab;
my @hgt_pos;
my @hgt_pos2;
my $mode;
my $total_hgt;
my $total_trivial;
my $val_retour=0;   #= nombre de hgt trouve
my @hgt_tab;
my $nbExec=0;
my $reset="false";

my $translationsFile="_translations.txt";
my $outputFile="output.txt";
my $minInternalNodes=3;
my $minExternalNodes=2;

foreach my $elt (@ARGV){
  if($elt =~ "inputfile"){
    @tmp_tab = split("=",$elt);
    $inputFile = $tmp_tab[1];
  }
  if($elt =~ "minexternalnodes"){
    @tmp_tab = split("=",$elt);
    $minExternalNodes = $tmp_tab[1];
  }
  if($elt =~ "mininternalnodes"){
    @tmp_tab = split("=",$elt);
    $minInternalNodes = $tmp_tab[1];
  }
  if($elt =~ "path"){
    @tmp_tab = split("=",$elt);
    $path = $tmp_tab[1];
  }
}

my $results    = "$path" . "results.txt";
my $hgtplus    = "$path" . "hgtplus.txt";
my $all_hgt    = "$path" . "all_hgt.txt";
my $tmp_input  = "$path" . "tmp_input.txt";
my $returnFile = "$path" . "return.txt";
my $logFile    = "$path" . "wbe.log";
my $filteredTree = "$path" . "_filteredLangueTree.new";


unlink($path.$outputFile);
unlink($hgtplus);	

my($initial_langue_tree,$initial_word_tree) = getTrees($path.$inputFile);
my($filtered_langue_tree,$filtered_word_tree) = filterTrees($initial_langue_tree,$initial_word_tree);
my $content = $filtered_langue_tree->as_text('newick') ."\n".$filtered_word_tree->as_text('newick') ;
save_to_file($content, $tmp_input);

#print STDOUT $initial_langue_tree->as_text('newick') . "\n" . $initial_word_tree->as_text('newick');
#print STDOUT "\n\n" . $filtered_langue_tree->as_text('newick') . "\n" . $filtered_word_tree->as_text('newick');

save_to_file(getTranslations($path.$inputFile),$path.$translationsFile);

#===========================================================================
#======================== EXECUTION DU PROGRAMME ===========================
#=========================================================================== 
$cmd .= "-inputfile=$tmp_input -translationsfile=$path$translationsFile"; # > $logFile";
    
#print STDERR "\nPERL : $cmd";
execute_hgt($cmd);
my $postTraitement = new PostTraitement($initial_langue_tree,$filtered_langue_tree,$minInternalNodes,$minExternalNodes,$results,$hgtplus);
$postTraitement->findAdditionnalsWBE();
saveResultats($hgtplus, $outputFile);

deleteTempFiles(qw/speciesRoot.txt input_.txt geneRoot.txt/);
exit_program($val_retour,$returnFile,"");
             
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
  my ($file) = @_;
  my $content  = "";
  open(IN,$file) or die("Cannot open $file");
  <IN>;
  <IN>;
  while(my $line =<IN>){
    chomp($line);
    my @tmp = split(" ",$line);
    if( scalar @tmp == 3){
      $content .= $tmp[1] . " " . $tmp[2] . "\n";
    }
    if( scalar @tmp == 2){
      $content .= $tmp[0] . " " . $tmp[1] . "\n";
    }
  }
  return $content;
}

sub save_to_file{
  my ($content,$file) = @_;
  open(OUT,">$file") or  die ("Cannot open $file");
  print OUT $content;
  close(OUT); 
}
