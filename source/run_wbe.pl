#!/usr/bin/perl

use strict;
use File::Copy qw/ copy /;


#======================================================
#= VERIFICATION DES ARGUMENTS DE LA LIGNE DE COMMANDE
#======================================================
if( scalar @ARGV < 0){
    print "\nErreur\nusage : $0";
    exit 0;
}

my $cmd = "./hgt ";
my $inputfile = "";
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

my $languesFile="langues.new";
my $wordsFile="words.new";
my $translationsFile="translations.txt";

foreach my $elt (@ARGV){
  $cmd .= $elt . " ";
  if($elt =~ "bootstrap"){
    @tmp_tab = split("=",$elt);
    $bootstrap = $tmp_tab[1];
    chomp($bootstrap);
    #print STDOUT "$bootstrap";
  }
  if($elt =~ "languesfile"){
    @tmp_tab = split("=",$elt);
    $languesFile = $tmp_tab[1];
  }
  
  if($elt =~ "wordsfile"){
    @tmp_tab = split("=",$elt);
    $wordsFile = $tmp_tab[1];
  }

  if($elt =~ "translationsfile"){
    @tmp_tab = split("=",$elt);
    $translationsFile = $tmp_tab[1];
  }

  if($elt =~ "path"){
    @tmp_tab = split("=",$elt);
    $path = $tmp_tab[1];
  }
  if($elt =~ "viewtree"){
    @tmp_tab = split("=",$elt);
    $viewtree = $tmp_tab[1];
  }
}

die("Cannot open $languesFile") if ( !-f $languesFile );
die("Cannot open $wordsFile") if ( !-f $wordsFile );
die("Cannot open $translationsFile") if ( !-f $translationsFile );

$inputfile    = "$path" . "$inputfile";
#$cmd          = "usagers/" . $cmd;
my $results   = "$path" . "results.txt";
my $hgtplus   = "$path" . "hgtplus.txt";
my $all_hgt   = "$path" . "all_hgt.txt";
my $tmp_input = "$path" . "tmp_input.txt";
my $input_no_space    = "$path" . "input_no_space.txt";
my $return_file = "$path" . "return.txt";
my $output;
my $output_tmp;
my $outputWeb = "$path" . "outputWeb.txt";
my $RF;
my $BD;
my $LS;

my $data_path = "data";
my $currentCpt;
my $currentInputFile;
my $currentOutputFile;
my $currentTransFile;
my $currentMotFile;
my $nbLtrans;

mkdir($data_path);
printToFile("",$all_hgt,">");

my $languesTree = get_content_file($languesFile);
my $wordsTree = get_content_file($wordsFile);

my $posLtrans = 1;
#===========================================================================
#======================== EXECUTION DU PROGRAMME ===========================
#=========================================================================== 
$cmd .= "-inputfile=$tmp_input -translationsfile=$translationsFile";
  
#== Fichier input
open(OUT,">$tmp_input") or  die ("Cannot open $tmp_input");
print OUT $languesTree;
print OUT $wordsTree;
close(OUT); 

unlink("hgtplus.txt");	
    
print STDERR "\nPERL : $cmd";
execute_hgt($cmd);
`perl postTraitement.pl $tmp_input`;
filtrerResultats($results, $hgtplus, $all_hgt);

system("rm speciesRoot.txt input_.txt geneRoot.txt");    
exit_program($val_retour,$return_file,"PERL : fin normale du programme");
              
#===============================================================================
#=============================== FUNCTIONS =====================================
#===============================================================================
sub par_num {return $a <=> $b}

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
  print STDOUT "\nexit=>$message";
  exit;
}

sub execute_hgt{
    my ($cmd) = @_;
    my $retour = 0;
    print STDOUT system($cmd);
    print STDOUT "============== $retour ================";
}

#
# Apres une exécution de hgt-detection, on filtre les résultats
# en considérant les tranferts ajoutés au niveau des noeuds internes
#
sub filtrerResultats {

  my ($results, $hgtplus, $all_hgt) = @_; 

  return if((! -f $results) or (! -f $hgtplus));
 
  printToFile("1 cognat\n", $all_hgt, ">>");

  open(IN, $results) or die($!);
  my $temoin=0;
  my $source = "";
  my $dest = "";
  my $br = "";

  open (IN,$hgtplus) or die($!);
  while( my $ligne =<IN>) {
    chomp($ligne);
    printToFile("$ligne\n", $all_hgt, ">>") if ($ligne !~ /^$/);
  }  
  close (IN);
}


sub findElements{
  
  my( $hgtplus, $source, $dest ) = @_;
 

#  print STDERR "\n$source -> $dest";
  open(IN2, $hgtplus) or die($!);

  my $delElt = 0;
  
  while( my $ligne = <IN2>){
    chomp($ligne);
    if($ligne !~ /^$/){
      my ($source2, $dest2, $facteur) = split("->", $ligne);
      # print STDERR "\n$source2 -> $dest2";

       my $nb_source=0;
       my $nb_dest=0;
      
      #= Tous les elements sources doivent être retrouves
      foreach my $elt ( split(" ", $source)){
        chomp($elt);
        if($source2 =~ /$elt/){
          $nb_source++;
        }
      }
      foreach my $elt ( split(" ", $dest)){
        chomp($elt);
        if($dest2 =~ /$elt/){
          $nb_dest++;
        }
      }
      #     print STDERR "\n1)nb_source=$nb_source nb_dest=$nb_dest";
      if($nb_source == scalar( split(" ",$source) ) and $nb_dest == scalar( split(" ",$dest) )  ){
        close(IN2);
        return 1;
      }

      $nb_source=0;
      $nb_dest=0;

      #= Tous les elements sources doivent être retrouves
      foreach my $elt ( split(" ", $source)){
        chomp($elt);
        #print STDERR "\n\t$source2 <> $elt";
        if($dest2 =~ /$elt/){
          $nb_source++;
        }
      }
      foreach my $elt ( split(" ", $dest)){
        chomp($elt);
        if($source2 =~ /$elt/){
          $nb_dest++;
        }
      }
      # print STDERR "\n2)nb_source=$nb_source nb_dest=$nb_dest";
      if($nb_source == scalar( split(" ",$source) ) and $nb_dest == scalar( split(" ",$dest) )  ){
        close(IN2);
        return 1;
      } 

    }
  }

  close(IN2);
  
  return 0;
}


sub printToFile{
  my ($content, $file, $mode) = @_;
  open(OUT, "$mode$file") or die($!);
  print OUT $content;
  close(OUT);
}


sub get_content_file{
  my $filename = $_[0];
  my $content  = "";
  open(IN,$filename) or die("Cannot open $filename");
  my $content = <IN>;
  close(IN);
  return $content;
}
