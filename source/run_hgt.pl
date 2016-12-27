#!/usr/bin/perl

use strict;
use File::Copy qw/ copy /;

require "postTraitement.pl";

#======================================================
#= VERIFICATION DES ARGUMENTS DE LA LIGNE DE COMMANDE
#======================================================
if( scalar @ARGV < 0){
    print "\nErreur\nusage : $0";
    exit 0;
}

my $cmd = "./hgt  ";
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

foreach my $elt (@ARGV){
  $cmd .= $elt . " ";
  if($elt =~ "bootstrap"){
    @tmp_tab = split("=",$elt);
    $bootstrap = $tmp_tab[1];
    chomp($bootstrap);
    #print STDOUT "$bootstrap";
  }
  if($elt =~ "inputfile"){
    @tmp_tab = split("=",$elt);
    $inputfile = $tmp_tab[1];
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

#=== LECTURE DE L'ARBRE D'ESPECES ===
open(IN,"$inputfile") or  die ("Cannot open $inputfile");
my @trees_tab = <IN>;
close(IN);
  
open(IN,"list_trans_custom.txt") or die ("Cannot open list_translations.txt");
my @lTrans = <IN>;
close(IN);
  
  my $posLtrans = 1;
  #===========================================================================
  #======================== EXECUTION DU PROGRAMME ===========================
  #===========================================================================
 
  $cmd .= "-inputfile=$tmp_input";
  
  my $nbTrees = scalar @trees_tab - 1;
  my $mot;
  if((scalar @trees_tab < 2)){
    exit_program(-1,$return_file,"PERL : nombre d'arbres invalide");
  }
  
  for (my $i=1;($i< scalar @trees_tab);$i++){
    #print STDOUT "\n==== $i\n";
    my $langue_tree = $trees_tab[0];
    #print IN $langue_tree;

    if($trees_tab[$i] =~ /^=>/){
      $currentCpt=0;

		  $mot = $trees_tab[$i]; chomp($mot);
      print STDOUT "\n>" . $trees_tab[$i++];
		  print STDERR "\n" . 	$mot;
      printToFile("$mot\n",$all_hgt,">>"); # if ($mot =~ "FRUIT") ;
      $mot =~ s/=> //g;
    }
    $currentCpt++;
    $currentInputFile  = "$data_path/" . $mot . "-input-$currentCpt"  . ".txt";
    $currentOutputFile = "$data_path/" . $mot . "-output-$currentCpt" . ".txt";
    $currentTransFile  = "$data_path/" . $mot . "-trans-$currentCpt"  . ".txt";
    $currentMotFile    = "$data_path/" . $mot . "-mot-$currentCpt"    . ".txt";
 
    if($reset eq "true"){
      unlink($currentInputFile);
      unlink($currentOutputFile);
      unlink($currentTransFile);
      unlink($currentMotFile);
    }

    my $trans_tree =  $trees_tab[$i];
    
    my $langue_tmp;
    
    #print STDERR "\n$trans_tree";
    for(my $it=1;$it<=84;$it++){    
	    for(my $iu=1;$iu<=5;$iu++){
			 
		    if ($it > 9){
			    $langue_tmp = $it . "-" . $iu ;
			  }
			  else{
	        $langue_tmp = "0" . $it . "-" . $iu;
			  }
			  if( index($trans_tree,$langue_tmp) >= 0){
			    #print STDERR "\nil est dans l'arbre de gene : $langue_tmp (" . index($trans_tree,$langue_tmp) . ")";
			  }
			  else{
			    #print STDERR "\nil n'est pas dans l'arbre de gene : $langue_tmp (" . index($trans_tree,$langue_tmp) . ")";
			    $langue_tree =~ s/$langue_tmp\:1\.0,//g;
			  }
	    }
    }
    
    #== Fichier input
    if(!-e $currentInputFile){
      open(OUT,">$currentInputFile") or  die ("Cannot open $currentInputFile");
      print OUT $langue_tree;
      print OUT $trees_tab[$i];
      close(OUT); 
	  }
    copy($currentInputFile,$tmp_input) or die($!);

    #== Fichier de traductions
  	open(OUT,">$currentTransFile") || die "Cannot open $currentTransFile";
	  $nbLtrans = $lTrans[$posLtrans++];
	  my $lTrans_ligne;
	  for(my $i1=0;$i1<$nbLtrans;$i1++){
		  $lTrans_ligne = $lTrans[$posLtrans++];
		  chomp($lTrans_ligne);
		  print OUT $lTrans_ligne;
		  if($i1 < ($nbLtrans-1)){
			  print OUT "\n";
		  }
	  }
	  close(OUT);
    copy($currentTransFile,"translations.txt") or die ($!);
  


  chomp($nbLtrans);
  
  # if( ($mot eq "MOTHER") ){
	  $nbExec++;
    unlink("hgtplus.txt");	
    
    unlink("mot.new");
    #== fichier de sortie
    print STDERR  "\n$currentOutputFile -- $currentOutputFile" ;
    if(!-e $currentOutputFile or !-e  $currentOutputFile ){
      print STDERR "\nPERL : $nbExec - [$mot - $nbLtrans] - $cmd";
      execute_hgt($cmd);
      copy($results,$currentOutputFile) or die($!);
      if(-e "mot.new"){
        copy("mot.new",$currentMotFile) or die($!);
      }
    }
    copy($currentOutputFile,$results) or die($!);

    if (-e $currentMotFile){
      copy($currentMotFile,"mot.new") or die($!);
    
      #print STDERR `cat mot.new`;
      # print STDERR "\n $currentOutputFile - $results";
      my $dibon =  `perl postTraitement.pl`;
    
      filtrerResultats($results, $hgtplus, $all_hgt);
    }
    #my $bidon = <>;
    #   }
}	

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


