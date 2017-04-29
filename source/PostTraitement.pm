package PostTraitement;

use Bio::TreeIO;
use Array::Utils qw(:all);
use Data::Dumper;
use strict;
use warnings;

use constant FALSE => 0;
use constant TRUE  => 1;
use constant DATE_GREEK  => 2400;
use constant DEBUG => 0; # 0,1,2

no warnings 'experimental::smartmatch';

my %tabGroup = ();
my $tree_langue;
my $tree_mot;

#= Constructor
sub new {
  my $class = shift;
  my $self = {
                LANGUE_TREE => shift,
                FILTERED_TREE => shift,
                TAB_GROUP => shift,
                MIN_INTERNAL_NODES => shift,
                MIN_EXTERNAL_NODES => shift,
                RESULTS_FILE => shift,
                wbePlusFile => shift,
             };
  return bless $self, $class;
}

sub findAdditionnalsWBE{
  my $self = shift;
  $tree_langue = $self->{LANGUE_TREE};
  $tree_mot    = $self->{FILTERED_TREE};
  %tabGroup    = %{$self->{TAB_GROUP}};
  $tree_langue->reroot(trouverRacine("Root",$tree_langue));
  $tree_mot->reroot(trouverRacine("Root",$tree_mot));

  my @results = ();
  @results = lectureTransfertsOriginaux($self->{RESULTS_FILE});
  @results = rechercheTransfertsSupplementaires($tree_mot,$self->{MIN_EXTERNAL_NODES},$self->{MIN_INTERNAL_NODES},@results);
  @results = ajustementTransfertsDates(@results);
  @results = ajustementTransferts2(@results);
  @results = ajustementTransferts(@results);
  @results = ajustementTransfertsDatesEtape2(@results);
  @results = ajustementTransferts05(@results);
  saveAdditionnalsWBE($self->{wbePlusFile},@results);
 }

sub saveAdditionnalsWBE{
  my ($filename,@results) = @_;
  open (OUT, ">$filename") or die($!);
  foreach my $item (@results){
    if ($item->{status} ne "del"){
      print OUT "\n" . join(" ", @{$item->{source}}) . " -> " . join(" ",@{$item->{destination}});
      print OUT " -> " . $item->{fact} if($item->{fact} ne "0");
    }
  }
  close(OUT);
}

sub printTransferts{
  my ($title,@results) = @_;
  print STDERR "$title";
  foreach my $toFind (@results) {
    print STDERR "\n" . join(" ",@{$toFind->{source}}) . " -> " . join(" ", @{$toFind->{destination}}) . " [". $toFind->{fact} ."]  (" . $toFind->{status} .")";
  }
}

#
# Lecture des résultats originaux
#
sub lectureTransfertsOriginaux{
  my $filename = $_[0];
  my @results = ();

  open(IN, $filename) or die($! . "($filename)");
  <IN>;<IN>;
  my $temoin=0;
  my $source = "";
  my $dest = "";
  my $br = "";

  while (my $ligne = <IN>){
    chomp($ligne);
    if($ligne =~ /^\d+$/){
      my $nbWbe = $ligne;
      for(my $i=0;$i<$nbWbe;$i++){
        <IN>;
        $source = <IN>;
        $dest = <IN>;
        chomp($source);
        chomp($dest);
        <IN>;<IN>;
        my @ids_fils1 = sort(split(" ", $source));
        my @ids_fils2 = sort(split(" ", $dest));
        push @results , {source => \@ids_fils1, destination => \@ids_fils2, fact=>"1.0", status=>"init"};
      }
    }
    <IN>;
  }
  close(IN);
  printTransferts("\n\nLecture des transferts initiaux",@results) if (DEBUG > 0);
  return @results;
}

#
# Ajustement des transferts en fonction des dates des sous-arbres
#
sub ajustementTransfertsDates{
  my @results = @_;
  my @results2 = ();
  foreach my $item (@results){

    next if($item->{status} eq "del");
    my @ids_fils1 = @{$item->{source}};
    my @ids_fils2 = @{$item->{destination}};
    my $status = $item->{status};
    my $fact = "0";

    if( (scalar(@ids_fils1) > 0) and (scalar(@ids_fils2) >0 ) ){
      #my @feuilles1 = map { my $v=$_; $v =~ s/-[0-9]+//g; $v } @ids_fils1;
      #my @feuilles2 = map { my $v=$_; $v =~ s/-[0-9]+//g; $v } @ids_fils2;

      my $node1 = lcaDateTree(@ids_fils1);
      my $node2 = lcaDateTree(@ids_fils2);

      my $isOlder =  1;
      $isOlder = subtreeIsOlder(join(" ",@ids_fils1),join(" ",@ids_fils2)) if(($status eq "plus"));
      print STDERR "\n" . join(" ",@ids_fils1) . " -> " . join(" ",@ids_fils2) . ":" . $isOlder if (DEBUG > 1);

      if( $isOlder == 1 ){
        $fact = "1.0" if($status eq "plus");
        push @results2 , {source => \@ids_fils1, destination => \@ids_fils2, fact =>$fact, status=>$status};
        #print STDERR "\n2:" . join(" ", @ids_fils1) if (DEBUG);
      }
      elsif( $isOlder == -1  ){
        $fact = "1.0" if($status eq "plus");
        push @results2 , {source => \@ids_fils2, destination => \@ids_fils1, fact =>$fact, status=>$status};
        #print STDERR "\n2:" . join(" ", @ids_fils1) if (DEBUG);
      }
      else{
        $fact = "0.5" if($status eq "plus");
        push @results2 , {source => \@ids_fils2, destination => \@ids_fils1, fact =>$fact, status=>$status};
        push @results2 , {source => \@ids_fils1, destination => \@ids_fils2, fact =>"0.5", status=>$status} if($status eq "plus");
        #print STDERR "\n3:" . join(" ", @ids_fils1) if (DEBUG);
      }
    }
  }
  printTransferts("\n\nGestion des dates",@results2) if (DEBUG > 0);
  return @results2;
}

#
# Ajustement des transferts en fonction des dates des sous-arbres
#
sub ajustementTransfertsDatesEtape2{
  my @results = @_;
  my @results2 = ();
  foreach my $item (@results){

    next if($item->{status} eq "del");
    my @ids_fils1 = @{$item->{source}};
    my @ids_fils2 = @{$item->{destination}};
    my $status = $item->{status};
    my $fact = $item->{fact};

    if( (scalar(@ids_fils1) > 0) and (scalar(@ids_fils2) >0 ) ){
      my @feuilles1 = @ids_fils1;
      my @feuilles2 = @ids_fils2;

      my $node1 = lcaDateTree(@feuilles1);
      my $node2 = lcaDateTree(@feuilles2);

      my $isOlder =  1;
      $isOlder = subtreeIsOlder(join(" ",@ids_fils1),join(" ",@ids_fils2)) if(($status eq "plus"));

      if ( $isOlder == 1){
        push @results2 , {source => \@ids_fils1, destination => \@ids_fils2, fact=>$fact, status=>$status};
      }
      else{
        push @results2 , {source => \@ids_fils2, destination => \@ids_fils1, fact=>$fact, status=>$status};
      }
    }
  }
  printTransferts("\n\nGestion des dates etape 2",@results2) if (DEBUG > 0);
  return @results2;
}


#
#
#
sub rechercheTransfertsSupplementaires{

  my ($tree_mot,$MIN_EXTERNAL_NODES,$MIN_INTERNAL_NODES,@results) = @_;
  my @mot_feuilles = getFeuilles($tree_mot->get_root_node());
  #
  # RECHERCHE DES TRANSFERTS SUPPLEMENTAIRES
  #

  foreach my $parent ( $tree_mot->get_nodes()){

    if( is_internal_node($parent) ){
      #== Recherche dans l'arbre du mot
      my @ids_parent = getFeuilles($parent);
      my ($fils1,$fils2) = getFils($parent);
      my @ids_fils1 = getFeuilles($fils1);
      my @ids_fils2 = getFeuilles($fils2);
      # print STDERR "\n=========================================";
      # print STDERR "\nParent  (" . $parent->internal_id . ") : " . join(",",@ids_parent);
      # print STDERR "\n\tFils 1  (" . $fils1->internal_id . ") : " . join(",",@ids_fils1);
      # print STDERR "\n\tFils 2  (" . $fils2->internal_id . ") : " . join(",",@ids_fils2);

      #== Recherche dans l'arbre de langues
      my $langue_parent = trouverNoeudCorrespondant( @ids_parent);
      # print STDERR "\n\nParent (Langue) = " . $langue_parent->internal_id . "->" .  join(",",getFeuilles($langue_parent));
      my $langue_fils1 = trouverNoeudCorrespondant( @ids_fils1);
      # print STDERR "\nFils (Langue) = " . $langue_fils1->internal_id . "->" .  join(",",getFeuilles($langue_fils1));
      my $langue_fils2 = trouverNoeudCorrespondant( @ids_fils2);
      # print STDERR "\nFils (Langue) = " . $langue_fils2->internal_id . "->" .  join(",",getFeuilles($langue_fils2));

      my $nbNoeud1 = nbNoeudIntermediaire($langue_parent,$langue_fils1);
      my $nbNoeud2 = nbNoeudIntermediaire($langue_parent,$langue_fils2);

      my @feuilles_fils1 = ();
      my @feuilles_fils2 = ();

      foreach my $f (@mot_feuilles){

        foreach my $f1 (getFeuilles($langue_fils1)){
          push ( @feuilles_fils1, $f1)  if ($f eq $f1);
        }
        foreach my $f2 (getFeuilles($langue_fils2)){
          push ( @feuilles_fils2, $f2)  if ($f eq $f2);
        }
      }

      my $sontDansLeMemeGroupe = &memeGroupe(\@ids_fils1,\@ids_fils2);


      print STDERR "\n" . join(" ",@ids_fils1) . "[" . $nbNoeud1 . "]" . " -> " . join(" ",@ids_fils2) . "[" . $nbNoeud2 . "]" . ":" . $sontDansLeMemeGroupe if (DEBUG > 1);

      if (
        ( ( $sontDansLeMemeGroupe == FALSE ) and  (($nbNoeud1 + $nbNoeud2) >= $MIN_EXTERNAL_NODES ) ) or
        ( ( $sontDansLeMemeGroupe == TRUE  ) and  (($nbNoeud1 + $nbNoeud2) >= $MIN_INTERNAL_NODES ) )
      ){
        print STDERR "\n" . join(" ",@ids_fils1) . " -> " . join(" ",@ids_fils2) if (DEBUG > 1);
        #if(($nbNoeud1 + $nbNoeud2) >= $MIN_EXTERNAL_NODES ){
        if( (scalar(@ids_fils1) > 0) and (scalar(@ids_fils2) > 0 ) ){

          my @feuilles1 =  do { my %seen; grep { !$seen{$_}++ }  getFeuillesUnique($fils1)};
          my @feuilles2 =  do { my %seen; grep { !$seen{$_}++ }  getFeuillesUnique($fils2)};

          my $node1 = lcaDateTree(@feuilles1);
          my $node2 = lcaDateTree(@feuilles2);

          my $isOlder =  1;
          $isOlder = subtreeIsOlder(join(" ",@ids_fils1),join(" ",@ids_fils2));

          if( $isOlder == 1 ){
            push @results , {source => \@ids_fils1, destination => \@ids_fils2, status=>"plus", fact=>""};
          }
          elsif( $isOlder == -1 ){
            push @results , {source => \@ids_fils2, destination => \@ids_fils1, status=>"plus", fact=>""};
          }
          else{
            push @results , {source => \@ids_fils1, destination => \@ids_fils2, status=>"plus", fact=>""};
          }
        }
      }
    }
  }
  printTransferts("\n\nAjout des transferts supplémentaires (plus)",@results) if (DEBUG > 0);
  return @results;
}

sub ajustementTransferts{
  my @results = @_;
  #
  # Suppression des sous transferts :
  # if (currentSource,currentDest) is included into anotherSource then
  #   remove currentDest from anotherSource
  # end if
  #
  my $event = 1;
  while( $event == 1 ){
    $event = 0;
    foo:{
      foreach my $toFind (@results) {
        my @tabToFind = (@{$toFind->{source}},@{$toFind->{destination}});

        foreach my $toCompare (@results) {
          #= (source + dest) est un sous-groupe d'une (source) ?
          if ( ($toFind->{fact} ne "0.5") && ($toCompare->{fact} ne "0.5") ) {

            my @minus = array_minus(@tabToFind,@{$toCompare->{source}});
            if ( @minus == 0){
              my @newSource =  array_minus(@{$toCompare->{source}},@{$toFind->{destination}});
              if (DEBUG > 1){
                print STDERR "\n\n>>" . join(" ",@tabToFind) . " is included into " . join(" ", @{$toCompare->{source}});
                print STDERR "\nOn supprime " . join(" ",@{$toFind->{destination}});
                print STDERR "\n==> New Source : " . join(" ",@newSource);
              }
              $toCompare->{source} = \@newSource;
              $event=1;
              last foo;
            }
            else{
              #= (source + dest) est un sous-groupe d'une (destionation) ?
              my @minus = array_minus(@tabToFind,@{$toCompare->{destination}});
              if ( @minus == 0){
                my @newDest =  array_minus(@{$toCompare->{destination}},@{$toFind->{destination}});
                if (DEBUG > 1){
                  print STDERR "\n2>>" . join(" ",@tabToFind) . " is included into " . join(" ", @{$toCompare->{destination}});
                  print STDERR "\nOn supprime " . join(" ",@{$toFind->{destination}});
                  print STDERR " ==> New dest : " . join(" ",@newDest);
                }
                $toCompare->{destination} = \@newDest;
                $event=1;
                last foo;
              }
            }
          }
        } #end - foreach
      } #end - foreach
    }
  }
  printTransferts("\n\nAjustement 1",@results) if (DEBUG > 0);
  return @results;
}


sub ajustementTransferts2{
  my @results = @_;
  #
  # Suppression des sous transferts :
  # if (currentSource,currentDest) is included into anotherSource then
  #   remove currentDest from anotherSource
  # end if
  #
  my $event = 1;
  while( $event == 1 ){
    $event = 0;
    foo:{
      foreach my $toFind (@results) {
        next if($toFind->{status} eq "del");

        foreach my $toCompare (@results) {
          next if($toCompare->{status} eq "del");

          if ( ($toFind->{status} eq "init") and ($toCompare->{status} eq "plus") and ($toCompare->{fact} ne "0.5") ){
            my @minus_source = array_minus(@{$toFind->{source}},@{$toCompare->{source}});
            my @minus_destination = array_minus(@{$toFind->{destination}},@{$toCompare->{destination}});
            my @minus_destination_inverse = array_minus(@{$toCompare->{destination}},@{$toFind->{destination}});
            #print STDERR "\n\n>>" . join(" ",@{$toFind->{source}}) . " => " . join(" ", @{$toCompare->{source}}) . " : " .join(" ", @minus_source)  ;
            #print STDERR "\n>>" . join(" ",@{$toFind->{destination}}) . " => " . join(" ", @{$toCompare->{destination}}) . " : " .join(" ", @minus_destination)  ;
            if( ((@minus_source == 0) and (@minus_destination==0)) or ((@minus_source==0) and (@minus_destination_inverse==0)) ){
              if (DEBUG > 1){
                print STDERR "\n>>" . join(" ",@{$toFind->{source}}) . " is included into " . join(" ", @{$toCompare->{source}});
                print STDERR "\n>>" . join(" ",@{$toFind->{destination}}) . " is included into " . join(" ", @{$toCompare->{destination}});
              }
              my @newDest = array_minus(@{$toCompare->{destination}},@{$toFind->{destination}});
              #print STDERR "\n>> nouvelle destination : " . join(" ", @newDest);
              if(@newDest == 0){
                $toCompare->{status} = "del"
              }
              else{
                $toCompare->{destination} = \@newDest;
              }
              #print STDERR "\nfini";
              $event=1;
              last foo;
            }
            else{
              my @minus_source = array_minus(@{$toFind->{source}},@{$toCompare->{destination}});
              my @minus_destination = array_minus(@{$toFind->{destination}},@{$toCompare->{source}});
              my @minus_destination_inverse = array_minus(@{$toCompare->{source}},@{$toFind->{destination}});
              #print STDERR "\n\n2>>" . join(" ",@{$toFind->{source}}) . " => " . join(" ", @{$toCompare->{destination}}) . " : " .join(" ", @minus_source)  ;
              #print STDERR "\n2>>" . join(" ",@{$toFind->{destination}}) . " => " . join(" ", @{$toCompare->{source}}) . " : " .join(" ", @minus_destination)  ;
              if( ((@minus_source == 0) and (@minus_destination==0)) or ((@minus_source==0) and (@minus_destination_inverse==0)) ){
                my @newSource = array_minus(@{$toCompare->{destination}},@{$toFind->{source}});
                if (DEBUG > 1){
                  print STDERR "\n2>>" . join(" ",@{$toFind->{source}}) . " is included into " . join(" ", @{$toCompare->{destination}});
                  print STDERR "\n2>>" . join(" ",@{$toFind->{destination}}) . " is included into " . join(" ", @{$toCompare->{source}});
                  print STDERR "\n2>>Nouvelle source : " . join(" ", @newSource);
                }
                #if ( (@{$toFind->{source}} == 1) and (@{$toFind->{destination}} == 1)) {
                my $tmp =  $toFind->{source};
                $toFind->{source} =  $toFind->{destination} ;
                $toFind->{destination} = $tmp;
                #}
                if(@newSource == 0){
                  $toCompare->{status} = "del"
                }
                else{
                  $toCompare->{destination} = \@newSource;
                }
                #print STDERR "\nfini";
                $event=1;
                last foo;
              }
            }
          }
        }
      }
    }
  }
  printTransferts("\n\nAjustement 2",@results) if (DEBUG > 0);
  return @results;
}

sub ajustementTransferts05{
  my @results = @_;
  #
  # Suppression des sous transferts :
  # if (currentSource,currentDest) is included into anotherSource then
  #   remove currentDest from anotherSource
  # end if
  #
  my $event = 1;
  while( $event == 1 ){
    $event = 0;
    foo:{
      foreach my $toFind (@results) {
        next if($toFind->{status} eq "del");

        foreach my $toCompare (@results) {
          next if($toCompare->{status} eq "del");

          if ( ($toFind->{fact} eq "0.5") and ($toCompare->{fact} ne "0.5") ){
            if (DEBUG > 1){
              print STDERR "\n>>" . join(" ",@{$toFind->{destination}}) . " ? " . join(" ", @{$toCompare->{destination}});
            }

            if ( join(" ",@{$toFind->{destination}})  eq  join(" ",@{$toCompare->{destination}}) ){
              if (DEBUG > 1){
                print STDERR "\n>>" . join(" ",@{$toFind->{destination}}) . " = " . join(" ", @{$toCompare->{destination}});
              }
              $toFind->{status} = "del";
              $event=1;
              last foo;
            }
          }
        }
      }
    }
  }
  printTransferts("\n\nAjustement 0.5",@results) if (DEBUG > 0);
  return @results;
}


sub nbNoeudIntermediaire{

  my ( $parent , $fils) = @_;
  my $tmp = $parent;
  my $cpt=0;

  #print "\n" . $parent->internal_id  . "-" .  $fils->internal_id;
  if( $parent->internal_id eq $fils->internal_id){
    #print "\nfils = parent";
    return 0;
  }
  foreach my $node ( reverse $tree_langue->get_lineage_nodes($fils)){
    #print "\n" . $parent->internal_id  . "-" .  $node->internal_id;
    if( $node->internal_id ne $parent->internal_id){
      $cpt++;
    }
    else{
      last;
    }
    # print STDOUT "\n" . $node->internal_id . "->" .  join(",",getFeuilles($node));
  }
  #print "\nresultat=$cpt";
  return $cpt;
}

sub trouverRacine{
  my ($racine,$tree) = @_;
  my @nodes = $tree->find_node(-id => "$racine");
  return $nodes[0];
}

sub trouverNoeudCorrespondant{

  my @ids_parent = @_;
  my @ids = ();
  foreach my $id (@ids_parent){
    my @nodes = $tree_langue->find_node(-id => "$id");
		push(@ids,$nodes[0]);
  }
  my $langue_parent = $ids[0];
  #print STDOUT  "\n nb nodes = " . scalar (@ids);

  $langue_parent = $tree_langue->get_lca(@ids) if ( scalar (@ids) > 1);

  return $langue_parent;
}

sub lcaDateTree{

  my @ids_parent = @_;
  my @ids = ();
  foreach my $id (@ids_parent){
    my @nodes = $tree_langue->find_node(-id => "$id");
		push(@ids,$nodes[0]);
  }
  my $langue_parent = $ids[0];
  my $id_output = $langue_parent->id_output;
  $langue_parent = $tree_langue->get_lca(@ids) if ( scalar (@ids) > 1);
  return $langue_parent;
}


sub getFeuilles{

  my $node = $_[0];
  my @ids = ();

  if( $node->is_Leaf){
    push(@ids,$node->id_output);
  }
  foreach my $descendant ( $node->get_all_Descendents){
      if( $descendant->is_Leaf){
        push(@ids,$descendant->id_output);
      }
    }
  return sort(@ids);
}

sub getFeuillesUnique{

  my $node = $_[0];
  my @ids = ();
  my $id = -1;
  if( $node->is_Leaf){
    $id = $node->id_output;
    #$id =~ s/-[0-9]+//g;
    push(@ids,$id);
  }
  foreach my $descendant ( $node->get_all_Descendents){
      if( $descendant->is_Leaf){
        $id = $descendant->id_output;
        #$id =~ s/-[0-9]+//g;
        push(@ids,$id);
      }
    }
  return sort(@ids);
}

sub getFils{

  my $node = $_[0];
  my $fils1 = "";
  my $fils2 = "";
  foreach my $descendant ($node->each_Descendent){
    if( $fils1 eq ""){
      $fils1 = $descendant;
    }
    else{
      $fils2 = $descendant;
    }
  }
  return $fils1,$fils2;
}


sub subtreeIsOlder{
  my($subtree1,$subtree2) = @_;

  my @feuilles1 = split(" ",$subtree1);
  my @feuilles2 = split(" ",$subtree2);

  my $node1 = lcaDateTree(@feuilles1);
  my $node2 = lcaDateTree(@feuilles2);

  print STDERR "\n" . join(" ",@feuilles1) . " -> " . join(" ",@feuilles2)  if (DEBUG > 1);

  my $height1 = $node1->height;
  my $height2 = $node2->height;

  print STDERR "\nsubtreeIsOlder : [" . $height1 . "," . $height2 . "]"  if (DEBUG > 1);
  return 1 if($height1 > $height2);
  return 0 if($height1 == $height2);
  return -1 if($height1 < $height2);
}

sub is_internal_node{

  my $node = $_[0];
  return 1 if (!$node->is_Leaf and (scalar($node->each_Descendent) > 1));
  return 0;
}

sub memeGroupe(\@\@){

  my ($tab1,$tab2) = @_;
  foreach my $elt1 (@$tab1){
    my $e1 = $elt1;
    foreach my $elt2 (@$tab2){
      my $e2 = $elt2;
      my $groupe_e1 = 0;
      $groupe_e1 = $tabGroup{$e1} if( exists $tabGroup{$e1} );
      my $groupe_e2 = 0;
      $groupe_e2 = $tabGroup{$e2} if( exists $tabGroup{$e2} );
      if (  $groupe_e1 ne  $groupe_e2  ){
        return FALSE;
      }
    }
  }
  return TRUE;
}

1;
