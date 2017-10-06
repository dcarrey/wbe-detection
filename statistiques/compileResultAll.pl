#!/usr/bin/perl

use strict;
use Data::Dumper;

#==============================================================================
#== GESTION DES PARAMETRES : le premier parametre est le fichier resultat de ==
#== HGT-Detection et le second						     ==	
#==============================================================================
if (scalar @ARGV != 3){
	print STDOUT "\nNombre de parametres insuffisant";
	print STDOUT "\nUsage : $0 fichier_resultat(hwt.txt) fichierDeMots(all.mots) fichierDarbres(arbres.txt)"; 
	exit();
}

print STDOUT `date`;

my $input     = $ARGV[0];
my $input2    = $ARGV[1];
my $dbnewick  = $ARGV[2];

open (RES, "$input") 	or die ($!);
open (IN2, "$input2") or die ($!);
my @motsCategorie = <IN2>;
close(IN2);

for(my $i=0;$i<scalar @motsCategorie;$i++){
	my $element = $motsCategorie[$i];
	chomp($element);
	$motsCategorie[$i] = $element;
}
#print STDOUT "\nNombre de mots dans cette categorie = " . scalar @motsCategorie ;
# print STDERR join("-",@motsCategorie);


#==============================================================================
#== CONSTANTES ET VARIABLES
#==============================================================================
my $NBLANGUES = 84;
my $NBGROUPES = 12;
my %nomGroup = ( '1'  => "Celtic",
				 '2'  => "Italic",
				 '3'  => "French/Iberian",
				 '4'  => "West Germanic",
				 '5'  => "North Germanic",
				 '6'  => "Baltic",
				 '7'  => "Slavic",
				 '8'  => "Indic",
				 '9'  => "Iranian",
				 '10' => "Albanian",
				 '11' => "Greek",
				 '12' => "Armenian",
			   );
			   
my %tabGroup = ("01"=>'1',"02"=>'1',"03"=>'1',"04"=>'1',"05"=>'1',"06"=>'1',"07"=>'1',
		"08"=>'2',"09"=>'2',"10"=>'2',"11"=>'2',"17"=>'2',"18"=>'2',"19"=>'2',
		"12"=>'3',"13"=>'3',"14"=>'3',"15"=>'3',"16"=>'3',"20"=>'3',"21"=>'3',"22"=>'3',"23"=>'3',
		"24"=>'4',"25"=>'4',"26"=>'4',"27"=>'4',"28"=>'4',"29"=>'4',"37"=>'4',"38"=>'4',
		"30"=>'5',"31"=>'5',"32"=>'5',"33"=>'5',"34"=>'5',"35"=>'5',"36"=>'5',
		"39"=>'6',"40"=>'6',"41"=>'6',
		"42"=>'7',"43"=>'7',"44"=>'7',"45"=>'7',"46"=>'7',"47"=>'7',"48"=>'7',"49"=>'7',"50"=>'7',"51"=>'7',"52"=>'7',"53"=>'7',"54"=>'7',
		"55"=>'8',"56"=>'8',"57"=>'8',"58"=>'8',"59"=>'8',"60"=>'8',"61"=>'8',"62"=>'8',"63"=>'8',"64"=>'8',"65"=>'8',
		"73"=>'9',"74"=>'9',"75"=>'9',"76"=>'9',"77"=>'9',"78"=>'9',"79"=>'9',
		"80"=>'10',"81"=>'10',"82"=>'10',"83"=>'10',"84"=>'10',
		"66"=>'11',"67"=>'11',"68"=>'11',"69"=>'11',"70"=>'11',
		"71"=>'12',"72"=>'12'
		);

my @age_groupe               = qw /3000 1800 1200 1400 800 1400 1300 4200 4000 1500 3500 2600/;
my @nombre_langue_par_groupe = qw /7 7 9 8 7 3 13 11 7 5 5 2/;
my %new_total_inter_groupes = ();
my %new_total_inter_affectees = ();
for(my $i=1;$i<=$NBGROUPES;$i++){
	for(my $j=1;$j<=$NBGROUPES;$j++){
		$new_total_inter_groupes   {$i,$j} = 0;
		$new_total_inter_affectees {$i,$j} = 0;
	}
}

my %tab_mots_impliques_intra = (); 	#== tableau indiquant les traductions impliquees dans un transfert
my %tab_mots_impliques_entrant = ();
my %tab_mots_impliques_sortant = ();
my %tab_mots_impliques_inter = ();
my %tab_mots_affectes = ();

my %nb_mots_par_groupe_intra = ();
my %nb_mots_par_groupe_inter = ();
my %nb_mots_par_groupe_entrant=();
my %nb_mots_par_groupe_sortant=();

my %nb_mots_affectes_par_groupe_intra = ();
my %nb_mots_affectes_par_groupe_inter = ();
my %nb_mots_affectes_par_groupe_entrant = ();
my %nb_mots_affectes_par_groupe_sortant = ();
my %nb_mots_non_affectes_par_groupe_inter = ();
my $mot_accepte;
my $mot_courant;
my @tab_mots;
my %nb_mots_affectes_sortant = ();
my %nb_mots_affectes = ();
my %nb_mots_total_sortant = ();
my %compteur_nb_transferts_inter = ();
my %nb_mots_moyen_par_groupe = ();
my %nb_mots_total_par_groupe2 = ();
my $nb_mots_moyen_total = 0;
my $nbMots=0;
my $nbCognats = 0;
my %list_traduction = ();
my @tableau_transferts_dest = ();
my %compteur_transferts = ();
my $nb_hwt_intra = 0;
my %nombre_total_transfert_decortique = ();
my %nombre_total_vieux_transferts = ();
my $nombre_total_transferts = 0;

#== pourcentage de mots affectes par langue (alix, novembre 2014)
my %nb_mots_affectes_par_langue = ();
my %nb_mots_total_par_langue = ();

for(my $i=1;$i<=$NBGROUPES;$i++){
	$nb_mots_moyen_par_groupe{$i} = 0;	
}
my $petit_cpt = 0;
my $list_trans = "";

#======================================================
#= CALCUL DU NOMBRE DE MOTS TOTAL PAR GROUPE
#======================================================
print STDERR "\nFiltrage des mots selectionnes : ";
my $i=0;
my $val=-1;
my $cpt_mot=0;

open (NEWICK,$dbnewick) or die ("impossible d'ouvrir $dbnewick ($!)");
#open (OUT,">tmp.txt") || die "impossible d'ouvrir tmp.txt";
my $print_in_file=0;
while(my $ligne = <NEWICK>){
	chomp($ligne);
	if( $ligne =~ /=>/ ){	
		(my $mot) =  ($ligne =~/=> ([_A-Z0-9]*)/); 
		#	print STDOUT "\n$mot : ";
		if ( estDansTableau($mot,@motsCategorie) == 1){
			$print_in_file=1;
			#print OUT "$ligne\n";
			#		print STDOUT "trouve";
			$cpt_mot ++;
		}
		else{
				#		print STDOUT "non trouve";
			exit;
			$print_in_file=0;
		}
	}
	else{
		if ($print_in_file==1){
				#		print OUT "$ligne\n" ;
			#print STDERR "$ligne\n" ;
			(my @tab_feuilles) = ($ligne =~ /([0-9][0-9])-/g);
			#print STDERR join("-",@tab_feuilles) . "\n";
			
			for(my $i=0;$i<scalar @tab_feuilles;$i++){
				my $elt = $tab_feuilles[$i];
				$nb_mots_total_par_langue{$elt} ++;
				if(exists($nb_mots_total_par_groupe2{$tabGroup{$elt}})){
				$nb_mots_total_par_groupe2{$tabGroup{$elt}} ++;
				}
				else{
					$nb_mots_total_par_groupe2{$tabGroup{$elt}} = 1;
				}
			}
		}
	}
}
close(NEWICK);
#close(OUT);
print STDERR "\nNombre de mots existants : $cpt_mot";



#======================================================
#== VARIABLES
#======================================================
my %nombre_elements_affectes 		= ();
my %nb_mots_total_par_groupe 		= ();
my %pourcentage_elements_affectes	= ();
my $nombre_total_cognats			= 0;
my $alix_tmp						= 0;
my %wbe_multigroupes				= ();
my %stat_cognats					= {"HWT"=>0 , "noHWT" => 0};
my $total_cognats					= 0;
#=====================
#== INITIALISATION
#=====================
for(my $i=1;$i<12;$i++){
	for(my $j=1;$j<12;$j++){
		$pourcentage_elements_affectes{$i}{$j} = 0;
	}
}


#===========================================================================================================================
#= TRAITEMENT D'UN NOUVEAU MOT
#===========================================================================================================================
my $ligne = <RES>;
my $current_status = "";
my $nbHWT;

open ARBRES , $ARGV[2] || die "$!";
my @tab_tmp_arbres = <ARBRES>;
close ARBRES;
	
#======================================================
#== LECTURE DES ARBRES
#======================================================
my $arbre = "";
my $nombre_mot_total = 0;

open (MOTS , $ARGV[1]) or die ($!);
while (my $mot = <MOTS>){
	chomp($mot);	
	for(my $i=0;$i<scalar @tab_tmp_arbres;$i++){
		$arbre = $tab_tmp_arbres[$i];
		chomp($arbre);
		if($arbre eq "=> $mot"){
			$i++;
			$arbre = $tab_tmp_arbres[$i];
			while ($arbre =~ /^\(/){
				chomp($arbre);
				(my @tab_feuilles) = ($arbre =~ /([0-9][0-9])-/g);
				for(my $j=0;$j<scalar @tab_feuilles;$j++){
					$nb_mots_total_par_groupe{$tabGroup{$tab_feuilles[$j]}} ++;
					$nombre_mot_total ++;
				}
				$nombre_total_cognats++;
				$i++;
				$arbre = $tab_tmp_arbres[$i];
			}
			$i = (scalar @tab_tmp_arbres) + 1;
		}	
	}
}

close MOTS;

	

open (MOTS , $ARGV[1]) or die ($!);
#======================================================
#== TANT QU'IL Y A DES MOTS A TRAITER
#======================================================
while (my $mot = <MOTS>){

	chomp($mot);
#	print STDOUT "\n$mot";
	
	my @tab_arbres      = ();	#= Tableau contenant les abres
	my $nbCognats       = 0;	#= Nombre de cognats pour le mot
	my %tab_scenario    = ();	#= Tableau des scenarios
	my $noHWT			= 0;	#= Nombre de hwt par scenario
	my %stat_scenario	= ();	#= Nombre de hwt par scenario
	
	#======================================================
	#== RECHERCHE DES ARBRES DU MOT EN COURS DE TRAITEMENT
	#======================================================
	open (ARBRES , $ARGV[2]) or die ($!);
	while (my $arbre = <ARBRES>){
		chomp($arbre);
		if($arbre eq "=> $mot"){
			$arbre = <ARBRES>;
			while ($arbre =~ /^\(/){
				chomp($arbre);
				(my @tab_feuilles) = ($arbre =~ /([0-9][0-9])-/g);
				#for(my $i=0;$i<scalar @tab_feuilles;$i++){
				#	$nb_mots_total_par_groupe{$tabGroup{$tab_feuilles[$i]}} ++;
				#}
				push(@tab_arbres,$arbre); 
				#$nombre_total_cognats++;
				$arbre = <ARBRES>;
			}
			last;
		}	
	}
	close ARBRES;
	#print STDOUT join("\n",@tab_arbres);
	
	#==============================================
	#== RECHERCHE DES SCENARIOS POUR CHAQUE ARBRE
	#==============================================	 
	open (HWT , $ARGV[0]) or die ($!);
	while(my $hwt = <HWT>){
		chomp($hwt);
		if($hwt =~ /=> $mot$/){
			$hwt = <HWT>;
			while($hwt =~ /^[0-9]/){
				chomp($hwt);
				#print STDOUT "\n" . $hwt;
				if($hwt =~ /cognat/){
					$alix_tmp++;
					#print STDOUT "\n$mot $alix_tmp $hwt : ";
					$nbCognats++;
					$stat_scenario{$nbCognats}{"nbHWT"} = 0;
					$noHWT=0;
				}
				else{
					$noHWT++;
					my @tmp = split("->",$hwt);
					my $status = $tmp[2];
					my $destination = $tmp[1];
					my $source = $tmp[0];
					$source = trim($source);
					$destination = trim($destination);
					chomp($status);
					$status = trim($status);
					$tab_scenario{$nbCognats}{$noHWT}{"source"} 		= $source;
					$tab_scenario{$nbCognats}{$noHWT}{"destination"} 	= $destination;
					$tab_scenario{$nbCognats}{$noHWT}{"status"} 		= $status;
          if($status eq "CASXX"){
          	$tab_scenario{$nbCognats}{$noHWT}{"status"} = "CAS250";
          }
					$stat_scenario{$nbCognats}{"nbHWT"}++; 
					
				}
				$hwt = <HWT>;
			}
			
		}		
	}
	
	close HWT;

	
	#=========================================
	#== TRAITEMENT DES TRANSFERTS PAR COGNAT
	#=========================================	 
	$total_cognats += $nbCognats;
	for(my $i=1;$i<=$nbCognats;$i++){
    #print STDOUT "\nCognat $i ";
		my $arbre = $tab_arbres[$i-1];
		#print STDOUT "\n$arbre";
		(my @tab_feuilles) = ($arbre =~ /([0-9][0-9])-/g);
		my %nb_mots_par_groupe = ();
		for(my $i=0;$i<scalar @tab_feuilles;$i++){
			$nb_mots_par_groupe{$tabGroup{$tab_feuilles[$i]}} ++;
		}
		#print STDOUT "\nCognat $i : " . $stat_scenario{$i}{"nbHWT"};
		
		if($stat_scenario{$i}{"nbHWT"} == 0){
			$stat_cognats{"noHWT"} ++;		
		}
		else{
			$stat_cognats{"HWT"} ++;
		}
		
    my %eltDejaComptabilise = ();

    for(my $j=1;$j<=$stat_scenario{$i}{"nbHWT"};$j++){
			my @tab_source 	= split(" ", $tab_scenario{$i}{$j}{"source"});
			my %groupe_source = ();
			
      foreach my $elt (@tab_source){
				chomp($elt);
        my $copy = $elt;

				$elt =~ s/-[0-9]//g ;
        
        if ( ! exists $eltDejaComptabilise{$copy} ){
          $eltDejaComptabilise{$copy} = 1;
          if ( exists $stat_scenario{$i}{"cptTraduction"}{$elt} ) { $stat_scenario{$i}{"cptTraduction"}{$elt} ++; }
          else { $stat_scenario{$i}{"cptTraduction"}{$elt} = 1; }
        }
			}
		
			my @tab_dest 	= split(" ", $tab_scenario{$i}{$j}{"destination"});
			my %groupe_destination = ();
			
			foreach my $elt (@tab_dest){				
        chomp($elt);
        my $copy = $elt;

			  $elt =~ s/-[0-9]//g ;

        if ( ! exists $eltDejaComptabilise{$copy} ){
           $eltDejaComptabilise{$copy} = 1;
          if ( exists $stat_scenario{$i}{"cptTraduction"}{$elt} ) { $stat_scenario{$i}{"cptTraduction"}{$elt} ++; }
          else { $stat_scenario{$i}{"cptTraduction"}{$elt} = 1; }
        }
			}
    }

    
		for(my $j=1;$j<=$stat_scenario{$i}{"nbHWT"};$j++){
      #print STDOUT "\n\t\t" . $tab_scenario{$i}{$j}{"source"} . " -> " . $tab_scenario{$i}{$j}{"destination"}; 
			my @tab_source 	= split(" ", $tab_scenario{$i}{$j}{"source"});
			my %groupe_source = ();
			#print STDOUT "\n";
			foreach my $elt (@tab_source){
				chomp($elt);

				$elt =~ s/-[0-9]//g ;
        
				$groupe_source{$tabGroup{$elt}} += (1/ $stat_scenario{$i}{"cptTraduction"}{$elt}) ;
        #	print STDOUT "\n\tSource -> $elt (" . 	$stat_scenario{$i}{"cptTraduction"}{$elt}  . ":" . $groupe_source{$tabGroup{$elt}}  .") ";
			}
		
			my @tab_dest 	= split(" ", $tab_scenario{$i}{$j}{"destination"});
			my %groupe_destination = ();
			
			#print STDOUT "->";
			foreach my $elt (@tab_dest){				
        chomp($elt);

			  $elt =~ s/-[0-9]//g ;

				$groupe_destination{$tabGroup{$elt}} += (1/$stat_scenario{$i}{"cptTraduction"}{$elt}) ;
        #	print STDOUT "\n\tDestination -> $elt (" . 	$stat_scenario{$i}{"cptTraduction"}{$elt}  . ":" . $groupe_destination{$tabGroup{$elt}}  .") ";
			}

=for comment

      for(my $j=1;$j<=$stat_scenario{$i}{"nbHWT"};$j++){
			  print STDOUT "\n>" . $tab_scenario{$i}{$j}{"source"} . "< -> >" . $tab_scenario{$i}{$j}{"destination"} . "< -> >" . $tab_scenario{$i}{$j}{"status"} . "<"; 
			  	my @tab_source 	= split(" ", $tab_scenario{$i}{$j}{"source"});
          foreach my $elt (@tab_source){
          $elt =~ s/-[0-9]//g ;
          print STDOUT "\n" . $elt . " : " .  $stat_scenario{$i}{"cptTraduction"}{$elt};
        }
      }
=cut


=for comment
			if( ((scalar keys %groupe_source) >= 1) || ((scalar keys %groupe_destination) >= 1) ){
				my $fact = 1;	
  			if($tab_scenario{$i}{$j}{"status"} eq "CAS50"){
					$fact = .5;
					$wbe_multigroupes{ $tab_scenario{$i}{$j}{"source"} . "->" . $tab_scenario{$i}{$j}{"destination"}}{"compteur"}++;
					$wbe_multigroupes{ $tab_scenario{$i}{$j}{"source"} . "->" . $tab_scenario{$i}{$j}{"destination"}}{"status"} += .5;
				}
				elsif($tab_scenario{$i}{$j}{"status"} eq "CAS250"){
					$fact = .25;
					$wbe_multigroupes{ $tab_scenario{$i}{$j}{"source"} 		. "->" . $tab_scenario{$i}{$j}{"destination"}}{"compteur"}++;
					$wbe_multigroupes{ $tab_scenario{$i}{$j}{"destination"} . "->" . $tab_scenario{$i}{$j}{"source"}}{"compteur"}++;
					$wbe_multigroupes{ $tab_scenario{$i}{$j}{"destination"} . "->" . $tab_scenario{$i}{$j}{"source"}}{"status"} += .25;
					$wbe_multigroupes{ $tab_scenario{$i}{$j}{"source"} 		. "->" . $tab_scenario{$i}{$j}{"destination"}}{"status"} += .25;
				}
				elsif($tab_scenario{$i}{$j}{"status"} eq "CAS2"){
					$fact = .5;
					$wbe_multigroupes{ $tab_scenario{$i}{$j}{"source"} 		. "->" . $tab_scenario{$i}{$j}{"destination"}}{"compteur"}++;
					$wbe_multigroupes{ $tab_scenario{$i}{$j}{"destination"} . "->" . $tab_scenario{$i}{$j}{"source"}}{"compteur"}++;
					$wbe_multigroupes{ $tab_scenario{$i}{$j}{"destination"} . "->" . $tab_scenario{$i}{$j}{"source"}}{"status"} += .5;
					$wbe_multigroupes{ $tab_scenario{$i}{$j}{"source"} 		. "->" . $tab_scenario{$i}{$j}{"destination"}}{"status"} += .5;
				}
        elsif($tab_scenario{$i}{$j}{"status"} ne ""){
          #print STDOUT "\n" . $tab_scenario{$i}{$j}{"status"} ;
          $fact = $tab_scenario{$i}{$j}{"status"} ;
          $wbe_multigroupes{ $tab_scenario{$i}{$j}{"source"} . "->" . $tab_scenario{$i}{$j}{"destination"}}{"compteur"}++;
				  $wbe_multigroupes{ $tab_scenario{$i}{$j}{"source"} . "->" . $tab_scenario{$i}{$j}{"destination"}}{"status"} += $fact;
        }
				else{
					$wbe_multigroupes{ $tab_scenario{$i}{$j}{"source"} . "->" . $tab_scenario{$i}{$j}{"destination"}}{"compteur"}++;
					$wbe_multigroupes{ $tab_scenario{$i}{$j}{"source"} . "->" . $tab_scenario{$i}{$j}{"destination"}}{"status"} += 1;
				}

			}
=cut		
			#=============================================================================
			#== Calcul des pourcentages d'éléments affectés
			#=============================================================================
			foreach my $elt_source (keys %groupe_source){
				foreach my $elt_dest (keys %groupe_destination){
			
					my $fact = 1;
					
					if(($tab_scenario{$i}{$j}{"status"} eq "CAS50") || ($tab_scenario{$i}{$j}{"status"} eq "CAS2")){
						$fact = 0.5;
					}
					elsif($tab_scenario{$i}{$j}{"status"} eq "CAS250"){
						$fact = 0.25;	
					}
          elsif($tab_scenario{$i}{$j}{"status"} ne ""){
						$fact = $tab_scenario{$i}{$j}{"status"};	
					}
          
					
					$pourcentage_elements_affectes{$elt_source}{$elt_dest} 		+= $fact * ($groupe_source{$elt_source}/ (scalar @tab_source)) * ($groupe_destination{$elt_dest}/$nb_mots_total_par_groupe{$elt_dest}) * 100;
					$nombre_total_transfert_decortique{$elt_source}{$elt_dest} 	+= $fact * ($groupe_source{$elt_source}/ (scalar @tab_source)) * ($groupe_destination{$elt_dest}/(scalar @tab_dest));
					
					if($tab_scenario{$i}{$j}{"status"} eq "CAS2"){
						$pourcentage_elements_affectes{$elt_dest}{$elt_source} 	   += $fact * ($groupe_destination{$elt_dest}/ (scalar @tab_dest)) * ($groupe_source{$elt_source}/$nb_mots_total_par_groupe{$elt_source}) * 100;
						$nombre_total_transfert_decortique{$elt_dest}{$elt_source} += $fact * ($groupe_destination{$elt_dest}/ (scalar @tab_dest)) * ($groupe_source{$elt_source}/(scalar @tab_source));
					}
					elsif($tab_scenario{$i}{$j}{"status"} eq "CAS250"){
						$pourcentage_elements_affectes{$elt_dest}{$elt_source} 	   += $fact * ($groupe_destination{$elt_dest}/ (scalar @tab_dest)) * ($groupe_source{$elt_source}/$nb_mots_total_par_groupe{$elt_source}) * 100;
						$nombre_total_transfert_decortique{$elt_dest}{$elt_source} += $fact * ($groupe_destination{$elt_dest}/ (scalar @tab_dest)) * ($groupe_source{$elt_source}/(scalar @tab_source));
					}	
				}	
			}
			
			#=============================================================================
			#== Calcul des pourcentages d'éléments affectés par langue
			#=============================================================================
			my @tab_source  = split(" ", $tab_scenario{$i}{$j}{"source"});
			my @tab_dest    = split(" ", $tab_scenario{$i}{$j}{"destination"});
			foreach my $elt_source (@tab_source){
				$elt_source =~ s/-[0-9]//g ;
				foreach my $elt_dest (@tab_dest){
					$elt_dest =~ s/-[0-9]//g ;
			
					my $fact = 1;
					
					if(($tab_scenario{$i}{$j}{"status"} eq "CAS50") || ($tab_scenario{$i}{$j}{"status"} eq "CAS2")){
						$fact = 0.5;
					}
					elsif($tab_scenario{$i}{$j}{"status"} eq "CAS250"){
						$fact = 0.25;	
					}
          elsif($tab_scenario{$i}{$j}{"status"} ne ""){
						$fact = $tab_scenario{$i}{$j}{"status"};	
					}


					$nb_mots_affectes_par_langue{$elt_source}{$elt_dest}  += ($fact / scalar(@tab_source));
					
					if($tab_scenario{$i}{$j}{"status"} eq "CAS2"){
						$nb_mots_affectes_par_langue{$elt_dest}{$elt_source}  += ($fact / scalar(@tab_dest));
					}
					elsif($tab_scenario{$i}{$j}{"status"} eq "CAS250"){
						$nb_mots_affectes_par_langue{$elt_dest}{$elt_source}  += ($fact / scalar(@tab_dest));
					}	
				}	
			}
	
			#=============================================================================
			#= Recherche des vieux transferts
			#=============================================================================
			if( (scalar @tab_source > 1) || (scalar keys @tab_dest > 1) ){
				my $liste_source = join("<>",keys %groupe_source);
				my $liste_dest   = join("<>",keys %groupe_destination);
				if($tab_scenario{$i}{$j}{"status"} eq "CAS50"){
					if($liste_source eq $liste_dest){
						$nombre_total_vieux_transferts{"intra"} += 0.5;
					}
					else{
						$nombre_total_vieux_transferts{"extra"} += 0.5;
					}
				}
				elsif($tab_scenario{$i}{$j}{"status"} eq "CAS2"){
					if($liste_source eq $liste_dest){
						$nombre_total_vieux_transferts{"intra"} += 2*0.5;
					}
					else{
						$nombre_total_vieux_transferts{"extra"} += 2*0.5;
					}
				}
				elsif($tab_scenario{$i}{$j}{"status"} eq "CAS250"){
					if($liste_source eq $liste_dest){
						$nombre_total_vieux_transferts{"intra"} += 2*0.25;
					}
					else{
						$nombre_total_vieux_transferts{"extra"} += 2*0.25;
					}
				}
        elsif($tab_scenario{$i}{$j}{"status"} ne ""){
					if($liste_source eq $liste_dest){
						$nombre_total_vieux_transferts{"intra"} += $tab_scenario{$i}{$j}{"status"};
					}
					else{
						$nombre_total_vieux_transferts{"extra"} += $tab_scenario{$i}{$j}{"status"};
					}
				}

				else{
					if($liste_source eq $liste_dest){
						$nombre_total_vieux_transferts{"intra"} += 1;
					}
					else{
						$nombre_total_vieux_transferts{"extra"} += 1;
					}
				}
			}
			
			#=============================================================================
			#= Calcul du nombre total de transferts
			#=============================================================================
			if($tab_scenario{$i}{$j}{"status"} eq "CAS50"){
				$nombre_total_transferts += 0.5;
			}
			elsif($tab_scenario{$i}{$j}{"status"} eq "CAS2"){
				$nombre_total_transferts += 2*0.5;
			}
			elsif($tab_scenario{$i}{$j}{"status"} eq "CAS250"){
				$nombre_total_transferts += 2*0.25;
			}
      elsif($tab_scenario{$i}{$j}{"status"} ne ""){
        #print STDOUT "\nNombre total de transferts=" . $tab_scenario{$i}{$j}{"status"};
				$nombre_total_transferts += $tab_scenario{$i}{$j}{"status"};
			}

			else{
				$nombre_total_transferts += 1;
			}
			
		}
	}
}

close MOTS;

foreach my $elt1 (keys %wbe_multigroupes){
	foreach my $elt2 (keys %wbe_multigroupes){
		#print STDOUT "\n" . $wbe_multigroupes{$elt} . " : $elt";
	}
}


#============================
#== AFFICHAGE DES RESULTATS
#============================
print STDOUT "\n\nNombre total de mots : $nombre_mot_total";

print STDOUT "\n\nNombre de mots par groupe\n";
for (my $i=1;$i<=$NBGROUPES;$i++){
	printf STDOUT "%4d\t", $i;
}
print STDOUT "\n";
for(my $i=1;$i<=$NBGROUPES;$i++){
	printf STDOUT "%4d\t" , $nb_mots_total_par_groupe{$i};
}


print STDOUT "\n\nNombre moyen de mots par groupe par cognat: \n";
for (my $i=1;$i<=$NBGROUPES;$i++){
	printf STDOUT "%4d\t", $i;
}
print STDOUT "\n";
for(my $i=1;$i<=$NBGROUPES;$i++){
	$nb_mots_moyen_par_groupe{$i} = $nb_mots_total_par_groupe{$i} /  $nombre_total_cognats; #$nbMots;
	$nb_mots_moyen_total += $nb_mots_moyen_par_groupe{$i};
	printf(STDOUT "%4.2lf\t",$nb_mots_moyen_par_groupe{$i});
}

printf STDOUT "\n\nNombre de mots moyens total = %4.2lf\n" , $nb_mots_moyen_total;

my @colonnes  = ();
my @lignes    = ();
my @diagonale = ();

my @colonnes_dec  = ();
my @lignes_dec    = ();
my @diagonale_dec = ();

my $nombre_transferts_intra = 0;
my $nombre_transferts_extra = 0;

my $magic_number = 0;
my $magic_number2 = 0;
#my $nombre_mot_total = 0;
print STDOUT "\nNombre total de cognats : $nombre_total_cognats ($alix_tmp)";
print STDOUT "\nNombre total de cognats avec au moins un HWT : " . $stat_cognats{"HWT"};
print STDOUT "\nNombre total de cognats sans HWT : " . $stat_cognats{"noHWT"};
&newline;

print STDOUT "\n\nPourcentage des elements affectes\n";
for (my $i=1;$i<=$NBGROUPES;$i++){
	printf STDOUT "\t%4d", $i;
}

#=========================================================================================
#== Affichage de la matrice de pourcentage des elements affectés
#== Creation du fichier de données pour le heatmap
#=========================================================================================
my $somme_test=0;

open OUT, ">heatmap_biolinguistique.dat" || die $!;
for (my $i=1;$i<=$NBGROUPES;$i++){
	print STDOUT "\n$i";
	for (my $j=1;$j<=$NBGROUPES;$j++){
		if(exists($pourcentage_elements_affectes{$i}{$j})){
			printf STDOUT "\t%7.2lf", $pourcentage_elements_affectes{$i}{$j};	
			printf OUT "\n%d\t\"%s\"\t%d\t\"%s\"\t%lf" , $i-1, $nomGroup{$i} , $j-1, $nomGroup{$j}, $pourcentage_elements_affectes{$i}{$j};
		}
		else{
			printf STDOUT "\t%7.2lf", 0;
			printf OUT "\n%d\t\"%s\"\t%d\t\"%s\"\t0" , $i-1, $nomGroup{$i} , $j-1, $nomGroup{$j};
		} 
		$colonnes[$j]  += $pourcentage_elements_affectes{$i}{$j} if(exists($pourcentage_elements_affectes{$i}{$j}) && ($i != $j));
		$lignes[$i]    += $pourcentage_elements_affectes{$i}{$j} * $nb_mots_total_par_groupe{$j} if(exists($pourcentage_elements_affectes{$i}{$j}) && ($i != $j));
		$diagonale[$i] += $pourcentage_elements_affectes{$i}{$j} if(exists($pourcentage_elements_affectes{$i}{$j}) && ($i == $j));

	}
}
close OUT;

print STDOUT "\n\nNombre total de transferts décortiqués par groupe";
for (my $i=1;$i<=$NBGROUPES;$i++){
	print STDOUT "\n$i";
	for (my $j=1;$j<=$NBGROUPES;$j++){
		if(exists($nombre_total_transfert_decortique{$i}{$j})){
			printf STDOUT "\t%7.2lf", $nombre_total_transfert_decortique{$i}{$j};
		}
		else{
			printf STDOUT "\t%7.2lf", 0;
		} 
		$colonnes_dec[$j]  += $nombre_total_transfert_decortique{$i}{$j} if($i != $j);
		$lignes_dec[$i]    += $nombre_total_transfert_decortique{$i}{$j} if($i != $j);
		$diagonale_dec[$i] += $nombre_total_transfert_decortique{$i}{$j} if($i == $j);
		$nombre_transferts_intra += $nombre_total_transfert_decortique{$i}{$j} if($i != $j);
		$nombre_transferts_extra += $nombre_total_transfert_decortique{$i}{$j} if($i == $j);
	}
}

&newline;

printf STDOUT "\nNombre total de vieux transferts INTRA-GROUPE : %d/%d (%3.2lf%%)" , $nombre_total_vieux_transferts{"intra"} , $nombre_transferts_intra , ($nombre_total_vieux_transferts{"intra"}*100/$nombre_transferts_intra);
printf STDOUT "\nNombre total de vieux transferts EXTRA-GROUPE : %d/%d (%3.2lf%%)" , $nombre_total_vieux_transferts{"extra"} , $nombre_transferts_extra , ($nombre_total_vieux_transferts{"extra"}*100/$nombre_transferts_extra) ;
printf STDOUT "\nNombre total de vieux transferts : %d (%3.2lf%%)" , ($nombre_total_vieux_transferts{"extra"} + $nombre_total_vieux_transferts{"intra"}) , ( ($nombre_total_vieux_transferts{"extra"} + $nombre_total_vieux_transferts{"intra"})*100/$nombre_total_transferts) ;
print STDOUT "\nNombre total de transferts avec pondération   : " . $nombre_total_transferts;

print STDOUT "\n\nWilcoxon signed-rank test";
for (my $i=1;$i<=$NBGROUPES;$i++){
	printf STDOUT "\n$i\t%d\t%d\t%3.2lf\t%3.2lf\t%3.2lf" , $age_groupe[$i-1] , $nombre_langue_par_groupe[$i-1] , $diagonale_dec[$i] , $lignes_dec[$i] , $colonnes_dec[$i];
}


print STDOUT "\n\nPourcentage des elements affectés ENTRANTS\n";
for (my $i=1;$i<=$NBGROUPES;$i++){
	printf STDOUT "%4d\t", $i;
}
print STDOUT "\n";
for (my $i=1;$i<=$NBGROUPES;$i++){
	$magic_number += ($colonnes[$i] + $diagonale[$i]) * ($nb_mots_moyen_par_groupe{$i} / $nb_mots_moyen_total);
	$magic_number2 +=  $diagonale[$i] * ($nb_mots_moyen_par_groupe{$i} / $nb_mots_moyen_total);
	#$nombre_mot_total += $nb_mots_total_par_groupe{$i};
	printf STDOUT "%4.1lf\t", $colonnes[$i];
}
#$magic_number2 = $magic_number2/$nombre_mot_total;

print STDOUT "\n\nPourcentage des elements affectés SORTANTS\n";
for (my $i=1;$i<=$NBGROUPES;$i++){
	printf STDOUT "%4d\t", $i;
}
print STDOUT "\n";
for (my $i=1;$i<=12;$i++){
	printf STDOUT "%4.5lf\t", $lignes[$i] /($nombre_mot_total-$nb_mots_total_par_groupe{$i}),$lignes[$i],$nombre_mot_total,$nb_mots_total_par_groupe{$i};#/($NBGROUPES-1);
	#printf STDOUT "%4.5lf(%lf,%d,%d)\t", $lignes[$i] /($nombre_mot_total-$nb_mots_total_par_groupe{$i}),$lignes[$i],$nombre_mot_total,$nb_mots_total_par_groupe{$i};#/($NBGROUPES-1);
}

print STDOUT "\n\nPourcentage des elements affectés SORTANTS relatif\n";
for (my $i=1;$i<=$NBGROUPES;$i++){
	printf STDOUT "%4d\t", $i;
}
print STDOUT "\n";
for (my $i=1;$i<=12;$i++){
	printf STDOUT "%4.1lf\t", ($lignes[$i] /($nombre_mot_total-$nb_mots_total_par_groupe{$i})) / $nb_mots_moyen_par_groupe{$i},  $lignes[$i],$nombre_mot_total,$nb_mots_total_par_groupe{$i},$nb_mots_moyen_par_groupe{$i};#/($NBGROUPES-1);
	#printf STDOUT "%4.1lf(%lf,%lf,%lf,%lf)\t", ($lignes[$i] /($nombre_mot_total-$nb_mots_total_par_groupe{$i})) / $nb_mots_moyen_par_groupe{$i},  $lignes[$i],$nombre_mot_total,$nb_mots_total_par_groupe{$i},$nb_mots_moyen_par_groupe{$i};#/($NBGROUPES-1);
}

print STDOUT "\n\nPourcentage des elements affectés INTRAS\n";
for (my $i=1;$i<=12;$i++){
	printf STDOUT "%4d\t", $i;
}
print STDOUT "\n";
for (my $i=1;$i<=12;$i++){
	printf STDOUT "%4.1lf\t", $diagonale[$i];
}

print STDOUT "\n\nNombre magique (somme des colonnes pondérées) : ";
printf STDOUT "%3.2lf%% de chance qu'un mot choisi au hasard soit affecté par un HWT" , $magic_number;
print STDOUT "\n\nNombre magique 2 (somme des diagonales pondérées) : ";
printf STDOUT "%3.2lf%% de chance qu'un mot choisi au hasard soit affecté par un HWT prevenant du même groupe" , $magic_number2;

print STDOUT "\n\nPour insertion facile dans le fichier excel\n";
for (my $i=1;$i<=12;$i++){
	printf STDOUT "%4.1lf\t", $diagonale[$i];
}
print STDOUT "\n";
for (my $i=1;$i<=$NBGROUPES;$i++){
	printf STDOUT "%4.1lf\t", $colonnes[$i];
}
print STDOUT "\n";
for (my $i=1;$i<=12;$i++){
	printf STDOUT "%4.1lf\t", $lignes[$i] /($nombre_mot_total-$nb_mots_total_par_groupe{$i});#/($NBGROUPES-1);
}
print STDOUT "\n";
for (my $i=1;$i<=12;$i++){
	printf STDOUT "%4.1lf\t", ($lignes[$i] /($nombre_mot_total-$nb_mots_total_par_groupe{$i})) / $nb_mots_moyen_par_groupe{$i};#/($NBGROUPES-1);
}
print STDOUT "\n";
for (my $i=1;$i<=$NBGROUPES;$i++){
	print STDOUT "\n";
	for (my $j=1;$j<=$NBGROUPES;$j++){
		if(exists($pourcentage_elements_affectes{$i}{$j})){
			printf STDOUT "%4.2lf\t", $pourcentage_elements_affectes{$i}{$j};
		}
		else{
			printf STDOUT "%4.2lf\t", 0;
		} 
	}
}


print STDOUT "\n";
foreach my $elt (keys %wbe_multigroupes){
	#printf STDOUT "\n%2d : %5.2lf : %s" , $wbe_multigroupes{$elt}{"compteur"} , $wbe_multigroupes{$elt}{"status"} , $elt ; # if ($wbe_multigroupes{$elt}{"compteur"} > 0);
}

print STDOUT "\nNombre de mots par Langues";
foreach my $noLang (sort par_num keys %nb_mots_total_par_langue){
	print STDOUT "\n$noLang : $nb_mots_total_par_langue{$noLang}";
}

my %total_colonne = ();

open (OUT,">stat_langue.txt") or die($!);
foreach my $elt1 (sort par_num keys %tabGroup){
	print OUT "$elt1\t";
	foreach my $elt2 (sort par_num keys %tabGroup){
		printf (OUT "%3.1lf " , ($nb_mots_affectes_par_langue{$elt1}{$elt2}/$nb_mots_total_par_langue{$elt2})*100);
		$total_colonne{$elt2} += ($nb_mots_affectes_par_langue{$elt1}{$elt2}/$nb_mots_total_par_langue{$elt2})*100;
	}
	print OUT "\n";
}
close(OUT);

foreach my $elt1 (sort par_num keys %total_colonne){
	print STDOUT "\n$elt1:$total_colonne{$elt1}";
}

#$nb_mots_affectes_par_langue

print STDOUT "\n\nFin normale du programme\n\n";



#============================================================================
#============================= SOUS-ROUTINES ================================
#============================================================================

sub par_num { return $a <=> $b }

sub estDansTableau{

	(my $mot,my @tableau) = @_;
	
	for(my $i=0;$i< scalar @tableau;$i++){
		my $mot2 = $tableau[$i];
		chomp($mot2);
		return 1 if($mot eq $mot2);
	}
	return 0;
}
sub memeGroupe{
	
	my @liste_elts = @_;
	
	#print STDOUT "\nDans memeGroupe " . join(",",@liste_elts);	
	my $pas_bon = 0;
	my $premier = 0;
	my $indice;
	foreach my $elt1 (@liste_elts){
		my $elt = $elt1;
		$elt =~ s/-[0-9]//g ;
		my $tmp = $tabGroup{$elt};
		
		#print STDOUT "\n$elt->" . $tabGroup{$elt};	
		if($premier == 0){		
			$indice = $tmp;
			$premier = 1;
		}
		else{
			if ($tmp != $indice){
				$pas_bon = 1;
				$indice = 0;
			}
		}
	}
	return ($indice,$pas_bon);
}

sub listGroupe{
	
	my @liste_elts = @_;
	my @liste = ();
	
	foreach my $elt1 (@liste_elts){
		my $elt = $elt1;
		$elt =~ s/-[0-9]//g ;
		my $tmp = $tabGroup{$elt};
		if(!grep(/$tmp/,@liste)){
			push(@liste,$tmp); 
		}
	}
	return @liste;
}
sub nbEltGroupe{
	
	my @liste_elts = @_;
	my @liste = (0,0,0,0,0,0,0,0,0,0,0,0,0);
	
	foreach my $elt1 (@liste_elts){
		my $elt = $elt1;
		$elt =~ s/-[0-9]//g ;
		my $tmp = $tabGroup{$elt};
		$liste[$tmp] ++;
	}
	return @liste;
}

sub nbGroupes{
	
	my @liste_elts = @_;
	my $cpt = 0;
	my @tab = (0,0,0,0,0,0,0,0,0,0,0,0,0);
	
	#print STDOUT "\nDans nbGroupes " . join(",",@liste_elts);
	
	foreach my $elt1 (@liste_elts){
		my $elt = $elt1;
		$elt =~ s/-[0-9]//g ;
		$tab[$tabGroup{$elt}] = 1;
	#	print STDOUT " -> " . $tabGroup{$elt};
	}
	
	for(my $i=1;$i<=$NBGROUPES;$i++){
		$cpt += $tab[$i];
	}
	
	return $cpt;
}
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub newline{
	print STDOUT "\n";
}
