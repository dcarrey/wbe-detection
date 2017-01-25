//================================================================================================
//=  HGT-DETECTION v3.1b
//=  Authors : Alix Boc and Vladimir Makarenkov
//=  Date : June 2009
//=
//=  Description : This program detect horizontal gene transfer (HGT). As input it takes 2 
//=  matrices: a species matrix and a gene matrix. the goal is to transform the species matrix
//=  into the gene matrix following a transfer scenario. There are 3 criteria : the robinson and
//=  Foulds distance, the least-square criterion and the bipartition distance. We also use the
//=  subtree constraint. With this version we can now perform simulation.
//=
//=	 input   : file with species tree and gene tree as matrix in the phylip or newick format.
//=			   In case of simulation, the species tree and all the gene trees in the same file in
//=            the phylip format or newick string
//=  output  : a list of HGT and the criteria values for each one.
//=	 options :
//=
//================================================================================================


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <signal.h>

#pragma warning(disable:4996)
#define ARRAY_SIZE(array) (sizeof((array))/sizeof((array[0])))
#define MAX( a, b ) ( ( a > b) ? a : b ) 

#ifndef strsep
#include "strsep.cpp"
#endif

#include "structures.h"
#include "utils_tree.cpp"
#include "fonctions.cpp"
#include "bioling.cpp"

#define binaireSpecies 0 
#define binaireGene    1

void traiterSignal(int sig){
	printf("\nMESSAGE : SEGMENTATION FAULT #%d DETECTED",sig);
	printf("\nUse valgrind or gdb to fix the problem");
	printf("\n");
	exit(-1);
}

//========================================================================================================
//============================================ MAIN ======================================================
//========================================================================================================
int main(int nargc,char **argv){

	struct InputTree SpeciesTree;				    //== initial species tree
	struct InputTree SpeciesTreeCurrent;		//== initial species tree
	struct InputTree FirstTree;
	struct InputTree GeneTree,tmpGeneTree;					    //== initial gene tree
	struct InputTree SpeciesTreeRed;			  //== reduced species tree
	struct InputTree GeneTreeRed;				    //== reduced gene tree
	struct ReduceTrace aMap;					      //== mapping structure between species tree and recuded species tree
	struct InputTree geneTreeSave;
	struct HGT * bestHGTRed = NULL;				  //== list of HGT for the reduced tree
	struct HGT * bestHGT = NULL;				    //== list of HGT for the normale tree
	struct HGT * outHGT = NULL;
	struct HGT * bestHGTmulticheck = NULL;
	int nbHGT_boot;
	int first = 1,k,l;
	int cpt_hgt,i,j,tmp,nbTree=0;
	int bootstrap = 0;
	int multigene = 0;
	int nbHgtFound = 0;
	struct CRITERIA * multicheckTab=NULL;
	struct CRITERIA aCrit;						       //== struture of all the criteria
	struct DescTree *DTSpecies,					     //== structure of submatrices for the species tree
					*DTGene;					               //== structure of submatrices for the gene tree
	struct Parameters param;
	FILE *in,*out;
	int max_hgt,nbHGT;
	int ktSpecies;
	int trivial = 1;
	int *speciesLeaves = NULL;
	int RFref;
	int imc;	
	char *mot = (char*)malloc(100);
	struct Translation *trans;
  	
	initInputTree(&geneTreeSave);


	//== read parameters
	if(readParameters(&param,argv,nargc)==-1){
		help();
		//getchar();
		exit(-1);
	}

	//if(strcmp(param.version,"consol")==0)
		signal(SIGSEGV,traiterSignal);

	//if(strcmp(param.version,"consol")==0)
	//	printf(startMessage);

	//== open the output file
	if((out=fopen(param.outputfile,"w+"))==NULL){
		printf("Can't open output file (%s)",param.outputfile);
		exit(-1);
	}

	PrintHeader(out,param);
	
	//== open the input file
	if((in=fopen(param.inputfile,"r"))==NULL){
		printf("Can't open input file (%s)",param.inputfile);
		exit(-1);
	}

	//== open the bootstrapFile
	if(strcmp(param.bootstrap,"yes") == 0){
		bootstrap = 1;
	}
	if(strcmp(param.multigene,"yes") == 0){
		multigene = 1;
		if(strcmp(param.speciesroot,"file"))
			strcpy(param.speciesroot,"midpoint");
	}
	remove(param.hgtResultFile);

	initInputTree(&FirstTree);
	initInputTree(&SpeciesTreeCurrent);

	FILE * results; // = fopen(param.results,"w+");
	
	if((results = fopen(param.results,"w+"))==NULL){
		printf("PROBLEME AVEC RESULTS");
		exit(0);
	}

	
//==============================================================================
//============================= LECTURE DES ARBRES =============================
//==============================================================================
	tmp = readInputFile(in, param.input/*,&SpeciesTree,&GeneTree*/);

	if(tmp==-1) {
		printf("\nCannot read input data !!\n");
		exit(-1);
	}
  
	cpt_hgt = 0;

	initInputTree(&SpeciesTree);
	initInputTree(&GeneTree);
	initInputTree(&tmpGeneTree);
	initInputTree(&SpeciesTreeRed);
	initInputTree(&GeneTreeRed);

	//== lecture des matrices ou chaines newick en entree
	if(readInput(SPECIE,param.input,&SpeciesTree) == -1){ printf("\nError in species tree\n"); exit(-1);}
	if(readInput(GENE,param.input,&GeneTree) == -1){ printf("\nError in gene tree\n"); exit(-1);}

  sortMatrices(GeneTree.Input,GeneTree.SpeciesName,SpeciesTree.SpeciesName,SpeciesTree.size);
/*
	for(i=1;i<=SpeciesTree.size;i++){
			printf("\n%s\t",SpeciesTree.SpeciesName[i]);
			for(j=1;j<=SpeciesTree.size;j++){
				printf("%lf ",SpeciesTree.Input[i][j]);
			}
		}
	printf("\n");printf("\n");
	
	for(i=1;i<=GeneTree.size;i++){
			printf("\n%s\t",GeneTree.SpeciesName[i]);
			for(j=1;j<=GeneTree.size;j++){
				printf("%lf ",GeneTree.Input[i][j]);
			}
		}
	printf("\n");printf("\n");*/
	
	NJ(SpeciesTree.Input,SpeciesTree.ADD,SpeciesTree.size);
	NJ(GeneTree.Input,GeneTree.ADD,GeneTree.size);
/*	
	for(i=1;i<=SpeciesTree.size;i++){
			printf("\n%s\t",SpeciesTree.SpeciesName[i]);
			for(j=1;j<=SpeciesTree.size;j++){
				printf("%lf ",SpeciesTree.Input[i][j]);
			}
		}
	printf("\n");printf("\n");
	
	for(i=1;i<=GeneTree.size;i++){
			printf("\n%s\t",GeneTree.SpeciesName[i]);
			for(j=1;j<=GeneTree.size;j++){
				printf("%lf ",GeneTree.Input[i][j]);
			}
		}
	printf("\n");printf("\n");*/
	
	//if(SpeciesTree.size > 50 ) exit(0);
	//== ajouter le 19 janvier : initialisation d'une colonne suplémentaire à 0
	/*for(i=1;i<=GeneTree.size+1;i++){
			GeneTree.ADD[GeneTree.size+1][i] = GeneTree.ADD[i][GeneTree.size+1] = epsilon;
	}*/
	//== fin ajout
	
	
	//== construction des differentes représentation des arbres (adjacence,aretes,longueur,degre)
	CreateSubStructures(&SpeciesTree,1,binaireSpecies);
	CreateSubStructures(&GeneTree,1,binaireGene);
	//copyInputTree(&tmpGeneTree,GeneTree,1,1);
	//for(i=1;i<=2*SpeciesTree.size-3-SpeciesTree.kt;i++){
	//	printf("\n%d -> %ld--%ld : %lf",i,SpeciesTree.ARETE[2*i-2],SpeciesTree.ARETE[2*i-1],SpeciesTree.LONGUEUR[i-1]);
	//}
	
/*
  printf("\nAffichage du groupe par numero de langue :");
  for(int i=1;i<= (sizeof(groupLang)/sizeof(*groupLang)); i++){
    printf("\n%d => %d", i,groupLang[i]);
  }

  exit(0);
*/
//==============================================================================
//============================ GESTION DES RACINES =============================
//==============================================================================
  //== selection de la racine
  for(int i=1;i<=SpeciesTree.size;i++){
    printf("\n%d : %s", i,SpeciesTree.SpeciesName[i]);
    if(strcmp(SpeciesTree.SpeciesName[i],"root") == 0)
       SpeciesTree.Root=i;   
  }
 
  if(SpeciesTree.Root == -1) {
	  printf("\nWBE-DETECTION : add root to languages tree");
		addRoot(&SpeciesTree,NULL,SpeciesBranch,param.speciesroot,param.speciesRootfile,NULL);
	}
  for(int i=1;i<=GeneTree.size;i++){
    printf("\n%d : %s", i,GeneTree.SpeciesName[i]);
    if(strcmp(GeneTree.SpeciesName[i],"root") == 0)
       GeneTree.Root=i;   
  }
  if(GeneTree.Root == -1) {
	  printf("\nWBE-DETECTION : add root to words tree");
	  addRoot(&GeneTree,NULL,GeneBranch,param.generoot,param.geneRootfile,NULL); //bestbipartition
	}
  
  //SAVEASNewick(SpeciesTree.LONGUEUR, SpeciesTree.ARETE, SpeciesTree.SpeciesName, SpeciesTree.size, SpeciesTree.kt, param.filteredLanguageTree) ;

	nbTree++;
	
	if(first ==1){
		FirstTree.ADD=NULL;
		FirstTree.ARETE=NULL;
		copyInputTree(&FirstTree,SpeciesTree,1,1);
		AdjustBranchLength(&FirstTree,GeneTree,binaireSpecies,1);
	}
	SpeciesTreeCurrent.ADD=NULL;
	SpeciesTreeCurrent.ARETE=NULL;
	copyInputTree(&SpeciesTreeCurrent,SpeciesTree,1,0);

	AdjustBranchLength(&SpeciesTreeCurrent,GeneTree,binaireSpecies,1);
	AdjustBranchLength(&SpeciesTree,GeneTree,binaireSpecies,1);

	InitCriteria(&aCrit,SpeciesTree.size);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	//printf("\nWBE-DETECTION : rf=%d n=%d rf/(2n-6)=%lf\n",aCrit.RF,SpeciesTree.size,(double)aCrit.RF/(2*(double)SpeciesTree.size-6));
	
	if(SpeciesTree.size > GeneTree.size) max_hgt = 4*GeneTree.size * GeneTree.size;
	else max_hgt = 4*SpeciesTree.size * SpeciesTree.size;

	bestHGTRed = (struct HGT*)malloc(max_hgt*sizeof(struct HGT)); 
	bestHGT = (struct HGT*)malloc(max_hgt*sizeof(struct HGT));

	for(i=0;i<max_hgt;i++){
		bestHGTRed[i].listSource = NULL;
		bestHGTRed[i].listDestination = NULL;
		bestHGT[i].listSource = NULL;
		bestHGT[i].listDestination = NULL;
	} 

	if(first==1 && strcmp(param.scenario,"multiple")!=0 && strcmp(param.mode,"multicheck")==0)
		multicheckTab = (struct CRITERIA *)malloc(max_hgt*sizeof(struct CRITERIA)); 

	
	if((bootstrap != 1) || (bootstrap==1 && first==1)){

	//	if(strcmp(param.version,"consol")==0){
  //		fprintf(out,"\n\n#################################",nbTree);
	//		fprintf(out,"\n## Detection #%d                 #",nbTree);
	//		fprintf(out,"\n#################################",nbTree);
	//	}
	//	printBranches(out,SpeciesTree,SpeciesBranchNewick,NULL,0);
	//  fprintf(out,"\n\nafficherlesarbresici");
	//	fprintf(out,"\n\n=========================================");
	//	fprintf(out,"\n= Criteria values before the computation ");
	//	fprintf(out,"\n=========================================");
	//	fprintf(out,"\nRobinson and Foulds distance (RF) = %d",aCrit.RF);
	//	fprintf(out,"\nLeast-squares coefficient(LS)     = %lf",aCrit.LS);
	//	fprintf(out,"\nBipartition dissimilarity         = %lf",aCrit.BD);
	
    fprintf(results,"%d,%lf,%lf\n",aCrit.RF,aCrit.LS,aCrit.BD);
		RFref=aCrit.RF;
	}

	DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);

	DTGene = (struct DescTree*)malloc((2*GeneTree.size-2-GeneTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartition(GeneTree.ARETE,GeneTree.ADD,GeneTree.Root,GeneTree.Adjacence,DTGene,GeneTree.size,GeneTree.kt);

	
	//============= lecture des mots ================= 
	int nbGroupesBloques=0;
	int ** groupesBloques = (int **)malloc(2*SpeciesTree.size*(sizeof(int*)));
	for(i=0;i<2*SpeciesTree.size;i++){
		groupesBloques[i] = (int *)malloc(2*SpeciesTree.size*(sizeof(int*)));
		for(j=0;j<2*SpeciesTree.size;j++){
			groupesBloques[i][j] = 0;
		}
	}
	
	trans = (struct Translation*) malloc (SpeciesTree.size * sizeof(struct Translation));

  char ** trans_langNum   = (char**) malloc (SpeciesTree.size * sizeof(char*));
  char ** trans_langTrans = (char**) malloc (SpeciesTree.size * sizeof(char*));
  
 // for(j=1;j<=SpeciesTree.size-1;j++){
 //   trans_langNum[j]   = (char*) malloc (10 * sizeof(char));
 //   trans_langTrans[j] = (char*) malloc (10 * sizeof(char));
 // }

// exit(0);
//  printf("\nNombre de feuilles : %d",SpeciesTree.size);
  FILE * transFp ;
	char * langNum;
	char * langTrans;
	
  if((transFp = fopen(param.translationsfile,"r"))==NULL){
		printf("\nError : Cannot open %s", param.translationsfile);
		exit(0);
	}
  
 // for(j=1;j<=SpeciesTree.size-1;j++){
  //  printf("\n%d",j);
  //  trans[j].langNum   = (char *)malloc(10  * sizeof(char));
	//	trans[j].langTrans = (char *)malloc(100 * sizeof(char));
  //}
	langNum = (char *)malloc(10);
	langTrans = (char *)malloc(100);
	
  for(i=1;i<=SpeciesTree.size-1;i++){
		fscanf(transFp,"%s %s",langNum,langTrans);
		printf("\n%s - %s",langNum,langTrans);
		for(j=1;j<=SpeciesTree.size-1;j++){
			if(strcmp(langNum,SpeciesTree.SpeciesName[j]) == 0){
        trans[j].langNum = (char *)malloc(10);
        strcpy(trans[j].langNum,langNum);
				trans[j].langTrans = (char *)malloc(100);
        strcpy(trans[j].langTrans,langTrans);
			}
		}
	}
  free(langNum);
  free(langTrans);
	

/*	for(i=1;i<=SpeciesTree.size-1;i++){
	  printf("\n%d : %s -> %s -> %s",i,trans[i].langNum,trans[i].langTrans,SpeciesTree.SpeciesName[i]);
	}
*/
  
  //printf("\navgdiffblock = %lf", param.avgdiffblock);
	double distMots1,distMots2;
	double totalDistMots;
	int nbLang;
	char *chaine_res1 = (char*)malloc(100); 
	char *chaine_res2 = (char*)malloc(100); 
	char *tabDoublon1 = (char*)malloc(100); 
	char *tabDoublon2 = (char*)malloc(100); 
	
	int nbD1,nbD2;
	
	
	for(i=SpeciesTree.size+1;i<=2*SpeciesTree.size-2-SpeciesTree.kt;i++){
		//printf("\n=== %d : ",i);
		distMots1 = 0;
		distMots2 = INFINI;
		
		nbLang = 0;
		totalDistMots = 0;
		
		for(j=1;j<=DTSpecies[i].nbSommet;j++){
			//printf("\n(%d) ",DTSpecies[i].Tableau[j]);		
			for(k=j+1;k<=DTSpecies[i].nbSommet;k++){
        //printf("->(%d) ",DTSpecies[i].Tableau[k]);	
				nbLang++;
        //printf(" %s - %s", trans[DTSpecies[i].Tableau[j]].langTrans,trans[DTSpecies[i].Tableau[k]].langTrans);
				distMots1 = levenshtein_distance(trans[DTSpecies[i].Tableau[j]].langTrans,trans[DTSpecies[i].Tableau[k]].langTrans);
				nbD1 = delDoublons(trans[DTSpecies[i].Tableau[j]].langTrans,chaine_res1,tabDoublon1);
				nbD2 = delDoublons(trans[DTSpecies[i].Tableau[k]].langTrans,chaine_res2,tabDoublon2);
				//Changement de distance
				if(distMots1 > 0){
					if((nbD1 > 0) || (nbD2> 0) ){
						if(tabDoublon1[0] == tabDoublon2[0]){
							distMots2 = (levenshtein_distance(chaine_res1,chaine_res2))/MAX(strlen(trans[DTSpecies[i].Tableau[j]].langTrans),strlen(trans[DTSpecies[i].Tableau[k]].langTrans));
						}
						else{
							distMots2 = (levenshtein_distance(chaine_res1,chaine_res2)+0.10)/MAX(strlen(trans[DTSpecies[i].Tableau[j]].langTrans),strlen(trans[DTSpecies[i].Tableau[k]].langTrans));
						}
					}
					else{
						distMots2 = (levenshtein_distance(chaine_res1,chaine_res2))/MAX(strlen(trans[DTSpecies[i].Tableau[j]].langTrans),strlen(trans[DTSpecies[i].Tableau[k]].langTrans));
					}
				}
				
				totalDistMots += (distMots1<distMots2)?distMots1:distMots2;
				//printf(" [%s-%s] = %3.1lf\n",trans[DTSpecies[i].Tableau[j]].langTrans,trans[DTSpecies[i].Tableau[k]].langTrans,(distMots1<distMots2)?distMots1:distMots2);
			}
		}
    //printf("\n\nAverage wordsdist/lang: %lf/%d = %lf\n\n",totalDistMots,nbLang,totalDistMots/nbLang);
		
		if((totalDistMots/nbLang) <= param.avgdiffblock){ //1.5
			groupesBloques[nbGroupesBloques][0] = DTSpecies[i].nbSommet;
     // printf("\n%lf=>%d", (totalDistMots/nbLang),  DTSpecies[i].nbSommet);
			for(j=1;j<=DTSpecies[i].nbSommet;j++){
				groupesBloques[nbGroupesBloques][j] = DTSpecies[i].Tableau[j];
     //   printf(" %d", DTSpecies[i].Tableau[j]);
			}
			nbGroupesBloques++;
		}
	}

	//== verifions les sous-groupes ===
	int * tabGroupeA =(int *)malloc(nbGroupesBloques*sizeof(int));
	
	for(i=0;i<nbGroupesBloques;i++){
		for(j=1;j<nbGroupesBloques;j++){
			if(estUnSousGroupe(groupesBloques[i],groupesBloques[j])){
				tabGroupeA[i] = 1;
			}
		}
	}
	
  printf("\n\nNombre de groupes a ne pas briser : %d\n",nbGroupesBloques);
	for(i=0;i<nbGroupesBloques;i++){
		printf("\n(%d) - ",groupesBloques[i][0]);
		for(j=1;j<=groupesBloques[i][0];j++){
			printf("%d ",groupesBloques[i][j]);
		}
		printf(" - (%d)",tabGroupeA[i]);
	}
	
	if((param.constraints == 0) || (param.constraints == 2)){
		//nbGroupesBloques=0;
		printf("\nNE PAS FAIRE LA CONTRAINTE 1 (%d)",param.constraints);
	}
	
			
		printf("\nPre-traitement .......");
		//== date : 26 janvier 2008
		//== l'appel a cette fonction permet de recreer le meme sous arbre dans l'arbre de gene et l'arbre d'especes
		//== si les distances entre les especes sont nulles.
		int tousLesCasSontTraitees = FALSE;
		int *tab_tous_les_sommets = (int*)malloc((SpeciesTree.size+1) * sizeof(int));
		int *tab_sommets_selectionnes = (int*)malloc((SpeciesTree.size+1) * sizeof(int));
		int *tab_branches = (int*)malloc( 4*(GeneTree.size+1) * sizeof(int));
		int nb_branches;
		int temoin_nouveau_cas;
		struct HGT aHGT;
		
		
		
		for(i=1;i<=SpeciesTree.size;i++) tab_tous_les_sommets[i] = 0;
		tab_tous_les_sommets[0] = FALSE;

	/*	for(i=1;i<=GeneTree.size;i++){
			printf("\n");
			for(j=1;j<=GeneTree.size;j++){
				printf("%lf ",GeneTree.ADD[i][j]);
			}
		}
		printf("\n");printf("\n");*/
		while(tousLesCasSontTraitees == FALSE){
			
			ListeSommets_taille_0(GeneTree.Input,tab_tous_les_sommets,GeneTree.size); //= 25 janvier 2010 : GeneTree.size-1);
			tousLesCasSontTraitees = tab_tous_les_sommets[0];
			
			
			if(tousLesCasSontTraitees == FALSE){
				tab_sommets_selectionnes[0] = 0;
				temoin_nouveau_cas=0;	
				//printf("\nGroupe d'especes : ");
				for(i=1;i<GeneTree.size;i++){
					if(tab_tous_les_sommets[i] == 1){
						if(temoin_nouveau_cas == 0){
							temoin_nouveau_cas=1;
							//printf("\nTraitement : ");
						}
						tab_tous_les_sommets[i] = 2;
						tab_sommets_selectionnes[0] = tab_sommets_selectionnes[0] + 1;
						tab_sommets_selectionnes[tab_sommets_selectionnes[0]] = i;
						//printf("%d(%s) ",i,GeneTree.SpeciesName[i]);
					}
				}
				
				//GeneTree.Root = GeneTree.size; //findRoot(tab_sommets_selectionnes,GeneTree.size);
				//printf("\nRoot = %d",GeneTree.Root);
				
				
				if(tab_sommets_selectionnes[0] > 0){
					
					ListesBranchesPourHGT(tab_sommets_selectionnes,GeneTree.ARETE,GeneTree.size,DTGene,tab_branches,&nb_branches);
					
					/*printf("\n->Branches :");
					for(i=1;i<=nb_branches;i++){
						printf("%d--%d ",GeneTree.ARETE[2*tab_branches[i]-1],GeneTree.ARETE[2*tab_branches[i]-2]);
					}*/
					
					//exit(0);
					
					while(findBestHGT_nombreLimite(DTGene,DTSpecies,tab_branches,nb_branches,GeneTree,SpeciesTree,param,&aHGT) > 0){
						printf("\n=====>HGT : %ld-%ld -> %ld--%ld",GeneTree.ARETE[2*aHGT.source-1],GeneTree.ARETE[2*aHGT.source-2],GeneTree.ARETE[2*aHGT.destination-1],GeneTree.ARETE[2*aHGT.destination-2]);
						
						applyHGT(SpeciesTree.ADD,&GeneTree,aHGT.source,aHGT.destination);
						
						AdjustBranchLength(&GeneTree,SpeciesTree,0,1);
						
						deleteBipartition(DTGene,GeneTree);
						DTGene = (struct DescTree*)malloc((2*GeneTree.size-2-GeneTree.kt+1)*sizeof(struct DescTree));
						RechercherBipartition(GeneTree.ARETE,GeneTree.ADD,GeneTree.Root,GeneTree.Adjacence,DTGene,GeneTree.size,GeneTree.kt);
						
						ListesBranchesPourHGT(tab_sommets_selectionnes,GeneTree.ARETE,GeneTree.size,DTGene,tab_branches,&nb_branches);
					/*	printf("\n->Branches :");
						for(i=1;i<=nb_branches;i++){
							printf("%ld--%ld ",GeneTree.ARETE[2*tab_branches[i]-1],GeneTree.ARETE[2*tab_branches[i]-2]);
						}
						*/
						
					}

				}
				
			}		
		}
	
		free(tab_tous_les_sommets);
		free(tab_sommets_selectionnes);
		free(tab_branches);
		
		deleteBipartition(DTGene,GeneTree);
		deleteBipartition(DTSpecies,SpeciesTree);
		
		//initInputTree(&GeneTree);
		//copyInputTree(&GeneTree,tmpGeneTree,1,1);
		
		GeneTree.size = GeneTree.size-1;
		CreateSubStructures(&GeneTree,1,binaireGene);
		for(i=1;i<=GeneTree.size+1;i++){
			GeneTree.ADD[GeneTree.size+1][i] = GeneTree.ADD[i][GeneTree.size+1] = epsilon;
		}
		//CreateSubStructures(&GeneTree,1,binaireGene);
		AdjustBranchLength(&SpeciesTree,GeneTree,binaireSpecies,1);
		
		/*
		printf("\n\nLangue Tree : ");
		for(i=1;i<=SpeciesTree.size;i++){
			printf("\n%s\t",SpeciesTree.SpeciesName[i]);
			for(j=1;j<=SpeciesTree.size;j++){
				printf("%lf ",SpeciesTree.ADD[i][j]);
			}
		}
		printf("\n\nWord Tree : ");
		for(i=1;i<=GeneTree.size;i++){
			printf("\n%s\t",GeneTree.SpeciesName[i]);
			for(j=1;j<=GeneTree.size;j++){
				printf("%lf ",GeneTree.ADD[i][j]);
			}
		}
		*/
		
		DTGene = (struct DescTree*)malloc(3*(2*GeneTree.size-2-GeneTree.kt+1)*sizeof(struct DescTree));
		RechercherBipartitionSansRacine(GeneTree.ARETE,GeneTree.ADD,GeneTree.Adjacence,DTGene,GeneTree.size,GeneTree.kt);	
		
		DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
		RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);
		
		InitCriteria(&aCrit,SpeciesTree.size);
		computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
		printf("\nWBE-DETECTION : rf=%d n=%d rf/(2n-6)=%lf\n",aCrit.RF,SpeciesTree.size,(double)aCrit.RF/(2*(double)SpeciesTree.size-6));
		
		
		
	ReduceTree(SpeciesTree,GeneTree,&SpeciesTreeRed,&GeneTreeRed,&aMap,DTSpecies,DTGene,binaireSpecies,binaireGene);

	tmp = 0;

	int * listRef = (int *)malloc(2*SpeciesTree.size*sizeof(int));
	int * listJ   = (int *)malloc(2*SpeciesTree.size*sizeof(int));
	
	listRef[0] = listJ[0] = 0;
	//=============================== DETECTION DES TRANSFERTS ====================================
	if(strcmp(param.version,"consol")==0){
		printf("\nDetection #%d in progress",nbTree);
    printf("\n==========================");	
    printf("\nRF distance  = %2d",aCrit.RF);
		printf("\nLS criterion = %2.1lf",aCrit.LS);
		printf("\nBD criterion = %2.1lf",aCrit.BD);
		printf("\nQD criterion = %d",aCrit.QD);
    printf("\n==========================\n");
  }

	if(strcmp(param.scenario,"multiple")==0){
		cpt_hgt = findAllHGT(SpeciesTreeRed,GeneTreeRed,param,bestHGTRed);
		for(i=1;i<=cpt_hgt;i++){
			expandBestHGT(bestHGTRed[i],&bestHGT[i],aMap,DTSpecies,SpeciesTree);
		}
		sortHGT(bestHGT,cpt_hgt,param);
		if(cpt_hgt > param.nbhgt) cpt_hgt = param.nbhgt;
	}
	else if(strcmp(param.scenario,"unique")==0){
		if( strcmp(param.subtree,"yes") == 0 && strcmp(param.mode,"multicheck")==0 ){

			if(first==1){
				multicheckTab[0].m = 0;
				imc = 1;
			}
//			printf("\ndebut de la recherche");
      
			//==============================================================================
			//============ RECHERCHE DE TRANSFERTS - PLUSIEURS PAR TOUR ====================
			//==============================================================================

			trivial = (SpeciesTree.kt == 0)?0:1;	
			while( findBestHGTtab(SpeciesTreeRed,GeneTreeRed,param,bestHGTRed,&nbHgtFound,&trivial,listRef,listJ) > 0){
				int reelsHGT = 0;
				trivial = (SpeciesTree.kt == 0)?0:1;

//				if(strcmp(param.version,"consol")==0) 
				printf("\n\n[%d HGTs]", nbHgtFound);	
				
				if(first==1){
					multicheckTab[0].m ++;
					multicheckTab[imc].nbHgtFound = nbHgtFound;
				}
        
				int temoin_zero=-1;
				int cpt_hgt2 = cpt_hgt;
				for(i=0;i<nbHgtFound;i++){
					if(bestHGTRed[i].crit.RF == 0){
						temoin_zero = i;
					}
				}
				// printf("\ndans le main");
        
				for(i=0;i<nbHgtFound;i++){
					cpt_hgt++;		
					
				//	printf("\navant expand");
					expandBestHGT(bestHGTRed[i],&bestHGT[cpt_hgt],aMap,DTSpecies,SpeciesTree);
					bestHGT[cpt_hgt].trivial = 0;
					bestHGTRed[i].listSource = NULL;
					bestHGTRed[i].listDestination = NULL;
					
					printf("\n-------------------------------------------------------\n");
					double avgDiff = calculDifferenceMoyenne(bestHGT[cpt_hgt],GeneTree.ADD,GeneTree.size);

         /* if( (avgDiff > param.avgdiff) && ((param.constraints == 3) || (param.constraints == 2))){ //0.35
            printf("\nCONTRAINTE 2 APLLIQUEE");
						bestHGT[cpt_hgt].valide = 0; 
						continue; 
          }
*/

					char cmd[1000];
          char buffer[100];

          sprintf(cmd,"perl ageWBE.pl %s _src_ ", param.inputfile);
          for(int m=1;m<=bestHGT[cpt_hgt].listSource[0];m++){
						//printf("%s ",SpeciesTree.SpeciesName[bestHGT[cpt_hgt].listSource[m]]);
						sprintf(cmd + strlen(cmd), "%s ",SpeciesTree.SpeciesName[bestHGT[cpt_hgt].listSource[m]]);
					}
					//printf("==> ");
					sprintf(cmd + strlen(cmd), "_dest_ ");
					for(int m=1;m<=bestHGT[cpt_hgt].listDestination[0];m++){
						//printf("%s ",SpeciesTree.SpeciesName[bestHGT[cpt_hgt].listDestination[m]]);
						sprintf(cmd + strlen(cmd),"%s ",SpeciesTree.SpeciesName[bestHGT[cpt_hgt].listDestination[m]]);
					}
          
          //== Calcul de l'age du transfert => buffer
          FILE* file = popen(cmd, "r");
          fgets(buffer, 100, file);
          
          //printf("\ncmd= %s , age=%s", cmd, buffer);

          //double new_avg = param.avgdiffblock - param.avgdiff * log( atof(buffer) );
          
          //== Distance de levenstein moyenne => myAvgDist
          double myDist;
          double myTotalDist=0;
          int cptDist=0;
          for(int m=1;m<=bestHGT[cpt_hgt].listSource[0];m++){
             for(int n=1;n<=bestHGT[cpt_hgt].listDestination[0];n++){
                myDist = levenshtein_distance(trans[bestHGT[cpt_hgt].listSource[m]].langTrans,trans[bestHGT[cpt_hgt].listDestination[n]].langTrans);
                myTotalDist += myDist / fmax( (double)strlen(trans[bestHGT[cpt_hgt].listSource[m]].langTrans),(double)strlen(trans[bestHGT[cpt_hgt].listDestination[n]].langTrans) );
                cptDist++;
             }
					}
          
          
          //=========================================================================================
          //== Calcul de la probabilité d'avoir un transfert
          //=========================================================================================
          double D_l = myTotalDist / cptDist;
          double Age = atof(buffer);
          double p = (1-exp(-Age/param.c1)) * pow((1-D_l),param.c2);

          printf("\np=%lf, Age=%lf, Distance moyenne entre les 2 groupes de mots=%lf ", p, Age , D_l );
          

          if(( p < 0.5) && ((param.constraints == 3) || (param.constraints == 2))){ //0.35
            //printf("\nCONTRAINTE 2 APLLIQUEE");
						bestHGT[cpt_hgt].valide = 0; 
						continue; 
          }
          
					if((temoin_zero != -1)&&(i!=temoin_zero)){
						bestHGT[cpt_hgt].valide = 0; 
						//printf("\ndans temoin_zero");
						continue;  
					}
          	  
					if((bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_A) || 
						(bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_B) || 
						(bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_A) || 
						(bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_B)){
						bestHGT[cpt_hgt].trivial = 1;
					}
					
					
					
					if((bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_A && SpeciesTree.degre[bestHGT[cpt_hgt].source_A] == 3) ||
					   (bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_B && SpeciesTree.degre[bestHGT[cpt_hgt].source_A] == 3) ||
					   (bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_A && SpeciesTree.degre[bestHGT[cpt_hgt].source_B] == 3) ||
					   (bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_B && SpeciesTree.degre[bestHGT[cpt_hgt].source_B] == 3) ){

						bestHGT[cpt_hgt].valide = 0;
					//	if(strcmp(param.version,"consol")==0){
							printf("\ninutile=> HGT #%d : [%2d--%2d] -> [%2d--%2d]",cpt_hgt,bestHGT[cpt_hgt].source_A,bestHGT[cpt_hgt].source_B,bestHGT[cpt_hgt].dest_A,bestHGT[cpt_hgt].dest_B);
					//	}
						//printf("\nTransfert %d inutile ",cpt_hgt);
						continue;
					}
					
					printf("\n=====> ");
					for(int m=1;m<=bestHGT[cpt_hgt].listDestination[0];m++)
						printf("%d ",bestHGT[cpt_hgt].listDestination[m]);
						
					printf("\n=====> nbGroupesBloques=%d",nbGroupesBloques);
					int temoinZ = 0;  //= 1 - On ne casse pas le groupe
					int temoinY = 0;  //= 
					int nbSimilaire = 0;
					int nbSimilaire1 = 0;
					for(int n=0;n<nbGroupesBloques;n++){
						nbSimilaire = nbSimilaire1 = 0;
						for(int k=1;k<=groupesBloques[n][0];k++){
							for(int m=1;m<=bestHGT[cpt_hgt].listDestination[0];m++){
								if(groupesBloques[n][k] == bestHGT[cpt_hgt].listDestination[m]){
									nbSimilaire++;
								}
							}
							for(int m1=1;m1<=bestHGT[cpt_hgt].listSource[0];m1++){
								if(groupesBloques[n][k] == bestHGT[cpt_hgt].listSource[m1]){
									nbSimilaire1++;
								}
							}
						}
            printf("nbSimilaire = %d",nbSimilaire);
						if((nbSimilaire < groupesBloques[n][0]) && (nbSimilaire > 0))
							temoinZ = 1;
						if((nbSimilaire > 0) &&(nbSimilaire1 > 0))
							temoinY = 1; // 1 = il y a des éléments du même groupe à la source et a la destination 
					}
          printf("\n==>temoinZ=%d, temoinY=%d", temoinZ, temoinY);  
					if((temoinZ == 1) && (temoinY == 0)){
						printf("\n====---------=====");
						bestHGT[cpt_hgt].valide = 0;
						continue;
					}
					
					
					
					
			
					applyHGT2(GeneTree.ADD,&SpeciesTree,bestHGT[cpt_hgt].source,bestHGT[cpt_hgt].destination);
					reelsHGT++;
					
					if(strcmp(param.version,"consol")==0){
						printf("\nHGT #%d : [%2d--%2d] -> [%2d--%2d]",cpt_hgt,bestHGT[cpt_hgt].source_A,bestHGT[cpt_hgt].source_B,bestHGT[cpt_hgt].dest_A,bestHGT[cpt_hgt].dest_B);
					}
				}
				printf("\n==> %d",reelsHGT);
				if (reelsHGT == 0){
					//strcpy(param.subtree,"no");
					for(int n=0;n<=listRef[0];n++){
						listJ[n] = listRef[n];
					}
				}
				         
//			AdjustBranchLength(&SpeciesTree,GeneTree,0,1);
				//exit(0);
        computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
				
				//if(strcmp(param.version,"consol")==0){
				  //printf("\nCriteres apres ajout des transferts trouves : ");
					printf("\nRF = %d | LS = %1.2lf | BD = %1.2lf | QD = %d\n",aCrit.RF,aCrit.LS,aCrit.BD,aCrit.QD);	 
				//}
				if(first==1){
					multicheckTab[imc].LS = aCrit.LS;
					multicheckTab[imc].RF = aCrit.RF;
					multicheckTab[imc].BD = aCrit.BD;
					multicheckTab[imc].QD = aCrit.QD;
					imc++;
				}	
				
/*				if((bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_A && SpeciesTree.degre[bestHGT[cpt_hgt].source_A] == 3) ||
					   (bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_B && SpeciesTree.degre[bestHGT[cpt_hgt].source_A] == 3) ||
					   (bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_A && SpeciesTree.degre[bestHGT[cpt_hgt].source_B] == 3) ||
					   (bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_B && SpeciesTree.degre[bestHGT[cpt_hgt].source_B] == 3) ){

						bestHGT[cpt_hgt].valide = 0;
						if(strcmp(param.version,"consol")==0){
							printf("\ninutile=> HGT #%d : [%2d--%2d] -> [%2d--%2d]",cpt_hgt,bestHGT[cpt_hgt].source_A,bestHGT[cpt_hgt].source_B,bestHGT[cpt_hgt].dest_A,bestHGT[cpt_hgt].dest_B);
						}
			}*/	
 			  if((cpt_hgt >= param.nbhgt)||(aCrit.RF == 0)) break;
				
				
				//== A REGLER (probleme : segmentation fault) 
				
				deleteBipartition(DTSpecies,SpeciesTreeCurrent);			
				
				copyInputTree(&SpeciesTreeCurrent,SpeciesTree,1,1);
				
				DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
				RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);
				FreeMemory_InputTreeReduced(&SpeciesTreeRed,SpeciesTreeRed.size);
				FreeMemory_InputTreeReduced(&GeneTreeRed,GeneTreeRed.size);
				initInputTree(&SpeciesTreeRed);
				initInputTree(&GeneTreeRed);
				free(aMap.map);
				free(aMap.gene);
				free(aMap.species);
				
				ReduceTree(SpeciesTree,GeneTree,&SpeciesTreeRed,&GeneTreeRed,&aMap,DTSpecies,DTGene,binaireSpecies,binaireGene);
				
		}
		free(aMap.map);
		free(aMap.gene);
		free(aMap.species);
		deleteBipartition(DTSpecies,SpeciesTreeCurrent);
	
		int retour=0;
		int cpt,nbTours,j;
	
//  do{						
			printf("\nHGT-DETECTION : recherche des transferts inutiles");
			retour = DeleteUseLessHGT(cpt_hgt,bestHGT,SpeciesTree,FirstTree); 
			//if(strcmp(param.version,"consol")==0){
			   printf("\nNombre de transferts inutile = %d",retour);
     // }
	//	}while(retour > 0);
		cpt=1;
		if(first==1){
			for(i=1;i<=multicheckTab[0].m;i++){
				nbTours = multicheckTab[i].nbHgtFound;
				for(j=1;j<=nbTours;j++){
					if(bestHGT[cpt].valide == 0){
						multicheckTab[i].nbHgtFound --;
					}
					cpt++;
				}
			}
		}
	}
	else
		{	
			int initial=1;
     
				
			while(findBestHGT(initial,SpeciesTreeRed,GeneTreeRed,param,&bestHGTRed[cpt_hgt+1]) > 0){
				cpt_hgt++;		
				expandBestHGT(bestHGTRed[cpt_hgt],&bestHGT[cpt_hgt],aMap,DTSpecies,SpeciesTree);
              	bestHGT[cpt_hgt].trivial = 0;
          		if((bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_A) || 
             		(bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_B) || 
             		(bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_A) || 
             		(bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_B)){
              		bestHGT[cpt_hgt].trivial = 1;
          		}	
				if((bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_A && SpeciesTree.degre[bestHGT[cpt_hgt].source_A] == 3) ||
					   (bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_B && SpeciesTree.degre[bestHGT[cpt_hgt].source_A] == 3) ||
					   (bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_A && SpeciesTree.degre[bestHGT[cpt_hgt].source_B] == 3) ||
					   (bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_B && SpeciesTree.degre[bestHGT[cpt_hgt].source_B] == 3) ){

						bestHGT[cpt_hgt].valide = 0;
						if(strcmp(param.version,"consol")==0){
							printf("\ninutile=> HGT #%d : [%2d--%2d] -> [%2d--%2d]",cpt_hgt,bestHGT[cpt_hgt].source_A,bestHGT[cpt_hgt].source_B,bestHGT[cpt_hgt].dest_A,bestHGT[cpt_hgt].dest_B);
						}
					continue;
				}

				bestHGTRed[i].listSource = NULL;
				bestHGTRed[i].listDestination = NULL;
				
				applyHGT(GeneTree.ADD,&SpeciesTree,bestHGT[cpt_hgt].source,bestHGT[cpt_hgt].destination);
				AdjustBranchLength(&SpeciesTree,GeneTree,0,1);
				
				computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
				loadCriteria(aCrit,&(bestHGT[cpt_hgt]));
				
				if(strcmp(param.version,"consol")==0){
					printf("\nHGT #%d %d--%d -> %d--%d",cpt_hgt,bestHGT[cpt_hgt].source_A,bestHGT[cpt_hgt].source_B,bestHGT[cpt_hgt].dest_A,bestHGT[cpt_hgt].dest_B);
					printf("\nRF = %d, LS = %lf, BD = %lf QD = %d\n",aCrit.RF,aCrit.LS,aCrit.BD,aCrit.QD);	 
				} 
				if(bestHGT[cpt_hgt].crit.RF == 0) break;// || bestHGT[cpt_hgt].crit.LS < epsilon) break;

				if(cpt_hgt >= param.nbhgt) break;
			
				deleteBipartition(DTSpecies,SpeciesTreeCurrent);
				copyInputTree(&SpeciesTreeCurrent,SpeciesTree,1,1);
				DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
				RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);
				
				free(aMap.map);
				free(aMap.gene);
				free(aMap.species);
				FreeMemory_InputTreeReduced(&SpeciesTreeRed,SpeciesTreeRed.size);
				FreeMemory_InputTreeReduced(&GeneTreeRed,GeneTreeRed.size);
				initInputTree(&SpeciesTreeRed);
				initInputTree(&GeneTreeRed);
					
				ReduceTree(SpeciesTree,GeneTree,&SpeciesTreeRed,&GeneTreeRed,&aMap,DTSpecies,DTGene,binaireSpecies,binaireGene);
				//printf("<br> prochlorococcus = %d", aMap.species[15]);
			}
			
			initial=0; 
			while( findBestHGT(initial,SpeciesTreeRed,GeneTreeRed,param,&bestHGTRed[cpt_hgt+1]) > 0){
			printf("\nicit");	
				cpt_hgt++;		
				expandBestHGT(bestHGTRed[cpt_hgt],&bestHGT[cpt_hgt],aMap,DTSpecies,SpeciesTree);
              	bestHGT[cpt_hgt].trivial = 0;
          		if((bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_A) || 
             		(bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_B) || 
             		(bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_A) || 
             		(bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_B)){
              		bestHGT[cpt_hgt].trivial = 1;
          		}	
					if((bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_A && SpeciesTree.degre[bestHGT[cpt_hgt].source_A] == 3) ||
					   (bestHGT[cpt_hgt].source_A == bestHGT[cpt_hgt].dest_B && SpeciesTree.degre[bestHGT[cpt_hgt].source_A] == 3) ||
					   (bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_A && SpeciesTree.degre[bestHGT[cpt_hgt].source_B] == 3) ||
					   (bestHGT[cpt_hgt].source_B == bestHGT[cpt_hgt].dest_B && SpeciesTree.degre[bestHGT[cpt_hgt].source_B] == 3) ){

						bestHGT[cpt_hgt].valide = 0;
						if(strcmp(param.version,"consol")==0){
							printf("\ninutile=> HGT #%d : [%2d--%2d] -> [%2d--%2d]",cpt_hgt,bestHGT[cpt_hgt].source_A,bestHGT[cpt_hgt].source_B,bestHGT[cpt_hgt].dest_A,bestHGT[cpt_hgt].dest_B);
						}
						//printf("\nTransfert %d inutile ",cpt_hgt);
						continue;
					}

				bestHGTRed[i].listSource = NULL;
				bestHGTRed[i].listDestination = NULL;
					
			//printf("\navant apply dans findBestHGT");
				applyHGT(GeneTree.ADD,&SpeciesTree,bestHGT[cpt_hgt].source,bestHGT[cpt_hgt].destination);
			//printf("\napres apply dans findBestHGT");
			
			//printf("\ntoto");
				AdjustBranchLength(&SpeciesTree,GeneTree,0,1);
			
			/*for(i=1;i<=2*GeneTree.size-3;i++)
				printf("\n%d-%d -->%lf",GeneTree.ARETE[2*i-1],GeneTree.ARETE[2*i-2],GeneTree.LONGUEUR[i-1]);*/
				
			
			computeCriteria(SpeciesTree.ADD,GeneTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
			loadCriteria(aCrit,&(bestHGT[cpt_hgt]));
			//printf("\navant Delete");
				//if (cpt_hgt > 1)
				//	DeleteUseLessHGT(cpt_hgt,bestHGT,SpeciesTree,FirstTree);
				/*	
					int retour;
					do{						
						printf("\ntest des transferts inutiles");
						retour = DeleteUseLessHGT(cpt_hgt,bestHGT,SpeciesTree,FirstTree); 
						printf("\nretour = %d",retour);
					}while(retour > 0);
				}*/
			//printf("\napres Delete");
			//printf("ici = %s, %d",param.version,cpt_hgt);
				if(strcmp(param.version,"consol")==0){
					printf("\nHGT #%d %d--%d -> %d--%d",cpt_hgt,bestHGT[cpt_hgt].source_A,bestHGT[cpt_hgt].source_B,bestHGT[cpt_hgt].dest_A,bestHGT[cpt_hgt].dest_B);
					printf("\nRF = %d, LS = %lf, BD = %lf QD = %d\n",aCrit.RF,aCrit.LS,aCrit.BD,aCrit.QD);	 
				}
			//printf("icitttt = %d , %d",bestHGT[cpt_hgt].crit.RF,param.nbhgt);
				if(bestHGT[cpt_hgt].crit.RF == 0) {
					int retour;
					do{						
					//	printf("\ntest des transferts inutiles");
						retour = DeleteUseLessHGT(cpt_hgt,bestHGT,SpeciesTree,FirstTree); 
					//	printf("\nretour = %d",retour);
					}while(retour > 0);
				
					break;// || bestHGT[cpt_hgt].crit.LS < epsilon) break;
				}
				if(cpt_hgt >= param.nbhgt) break;
			
				deleteBipartition(DTSpecies,SpeciesTreeCurrent);
				copyInputTree(&SpeciesTreeCurrent,SpeciesTree,1,1);
				free(aMap.map);
				free(aMap.gene);
				free(aMap.species);
				DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
				RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);
				FreeMemory_InputTreeReduced(&SpeciesTreeRed,SpeciesTreeRed.size);
				FreeMemory_InputTreeReduced(&GeneTreeRed,GeneTreeRed.size);
				initInputTree(&SpeciesTreeRed);
				initInputTree(&GeneTreeRed);
			
				ReduceTree(SpeciesTree,GeneTree,&SpeciesTreeRed,&GeneTreeRed,&aMap,DTSpecies,DTGene,binaireSpecies,binaireGene);
			//printf("==>%d",SpeciesTreeRed.size);
				//printf("<br> prochlorococcus = %d", aMap.map[6]);
		}
		free(aMap.map);
		free(aMap.gene);
		free(aMap.species);
		deleteBipartition(DTSpecies,SpeciesTreeCurrent);
		//FreeMemory_InputTreeReduced(&SpeciesTreeRed,SpeciesTreeRed.size);
		//FreeMemory_InputTreeReduced(&GeneTreeRed,GeneTreeRed.size);
	   }
		
	}
	else if(strcmp(param.scenario,"all")==0){
		printf("\nTous les scenarios de taille minimale\n");
		//findAllMinimalScenario(SpeciesTree,GeneTree);
		exit(0);
	}
	
//==============================================================================
//============================ TRAITEMENT DES RESULTATS ========================
//==============================================================================
  outHGT = (struct HGT*)malloc(2*param.nbhgt*sizeof(struct HGT)); 
  nbHGT = formatResult(bestHGT,cpt_hgt,outHGT,FirstTree);
  printHGT(results,multicheckTab,param.mode,RFref,out,FirstTree,outHGT,nbHGT,NULL,param.subtree,param.bootmin);

	printf("\nHGT-DETECTION : apres printHGT");
	if((strcmp(param.version,"web")==0) && (strcmp(param.printWeb,"yes")==0)){
		saveTree(param.outputWeb,FirstTree,outHGT,1,nbHGT,param.subtree,param.scenario,NULL);
	}

//==============================================================================
//============================ LIBERATION DE LA MEMOIRE ========================
//==============================================================================		
	deleteBipartitionSansRacine(DTGene,GeneTree.size);
	FreeMemory_InputTreeReduced(&SpeciesTreeRed,SpeciesTreeRed.size);
	FreeMemory_InputTreeReduced(&GeneTreeRed,GeneTreeRed.size);
	FreeMemory_InputTree(&SpeciesTreeCurrent,SpeciesTreeCurrent.size);
	FreeMemory_InputTree(&GeneTree,GeneTree.size);
	FreeMemory_InputTree(&SpeciesTree,SpeciesTree.size);
	if(bootstrap !=1)
		FreeMemory_InputTree(&FirstTree,FirstTree.size);
	FreeCriteria(&aCrit,SpeciesTree.size);
	
	for(i=1;i<=cpt_hgt;i++){
		free(bestHGT[i].listSource);
		free(bestHGT[i].listDestination);
	}
	free(bestHGTRed);
	free(bestHGT); 
	
  fclose(results);
 
//==============================================================================
//============================ SUPPRESSION DES FICHIERS ========================
//==============================================================================
/*  remove(param.input);
	remove("t1");
	remove("t2");
	//remove("results.txt");
	//remove(param.inputfile);
	//remove(param.outputfile);	
  if(strcmp(param.version,"web")!=0){
	   remove(param.hgtResultFile);
	   remove(param.speciesTree);
	   remove(param.geneTree);
	   remove(param.speciesRootfile);
	   remove(param.geneRootfile);
	   remove(param.speciesTreeWeb);
	   remove(param.geneTreeWeb);
	}*/
	//remove(param.outputWeb);
	
  //if(strcmp(param.version,"consol") == 0)
		printf("\nvaleur retournee : %d \n",nbHGT);
  
  exit(nbHGT);
}
