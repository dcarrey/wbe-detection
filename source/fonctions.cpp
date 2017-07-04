//===================================================================================
//============================ DEFINITION DES FONCTIONS =============================
//===================================================================================

void printHGTandGroups(struct HGT *aHGT, struct InputTree SpeciesTree){

  printf("\nAffichage du transfert : \n");
  for(int m=1;m<=aHGT->listSource[0];m++)
	  printf("%s ",SpeciesTree.SpeciesName[aHGT->listSource[m]]);
	printf("=>");
  for(int m=1;m<=aHGT->listDestination[0];m++)
	  printf("%s ",SpeciesTree.SpeciesName[aHGT->listDestination[m]]);

}


int estUnSousGroupe(int * tab1,int * tab2){

	int temoin=1;
	int trouve;
	for(int i=1;i<=tab1[0];i++){
		trouve=0;
		for(int j=1;j<=tab2[0];j++){
			if(tab1[i] == tab2[j]){
				trouve=1;
			}
		}
		if(trouve==0){
			temoin = 0;
		}
	}

	return temoin;
}

/****************************************/
/*Implementation of Levenshtein distance*/
/****************************************/

double minimum(double a,double b,double c)
/*Gets the minimum of three values*/
{
  double min=a;
  if(b<min)
    min=b;
  if(c<min)
    min=c;
  return min;
}

int delDoublons(const char * chaine, char * chaine_res,char *tabDoublons)
{
	chaine_res[0] = chaine[0];
	int j=1;
	int i;
	int temoin = 1;
	int nbDoublons = 0;
	for(i=1;i<strlen(chaine);i++){
		if(chaine[i] != chaine_res[j-1]){
			chaine_res[j++] = chaine[i];
		}
		else{
			if(temoin == 1){
				tabDoublons[nbDoublons] = chaine[i];
				nbDoublons++;
				temoin = 0;
			}
		}
	}
	chaine_res[j] = '\0';
	tabDoublons[nbDoublons] = '\0';

	return (nbDoublons);
}

double levenshtein_distance(char *s,char*t)
/*Compute levenshtein distance between s and t*/
{
  //Step 1
  int k,i,j,n,m;
  double cost,*d,distance=0;
  n=strlen(s);
  m=strlen(t);
  if(n!=0&&m!=0)
  {
    d=(double*)malloc((sizeof(double))*(m+1)*(n+1));
    m++;
    n++;
    //Step 2
    for(k=0;k<n;k++)
	d[k]=k;
    for(k=0;k<m;k++)
      d[k*n]=k;
    //Step 3 and 4
    for(i=1;i<n;i++)
      for(j=1;j<m;j++)
	{
        //Step 5
        if(s[i-1]==t[j-1])
          cost=0;
        else{
				cost = 1;
				//QG, VF, QW, WB 0.800
				if(((s[i-1] == 'Q') && (t[j-1] == 'G')) || ((s[i-1] == 'G') && (t[j-1] == 'Q'))){
					cost = 0.800;
				}
				if(((s[i-1] == 'V') && (t[j-1] == 'F')) || ((s[i-1] == 'F') && (t[j-1] == 'V'))){
					cost = 0.800;
				}
				if(((s[i-1] == 'O') && (t[j-1] == 'W')) || ((s[i-1] == 'W') && (t[j-1] == 'O'))){
					cost = 0.800;
				}
				if(((s[i-1] == 'W') && (t[j-1] == 'B')) || ((s[i-1] == 'B') && (t[j-1] == 'W'))){
					cost = 0.800;
				}
				//CT, QC,KQ, KG, IW, WJ, OU, UW, VW 0.600
				if(((s[i-1] == 'C') && (t[j-1] == 'T')) || ((s[i-1] == 'T') && (t[j-1] == 'C'))){
					cost = 0.600;
				}
				if(((s[i-1] == 'Q') && (t[j-1] == 'C')) || ((s[i-1] == 'C') && (t[j-1] == 'Q'))){
					cost = 0.600;
				}
				if(((s[i-1] == 'K') && (t[j-1] == 'C')) || ((s[i-1] == 'C') && (t[j-1] == 'K'))){
					cost = 0.600;
				}
				if(((s[i-1] == 'K') && (t[j-1] == 'Q')) || ((s[i-1] == 'Q') && (t[j-1] == 'K'))){
					cost = 0.600;
				}
				if(((s[i-1] == 'K') && (t[j-1] == 'G')) || ((s[i-1] == 'G') && (t[j-1] == 'K'))){
					cost = 0.600;
				}
				if(((s[i-1] == 'I') && (t[j-1] == 'Y')) || ((s[i-1] == 'Y') && (t[j-1] == 'I'))){
					cost = 0.600;
				}
				if(((s[i-1] == 'Y') && (t[j-1] == 'J')) || ((s[i-1] == 'J') && (t[j-1] == 'Y'))){
					cost = 0.600;
				}
				if(((s[i-1] == 'O') && (t[j-1] == 'U')) || ((s[i-1] == 'U') && (t[j-1] == 'O'))){
					cost = 0.600;
				}
				if(((s[i-1] == 'U') && (t[j-1] == 'W')) || ((s[i-1] == 'W') && (t[j-1] == 'U'))){
					cost = 0.600;
				}
				if(((s[i-1] == 'V') && (t[j-1] == 'W')) || ((s[i-1] == 'W') && (t[j-1] == 'U'))){
					cost = 0.600;
				}
				//NM, TD, IJ, VB 0.400
				if(((s[i-1] == 'N') && (t[j-1] == 'M')) || ((s[i-1] == 'M') && (t[j-1] == 'N'))){
					cost = 0.400;
				}
				if(((s[i-1] == 'T') && (t[j-1] == 'D')) || ((s[i-1] == 'D') && (t[j-1] == 'T'))){
					cost = 0.400;
				}
				if(((s[i-1] == 'I') && (t[j-1] == 'J')) || ((s[i-1] == 'J') && (t[j-1] == 'I'))){
					cost = 0.400;
				}
				if(((s[i-1] == 'V') && (t[j-1] == 'B')) || ((s[i-1] == 'B') && (t[j-1] == 'V'))){
					cost = 0.400;
				}
				// XH, PB, ZS H- 0.2
				if(((s[i-1] == 'X') && (t[j-1] == 'H')) || ((s[i-1] == 'H') && (t[j-1] == 'X'))){
					cost = 0.200;
				}
				if(((s[i-1] == 'P') && (t[j-1] == 'B')) || ((s[i-1] == 'B') && (t[j-1] == 'P'))){
					cost = 0.200;
				}
				if(((s[i-1] == 'Z') && (t[j-1] == 'S')) || ((s[i-1] == 'S') && (t[j-1] == 'Z'))){
					cost = 0.200;
				}
				// Dash
				if(((s[i-1] == '*') && (t[j-1] != '*')) || ((s[i-1] != '*') && (t[j-1] == '*'))){
					cost = 0.100;
				}
				if(((s[i-1] == '\'') && (t[j-1] != '\'')) || ((s[i-1] != '\'') && (t[j-1] == '\''))){
					cost = 0.100;
				}
				if(((s[i-1] == '-') && (t[j-1] != '-')) || ((s[i-1] != '-') && (t[j-1] == '-'))){
					cost = 0.100;
				}


			}
        //Step 6


			// # Cell d[j*n-i] equals the minimum of:
            //
            // # - The cell immediately above plus 1
            // # - The cell immediately to the left plus 1
            // # - The cell diagonally above and to the left plus the cost
            // #
            // # We can either insert a new char, delete a char or
            // # substitute an existing char (with an associated cost)
            // #
			// # Handle the case with H and dash (-)

			double cost_i=1.0;
			double cost_j=1.0;

			d[j*n+i]=minimum(d[(j-1)*n+i]+cost_j,d[j*n+i-1]+cost_i,d[(j-1)*n+i-1]+cost);

			if (d[j*n+i]==(d[j*n+i-1]+cost_i)&&(s[i-1] == 'H')) {
				 cost_i=0.2;
			}
			if (d[j*n+i]==(d[(j-1)*n+i]+cost_j)&&(t[j-1] == 'H')) {
				 cost_j=0.2;
			}
			if (cost_i != 1.0 || cost_j != 1) {
				 d[j*n+i] = minimum(d[(j-1)*n+i]+cost_j,d[j*n+i-1]+cost_i,d[(j-1)*n+i-1]+cost);
			}

      }
       distance=d[n*m-1];
    free(d);
    return distance;
  }
  else
    return -1; //a negative return value means that one or both strings are empty.
}

double levenshtein_distance_old(char *s,char*t)
/*Compute levenshtein distance between s and t*/
{
  //Step 1
  int k,i,j,n,m;
  double cost,*d,distance;
  n=strlen(s);
  m=strlen(t);
  if(n!=0&&m!=0)
  {
    d=(double*)malloc((sizeof(double))*(m+1)*(n+1));
    m++;
    n++;
    //Step 2
    for(k=0;k<n;k++)
	d[k]=k;
    for(k=0;k<m;k++)
      d[k*n]=k;
    //Step 3 and 4
    for(i=1;i<n;i++)
      for(j=1;j<m;j++)
	{
        //Step 5
        if(s[i-1]==t[j-1])
          cost=0;
        else{
				cost = 1;
				if(((s[i-1] == 'D') && (t[j-1] == 'T')) || ((s[i-1] == 'T') && (t[j-1] == 'D'))){
					cost = 0.5;
				}
				if(((s[i-1] == 'V') && (t[j-1] == 'F')) || ((s[i-1] == 'F') && (t[j-1] == 'V'))){
					cost = 0.5;
				}
				if(((s[i-1] == 'I') && (t[j-1] == 'Y')) || ((s[i-1] == 'Y') && (t[j-1] == 'I'))){
					cost = 0.5;
				}
				if(((s[i-1] == 'Q') && (t[j-1] == 'K')) || ((s[i-1] == 'K') && (t[j-1] == 'Q'))){
					cost = 0.5;
				}
				if(((s[i-1] == 'W') && (t[j-1] == 'V')) || ((s[i-1] == 'V') && (t[j-1] == 'W'))){
					cost = 0.5;
				}
				if(((s[i-1] == 'K') && (t[j-1] == 'C')) || ((s[i-1] == 'C') && (t[j-1] == 'K'))){
					cost = 0.5;
				}
				if(((s[i-1] == 'S') && (t[j-1] == 'Z')) || ((s[i-1] == 'Z') && (t[j-1] == 'S'))){
					cost = 0.5;
				}
				if(((s[i-1] == 'C') && (t[j-1] == 'Q')) || ((s[i-1] == 'Q') && (t[j-1] == 'C'))){
					cost = 0.5;
				}
				if(((s[i-1] == '*') && (t[j-1] != '*')) || ((s[i-1] != '*') && (t[j-1] == '*'))){
					cost = 0.25;
				}
			}
        //Step 6
        d[j*n+i]=minimum(d[(j-1)*n+i]+1,d[j*n+i-1]+1,d[(j-1)*n+i-1]+cost);
      }
    distance=d[n*m-1];
    free(d);
    return distance;
  }
  else
    return -1; //a negative return value means that one or both strings are empty.
}
double levenshtein_distance1(char *s,char*t)
/*Compute levenshtein distance between s and t*/
{
  //Step 1
   int k,i,j,n,m;
  double cost,*d,distance;
  n=strlen(s);
  m=strlen(t);
  if(n!=0&&m!=0)
  {
    d=(double*)malloc((sizeof(double))*(m+1)*(n+1));
    m++;
    n++;
    //Step 2
    for(k=0;k<n;k++)
	d[k]=k;
    for(k=0;k<m;k++)
      d[k*n]=k;
    //Step 3 and 4
    for(i=1;i<n;i++)
      for(j=1;j<m;j++)
	{
        //Step 5
        if(s[i-1]==t[j-1])
          cost=0;
        else
          cost=1;
        //Step 6
        d[j*n+i]=minimum(d[(j-1)*n+i]+1,d[j*n+i-1]+1,d[(j-1)*n+i-1]+cost);
      }
    distance=d[n*m-1];
    free(d);
    return distance;
  }
  else
    return -1; //a negative return value means that one or both strings are empty.
}


void deleteFiles(){

}

//=================================================================
//==
//=================================================================
void addLeafAndUpdate(struct InputTree *aTree, int choix){

	int n = aTree->size;
	int i,j;

/*	for(i=1;i<=2*n-3-aTree->kt;i++)
		printf("\n%d-%d -> %lf",aTree->ARETE[2*i-1],aTree->ARETE[2*i-2],aTree->LONGUEUR[i-1]);
	printf("\n");
*/
	for(i=1;i<=2*n-3-aTree->kt;i++){
		if(aTree->ARETE[2*i-1] > n) aTree->ARETE[2*i-1] += 1 ;
		if(aTree->ARETE[2*i-2] > n) aTree->ARETE[2*i-2] += 1 ;
	}


	//= ajout des 2 nouvelles branches
	aTree->ARETE[2*(2*n-2-aTree->kt)-1] = aTree->ARETE[2*choix-1];
	aTree->ARETE[2*(2*n-2-aTree->kt)-2] = 2*n-aTree->kt;
	aTree->LONGUEUR[(2*n-2-aTree->kt)-1] = ((aTree->LONGUEUR[choix-1]/2.0)>5*epsilon)?(aTree->LONGUEUR[choix-1]/2.0):5*epsilon;

	aTree->ARETE[2*(2*n-1-aTree->kt)-1] = aTree->ARETE[2*choix-2];
	aTree->ARETE[2*(2*n-1-aTree->kt)-2] = 2*n-aTree->kt;
	aTree->LONGUEUR[(2*n-1-aTree->kt)-1] = ((aTree->LONGUEUR[choix-1]/2.0)>5*epsilon)?(aTree->LONGUEUR[choix-1]/2.0):5*epsilon;

	aTree->ARETE[2*choix-1] = n+1;
	aTree->ARETE[2*choix-2] = 2*n-aTree->kt;
	aTree->LONGUEUR[choix-1]= 5*epsilon + 1.0;

	aTree->size = aTree->size+1;

	/*for(i=1;i<=2*aTree->size-3-aTree->kt;i++)
		printf("\n%d-%d -> %lf",aTree->ARETE[2*i-1],aTree->ARETE[2*i-2],aTree->LONGUEUR[i-1]);*/

	loadAdjacenceMatrix(aTree->Adjacence,aTree->ARETE, aTree->LONGUEUR,n+1,aTree->kt);
	/*for(i=1;i<=2*aTree->size-2-aTree->kt;i++){
		printf("\n");
		for(j=1;j<=2*aTree->size-2-aTree->kt;j++)
			printf("%lf ",aTree->Adjacence[i][j]);
	}
	printf("\n");*/
	Floyd(aTree->Adjacence,aTree->ADD,aTree->Input,n+1,aTree->kt);//5eme fois
//  global_cpt5++;

	strcpy(aTree->SpeciesName[aTree->size],"Root");

}

int findRoot(int * tab, int size){
	int i,j;
	int trouve;
	for(i=1;i<=size;i++){
		trouve=0;
		for(j=1;j<=tab[0];j++){
			if(i==tab[j]) trouve=1;
		}
		if(trouve==0){
			return i;
		}
	}
	return 0;
}

void ListeSommets_taille_0(double ** matrix,int * tab_sommet,int size){

	int cpt=1;
	int cpt_init;
	int nouveau_cas=0;

	//== recherche du premier element admissible
	while(cpt <= size){
		if(tab_sommet[cpt] == 0)
			break;
		else
			cpt++;
	}

	if(cpt < size){
		//== recherche des distances nulles
		cpt_init = cpt++;
		while(cpt <= size){
			if(matrix[cpt_init][cpt] < 100*epsilon ){
				tab_sommet[cpt] = 1;
				nouveau_cas = 1;
			}
			cpt++;
		}
		if(nouveau_cas == 1){
			tab_sommet[cpt_init] = 1;
		}
		else
			tab_sommet[cpt_init] = 2;
	}
	else
		tab_sommet[0] = TRUE; //== aucun elements trouves
}

int Est_un_sous_ensemble_exact(int * sous_ensemble,int * ensemble){

	int temoin;

	if(sous_ensemble[0] > ensemble[0]){
		return 0;
	}
	else{
		temoin=0;
		for(int i=1;i<=sous_ensemble[0];i++){
			for(int j=1;j<=ensemble[0];j++){
				if(sous_ensemble[i] == ensemble[j]){
					temoin++;
				}
			}
		}
		if(temoin == sous_ensemble[0]) return 1;
	}
	return 0;
}

void ListesBranchesPourHGT(int *tab_sommets,long int * ARETE, int taille,struct DescTree *DT,int *tab_branches,int *nb_branches){

	(*nb_branches) = 0;
	for(int i=1;i<=2*taille-3;i++){
		DT[ARETE[2*i-1]].Tableau[0] = DT[ARETE[2*i-1]].nbSommet;
		DT[ARETE[2*i-2]].Tableau[0] = DT[ARETE[2*i-2]].nbSommet;

		if( ((Est_un_sous_ensemble_exact(DT[ARETE[2*i-1]].Tableau,tab_sommets)==1) || (Est_un_sous_ensemble_exact(DT[ARETE[2*i-2]].Tableau,tab_sommets)==1)) ||
		    ((Est_un_sous_ensemble_exact(DT[ARETE[2*i-2]].Tableau,tab_sommets)==1) || (Est_un_sous_ensemble_exact(DT[ARETE[2*i-1]].Tableau,tab_sommets)==1)) ){
			(*nb_branches) ++;

			tab_branches[(*nb_branches)] = i;

		}
	}
}

void ListesBranchesPourHGTSansRacine(int *tab_sommets,long int * ARETE, int taille,struct DescTree *DT,int *tab_branches,int *nb_branches){

	(*nb_branches) = 0;
	int A,B;
	for(int i=taille+1;i<=2*taille-3;i++){
		A = ARETE[2*i-1];
		B = ARETE[2*i-2];

		if((Est_un_sous_ensemble_exact(DT[A].Tableau,tab_sommets)==1) || (Est_un_sous_ensemble_exact(DT[B].Tableau,tab_sommets)==1)){
			(*nb_branches) ++;

			tab_branches[(*nb_branches)] = i;

		}
	}
}


void SAVEASNewick(double *LONGUEUR, long int *ARETE,int nn,const char* t)
{
	int n,root,a;
	int Ns;
	int i, j, sd, sf, *Suc, *Fre, *Tree, *degre, *Mark;
	double *Long;
	int *boot;
  printf("\nOn et vraiment au debut !!");
	char *string = (char*)malloc(10000);
	n = nn;
	Ns=2*n-3;

	double * bootStrap= NULL;

	Suc =(int*) malloc((2*n) * sizeof(int));
	Fre =(int*) malloc((2*n) * sizeof(int));
	degre =(int*) malloc((2*n) * sizeof(int));
	Long = (double*) malloc((2*n) * sizeof(double));
	boot = (int*) malloc((2*n) * sizeof(int));
	Tree = (int*) malloc((2*n) * sizeof(int));
	Mark =(int*) malloc((2*n) * sizeof(int));

  printf("\nOn et ici au debut !!");
	if ((degre==NULL)||(Mark==NULL)||(string==NULL)||(Suc==NULL)||(Fre==NULL)||(Long==NULL)||(Tree==NULL)||(ARETE==NULL)||(LONGUEUR==NULL))
		{ printf("Tree is too large to be saved"); return;}

	for (i=1;i<=2*n-3;i++)
	{

		if (i<=n) degre[i]=1;
		else degre[i]=3;
	} degre[2*n-2]=3;

	root=Ns+1;

	int cpt=0;

	for (;;)
	{
		cpt++;
		if(cpt > 1000000) exit(1);
		a=0; a++;
		for (j=1;j<=2*n-2;j++)
			Mark[j]=0;

		for (i=1;i<=2*n-3;i++)
		{
			if ((degre[ARETE[2*i-2]]==1)&&(degre[ARETE[2*i-1]]>1)&&(Mark[ARETE[2*i-1]]==0)&&(Mark[ARETE[2*i-2]]==0))
			{
				Tree[ARETE[2*i-2]]=ARETE[2*i-1]; degre[ARETE[2*i-1]]--; degre[ARETE[2*i-2]]--; Mark[ARETE[2*i-1]]=1; Mark[ARETE[2*i-2]]=1;
				Long[ARETE[2*i-2]]=LONGUEUR[i-1];
				if(bootStrap != NULL) boot[ARETE[2*i-2]] = (int) bootStrap[i-1];

			}
			else if ((degre[ARETE[2*i-1]]==1)&&(degre[ARETE[2*i-2]]>1)&&(Mark[ARETE[2*i-1]]==0)&&(Mark[ARETE[2*i-2]]==0))
			{
				Tree[ARETE[2*i-1]]=ARETE[2*i-2]; degre[ARETE[2*i-1]]--; degre[ARETE[2*i-2]]--; Mark[ARETE[2*i-1]]=1; Mark[ARETE[2*i-2]]=1;
				Long[ARETE[2*i-1]]=LONGUEUR[i-1];
				if(bootStrap != NULL) boot[ARETE[2*i-1]] = (int) bootStrap[i-1];

			}
			else if ((degre[ARETE[2*i-2]]==1)&&(degre[ARETE[2*i-1]]==1)&&(Mark[ARETE[2*i-2]]==0)&&(Mark[ARETE[2*i-1]]==0))
			{
				Tree[ARETE[2*i-2]]=ARETE[2*i-1]; root=ARETE[2*i-1]; degre[ARETE[2*i-1]]--; degre[ARETE[2*i-2]]--; a=-1;
				Long[ARETE[2*i-2]]=LONGUEUR[i-1];
				if(bootStrap != NULL) boot[ARETE[2*i-2]] = (int) bootStrap[i-1];
			}
			if (a==-1) break;
		}
		if (a==-1) break;
	}


	/*  On decale et on complete la structure d'arbre avec Successeurs et Freres  */
	for (i=Ns+1;i>0;i--)
	{ 	Fre[i]=0; Suc[i]=0;
		//Tree[i]=Tree[i-1]+1; Long[i]=Long[i-1];
	}	Tree[root]=0;/*Tree[Ns+1]=0;*/

	for (i=1;i<=Ns+1/*Ns*/;i++)
	{
		if (i!=root)
		{
			sd=i; sf=Tree[i];
			if (Suc[sf]==0) Suc[sf]=sd;
			else {
				sf=Suc[sf];
				while (Fre[sf]>0) sf=Fre[sf];
				Fre[sf]=sd;
			}
		}
	}


	/* On compose la chaine parenthesee */
	string[0]=0; i=root;/*i=Ns+1;*/
	cpt=0;
	for (;;)
	{
		if(cpt > 1000000) exit(1);

		if (Suc[i]>0)
		{	sprintf(string,"%s(",string);
		Suc[i]=-Suc[i]; i=-Suc[i]; }
		else if (Fre[i]!=0)
		{	if (Suc[i]==0) sprintf(string,"%s%d:%.5f,",string,i,Long[i]);
			else {
				if(bootStrap != NULL)
					sprintf(string,"%s%d:%.5f,",string,boot[i],Long[i]);
				else
					sprintf(string,"%s:%.5f,",string,Long[i]);
			}
		i=Fre[i]; }
		else if (Tree[i]!=0)
		{	if (Suc[i]==0) sprintf(string,"%s%d:%.5f)",string,i,Long[i]);
		else {
			if(bootStrap != NULL)
				sprintf(string,"%s%d:%.5f)",string,boot[i],Long[i]);
			else
				sprintf(string,"%s:%.5f)",string,Long[i]);
			}
		i=Tree[i]; }
		else break;
	}
	strcat(string,";");

	FILE *pt_t = fopen(t,"w+");
	fprintf(pt_t,"%s",string);
	fclose(pt_t);

	free(Suc); free(Fre); free(Tree); free(Long); free(degre); free(Mark);	free(string);
}


void deleteBipartition(DescTree *DT,InputTree aTree){
	int i,j,k;
	//printf("\ndeleteBipartition, suppression des noeuds (size=%d - kt=%d):",aTree.size,aTree.kt);
	for(i=1;i<=2*aTree.size-2-aTree.kt;i++){
		//printf("%d ",i);
		if(i!=aTree.size){
			for(j=0;j<=DT[i].nbSommet+1;j++)
				free(DT[i].Matrice[j]);
			free(DT[i].Matrice);
		}
		free(DT[i].Tableau);
	}
	free(DT);
}

void PrintHeader(FILE *out,Parameters param){
	char criterion[50];
	if(strcmp(param.criterion,"rf")==0)
		strcpy(criterion,"Robinson and Foulds distance");
	if(strcmp(param.criterion,"ls")==0)
		strcpy(criterion,"least-squares");
	if(strcmp(param.criterion,"bd")==0)
		strcpy(criterion,"bipartition distance");
	if(strcmp(param.criterion,"qd")==0)
		strcpy(criterion,"quartet distance");

	fprintf(out,"%s",description);
	fprintf(out,"\n");
	//fprintf(out,"\nScenario             :%s",param.scenario);
	fprintf(out,"\nSubtree constraint   :%s",param.subtree);
	fprintf(out,"\nCriterion            :%s",criterion);
	fprintf(out,"\nBootstrap            :%s",param.bootstrap);
	fprintf(out,"\n");
}

//=================================================================
//==
//=================================================================
void initInputTree(struct InputTree *aTree){
	aTree->Adjacence = NULL;
	aTree->ARETE = NULL;
	aTree->LONGUEUR = NULL;
	aTree->ADD = NULL;
	aTree->Root = -1;
	aTree->size = -1;
	aTree->SpeciesName = NULL;
	aTree->Input = NULL;
	aTree->kt = 0;
}

//============================================
//=
//============================================
void printRoot(char *fichier, int R1, int R2){

	FILE *out;

	if((out=fopen(fichier,"w+"))==NULL){
		printf("\nCan't open root file (%s)",fichier);
		exit(-1);
	}
	fprintf(out,"%d %d",R1,R2);
	fclose(out);
}

void printRootByLeaves(char *fichier,int choix,struct InputTree *aTree){

	FILE *out;
	int i;
  printf("\nprintRootByLeaves (%s)",fichier);
	if((out=fopen(fichier,"w+"))==NULL){
		printf("\nCan't open root file (%s)",fichier);
		exit(-1);
	}

	for(i=1;i<=aTree->size;i++)
		if(aTree->ADD[i][aTree->ARETE[2*choix-1]] < aTree->ADD[i][aTree->ARETE[2*choix-2]])
			fprintf(out,"%s ",aTree->SpeciesName[i]);

	fprintf(out,"\n<>\n");

	for(i=1;i<=aTree->size;i++)
		if(aTree->ADD[i][aTree->ARETE[2*choix-1]] > aTree->ADD[i][aTree->ARETE[2*choix-2]])
			fprintf(out,"%s ",aTree->SpeciesName[i]);

	fclose(out);

}

//=====================================================
//
//=====================================================
void allocMemmory(struct InputTree *aTree, int n){

	int i;

	if(aTree->ARETE == NULL){

		aTree->degre = (int*)malloc(2*n*sizeof(int));
		aTree->ADD = (double**)malloc(2*n*sizeof(double*));
		aTree->Adjacence = (double**)malloc(2*n*sizeof(double*));
		aTree->Input = (double**)malloc(2*n*sizeof(double*));
		aTree->W = (double**)malloc((n+1)*sizeof(double*));

		for(i=0;i<2*n;i++){
			aTree->ADD[i] = (double*)malloc(2*n*sizeof(double));
			aTree->Adjacence[i] = (double*)malloc(2*n*sizeof(double));
			aTree->Input[i] = (double*)malloc(2*n*sizeof(double));
			if(i<=n)
				aTree->W[i] = (double*)malloc(2*n*sizeof(double));
		}

		aTree->ARETE    =(long int*)malloc(4*(2*(n))*sizeof(long int));
		aTree->LONGUEUR	=(double*)malloc((4*(n))*sizeof(double));
		aTree->SpeciesName = (char**)malloc(2*n*sizeof(char*));

		for(i=0;i<=n;i++)
			aTree->SpeciesName[i] = (char*)malloc(50);
	}

}

void freeReducedTree(struct InputTree *aTree,int n){
	int i,j;

	for(i=0;i<2*n;i++){
		free(aTree->ADD[i]);
		free(aTree->Adjacence[i]);
		free(aTree->Input[i]);
		free(aTree->W[i]);
		free(aTree->SpeciesName[i]);
	}
	free(aTree->ADD);
	free(aTree->Adjacence);
	free(aTree->Input);
	free(aTree->W);
	free(aTree->ARETE);
	free(aTree->SpeciesName);

}
//================================================================================
//==
//================================================================================

int TestSubTreeConstraint(struct InputTree aTree,int source, int dest,struct DescTree *DTSpecies, struct DescTree *DTGene){

	int basSource,basDest;
	int etape1=0,etape2=0;
	int LGTpossible = 0;
	int pv,i,y,z,tk,tl;
	int n = aTree.size;
	int mk1,mk2,nbSom;
	struct DescTree DTemp;
	double ** DISTemp =(double **)malloc((2*n+1)*sizeof(double*));

	int * PLACEk1=(int *) malloc((2*n-2)*sizeof(int));
	int * PLACEk2=(int *) malloc((2*n-2)*sizeof(int));

	int ** Bk1=(int **) malloc((2*n-2)*sizeof(int*));
	int ** Bk2=(int **) malloc((2*n-2)*sizeof(int*));

	for (i=0;i<2*n-2;i++)
	{
		Bk1[i]=(int *) malloc((n)*sizeof(int));
		Bk2[i]=(int *) malloc((n)*sizeof(int));
		DISTemp[i]=(double *)malloc((2*n+1)*sizeof(double));
	}

	//== find the node inf for the source branch
	if(aTree.ADD[aTree.ARETE[2*source-1]][aTree.Root] <  aTree.ADD[aTree.ARETE[2*source-2]][aTree.Root])
		basSource = aTree.ARETE[2*source-2];
	else
		basSource = aTree.ARETE[2*source-1];

	//== find the node inf for the dest branch
	if(aTree.ADD[aTree.ARETE[2*dest-1]][aTree.Root] <  aTree.ADD[aTree.ARETE[2*dest-2]][aTree.Root])
		basDest = aTree.ARETE[2*dest-2];
	else
		basDest = aTree.ARETE[2*dest-1];


	if(basSource<n)
		etape1=1;
	else{
		for(pv=3*(n-1+1)-2;pv<=3*(2*(n-1)-2);pv++){
			if(vecteursEgaux(DTSpecies[basSource],DTGene[pv]) == 1){
				/*if(DTSpecies[basSource].nbSommet == 3)
					etape1= 0;
				else*/{
					mk1=Bipartition_Table(DTGene[pv].Matrice,Bk1,PLACEk1,DTGene[pv].nbSommet/*+1*/);
					mk2=Bipartition_Table(DTSpecies[basSource].Matrice,Bk2,PLACEk2,DTSpecies[basSource].nbSommet/*+1*/);
					if(Table_Comparaison(Bk1,Bk2,PLACEk1,PLACEk2,mk1,mk2,DTSpecies[basSource].nbSommet/*+1*/) == 0)
						etape1 = 1;
				}
				pv = (int)INFINI;
			}
		}
	}

	if(etape1){/*/== comparer les sous arbres au niveau destination*/
		if(basDest<n)
			etape2 = 1;
		else{
			for(pv=3*(n-1+1)-2;pv<=3*(2*(n-1)-2);pv++){
				//printf("[(etape1) basDest=%d,pv = %d] ",basDest,pv);
				if(vecteursEgaux(DTSpecies[basDest],DTGene[pv]) == 1){
					/*if(DTSpecies[basDest].nbSommet == 3)
						etape2 = 0;
					else*/{
						mk1=Bipartition_Table(DTGene[pv].Matrice,Bk1,PLACEk1,DTGene[pv].nbSommet/*+1*/);
						mk2=Bipartition_Table(DTSpecies[basDest].Matrice,Bk2,PLACEk2,DTSpecies[basDest].nbSommet/*+1*/);
						if(Table_Comparaison(Bk1,Bk2,PLACEk1,PLACEk2,mk1,mk2,DTGene[pv].nbSommet/*+1*/) == 0)
							etape2 = 1;
					}
					pv = (int)INFINI;
				}
			}
		}
	}

	if(etape1==1 && etape2==1){

		nbSom = DTemp.nbSommet = DTSpecies[basDest].nbSommet + DTSpecies[basSource].nbSommet;

		DTemp.Tableau = (int*)malloc((DTemp.nbSommet+1)*sizeof(int));

		for(y=1;y<=DTSpecies[basDest].nbSommet;y++)
			DTemp.Tableau[y] = DTSpecies[basDest].Tableau[y];
		y--;
		for(z=1;z<=DTSpecies[basSource].nbSommet;z++)
			DTemp.Tableau[z+y] = DTSpecies[basSource].Tableau[z];

		TrierTableau(DTemp.Tableau,DTemp.nbSommet);

		for(pv=3*(n-1+1)-2;pv<=3*(2*(n-1)-2);pv++){
			if(vecteursEgaux(DTemp,DTGene[pv])){
				//printf("[(etape2) basSource=%d,pv = %d] ",basSource,pv);
				for(tk=1;tk<=DTGene[pv].nbSommet;tk++){
					for(tl=1;tl<=DTGene[pv].nbSommet;tl++)
						DISTemp[tk][tl] = aTree.ADD[DTemp.Tableau[tk]][DTemp.Tableau[tl]];
				}

				//for(tk=1;tk<=DTGene[pv].nbSommet;tk++)
				//	DISTemp[DTGene[pv].nbSommet+1][tk] =  DISTemp[tk][DTGene[pv].nbSommet+1] = ;

				mk1=Bipartition_Table(DTGene[pv].Matrice,Bk1,PLACEk1,DTGene[pv].nbSommet);
				mk2=Bipartition_Table(DISTemp,Bk2,PLACEk2,DTemp.nbSommet);
				if(Table_Comparaison(Bk1,Bk2,PLACEk1,PLACEk2,mk1,mk2,DTGene[pv].nbSommet) == 0)
					LGTpossible = 1;
			}
		}
		free(DTemp.Tableau);
	}

	free(PLACEk1);
	free(PLACEk2);

	for (i=0;i<2*n-2;i++)
	{
		free(Bk1[i]);
		free(Bk2[i]);
		free(DISTemp[i]);
	}
	free(DISTemp);
	free(Bk1);
	free(Bk2);

	//printf("\nLGTpossible=%d",LGTpossible);
	return LGTpossible;

}

//=================================================================
//==
//=================================================================
void loadCriteria(struct CRITERIA aCrit, struct HGT *aHGT){

	aHGT->crit.LS = aCrit.LS;
	aHGT->crit.RF = aCrit.RF;
	aHGT->crit.BD = aCrit.BD;
	aHGT->crit.QD = aCrit.QD;
	aHGT->valide = 1;
}

//=================================================================
//==
//=================================================================
void computeCriteria(double ** Matrix1, double ** Matrix2, int taille,struct CRITERIA *aCrit,double *L1, long int *A1,double *L2, long int *A2){

	int mI,m,i,j,RF,QD;
	double LS,BD=0;
	int size = taille-1;

	//= robinson and foulds
	m=Bipartition_Table(Matrix1,aCrit->B,aCrit->PLACE,size);
	mI=Bipartition_Table(Matrix2,aCrit->BI,aCrit->PLACEI,size);
	RF = Table_Comparaison(aCrit->B,aCrit->BI,aCrit->PLACE,aCrit->PLACEI,m,mI,size);

	//= least-squares
	LS = 0.0;
	for (i=1;i<=size-1;i++)
	{
		for (j=i+1;j<=size;j++){
			LS=LS + (Matrix1[i][j]-Matrix2[i][j])*(Matrix1[i][j]-Matrix2[i][j]);
			if(LS > INFINI){
					//printf("\n\n %lf,%lf (i=%d,j=%d,size=%d)",Matrix1[i][j],Matrix2[i][j],i,j,size);
					//getchar();
			}
		}
	}
	//LS = sqrt(LS);

	//= Bipartition Distance
	BD = BipartitionDistance(aCrit->B,aCrit->BI,size);

	//= Quartet Distance
	//= le calcul de QD necessite la creation de 2 fichiers d'arbres t1 et t2
	//= sinon la distance sera de -1 (distance non calculee)
/*	if(L1 != NULL){
		if(size > 3){
			SAVEASNewick(L1,A1,size,"t1");
			SAVEASNewick(L2,A2,size,"t2");
			system("./quartetDistance.sh");
			FILE *out = fopen("t3","r");
			fscanf(out,"%d",&QD);
			fclose(out);
		}
		else{
			QD=0;
		}
	}
	else{
		QD = (int)INFINI;
	}
*/
	aCrit->LS = LS;
	aCrit->BD = BD;
	aCrit->RF = RF;
	aCrit->QD = QD;

}

double MIN2 (double A,double B){
	if ( A > B) return B;
	return A;
}

//=================================================================
//==
//=================================================================
void applyHGT(double**ref,struct InputTree * aTree,int i,int j){

	int nodeToDel,otherNode,neighbor1=0,neighbor2=0,p,q,branch1=0,branch2=0,newNode;
	int nouveauNoeud;
	double tailleBrancheSource,tailleTransfert;
	int Ssup,Sinf,Sdestination;
	int kt = aTree->kt;
	int iternumber=100;

	if(i==0 || j==0) { printf("\ntransfert impossible"); exit(-1);}
	double ** DIST = (double**)malloc(2*aTree->size*(sizeof(double*)));
	for(p=0;p<2*aTree->size;p++)
		DIST[p] = (double *)malloc(2*aTree->size*(sizeof(double)));

  //= copie du contenue de ADD dans DIST
	for(p=1;p<=2*aTree->size-2-kt;p++)
		for(q=1;q<=2*aTree->size-2-kt;q++)
			DIST[p][q] = aTree->ADD[p][q];

	//= recherche des sommets superieurs et inferieurs sources
	tailleBrancheSource = aTree->LONGUEUR[i];
	if ( aTree->ADD[aTree->ARETE[2*i-1]][aTree->Root] < aTree->ADD[aTree->ARETE[2*i-2]][aTree->Root]){
			Ssup = aTree->ARETE[2*i-1];
			Sinf = aTree->ARETE[2*i-2];
	}
	else{
			Ssup = aTree->ARETE[2*i-2];
			Sinf = aTree->ARETE[2*i-1];
	}

	//= recherche des sommets superieurs et inferieurs destinations
	if(aTree->ADD[aTree->ARETE[2*j-1]][aTree->Root] < aTree->ADD[aTree->ARETE[2*j-2]][aTree->Root])
	{nodeToDel = aTree->ARETE[2*j-1]; otherNode = Sdestination = aTree->ARETE[2*j-2];}
	else
	{nodeToDel = aTree->ARETE[2*j-2]; otherNode = Sdestination = aTree->ARETE[2*j-1];}



	//= cas d'un noeud binaire, on va deplacer le noeud superieur (nodeToDel)
	if(aTree->degre[nodeToDel] == 3){
		//== 1. recherche des voisins du noeud a supprimer
		for(p=1;p<=2*(aTree->size)-3-aTree->kt;p++){
			if(aTree->ARETE[2*p-1] == nodeToDel && aTree->ARETE[2*p-2] != otherNode){
				if(neighbor1==0){neighbor1 = aTree->ARETE[2*p-2]; branch1=p;}
				else {neighbor2 = aTree->ARETE[2*p-2]; branch2=p;}
			}
			if(aTree->ARETE[2*p-2] == nodeToDel && aTree->ARETE[2*p-1] != otherNode){
				if(neighbor1==0){neighbor1 = aTree->ARETE[2*p-1]; branch1=p;}
				else {neighbor2 = aTree->ARETE[2*p-1]; branch2=p;}
			}
		}
		if(aTree->LONGUEUR[i-1] <= 2*5*epsilon){
		//printf("icit");
			for(p=1;p<=2*aTree->size-2-aTree->kt;p++)
				for(q=p+1;q<=2*aTree->size-2-aTree->kt;q++){
					if((fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-1]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-2]][q]) < 2*epsilon) ||
					   (fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-2]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-1]][q]) < 2*epsilon) )
					{
					  if((DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q]) <(DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q])){
						  aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q] + 2*5*epsilon;
						}
						else{
               aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q] + 2*5*epsilon;
            }
					}
				}
		}

		//== 3. branch1 is used to connect the two neighbors
		aTree->ARETE[2*branch1-1] = neighbor1;
		aTree->ARETE[2*branch1-2] = neighbor2;
		aTree->LONGUEUR[branch1-1]= aTree->LONGUEUR[branch1-1] + aTree->LONGUEUR[branch2-1];

		//== 4. branch2 is used to connect one of the source node (i)
		aTree->ARETE[2*branch2-1] = aTree->ARETE[2*i-1];
		aTree->ARETE[2*branch2-2] = nodeToDel;
		aTree->LONGUEUR[branch2-1]= (aTree->LONGUEUR[i-1]/2.0 > 5*epsilon)?aTree->LONGUEUR[i-1]/2.0:5*epsilon;

		//== 5. branch i is used to connect the other source node (i)
		aTree->ARETE[2*i-1] = aTree->ARETE[2*i-2];
		aTree->ARETE[2*i-2] = nodeToDel;
		aTree->LONGUEUR[i-1]= (aTree->LONGUEUR[i-1]/2.0 > 5*epsilon)?aTree->LONGUEUR[i-1]/2.0:5*epsilon;
		nouveauNoeud = nodeToDel;

		for(p=1;p<=2*aTree->size-2-aTree->kt;p++){

			if (fabs(DIST[p][aTree->Root] - DIST[p][Sdestination] - DIST[Sdestination][aTree->Root]) < 3*epsilon	){
			 //printf("\n==%d,%d,%d",p,Sdestination,nodeToDel);
      	aTree->ADD[nodeToDel][p] = aTree->ADD[p][nodeToDel] = 1.0 + DIST[Sdestination][p];
			}
			else	{
			  if(DIST[p][Ssup] < DIST[p][Sinf]){
				  aTree->ADD[nodeToDel][p] = aTree->ADD[p][nodeToDel] = DIST[p][Ssup] + aTree->LONGUEUR[i-1];
				}
				else{
          aTree->ADD[nodeToDel][p] = aTree->ADD[p][nodeToDel] = DIST[p][Sinf] + aTree->LONGUEUR[i-1];
        }
			}
		}

		aTree->ADD[nodeToDel][nodeToDel] = 0.0;
		tailleBrancheSource = aTree->LONGUEUR[i-1];
		tailleTransfert = 1.0;
		aTree->LONGUEUR[j-1]=1.0;
		nouveauNoeud=nodeToDel;
	}
	//= noeud de degre > a 3
	else{
		if(aTree->LONGUEUR[i-1] <= 2*5*epsilon){
			for(p=1;p<=2*aTree->size-2-aTree->kt;p++)
				for(q=p+1;q<=2*aTree->size-2-aTree->kt;q++){
					if((fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-1]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-2]][q]) < 2*epsilon) ||
					   (fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-2]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-1]][q]) < 2*epsilon) )
					{
					   if((DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q]) < (DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q])){
						    aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q] + 2*5*epsilon;
						 }
						 else{
                aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q] + 2*5*epsilon;
             }
					}
				}
		}

		//== 2. Ajouter un nouveau noeud sur la branche i (ajout d'une nouvelle arete)
		newNode = nouveauNoeud = 2*aTree->size-2-aTree->kt+1;

		double lt = (aTree->LONGUEUR[i-1]/2.0 >= 5*epsilon)?aTree->LONGUEUR[i-1]/2.0:5*epsilon;
		//== 4. connecter le nouveau noeud
		aTree->ARETE[2*j-1]=newNode;
		aTree->ARETE[2*j-2]=otherNode;
		aTree->LONGUEUR[j-1]= 1.0;

		//== 5. connecter la nouvelle arete
		aTree->ARETE[2*(2*aTree->size-3-aTree->kt+1)-1] = aTree->ARETE[2*i-1];
		aTree->ARETE[2*(2*aTree->size-3-aTree->kt+1)-2] = newNode;
		aTree->LONGUEUR[(2*aTree->size-3-aTree->kt+1)-1]= lt; //1.0;

		aTree->ARETE[2*i-1] = aTree->ARETE[2*i-2];
		aTree->ARETE[2*i-2] = newNode;
		aTree->LONGUEUR[i-1]= lt; //1.0;

		aTree->degre[newNode] = 3;
		aTree->degre[nodeToDel]--;
		aTree->kt--;

		for(p=1;p<=newNode;p++){
			if (fabs(DIST[p][aTree->Root] - DIST[p][Sdestination] - DIST[Sdestination][aTree->Root]) < 3*epsilon	){
				aTree->ADD[newNode][p] = aTree->ADD[p][newNode] = 1.0 + DIST[Sdestination][p];
			}
			else	{
				if((DIST[p][Ssup] + lt) < (DIST[p][Sinf] + lt)){
          aTree->ADD[newNode][p] = aTree->ADD[p][newNode] = DIST[p][Ssup] + lt;
        }
        else{
          aTree->ADD[newNode][p] = aTree->ADD[p][newNode] = DIST[p][Sinf] + lt;
        }
			}
		}



		aTree->ADD[newNode][newNode] = 0.0;
		tailleBrancheSource = lt;
		tailleTransfert=1.0;
	}

	//== Mise a jour des distances
	int tab_noeud1[2*aTree->size];
	int nbNoeud1=0;
	int tab_noeud2[2*aTree->size];
	int nbNoeud2=0;

  //== liste des noeuds sous sDestination

	for(int i=1;i<2*aTree->size-2-aTree->kt;i++){
      if( fabs(DIST[i][aTree->size]-DIST[i][Sdestination] - DIST[Sdestination][aTree->size]) < 5*epsilon ){
        tab_noeud1[nbNoeud1++] = i;
        //printf("\n1--%d",i);
      }
      else{
        tab_noeud2[nbNoeud2++] = i;
       // printf("\n2--%d",i);
      }
  }
    //== mise a jour des distances
  	for(int i=0;i<nbNoeud1;i++){
      for(int j=0;j<nbNoeud2;j++){
        aTree->ADD[tab_noeud1[i]][tab_noeud2[j]] =  aTree->ADD[tab_noeud2[j]][tab_noeud1[i]] = aTree->ADD[tab_noeud2[j]][nouveauNoeud] + tailleTransfert + aTree->ADD[Sdestination ][tab_noeud1[i]];
      }
    }


	//printf("\ndemi-branche = %lf",tailleBrancheSource);
	kt = aTree->kt;
	//printf("\navant load");

  //ComputeNewDistances(DIST,aTree->ADD,aTree->size,aTree->kt,nouveauNoeud,tailleBrancheSource,tailleTransfert,Ssup,Sinf,Sdestination);
  //approx_arb(ref,aTree->ADD,aTree->ADD,aTree->W,&iternumber,aTree->ARETE,aTree->LONGUEUR,1,&kt,0,aTree->size);
	loadAdjacenceMatrix(aTree->Adjacence,aTree->ARETE, aTree->LONGUEUR,aTree->size, aTree->kt);
	/*for(int i=1;i<=2*(aTree->size)-3-aTree->kt;i++){
		aTree->Adjacence[aTree->ARETE[2*i-1]][aTree->ARETE[2*i-2]] = aTree->Adjacence[aTree->ARETE[2*i-2]][aTree->ARETE[2*i-1]] = aTree->LONGUEUR[i-1];
	}*/
	//printf("\navant approx");
 //Floyd(aTree->Adjacence,aTree->ADD,aTree->size,aTree->kt);

/*
printf("\n\nFin de ApplyHGT : %d,kt=%d\n",aTree->size,aTree->kt);
	for(i=1;i<=2*aTree->size-2-aTree->kt;i++){
		printf("\n%d ",i);
		for(j=1;j<=2*aTree->size-2-aTree->kt;j++){
			printf("%.2lf ",aTree->ADD[i][j]);
		}
	}*/
	//printf("\nfin apply");


  //== Liste des noeuds sous





	for(p=0;p<2*aTree->size;p++)
		free(DIST[p]);
	free(DIST);
}

//=================================================================
//==
//=================================================================
void applyHGT2(double**ref,struct InputTree * aTree,int i,int j){

	int nodeToDel,otherNode,neighbor1=0,neighbor2=0,p,q,branch1=0,branch2=0,newNode;
	int nouveauNoeud;
	double tailleBrancheSource,tailleTransfert;
	int Ssup,Sinf,Sdestination;
	int kt = aTree->kt;
	int iternumber=100;

	if(i==0 || j==0) { printf("\ntransfert impossible"); exit(-1);}
	double ** DIST = (double**)malloc(2*aTree->size*(sizeof(double*)));
	for(p=0;p<2*aTree->size;p++)
		DIST[p] = (double *)malloc(2*aTree->size*(sizeof(double)));

  //= copie du contenue de ADD dans DIST
	for(p=1;p<=2*aTree->size-2-kt;p++)
		for(q=1;q<=2*aTree->size-2-kt;q++)
			DIST[p][q] = aTree->ADD[p][q];

	//= recherche des sommets superieurs et inferieurs sources
	tailleBrancheSource = aTree->LONGUEUR[i];
	if ( aTree->ADD[aTree->ARETE[2*i-1]][aTree->Root] < aTree->ADD[aTree->ARETE[2*i-2]][aTree->Root]){
			Ssup = aTree->ARETE[2*i-1];
			Sinf = aTree->ARETE[2*i-2];
	}
	else{
			Ssup = aTree->ARETE[2*i-2];
			Sinf = aTree->ARETE[2*i-1];
	}

	//= recherche des sommets superieurs et inferieurs destinations
	if(aTree->ADD[aTree->ARETE[2*j-1]][aTree->Root] < aTree->ADD[aTree->ARETE[2*j-2]][aTree->Root])
	{nodeToDel = aTree->ARETE[2*j-1]; otherNode = Sdestination = aTree->ARETE[2*j-2];}
	else
	{nodeToDel = aTree->ARETE[2*j-2]; otherNode = Sdestination = aTree->ARETE[2*j-1];}

	//= cas d'un noeud binaire, on va deplacer le noeud superieur (nodeToDel)
	if(aTree->degre[nodeToDel] == 3){
		//== 1. recherche des voisins du noeud a supprimer
		for(p=1;p<=2*(aTree->size)-3-aTree->kt;p++){
			if(aTree->ARETE[2*p-1] == nodeToDel && aTree->ARETE[2*p-2] != otherNode){
				if(neighbor1==0){neighbor1 = aTree->ARETE[2*p-2]; branch1=p;}
				else {neighbor2 = aTree->ARETE[2*p-2]; branch2=p;}
			}
			if(aTree->ARETE[2*p-2] == nodeToDel && aTree->ARETE[2*p-1] != otherNode){
				if(neighbor1==0){neighbor1 = aTree->ARETE[2*p-1]; branch1=p;}
				else {neighbor2 = aTree->ARETE[2*p-1]; branch2=p;}
			}
		}
		if(aTree->LONGUEUR[i-1] <= 2*5*epsilon){
		printf("icit");
			for(p=1;p<=2*aTree->size-2-aTree->kt;p++)
				for(q=p+1;q<=2*aTree->size-2-aTree->kt;q++){
					if((fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-1]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-2]][q]) < 2*epsilon) ||
					   (fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-2]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-1]][q]) < 2*epsilon) )
					{
					  if((DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q]) <(DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q])){
						  aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q] + 2*5*epsilon;
						}
						else{
               aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q] + 2*5*epsilon;
            }
					}
				}
		}

		//== 3. branch1 is used to connect the two neighbors
		aTree->ARETE[2*branch1-1] = neighbor1;
		aTree->ARETE[2*branch1-2] = neighbor2;
		aTree->LONGUEUR[branch1-1]= aTree->LONGUEUR[branch1-1] + aTree->LONGUEUR[branch2-1];

		//== 4. branch2 is used to connect one of the source node (i)
		aTree->ARETE[2*branch2-1] = aTree->ARETE[2*i-1];
		aTree->ARETE[2*branch2-2] = nodeToDel;
		aTree->LONGUEUR[branch2-1]= (aTree->LONGUEUR[i-1]/2.0 > 5*epsilon)?aTree->LONGUEUR[i-1]/2.0:5*epsilon;

		//== 5. branch i is used to connect the other source node (i)
		aTree->ARETE[2*i-1] = aTree->ARETE[2*i-2];
		aTree->ARETE[2*i-2] = nodeToDel;
		aTree->LONGUEUR[i-1]= (aTree->LONGUEUR[i-1]/2.0 > 5*epsilon)?aTree->LONGUEUR[i-1]/2.0:5*epsilon;
		nouveauNoeud = nodeToDel;

		for(p=1;p<=2*aTree->size-2-aTree->kt;p++){

			if (fabs(DIST[p][aTree->Root] - DIST[p][Sdestination] - DIST[Sdestination][aTree->Root]) < 3*epsilon	){
			 //printf("\n==%d,%d,%d",p,Sdestination,nodeToDel);
      	aTree->ADD[nodeToDel][p] = aTree->ADD[p][nodeToDel] = 1.0 + DIST[Sdestination][p];
			}
			else	{
			  if(DIST[p][Ssup] < DIST[p][Sinf]){
				  aTree->ADD[nodeToDel][p] = aTree->ADD[p][nodeToDel] = DIST[p][Ssup] + aTree->LONGUEUR[i-1];
				}
				else{
          aTree->ADD[nodeToDel][p] = aTree->ADD[p][nodeToDel] = DIST[p][Sinf] + aTree->LONGUEUR[i-1];
        }
			}
		}

		aTree->ADD[nodeToDel][nodeToDel] = 0.0;
		tailleBrancheSource = aTree->LONGUEUR[i-1];
		tailleTransfert = 1.0;
		aTree->LONGUEUR[j-1]=1.0;
		nouveauNoeud=nodeToDel;
	}
	//= noeud de degre > a 3
	else{
		if(aTree->LONGUEUR[i-1] <= 2*5*epsilon){
			for(p=1;p<=2*aTree->size-2-aTree->kt;p++)
				for(q=p+1;q<=2*aTree->size-2-aTree->kt;q++){
					if((fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-1]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-2]][q]) < 2*epsilon) ||
					   (fabs(DIST[p][q] - DIST[p][aTree->ARETE[2*i-2]] - aTree->LONGUEUR[i-1] - DIST[aTree->ARETE[2*i-1]][q]) < 2*epsilon) )
					{
					   if((DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q]) < (DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q])){
						    aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-1]]+DIST[aTree->ARETE[2*i-2]][q] + 2*5*epsilon;
						 }
						 else{
                aTree->ADD[p][q] = aTree->ADD[q][p] = DIST[p][q] = DIST[q][p] = DIST[p][aTree->ARETE[2*i-2]]+DIST[aTree->ARETE[2*i-1]][q] + 2*5*epsilon;
             }
					}
				}
		}

		//== 2. Ajouter un nouveau noeud sur la branche i (ajout d'une nouvelle arete)
		newNode = nouveauNoeud = 2*aTree->size-2-aTree->kt+1;

		double lt = (aTree->LONGUEUR[i-1]/2.0 >= 5*epsilon)?aTree->LONGUEUR[i-1]/2.0:5*epsilon;
		//== 4. connecter le nouveau noeud
		aTree->ARETE[2*j-1]=newNode;
		aTree->ARETE[2*j-2]=otherNode;
		aTree->LONGUEUR[j-1]= 1.0;

		//== 5. connecter la nouvelle arete
		aTree->ARETE[2*(2*aTree->size-3-aTree->kt+1)-1] = aTree->ARETE[2*i-1];
		aTree->ARETE[2*(2*aTree->size-3-aTree->kt+1)-2] = newNode;
		aTree->LONGUEUR[(2*aTree->size-3-aTree->kt+1)-1]= lt; //1.0;

		aTree->ARETE[2*i-1] = aTree->ARETE[2*i-2];
		aTree->ARETE[2*i-2] = newNode;
		aTree->LONGUEUR[i-1]= lt; //1.0;

		aTree->degre[newNode] = 3;
		aTree->degre[nodeToDel]--;
		aTree->kt--;

		for(p=1;p<=newNode;p++){
			if (fabs(DIST[p][aTree->Root] - DIST[p][Sdestination] - DIST[Sdestination][aTree->Root]) < 3*epsilon	){
				aTree->ADD[newNode][p] = aTree->ADD[p][newNode] = 1.0 + DIST[Sdestination][p];
			}
			else	{
				if((DIST[p][Ssup] + lt) < (DIST[p][Sinf] + lt)){
          aTree->ADD[newNode][p] = aTree->ADD[p][newNode] = DIST[p][Ssup] + lt;
        }
        else{
          aTree->ADD[newNode][p] = aTree->ADD[p][newNode] = DIST[p][Sinf] + lt;
        }
			}
		}

		aTree->ADD[newNode][newNode] = 0.0;
		tailleBrancheSource = lt;
		tailleTransfert=1.0;
	}

	//== Mise a jour des distances
	int tab_noeud1[2*aTree->size];
	int nbNoeud1=0;
	int tab_noeud2[2*aTree->size];
	int nbNoeud2=0;

  //== liste des noeuds sous sDestination

	for(int i=1;i<2*aTree->size-2-aTree->kt;i++){
      if( fabs(DIST[i][aTree->size]-DIST[i][Sdestination] - DIST[Sdestination][aTree->size]) < 5*epsilon ){
        tab_noeud1[nbNoeud1++] = i;
        //printf("\n1--%d",i);
      }
      else{
        tab_noeud2[nbNoeud2++] = i;
       // printf("\n2--%d",i);
      }
  }
    //== mise a jour des distances
  	for(int i=0;i<nbNoeud1;i++){
      for(int j=0;j<nbNoeud2;j++){
        aTree->ADD[tab_noeud1[i]][tab_noeud2[j]] =  aTree->ADD[tab_noeud2[j]][tab_noeud1[i]] = aTree->ADD[tab_noeud2[j]][nouveauNoeud] + tailleTransfert + aTree->ADD[Sdestination ][tab_noeud1[i]];
      }
    }


	//printf("\ndemi-branche = %lf",tailleBrancheSource);
	kt = aTree->kt;
	//printf("\navant load");

  //ComputeNewDistances(DIST,aTree->ADD,aTree->size,aTree->kt,nouveauNoeud,tailleBrancheSource,tailleTransfert,Ssup,Sinf,Sdestination);
  approx_arb(ref,aTree->ADD,aTree->ADD,aTree->W,&iternumber,aTree->ARETE,aTree->LONGUEUR,1,&kt,0,aTree->size);
  loadAdjacenceMatrix(aTree->Adjacence,aTree->ARETE, aTree->LONGUEUR,aTree->size, aTree->kt);
  /*
	for(int i=1;i<=2*(aTree->size)-3-aTree->kt;i++){
		aTree->Adjacence[aTree->ARETE[2*i-1]][aTree->ARETE[2*i-2]] = aTree->Adjacence[aTree->ARETE[2*i-2]][aTree->ARETE[2*i-1]] = aTree->LONGUEUR[i-1];
	}*/
	//printf("\navant approx");
  Floyd(aTree->Adjacence,aTree->ADD,aTree->size,aTree->kt);

/*
printf("\n\nFin de ApplyHGT : %d,kt=%d\n",aTree->size,aTree->kt);
	for(i=1;i<=2*aTree->size-2-aTree->kt;i++){
		printf("\n%d ",i);
		for(j=1;j<=2*aTree->size-2-aTree->kt;j++){
			printf("%.2lf ",aTree->ADD[i][j]);
		}
	}*/
	//printf("\nfin apply");


  //== Liste des noeuds sous





	for(p=0;p<2*aTree->size;p++)
		free(DIST[p]);
	free(DIST);
}
//=================================================================
//==
//=================================================================
void findListSpecies(struct HGT *bestHGT, struct DescTree *DTSpecies,struct InputTree aTree){

	int source,dest,i;
//	exit(0);

  //printf("\n%d-%d",bestHGT->source,aTree.Root);
	//== find the node inf for the source branch
	if(aTree.ADD[aTree.ARETE[2*(bestHGT->source)-1]][aTree.Root] <  aTree.ADD[aTree.ARETE[2*(bestHGT->source)-2]][aTree.Root])
		source = aTree.ARETE[2*(bestHGT->source)-2];
	else
		source = aTree.ARETE[2*(bestHGT->source)-1];

	//== find the node inf for the dest branch
	if(aTree.ADD[aTree.ARETE[2*(bestHGT->destination)-1]][aTree.Root] <  aTree.ADD[aTree.ARETE[2*(bestHGT->destination)-2]][aTree.Root])
		dest = aTree.ARETE[2*(bestHGT->destination)-2];
	else
		dest = aTree.ARETE[2*(bestHGT->destination)-1];

	if(bestHGT->listSource != NULL)
		free(bestHGT->listSource);
	bestHGT->listSource = (int *)malloc((DTSpecies[source].nbSommet+1)*sizeof(int));
	bestHGT->listSource[0] = DTSpecies[source].nbSommet;
	for(i=1;i<=DTSpecies[source].nbSommet;i++) bestHGT->listSource[i] = DTSpecies[source].Tableau[i];

	if(bestHGT->listDestination != NULL)
		free(bestHGT->listDestination);
	bestHGT->listDestination = (int *)malloc((DTSpecies[dest].nbSommet+1)*sizeof(int));
	bestHGT->listDestination[0] = DTSpecies[dest].nbSommet;
	for(i=1;i<=DTSpecies[dest].nbSommet;i++) bestHGT->listDestination[i] = DTSpecies[dest].Tableau[i];
}

//=================================================================
//==
//=================================================================
int TestCriterionAndUpdate(int *first,const char *CRITERION, struct CRITERIA aCrit, struct HGT * aHGT,int i,int j,int flag){

	int best = 0;

	if (strcmp(CRITERION,"rf")==0){
		if (aHGT->crit.RF > aCrit.RF) {
			best = 1;
			(*first)=0;
		}
		if (aHGT->crit.RF == aCrit.RF && aHGT->crit.LS > aCrit.LS && (*first) == 0)
			best = 1;
	}
	else if (strcmp(CRITERION,"ls")==0){
		if (aHGT->crit.LS > aCrit.LS) best = 1;
	}
	else if (strcmp(CRITERION,"bd")==0){
	  if (aHGT->crit.BD > aCrit.BD) best = 1;
    if (flag == 2){
        printf("LALA");
        best=1;
    }

	}
	else if (strcmp(CRITERION,"qd")==0){
		if (aHGT->crit.QD > aCrit.QD) best = 1;
	}
	else{
		printf("\nunknown criterion %s, try again...\n",CRITERION);
		exit(-1);
	}

	if(best == 1){
  	//printf("FIN TCAU");
		aHGT->destination = j;
		aHGT->source = i;
		aHGT->crit.BD = aCrit.BD;
		aHGT->crit.LS = aCrit.LS;
		aHGT->crit.RF = aCrit.RF;
		aHGT->crit.QD = aCrit.QD;
	}

	return best;
}

//=================================================================
//==
//=================================================================
void UpdateCriterion(int *first,const char *CRITERION, struct CRITERIA aCrit, struct HGT * aHGT,int i,int j){

		aHGT->destination = j;
		aHGT->source = i;
		aHGT->crit.BD = aCrit.BD;
		aHGT->crit.LS = aCrit.LS;
		aHGT->crit.RF = aCrit.RF;
		aHGT->crit.QD = aCrit.QD;
}

//=================================================================
//==
//=================================================================
int isAValidHGT(struct InputTree SpeciesTree,int i, int j){

	int d1,d2,d3,d4,root;

	/*if(SpeciesTree.ARETE[2*i-1]==2 && SpeciesTree.ARETE[2*i-2] ==12 && SpeciesTree.ARETE[2*j-1] == 11 && SpeciesTree.ARETE[2*j-2] == 13)
	{
		printEdges(0,SpeciesTree.ARETE,SpeciesTree.LONGUEUR,SpeciesTree.SpeciesName,SpeciesTree.size);
		printMatrix(SpeciesTree.SpeciesName,SpeciesTree.Matrix,2*SpeciesTree.size-2);
	}*/
	//== if is the root branch
	if(SpeciesTree.ARETE[2*i-1] == SpeciesTree.Root || SpeciesTree.ARETE[2*i-2] == SpeciesTree.Root ||
	   SpeciesTree.ARETE[2*j-1] == SpeciesTree.Root || SpeciesTree.ARETE[2*j-2] == SpeciesTree.Root )
		return 0;

	//== if there is a direct common ancestor
	if(((SpeciesTree.ARETE[2*i-1] == SpeciesTree.ARETE[2*j-1])&&(SpeciesTree.degre[SpeciesTree.ARETE[2*i-1]]==3)) ||
	   ((SpeciesTree.ARETE[2*i-1] == SpeciesTree.ARETE[2*j-2])&&(SpeciesTree.degre[SpeciesTree.ARETE[2*i-1]]==3)) ||
       ((SpeciesTree.ARETE[2*i-2] == SpeciesTree.ARETE[2*j-1])&&(SpeciesTree.degre[SpeciesTree.ARETE[2*i-2]]==3)) ||
	   ((SpeciesTree.ARETE[2*i-2] == SpeciesTree.ARETE[2*j-2])&&(SpeciesTree.degre[SpeciesTree.ARETE[2*i-2]]==3)) )
		return 0;

	//== if j is a branch connected to the root branch
	//if(SpeciesTree.ARETE[2*j-1] == 2*(SpeciesTree.size)-2 || SpeciesTree.ARETE[2*j-2] == 2*(SpeciesTree.size)-2)
	//	return 0;

	//== if i and j are on the same lineage

	root = SpeciesTree.Root;

	if(SpeciesTree.ADD[root][SpeciesTree.ARETE[2*i-1]] < SpeciesTree.ADD[root][SpeciesTree.ARETE[2*i-2]]){
		d1 = SpeciesTree.ARETE[2*i-1];
		d2 = SpeciesTree.ARETE[2*i-2];
	}
	else{
		d1 = SpeciesTree.ARETE[2*i-2];
		d2 = SpeciesTree.ARETE[2*i-1];
	}

	if(SpeciesTree.ADD[root][SpeciesTree.ARETE[2*j-1]] < SpeciesTree.ADD[root][SpeciesTree.ARETE[2*j-2]]){
		d3 = SpeciesTree.ARETE[2*j-1];
		d4 = SpeciesTree.ARETE[2*j-2];
	}
	else{
		d3 = SpeciesTree.ARETE[2*j-2];
		d4 = SpeciesTree.ARETE[2*j-1];
	}

	if(fabs(SpeciesTree.ADD[root][d4]-SpeciesTree.ADD[root][d1]-SpeciesTree.ADD[d1][d2]-SpeciesTree.ADD[d2][d3]-SpeciesTree.ADD[d3][d4]) < epsilon)
		return 0;
	if(fabs(SpeciesTree.ADD[root][d2]-SpeciesTree.ADD[root][d3]-SpeciesTree.ADD[d3][d4]-SpeciesTree.ADD[d4][d1]-SpeciesTree.ADD[d1][d2]) < epsilon)
		return 0;

	return 1;
}



//=================================================================
//==
//=================================================================
void copyInputTree(struct InputTree *tmpTree, struct InputTree aTree,int all,int allW){

	int n= aTree.size;
	int i,j;

	if(tmpTree->ADD == NULL){
		allocMemmory(tmpTree,n);
		/*tmpTree->SpeciesName = (char**)malloc(n*sizeof(char*));
		tmpTree->ARETE    =(long int*)malloc(4*(2*n-3)*sizeof(long int));
		tmpTree->LONGUEUR	=(double*)malloc(2*(2*n-3)*sizeof(double));
		tmpTree->Adjacence=(double**)malloc((2*n)*sizeof(double*));
		tmpTree->ADD = (double**)malloc((2*n)*sizeof(double*));
		tmpTree->Input = (double**)malloc((2*n)*sizeof(double*));
		tmpTree->W = (double**)malloc((2*n)*sizeof(double*));
		for(i=0;i<2*n;i++){
			tmpTree->Adjacence[i]=(double*)malloc((2*n)*sizeof(double));
			tmpTree->ADD[i] = (double*)malloc((2*n)*sizeof(double));
			if(all==1)
				tmpTree->Input[i] = (double*)malloc((2*n)*sizeof(double));
			tmpTree->W[i] = (double*)malloc((2*n)*sizeof(double));
		}*/
	}

	for(i=1;i<=2*n-1;i++)
		tmpTree->degre[i] = aTree.degre[i];

	for(i=1;i<=2*n-2;i++){
		for(j=1;j<=2*n-2;j++){
			tmpTree->Adjacence[i][j] = aTree.Adjacence[i][j];
			tmpTree->ADD[i][j] = aTree.ADD[i][j];
			if(all==1)
				tmpTree->Input[i][j] = aTree.Input[i][j];
		}
	}
  if(all==1){
  	 for(i=1;i<=2*n-2;i++){
		  for(j=1;j<=2*n-2;j++){
				tmpTree->Input[i][j] = aTree.Input[i][j];
		  }
	   }

    for(i=1;i<=n;i++){
			strcpy(tmpTree->SpeciesName[i],(const char*)aTree.SpeciesName[i]);
		}
	}

	if(allW==1){
    for(i=1;i<=n;i++){
		  for(j=1;j<=n;j++){
			 tmpTree->W[i][j] = aTree.W[i][j];
		  }
	 }
  }
	for(i=1;i<=2*n-3;i++){
		tmpTree->ARETE[2*i-1]  = aTree.ARETE[2*i-1];
		tmpTree->ARETE[2*i-2]  = aTree.ARETE[2*i-2];
		tmpTree->LONGUEUR[i-1] = aTree.LONGUEUR[i-1];
	}

	tmpTree->Root = aTree.Root;
	tmpTree->size = aTree.size;
	tmpTree->kt = aTree.kt;
}

void printLeaves(char *fichier,int nbHGT,struct HGT *outHGT,int noTree,char **NomsSpecies){

	FILE *out;
	int i,j;

	out = fopen(fichier,"a+");
/*	fprintf(out,"#detection %d",noTree);
	for(i=1;i<=nbHGT;i++){
		if(outHGT[i].valide == 1) {
			fprintf(out,"\n%d ",outHGT[i].listSource[0]);
			for(j=1;j<=outHGT[i].listSource[0];j++)
				fprintf(out,"%s ",NomsSpecies[outHGT[i].listSource[j]]);
			fprintf(out,"\n%d ",outHGT[i].listDestination[0]);
			for(j=1;j<=outHGT[i].listDestination[0];j++)
				fprintf(out,"%s ",NomsSpecies[outHGT[i].listDestination[j]]);
		}
	}
	fprintf(out,"\n");*/
	for(i=1;i<=nbHGT;i++){

		if(outHGT[i].valide == 1) {
			fprintf(out,"\ndetection %d : ",noTree);
			//fprintf(out,"\n%d ",outHGT[i].listSource[0]);
			for(j=1;j<=outHGT[i].listSource[0];j++){
				fprintf(out,"%s",NomsSpecies[outHGT[i].listSource[j]]);
				if(j<outHGT[i].listSource[0])
					fprintf(out,";");
				else
					fprintf(out," ");
			}
			//fprintf(out,"\n%d ",outHGT[i].listDestination[0]);
			for(j=1;j<=outHGT[i].listDestination[0];j++){
				fprintf(out,"%s",NomsSpecies[outHGT[i].listDestination[j]]);
				if(j<outHGT[i].listDestination[0])
					fprintf(out,";");
			}
		}
	}
	fprintf(out,"\n");
	fclose(out);
}
int sameHGT(struct HGT HGT1, struct HGT HGT2){
	int i;

	if(HGT1.listSource[0] == HGT2.listSource[0] && HGT1.listDestination[0] == HGT2.listDestination[0]){

		for(i=1;i<=HGT1.listSource[0];i++){
			if(HGT1.listSource[i] != HGT2.listSource[i])
				return 0;
		}
		for(i=1;i<=HGT1.listDestination[0];i++){
			if(HGT1.listDestination[i] != HGT2.listDestination[i])
				return 0;
		}
		return 1;
	}

	return 0;
}

int sameHGT2(struct HGT HGT1, struct HGT HGT2){
	int i,h1=0,h2=0;

	if(HGT1.listSource[0] == HGT2.listSource[0]){
		h1 = 1;
		for(i=1;i<=HGT1.listSource[0];i++){
			if(HGT1.listSource[i] != HGT2.listSource[i])
				h1=0;
		}
	}

	if(HGT1.listDestination[0] == HGT2.listDestination[0]){
		h2=1;
		for(i=1;i<=HGT1.listDestination[0];i++){
			if(HGT1.listDestination[i] != HGT2.listDestination[i])
				h2 = 0;
		}
	}

	if(h2 == 1 && h1 == 1) return 1;


	if(HGT1.listSource[0] == HGT2.listDestination[0]){
		h1=1;
		for(i=1;i<=HGT1.listSource[0];i++){
			if(HGT1.listSource[i] != HGT2.listDestination[i])
				h1=0;
		}
	}

	if(HGT1.listDestination[0] == HGT2.listSource[0]){
		h2=1;
		for(i=1;i<=HGT1.listDestination[0];i++){
			if(HGT1.listDestination[i] != HGT2.listSource[i])
				h2 = 0;
		}
	}

	if(h2 == 1 && h1 == 1) return 1;

	return 0;
}

void copyHGT(struct HGT source,struct HGT *dest){
	int j;
	dest->trivial = source.trivial;
	dest->valide = source.valide;
	dest->source_A = source.source_A;
	dest->source_B = source.source_B;
	dest->dest_A = source.dest_A;
	dest->dest_B = source.dest_B;
	dest->crit.LS = source.crit.LS;
	dest->crit.RF = source.crit.RF;
	dest->crit.BD = source.crit.BD;
	dest->crit.QD = source.crit.QD;
	dest->listSource = (int*) malloc((source.listSource[0]+1)*sizeof(int));
	dest->listDestination = (int*) malloc((source.listDestination[0]+1)*sizeof(int));
	dest->listSource[0] = source.listSource[0];
	for(j=1;j<=dest->listSource[0];j++){
		dest->listSource[j] = source.listSource[j];
	}
	dest->listDestination[0] = source.listDestination[0];
	for(j=1;j<=dest->listDestination[0];j++){
		dest->listDestination[j] = source.listDestination[j];
	}

}
//=============================================================================================================
//
//=============================================================================================================
void updateBootHGT(int first,struct HGT *bestHGT,int cpt_hgt,struct HGT *bootHGT, int *nbHGT_boot, double *tabBoot){

	int i,cpt=0,j;

	if(first==1){
		for(i=1;i<=cpt_hgt;i++){
			if(bestHGT[i].valide == 1){
				cpt++;
				tabBoot[cpt] = 1;
				copyHGT(bestHGT[i],&bootHGT[cpt]);
			}
		}
		(*nbHGT_boot) = cpt;
	}
	else{
		for(i=1;i<=(*nbHGT_boot);i++){
			for(j=1;j<=cpt_hgt;j++){
				if(bestHGT[j].valide == 1)
					if (sameHGT(bootHGT[i],bestHGT[j])==1)
						tabBoot[i]+=1;
			}
		}
	}



}
//=====================================================
//
//=====================================================
void sortHGT(struct HGT *tabHGT,int nbHGT,struct Parameters param){

	int i,j;
	int criterion;
	int a,b,c,d,rf,qd;
	double ls,bd;

	if(strcmp(param.criterion,"rf")==0)
		criterion = 1;
	else if(strcmp(param.criterion,"ls")==0)
		criterion = 2;
	else if(strcmp(param.criterion,"bd")==0)
		criterion = 3;
	else
		criterion = 4;

	for(i=1;i<=nbHGT;i++){
		for(j=1;j<=nbHGT-i;j++){
			if( ((criterion == 1) && (tabHGT[j].crit.RF > tabHGT[j+1].crit.RF)) ||
				((criterion == 2) && (tabHGT[j].crit.LS > tabHGT[j+1].crit.LS)) ||
				((criterion == 4) && (tabHGT[j].crit.QD > tabHGT[j+1].crit.QD)) ||
				((criterion == 3) && (tabHGT[j].crit.BD > tabHGT[j+1].crit.BD)) ){
				rf = tabHGT[j].crit.RF;
				ls = tabHGT[j].crit.LS;
				bd = tabHGT[j].crit.BD;
				qd = tabHGT[j].crit.QD;
				a = tabHGT[j].source_A;
				b = tabHGT[j].source_B;
				c = tabHGT[j].dest_A;
				d = tabHGT[j].dest_B;
				tabHGT[j].crit.RF = tabHGT[j+1].crit.RF;
				tabHGT[j].crit.LS = tabHGT[j+1].crit.LS;
				tabHGT[j].crit.BD = tabHGT[j+1].crit.BD;
				tabHGT[j].crit.QD = tabHGT[j+1].crit.QD;
				tabHGT[j].source_A = tabHGT[j+1].source_A;
				tabHGT[j].source_B = tabHGT[j+1].source_B;
				tabHGT[j].dest_A = tabHGT[j+1].dest_A;
				tabHGT[j].dest_B = tabHGT[j+1].dest_B;
				tabHGT[j+1].crit.RF = rf;
				tabHGT[j+1].crit.LS = ls;
				tabHGT[j+1].crit.BD = bd;
				tabHGT[j+1].crit.QD = qd;
				tabHGT[j+1].source_A = a;
				tabHGT[j+1].source_B = b;
				tabHGT[j+1].dest_A = c;
				tabHGT[j+1].dest_B = d;
			}
		}
	}
}
//=========================================================
//
//=========================================================
void FreeMemory_InputTreeReduced(struct InputTree *aTree,int size){

	int i;
	int n=size;

	if(aTree->ADD != NULL){

		free(aTree->LONGUEUR);
		free(aTree->ARETE);

		for(i=0;i<2*aTree->size-1;i++){
			free(aTree->ADD[i]);
			//free(aTree->Adjacence[i]);
			//free(aTree->Input[i]);
			free(aTree->W[i]);
		}
		for(i=0;i<2*n;i++)
			free(aTree->Adjacence[i]);

		if(aTree->SpeciesName != NULL){
			for(i=0;i<=n;i++){
				free(aTree->SpeciesName[i]);
			}
			free(aTree->SpeciesName);
		}

		free(aTree->ADD);
		free(aTree->Adjacence);
		//free(aTree->Input);
		free(aTree->W);
		free(aTree->degre);
	}
}

//=========================================================
//
//=========================================================
void FreeMemory_InputTree(struct InputTree *aTree,int size){

	int i;
	int n=size;

	if(aTree->ADD != NULL){

		free(aTree->LONGUEUR);
		free(aTree->ARETE);

		for(i=0;i<2*n;i++){
			free(aTree->ADD[i]);
			free(aTree->Input[i]);
			if(i<=n)
				free(aTree->W[i]);
		}
		for(i=0;i<2*n;i++){
			free(aTree->Adjacence[i]);
		}

		if(aTree->SpeciesName != NULL){
			for(i=0;i<=n;i++){
				free(aTree->SpeciesName[i]);
			}
			free(aTree->SpeciesName);
		}

		free(aTree->ADD);
		free(aTree->Adjacence);
		free(aTree->Input);
		free(aTree->W);
		free(aTree->degre);
	}
}

//==============================
//
//==============================
int nextTreeIsNewick(FILE *in){

	char c;

	while(fscanf(in,"%c",&c)!= EOF){
		if(c != ' ' && c !='\t' && c != '\n' && c !='\r'){
			fseek(in,-1,SEEK_CUR);
			if(c=='(') return 1;
			return 0;
		}
	}

	return -1;
}

//==============================
//
//==============================
char * readNewick(FILE *in){

	int cpt=0;
	char c;
	char * newick = (char*)malloc(100000);

	do{
		c=(char)fgetc(in);
		newick[cpt] = c;
		cpt++;
	}while(c!=';');

	newick[cpt] = '\0';

	return newick;
}

//=======================================
//
//=======================================
int nbSpeciesNewick(const char * newick){
	int i=0;
	int n = 0;
	char symbol,parent=' ';
	char symbolOld =' ';
	int temoin =0;

	do{
		symbol = newick[i];
		i++;
		//if (symbol == ':' && symbolOld != ')') n++;

		if(symbol == ':' && symbolOld !=')' && temoin != 1) n++;
		if(symbol >= '0' && symbol <= '9' && symbolOld==')') temoin=1;
		if(symbol==':' && temoin==1) temoin=0;
		symbolOld = symbol;
	}while(symbol != ';');

	return n;
}


//==================================================================
//
//==================================================================
void newickToMatrix(const char *newick,struct InputTree *aTree){

	int i,j;
	int pasBinaire=0;
	int pos_racine=-1;
//	double**DIST;


	aTree->size = nbSpeciesNewick(newick);
	allocMemmory(aTree,aTree->size+1);
	pos_racine = lectureNewick(newick,aTree->ARETE,aTree->LONGUEUR,aTree->SpeciesName,&aTree->kt);
	//printf("\nposition racine =%d",pos_racine);
	/*for(i=1;i<=2*aTree->size-3-aTree->kt;i++){
		printf("\n%d -> %d--%d : %lf",i,aTree->ARETE[2*i-2],aTree->ARETE[2*i-1],aTree->LONGUEUR[i-1]);
	}*/


 /* if(pos_racine != -1){
    addLeafAndUpdate(aTree,pos_racine);
	  aTree->Root = aTree->size;
  }
  else*/
  {
    loadAdjacenceMatrix(aTree->Adjacence,aTree->ARETE, aTree->LONGUEUR,aTree->size,aTree->kt);
    Floyd(aTree->Adjacence,aTree->ADD,aTree->Input,aTree->size,aTree->kt); //1ere fois
  }
  /*for(int i=1;i<=aTree->size;i++){
		printf("\n%s\t",aTree->SpeciesName[i]);
		for(int j=1;j<=aTree->size;j++){
			printf("%lf ",aTree->Input[i][j]);
		}
	}
	printf("\n");*/
}

//===============================================
//
//===============================================
void readMatrix(FILE *in,struct InputTree *aTree){

	int n,i,j;
	char name[50];
	double value;

	fscanf(in,"%d",&n);
	aTree->size = n;
	allocMemmory(aTree,aTree->size+1);

	for(i=1;i<=n;i++){
		fscanf(in,"%s",name); strcpy(aTree->SpeciesName[i],name);
		for(j=1;j<=n;j++){
			fscanf(in,"%lf",&value);
			aTree->ADD[i][j] = aTree->Input[i][j] = value;
		}
	}
}

//===================================================================================================================
//
//===================================================================================================================
int readInputFile(FILE *in, const char *tmpFile,bool addroot){

	char *newick;
	struct InputTree speciesTree_t;
	struct InputTree geneTree_t;
	int ret;

	initInputTree(&speciesTree_t);
	initInputTree(&geneTree_t);

	//== read the species tree
  ret = nextTreeIsNewick(in);
	if(ret==-1)
		return -1;
	else if(ret==1){
		newick = readNewick(in);
		//printf("newick=%s",newick);
		//ret = (int)strlen(newick);
		newickToMatrix(newick,&speciesTree_t);
		free(newick);
	}
	else
		readMatrix(in,&speciesTree_t);

	//== read the gene tree
	ret = nextTreeIsNewick(in);
	if(ret==-1)
		return -1;
	else if(ret==1){
		newick = readNewick(in);
    //printf("newick=%s",newick);
		newickToMatrix(newick,&geneTree_t);
		free(newick);
	}
	else
		readMatrix(in,&geneTree_t);

	if (!addroot){
    printf("not addroot");
    filtrerMatrice(speciesTree_t.Input,geneTree_t.Input,speciesTree_t.SpeciesName,geneTree_t.SpeciesName,speciesTree_t.size,geneTree_t.size);
  }
  //printf("%d--%d",speciesTree_t.size,geneTree_t.size);
	if(ecrireMatrice(speciesTree_t.Input,tmpFile,speciesTree_t.size,speciesTree_t.SpeciesName) == -1)
		return -2;
	ajouterMatriceGene(geneTree_t.Input,tmpFile,geneTree_t.size,geneTree_t.SpeciesName);
	/*
	TrierMatrices(geneTree_t.Input,geneTree_t.SpeciesName,speciesTree_t.SpeciesName,speciesTree_t.size);

  (*speciesTree) = speciesTree_t;
  (*geneTree)    = geneTree_t;
  */
	return 0;
}

//=================================================================================
//
//=================================================================================
int midPoint(long int *ARETE,double **DIST,int n,int kt){

	double max;
	int i,j,i1,j1,P;


	/*/== recherche la plus grande distance*/
	max = 0;
	for(i=1;i<=n;i++){
	//	printf("\n");
		for(j=1;j<=n;j++){
	//		printf("%lf ",DIST[i][j]);
			if(DIST[i][j] > max){
				max = DIST[i][j];
				i1=i; j1=j;
			}
		}
	}
	//	printf("\nnouvel arbre :  ");
	//	printf(" %d-%d, max =%lf,kt=%d",i1,j1,max,kt);
	P = -1;
	for(i=1;i<=2*n-3-kt;i++){
	//	printf("\n%ld--%ld : %lf",ARETE[2*i-1],ARETE[2*i-2],DIST[ARETE[2*i-1]][ARETE[2*i-2]]);
		if(ARETE[2*i-1] == i1 || ARETE[2*i-1] == j1 || ARETE[2*i-2] == i1 || ARETE[2*i-2] == j1 ){
			if(DIST[ARETE[2*i-1]][ARETE[2*i-2]] >= max/2.0){
				P = i;
			//	printf("icit meme ");
				break;
			}
		}
		else{
			if(ARETE[2*i-1] > n && ARETE[2*i-2]> n){
		       /* printf("\nmax/2.0 = %lf",max/2.0);
		        printf("\nDIST[%d][%ld] = %lf",i1,ARETE[2*i-1],DIST[i1][ARETE[2*i-1]]);
		        printf("\nDIST[%d][%ld] = %lf",j1,ARETE[2*i-2],DIST[j1][ARETE[2*i-2]]);
		        printf("\nDIST[%d][%ld] = %lf",i1,ARETE[2*i-2],DIST[i1][ARETE[2*i-2]]);
		        printf("\nDIST[%d][%ld] = %lf",j1,ARETE[2*i-1],DIST[j1][ARETE[2*i-1]]);
*/

			if( ( ((DIST[i1][ARETE[2*i-1]]>(max/2.0)) || (fabs(DIST[i1][ARETE[2*i-1]]-(max/2.0))<0.00001)) &&
				  ((DIST[j1][ARETE[2*i-2]]>(max/2.0)) || (fabs(DIST[j1][ARETE[2*i-2]]-(max/2.0))<0.00001))) ||
		        ( ((DIST[i1][ARETE[2*i-2]]>(max/2.0)) || (fabs(DIST[i1][ARETE[2*i-2]]-(max/2.0))<0.00001)) &&
				  ((DIST[j1][ARETE[2*i-1]]>(max/2.0)) || (fabs(DIST[j1][ARETE[2*i-1]]-(max/2.0))<0.00001)))
		     ){
					P = i;
				//	printf("la meme ");

					if( fabs(DIST[i1][j1] - DIST[i1][ARETE[2*i-1]] - DIST[ARETE[2*i-1]][j1]) < 0.00001 )
						break;
			}
			}
		}
	}
//	printf(" ,arete : %d--%d (%d) \n",ARETE[2*P-1],ARETE[2*P-2],P);

	return P;
}

//=================================================================
//==
//=================================================================
void AdjustBranchLength(struct InputTree *aTree1, struct InputTree aTree2,int binaire,int useFloyd){

	int i,j,iternumber=100;
	int kt = aTree1->kt;

	approx_arb(aTree2.ADD,aTree1->ADD,aTree1->ADD,aTree1->W,&iternumber,aTree1->ARETE,aTree1->LONGUEUR,1,&kt,binaire,aTree1->size);
  loadAdjacenceMatrix(aTree1->Adjacence,aTree1->ARETE, aTree1->LONGUEUR,aTree1->size,aTree1->kt);

	if(useFloyd == 1){
		Floyd(aTree1->Adjacence,aTree1->ADD,aTree1->size,aTree1->kt); // 3eme fois
	//	global_cpt3 ++;
  }

}

//=================================================================
//==
//=================================================================
void InitCriteria(struct CRITERIA * oldCrit, int size){

	int i;

	oldCrit->PLACE=(int *) malloc((2*size-3+1)*sizeof(int));
	oldCrit->PLACEI=(int *) malloc((2*size-3+1)*sizeof(int));
	oldCrit->B=(int **) malloc((2*size-3+1)*sizeof(int*));
	oldCrit->BI=(int **) malloc((2*size-3+1)*sizeof(int*));

	for (i=0;i<=2*size-3;i++)
	{
		oldCrit->B[i]=(int *) malloc((size+1)*sizeof(int));
		oldCrit->BI[i]=(int *) malloc((size+1)*sizeof(int));
	}

	oldCrit->BD=0;
	oldCrit->LS=0.0;
	oldCrit->RF=0;
}

void FreeCriteria(struct CRITERIA * Crit,int size){

	int i;

	free(Crit->PLACE);
	free(Crit->PLACEI);

	for(i=0;i<=2*size-3;i++){
		free(Crit->B[i]);
		free(Crit->BI[i]);
	}
	free(Crit->B);
	free(Crit->BI);
}

//====================================
//
//====================================

//================================================================================================================
//=
//================================================================================================================
int findAllMinimalScenario(struct InputTree SpeciesTree , struct InputTree GeneTree,int binaireSpecies,int binaireGene){
/*
	int encore=1;
	struct TreeHGT *aTree = (struct TreeHGT*)malloc(sizeof(struct TreeHGT));  //= sommet de l'arbre
	struct TreeHGT *courant;												  //= pointeur utilis� pour parcourrir l'arbre
	aTree->parent = NULL;													  //= le sommet n'a pas de parent

	struct DescTree * DTGene = (struct DescTree*)malloc((2*GeneTree.size-2-GeneTree.kt+1)*sizeof(struct DescTree));
	struct DescTree * DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
	struct ReduceTrace aMap;					//== mapping structure between species tree and recuded species tree
	struct InputTree SpeciesTreeRed;			//== reduced species tree
	struct InputTree GeneTreeRed;				//== reduced gene tree

	initInputTree(&SpeciesTreeRed);
	initInputTree(&GeneTreeRed);

	RechercherBipartition(GeneTree.ARETE,GeneTree.ADD,GeneTree.Root,GeneTree.Adjacence,DTGene,GeneTree.size,GeneTree.kt);
	RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);
	ReduceTree(SpeciesTree,GeneTree,&SpeciesTreeRed,&GeneTreeRed,&aMap,DTSpecies,DTGene,binaireSpecies,binaireGene);

	do{
		findBestHGTtab(SpeciesTreeRed,GeneTreeRed,param,bestHGTRed,&nbHgtFound);
		for(int i=0;i<nbHgtFound;i++){
			expandBestHGT(bestHGTRed[i],&bestHGT[cpt_hgt],aMap,DTSpecies,SpeciesTree);
			bestHGTRed[i].listSource = NULL;
			bestHGTRed[i].listDestination = NULL;
		}

	}while(encore);
*/
	return 0;
}

//================================================================================================================
//=
//================================================================================================================
int findAllHGT_no_criterion(struct InputTree SpeciesTree, struct InputTree GeneTree,struct Parameters param,struct HGT *tabHGT){

	int nbHGT=0;
	struct InputTree tmpTree;
	int i,j,ktSpecies;
	int size = SpeciesTree.size;
	struct CRITERIA aCrit;
	struct DescTree *DTSpecies,*DTGene;

	initInputTree(&tmpTree);

	//== compute criteria
	InitCriteria(&aCrit,size);

	DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size)*sizeof(struct DescTree));
	RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);

	for(i=1;i<2*size-3-SpeciesTree.kt;i++)
		for(j=i+1;j<2*size-3-SpeciesTree.kt;j++){

			//== is it a valid hgt ?
			if(isAValidHGT(SpeciesTree,i,j)==1){

				nbHGT++;
				tabHGT[nbHGT].source = i;
				tabHGT[nbHGT].destination = j;
				tabHGT[nbHGT].valide = 1;

			}
		}

	FreeCriteria(&aCrit,size);

	for(i=1;i<=nbHGT;i++)
		findListSpecies(&tabHGT[i],DTSpecies,SpeciesTree);


	//printf("\nfin de la fonction");

	return nbHGT;
}

//================================================================================================================
//=
//================================================================================================================
int findAllHGT(struct InputTree SpeciesTree, struct InputTree GeneTree,struct Parameters param,struct HGT *tabHGT){

	int nbHGT=0;
	struct InputTree tmpTree;
	int i,j,ktSpecies;
	int size = SpeciesTree.size;
	struct CRITERIA aCrit;
	struct DescTree *DTSpecies,*DTGene;

	initInputTree(&tmpTree);

	//== compute criteria
	InitCriteria(&aCrit,size);

	DTGene = (struct DescTree*)malloc((2*GeneTree.size)*sizeof(struct DescTree));
	RechercherBipartition(GeneTree.ARETE,GeneTree.ADD,GeneTree.Root,GeneTree.Adjacence,DTGene,GeneTree.size,GeneTree.kt);

	DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size)*sizeof(struct DescTree));
	RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);
	//printf("\n\n RF = %d \n LS = %lf",aCrit.RF,aCrit.LS);

	//printf("<br>size = %d<br>",size);
	//== Test all possible HGT
	//printEdges(0,SpeciesTree.ARETE,SpeciesTree.LONGUEUR,SpeciesTree.SpeciesName,SpeciesTree.size);


	for(i=1;i<2*size-3-SpeciesTree.kt;i++)
		for(j=1;j<2*size-3-SpeciesTree.kt;j++){

			//== is it a valid hgt ?
			if(isAValidHGT(SpeciesTree,i,j)==1){

			//	printf(".");

				copyInputTree(&tmpTree,SpeciesTree,0,0);

				//printf("\nTransfert : %d--%d -> %d--%d",tmpTree.ARETE[2*i-1],tmpTree.ARETE[2*i-2],tmpTree.ARETE[2*j-1],tmpTree.ARETE[2*j-2]);
				//getchar();

				if(strcmp(param.subtree,"yes") == 0)
					if(TestSubTreeConstraint(SpeciesTree,i,j,DTSpecies,DTGene) == 0) continue;

				//printf("\n\n");
				//printEdges(0,tmpTree.ARETE,tmpTree.LONGUEUR,SpeciesTree.SpeciesName,tmpTree.size);

				applyHGT(GeneTree.ADD,&tmpTree,i,j);
				//printf("\nTransfert : %d--%d -> %d--%d",tmpTree.ARETE[2*i-1],tmpTree.ARETE[2*i-2],tmpTree.ARETE[2*j-1],tmpTree.ARETE[2*j-2]);

				/*printEdges(0,tmpTree.ARETE,tmpTree.LONGUEUR,SpeciesTree.SpeciesName,tmpTree.size);
				printEdges(0,GeneTree.ARETE,GeneTree.LONGUEUR,SpeciesTree.SpeciesName,GeneTree.size);
				printMatrix(SpeciesTree.SpeciesName,tmpTree.ADD,tmpTree.size);
				printMatrix(SpeciesTree.SpeciesName,GeneTree.ADD,GeneTree.size);
				getchar();*/
				//printf("avant");
				AdjustBranchLength(&tmpTree,GeneTree,0,1);

				computeCriteria(tmpTree.ADD,GeneTree.ADD,size,&aCrit,tmpTree.LONGUEUR,tmpTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
				nbHGT++;
				loadCriteria(aCrit,&tabHGT[nbHGT]);
				tabHGT[nbHGT].source = i;
				tabHGT[nbHGT].destination = j;

				//printf("\nTransfert : %d--%d (%d)-> %d--%d (%d)",SpeciesTree.ARETE[2*i-1],SpeciesTree.ARETE[2*i-2],i,SpeciesTree.ARETE[2*j-1],SpeciesTree.ARETE[2*j-2],j);
				//printf("=====>RF = %d, LS = %lf",aCrit.RF,aCrit.LS);

			}
		}

	FreeCriteria(&aCrit,size);

	for(i=1;i<=nbHGT;i++)
		findListSpecies(&tabHGT[i],DTSpecies,SpeciesTree);

	return nbHGT;
}
//===============================================================================================================
//==
//===============================================================================================================
int findBestHGT(int initial,struct InputTree SpeciesTree,struct InputTree GeneTree,struct Parameters param,struct HGT *aHGT){

	struct InputTree tmpTree;
	int i,j,k,l,first=1,ret=0;
	int size = SpeciesTree.size;
	int ktSpecies;
	struct CRITERIA aCrit,aCritRef;
	struct DescTree *DTSpecies,*DTGene;

	initInputTree(&tmpTree);

	//printf("\ndebut findBestHgt");
	//== compute criteria
	InitCriteria(&aCrit,size);
	InitCriteria(&aCritRef,size);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,size,&aCritRef,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	loadCriteria(aCrit,aHGT);

	DTGene = (struct DescTree*)malloc((2*GeneTree.size-2-GeneTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartition(GeneTree.ARETE,GeneTree.ADD,GeneTree.Root,GeneTree.Adjacence,DTGene,GeneTree.size,GeneTree.kt);

	DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);
	//printf("\nRF = %d , LS = %lf , BD = %lf , QD = %d",aCrit.RF,aCrit.LS,aCrit.BD,aCrit.QD);

	/*for(i=1;i<=2*SpeciesTree.size-3-SpeciesTree.kt;i++)
		printf("\n%d-%d -->%lf",SpeciesTree.ARETE[2*i-1],SpeciesTree.ARETE[2*i-2],SpeciesTree.LONGUEUR[i-1]);
	printf("\n");*/

	//== Test all possible HGT
	for(i=1;i<2*size-3-SpeciesTree.kt;i++)
		for(j=1;j<2*size-3-SpeciesTree.kt;j++){
			//printf("\nicit : %d %d",i,j);
			//== is it a valid hgt ?
			if(isAValidHGT(SpeciesTree,i,j)==1 && i!=j ){

				//printf("\n%d-%d -> %d-%d",SpeciesTree.ARETE[2*i-1],SpeciesTree.ARETE[2*i-2],SpeciesTree.ARETE[2*j-1],SpeciesTree.ARETE[2*j-2]);
				//printf(".");

				copyInputTree(&tmpTree,SpeciesTree,0,0);

				//printf("\nTransfert : %d--%d -> %d--%d",tmpTree.ARETE[2*i-1],tmpTree.ARETE[2*i-2],tmpTree.ARETE[2*j-1],tmpTree.ARETE[2*j-2]);

				//if(tmpTree.ARETE[2*i-1]==6 && tmpTree.ARETE[2*i-2]==15 && tmpTree.ARETE[2*j-1]==7 && tmpTree.ARETE[2*j-2]==14)
				//	printf("");

				if(strcmp(param.subtree,"yes") == 0)
					if(TestSubTreeConstraint(SpeciesTree,i,j,DTSpecies,DTGene) == 0) {
						continue;
						//FreeMemory_InputTree(&tmpTree);
						tmpTree.ADD = NULL;
						tmpTree.ARETE = NULL;

					}

				//printf("\n\n");
				//printEdges(0,tmpTree.ARETE,tmpTree.LONGUEUR,SpeciesTree.SpeciesName,tmpTree.size);
				//printf("\navant apply");
				applyHGT(GeneTree.ADD,&tmpTree,i,j);
				//printf("\napres apply");

				//printEdges(0,tmpTree.ARETE,tmpTree.LONGUEUR,SpeciesTree.SpeciesName,tmpTree.size);
				//printEdges(0,GeneTree.ARETE,GeneTree.LONGUEUR,SpeciesTree.SpeciesName,GeneTree.size);
				/*for(k=1;k<=2*tmpTree.size-2-tmpTree.kt;k++){
					printf("\n");
					for(l=1;l<=2*tmpTree.size-2-tmpTree.kt;l++)ListeSommets_taille_0(SpeciesTree.Input,tab_sommet);
						printf("%lf ",tmpTree.ADD[k][l]);
				}
				for(k=1;k<=2*tmpTree.size-3-tmpTree.kt;k++)
					printf("\n%d-%d -->%lf",tmpTree.ARETE[2*k-1],tmpTree.ARETE[2*k-2],tmpTree.LONGUEUR[k-1]);
				printf("\n");*/
				//printMatrix(SpeciesTree.SpeciesName,GeneTree.ADD,GeneTree.size);

				AdjustBranchLength(&tmpTree,GeneTree,0,1);
				//printMatrix(SpeciesTree.SpeciesName,tmpTree.ADD,tmpTree.size);


			//	computeCriteria(tmpTree.ADD,GeneTree.ADD,size,&aCrit);
			//	printf("<br>Transfert : %d--%d (%d)-> %d--%d (%d)",SpeciesTree.ARETE[2*i-1],SpeciesTree.ARETE[2*i-2],i,SpeciesTree.ARETE[2*j-1],SpeciesTree.ARETE[2*j-2],j);
			//	printf("=====>RF = %d, LS = %lf",aCrit.RF,aCrit.LS);

			//	if(i==12 && j==7){
			//		aCrit.RF = 100;
			//	}
				if(initial == 1){
					computeCriteria(tmpTree.ADD,SpeciesTree.ADD,size,&aCrit,tmpTree.LONGUEUR,tmpTree.ARETE,SpeciesTree.LONGUEUR,SpeciesTree.ARETE);
					if(aCrit.RF == 1){
						ret += TestCriterionAndUpdate(&first,param.criterion,aCrit,aHGT,i,j,0);
						i=j=(int)INFINI;
						break;
					}
				}
				else{
					computeCriteria(tmpTree.ADD,GeneTree.ADD,size,&aCrit,tmpTree.LONGUEUR,tmpTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
			//		printf("<br>aCrit : RF = %d, LS = %lf",aCrit.RF,aCrit.LS);
			//		printf("<br>aHGT : RF = %d, LS = %lf",aHGT->crit.RF,aHGT->crit.LS);
		//			if(i==13 && j==3){
						//printf("Le gros changement");
		//				aCrit.RF=7;
		//			}
					ret += TestCriterionAndUpdate(&first,param.criterion,aCrit,aHGT,i,j,0);
					if(aHGT->crit.RF == 0 /*|| aHGT->crit.LS < epsilon */|| aHGT->crit.BD == 0 || aHGT->crit.QD == 0) { i=j=(int)INFINI;}
				}


				//FreeMemory_InputTree(&tmpTree);
				//tmpTree.ADD = NULL;
				//tmpTree.ARETE = NULL;

			}
			//printf("apres");
		}
	//printf("\n===> %d-%d -> %d-%d , RF = %d , LS = %lf QD = %d",SpeciesTree.ARETE[2*(aHGT->source)-1],SpeciesTree.ARETE[2*(aHGT->source)-2],SpeciesTree.ARETE[2*(aHGT->destination)-1],SpeciesTree.ARETE[2*(aHGT->destination)-2],aHGT->crit.RF,aHGT->crit.LS,aHGT->crit.QD);
	//printf("icit");
	if(ret > 0 )
		findListSpecies(aHGT,DTSpecies,SpeciesTree);

	deleteBipartition(DTSpecies,SpeciesTree);
	deleteBipartition(DTGene,GeneTree);
	FreeCriteria(&aCrit,size);
	FreeCriteria(&aCritRef,size);
	FreeMemory_InputTree(&tmpTree,tmpTree.size);

	//printf("\nfin findBestHgt");
	return ret;
}

//===========================================================================================================================
//==
//===========================================================================================================================
int findBestHGT_nombreLimite(struct DescTree *DTSpecies,struct DescTree *DTGene,int * tab_branches,int nb_branches,struct InputTree SpeciesTree,struct InputTree GeneTree,struct Parameters param,struct HGT *aHGT){

	struct InputTree tmpTree;
	int i,j,k,l,first=1,ret=0;
	int size = SpeciesTree.size;
	int ktSpecies;
	struct CRITERIA aCrit,aCritRef;

	initInputTree(&tmpTree);

	//printf("\ndebut findBestHgt");
	//== compute criteria
	InitCriteria(&aCrit,size);
	InitCriteria(&aCritRef,size);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,size,&aCritRef,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	loadCriteria(aCrit,aHGT);

	/*DTGene = (struct DescTree*)malloc((2*GeneTree.size-2-GeneTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartition(GeneTree.ARETE,GeneTree.ADD,GeneTree.Root,GeneTree.Adjacence,DTGene,GeneTree.size,GeneTree.kt);

	DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);*/
	//printf("\nRF = %d , LS = %lf , BD = %lf , QD = %d",aCrit.RF,aCrit.LS,aCrit.BD,aCrit.QD);

	//for(i=1;i<=2*SpeciesTree.size-3-SpeciesTree.kt;i++)
	//	printf("\n%d-%d -->%lf",SpeciesTree.ARETE[2*i-1],SpeciesTree.ARETE[2*i-2],SpeciesTree.LONGUEUR[i-1]);
	//printf("\n");

	//== Test all possible HGT
	for(i=1;i<=nb_branches;i++){
		for(j=1;j<=nb_branches;j++){
			//printf("\navant la validation : %d--%d -> %d--%d ",SpeciesTree.ARETE[2*tab_branches[i]-1],SpeciesTree.ARETE[2*tab_branches[i]-2],SpeciesTree.ARETE[2*tab_branches[j]-1],SpeciesTree.ARETE[2*tab_branches[j]-2]);
			if(isAValidHGT(SpeciesTree,tab_branches[i],tab_branches[j])==1 && i!=j ){
				//printf("\nc'est valide");
			//if(i!=j ){
				//printf("\n%d-%d -->%d--%d",SpeciesTree.ARETE[2*tab_branches[i]-1],SpeciesTree.ARETE[2*tab_branches[i]-2],SpeciesTree.ARETE[2*tab_branches[j]-1],SpeciesTree.ARETE[2*tab_branches[j]-2]);
				//printf("\nPOS 1");
				copyInputTree(&tmpTree,SpeciesTree,0,0);
				//printf("\nPOS 2");
				/*if(strcmp(param.subtree,"yes") == 0)
					if(TestSubTreeConstraint(SpeciesTree,tab_branches[i],tab_branches[j],DTSpecies,DTGene) == 0) {
						continue;
						//FreeMemory_InputTree(&tmpTree);
						tmpTree.ADD = NULL;
						tmpTree.ARETE = NULL;

					}*/
				//printf("\nPOS 3");
				applyHGT(GeneTree.ADD,&tmpTree,tab_branches[i],tab_branches[j]);
				//printf("\nPOS 4");
				AdjustBranchLength(&tmpTree,GeneTree,0,1);
				//printf("\nPOS 5");
				computeCriteria(tmpTree.ADD,GeneTree.ADD,size,&aCrit,tmpTree.LONGUEUR,tmpTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
				//printf("\nPOS 6");
				ret += TestCriterionAndUpdate(&first,param.criterion,aCrit,aHGT,tab_branches[i],tab_branches[j],0);
				//printf("\nRF = %d , LS = %lf , BD = %lf , QD = %d",aHGT->crit.RF,aHGT->crit.LS,aHGT->crit.BD,aHGT->crit.QD);
				if(aHGT->crit.RF == 0 /*|| aHGT->crit.LS < epsilon */|| aHGT->crit.BD == 0 || aHGT->crit.QD == 0) { i=j=(int)INFINI;}

			}
		}

	}
	/*if(ret > 0 )
		findListSpecies(aHGT,DTSpecies,SpeciesTree);*/

/*	deleteBipartition(DTSpecies,SpeciesTree);
	deleteBipartition(DTGene,GeneTree);*/
	FreeCriteria(&aCrit,size);
	FreeCriteria(&aCritRef,size);
	FreeMemory_InputTree(&tmpTree,tmpTree.size);

	printf("\nfin findBestHgt (%d)",ret);
	return ret;
}


//===============================================================================================================
//==
//===============================================================================================================
int findBestHGTtab(struct InputTree SpeciesTree,struct InputTree GeneTree,struct Parameters param,struct HGT *aHGT, int *nbHgtFound,int *initial, int *listRef,int * listJ){

	struct InputTree tmpTree;
	int i,j,k,l,first=1,ret=0,trouve;
	int size = SpeciesTree.size;
	int ktSpecies;
	struct CRITERIA aCrit,aCritRef,aCritRef2;
	struct DescTree *DTSpecies,*DTGene;
  int encore = 0;
  struct HGT *tmpHGT;

	initInputTree(&tmpTree);
	(*nbHgtFound) = 0;

	printf("\n== NOUVELLE RECHERCHE == [size=%d]\n",SpeciesTree.size);

	//== compute criteria
	InitCriteria(&aCrit,size);
	InitCriteria(&aCritRef,size);
	InitCriteria(&aCritRef2,size);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,size,&aCritRef,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	computeCriteria(SpeciesTree.ADD,GeneTree.ADD,size,&aCritRef2,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
	loadCriteria(aCrit,&aHGT[0]);

	DTGene = (struct DescTree*)malloc(3*(2*GeneTree.size-2-GeneTree.kt+1)*sizeof(struct DescTree));
	//RechercherBipartition(GeneTree.ARETE,GeneTree.ADD,GeneTree.Root,GeneTree.Adjacence,DTGene,GeneTree.size,GeneTree.kt);
	RechercherBipartitionSansRacine(GeneTree.ARETE,GeneTree.ADD,GeneTree.Adjacence,DTGene,GeneTree.size,GeneTree.kt);

	DTSpecies = (struct DescTree*)malloc((2*SpeciesTree.size-2-SpeciesTree.kt+1)*sizeof(struct DescTree));
	RechercherBipartition(SpeciesTree.ARETE,SpeciesTree.ADD,SpeciesTree.Root,SpeciesTree.Adjacence,DTSpecies,SpeciesTree.size,SpeciesTree.kt);
	//printf("\n\n RF = %d \n LS = %lf",aCrit.RF,aCrit.LS);
	//printf("\n%d-%d",GeneTree.size,SpeciesTree.size);
	/*for(i=1;i<=2*SpeciesTree.size-3-SpeciesTree.kt;i++)
		printf("\n%d-%d -->%lf",SpeciesTree.ARETE[2*i-1],SpeciesTree.ARETE[2*i-2],SpeciesTree.LONGUEUR[i-1]);
	printf("\n");*/
	//== Test all possible HGT
/*
	printf("\n\ndebut detection\n");
	for(i=1;i<=SpeciesTree.size;i++){
		printf("\n%d ",i);
		for(j=1;j<=SpeciesTree.size;j++){
			printf("%.2lf ",SpeciesTree.ADD[i][j]);
		}
	}

  printf("\n\ndebut detection\n");
	for(i=1;i<=2*SpeciesTree.size-3-SpeciesTree.kt;i++){
		printf("\n%d ",i);
		for(j=1;j<=2*SpeciesTree.size-3-SpeciesTree.kt;j++){
			printf("%.2lf ",SpeciesTree.ADD[i][j]);
		}
	}*/

	int flag2 = 1;
	int interdit1,interdit2;
	listRef[0]=0;

	do{
	for(i=1;i<2*size-3-SpeciesTree.kt;i++)
		  for(j=i+1;j<2*size-3-SpeciesTree.kt;j++){
			//== is it a valid hgt ?
			trouve=0;
			interdit1 = 0;
			interdit2 = 0;
			for(int k=1;k<=listJ[0];k++){
				if(listJ[k] == j) interdit1 = 1;
				if(listJ[k] == i) interdit2 = 1;
			}

			if(isAValidHGT(SpeciesTree,i,j)==1 && (i!=j) && (interdit1==0)){

				copyInputTree(&tmpTree,SpeciesTree,0,0);

				if(strcmp(param.subtree,"yes") == 0)
					if(TestSubTreeConstraint(SpeciesTree,i,j,DTSpecies,DTGene) == 0) {
						continue;
						tmpTree.ADD = NULL;
						tmpTree.ARETE = NULL;
					}
				  applyHGT(GeneTree.ADD,&tmpTree,i,j);
				  if(strcmp(param.criterion,"ls") == 0)
				    AdjustBranchLength(&tmpTree,GeneTree,0,1);
				  computeCriteria(tmpTree.ADD,GeneTree.ADD,size,&aCrit,tmpTree.LONGUEUR,tmpTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);
				  first=1;
				  loadCriteria(aCritRef,&aHGT[(*nbHgtFound)]);

				if(((aCritRef.RF-aCrit.RF) == 1) && (*initial==1) && (
				  (SpeciesTree.ARETE[2*i-1] == SpeciesTree.ARETE[2*j-1]) ||
				  (SpeciesTree.ARETE[2*i-1] == SpeciesTree.ARETE[2*j-2]) ||
				  (SpeciesTree.ARETE[2*i-2] == SpeciesTree.ARETE[2*j-2]) ||
				  (SpeciesTree.ARETE[2*i-2] == SpeciesTree.ARETE[2*j-2])
          )) {

				  UpdateCriterion(&first,param.criterion,aCrit,&aHGT[(*nbHgtFound)],i,j);

					aHGT[(*nbHgtFound)].source_A = SpeciesTree.ARETE[2*i-1];
					aHGT[(*nbHgtFound)].source_B = SpeciesTree.ARETE[2*i-2];
					aHGT[(*nbHgtFound)].dest_A = SpeciesTree.ARETE[2*j-1];
					aHGT[(*nbHgtFound)].dest_B = SpeciesTree.ARETE[2*j-2];

					findListSpecies(&aHGT[(*nbHgtFound)],DTSpecies,SpeciesTree);
					listRef[listRef[0]+1] = j;
					trouve=1;
				}
				else  if(*initial==0){
          if(TestCriterionAndUpdate(&first,param.criterion,aCrit,&aHGT[(*nbHgtFound)],i,j,0) == 1){
    					aHGT[(*nbHgtFound)].source_A = SpeciesTree.ARETE[2*i-1];
    					aHGT[(*nbHgtFound)].source_B = SpeciesTree.ARETE[2*i-2];
    					aHGT[(*nbHgtFound)].dest_A = SpeciesTree.ARETE[2*j-1];
    					aHGT[(*nbHgtFound)].dest_B = SpeciesTree.ARETE[2*j-2];
    					findListSpecies(&aHGT[(*nbHgtFound)],DTSpecies,SpeciesTree);
						  listRef[listRef[0]+1] = j;
    					trouve=1;

          }
        }
      }
			if(isAValidHGT(SpeciesTree,j,i)==1 && (i!=j) && (interdit2 == 0) ){
				copyInputTree(&tmpTree,SpeciesTree,0,0);

				if(strcmp(param.subtree,"yes") == 0)
					if(TestSubTreeConstraint(SpeciesTree,j,i,DTSpecies,DTGene) == 0) {
						continue;
						tmpTree.ADD = NULL;
						tmpTree.ARETE = NULL;
					}
				  applyHGT(GeneTree.ADD,&tmpTree,j,i);
				  if(strcmp(param.criterion,"ls") == 0)
				    AdjustBranchLength(&tmpTree,GeneTree,0,1);
				  computeCriteria(tmpTree.ADD,GeneTree.ADD,size,&aCrit,tmpTree.LONGUEUR,tmpTree.ARETE,GeneTree.LONGUEUR,GeneTree.ARETE);

			  	int flag=0;
				  first=1;
				  if(trouve==0)
					 loadCriteria(aCritRef,&aHGT[(*nbHgtFound)]);
				  else{
				    aHGT[(*nbHgtFound)].crit.diff_bd = fabs(aCrit.BD - aHGT[(*nbHgtFound)].crit.BD);
            if((flag2 == 1) && (fabs(aCrit.BD - aHGT[(*nbHgtFound)].crit.BD) <= 1)){
              if( (aHGT[(*nbHgtFound)].crit.RF - aCrit.RF) >= 2 ){
                  trouve = 0;
                  flag = 2;
               }
               else{
                  flag = 1;
               }
            }
            else{
				       printf("\n[%2d--%2d -> %2d--%2d]",aHGT[(*nbHgtFound)].source_A,aHGT[(*nbHgtFound)].source_B,aHGT[(*nbHgtFound)].dest_A,aHGT[(*nbHgtFound)].dest_B);
			         printf("\tRF1=%d | RF2=%d | BD1=%2.1lf | BD2=%2.1lf",aHGT[(*nbHgtFound)].crit.RF,aCrit.RF,aHGT[(*nbHgtFound)].crit.BD,aCrit.BD);
          }
        }
        if(flag != 1){
				if(((aCritRef.RF-aCrit.RF) == 1) && (*initial==1) && (
				  (SpeciesTree.ARETE[2*i-1] == SpeciesTree.ARETE[2*j-1]) ||
				  (SpeciesTree.ARETE[2*i-1] == SpeciesTree.ARETE[2*j-2]) ||
				  (SpeciesTree.ARETE[2*i-2] == SpeciesTree.ARETE[2*j-2]) ||
				  (SpeciesTree.ARETE[2*i-2] == SpeciesTree.ARETE[2*j-2])
          )) {

				   UpdateCriterion(&first,param.criterion,aCrit,&aHGT[(*nbHgtFound)],j,i);
					aHGT[(*nbHgtFound)].source_A = SpeciesTree.ARETE[2*j-1];
					aHGT[(*nbHgtFound)].source_B = SpeciesTree.ARETE[2*j-2];
					aHGT[(*nbHgtFound)].dest_A = SpeciesTree.ARETE[2*i-1];
					aHGT[(*nbHgtFound)].dest_B = SpeciesTree.ARETE[2*i-2];
					findListSpecies(&aHGT[(*nbHgtFound)],DTSpecies,SpeciesTree);
					listRef[listRef[0]+1] = i;
					trouve=1;
				}
				else if(*initial==0){
              if(TestCriterionAndUpdate(&first,param.criterion,aCrit,&aHGT[(*nbHgtFound)],j,i,flag) == 1){
    					aHGT[(*nbHgtFound)].source_A = SpeciesTree.ARETE[2*j-1];
    					aHGT[(*nbHgtFound)].source_B = SpeciesTree.ARETE[2*j-2];
    					aHGT[(*nbHgtFound)].dest_A = SpeciesTree.ARETE[2*i-1];
    					aHGT[(*nbHgtFound)].dest_B = SpeciesTree.ARETE[2*i-2];
    					findListSpecies(&aHGT[(*nbHgtFound)],DTSpecies,SpeciesTree);
						listRef[listRef[0]+1] = i;
    					trouve=1;
			     }
				}
			}

			//	if(aHGT->crit.RF == 0 /*|| aHGT->crit.LS < epsilon */|| aHGT->crit.BD == 0) { i=j=(int)INFINI;}
			}

			if (trouve == 1){

				listRef[0] += 1;
				(*nbHgtFound)++;
				printf("\n[%2d--%2d -> %2d--%2d]",aHGT[(*nbHgtFound)-1].source_A,aHGT[(*nbHgtFound)-1].source_B,aHGT[(*nbHgtFound)-1].dest_A,aHGT[(*nbHgtFound)-1].dest_B);
			}

		}
		encore = 0;

  	if((*nbHgtFound == 0) && (flag2==1)){
      flag2 = 0;
      encore = 1;
    }

    if((*nbHgtFound == 0) && (*initial==1)){
      (*initial) = 0;
      encore = 1;
      flag2=1;
    }

 //   printf("On est rendu icit");
	}while(encore == 1);

	deleteBipartition(DTSpecies,SpeciesTree);
	//deleteBipartition(DTGene,GeneTree);
	deleteBipartitionSansRacine(DTGene,GeneTree.size);
	FreeCriteria(&aCrit,size);
	FreeCriteria(&aCritRef,size);
	FreeMemory_InputTree(&tmpTree,tmpTree.size);

	//printf("\n3)RF=%d",aHGT[0].crit.RF);
	return (*nbHgtFound);
}



//=================================================================
//==
//=================================================================
void findBranch(struct InputTree aTree,int *branch, int * elt){

	double max = INFINI; //= distance between root and intersection
	int j;

	(*branch) = 0;

	if(elt[0]==1){
		for(j=1;j<=2*aTree.size-3-aTree.kt;j++)
			if(aTree.ARETE[2*j-1]==elt[1] || aTree.ARETE[2*j-2]==elt[1])
				(*branch) = j;
	}
	else{
		for(j=2;j<=elt[0];j++){
			if(max > ( (aTree.ADD[aTree.Root][elt[1]] + aTree.ADD[aTree.Root][elt[j]] - aTree.ADD[elt[1]][elt[j]]) / 2.0 ) )
				max = (aTree.ADD[aTree.Root][elt[1]] + aTree.ADD[aTree.Root][elt[j]] - aTree.ADD[elt[1]][elt[j]]) / 2.0 ;
		}
		//printf("\nmax=%lf",max);

		for(j=1;j<=2*aTree.size-3-aTree.kt;j++){
			if( (fabs(aTree.ADD[aTree.Root][aTree.ARETE[2*j-1]]- max) < 2*epsilon) && (aTree.ADD[aTree.Root][aTree.ARETE[2*j-2]] < aTree.ADD[aTree.Root][aTree.ARETE[2*j-1]]) && (fabs(aTree.ADD[aTree.Root][aTree.ARETE[2*j-1]]+aTree.ADD[aTree.ARETE[2*j-1]][elt[1]] - aTree.ADD[aTree.Root][elt[1]]) < 2*epsilon)  )
				(*branch) = j;
			if( (fabs(aTree.ADD[aTree.Root][aTree.ARETE[2*j-2]]- max) < 2*epsilon) && (aTree.ADD[aTree.Root][aTree.ARETE[2*j-1]] < aTree.ADD[aTree.Root][aTree.ARETE[2*j-2]]) && (fabs(aTree.ADD[aTree.Root][aTree.ARETE[2*j-2]]+aTree.ADD[aTree.ARETE[2*j-2]][elt[1]] - aTree.ADD[aTree.Root][elt[1]]) < 2*epsilon))
				(*branch) = j;
		}
	}
}

double calculDifferenceMoyenne(struct HGT aHGT, double ** matrice, int size){
	double diff = 0;
	int i,j;

	for(i=1;i<=aHGT.listSource[0];i++){
		for(j=1;j<=aHGT.listDestination[0];j++){
			diff += matrice[aHGT.listSource[i]][aHGT.listDestination[j]];
		}
	}
	return diff/(aHGT.listSource[0]*aHGT.listDestination[0]);
}

//======================================================
//==
//======================================================
void expandBestHGT(struct HGT bestHGTRed,struct HGT *bestHGT,struct ReduceTrace aMap,struct DescTree * DTSpecies,struct InputTree SpeciesTree){

	//bestHGT->source = (int*)malloc(

	int i,j,nbSource=0,nbDest=0,tmp;

	for(i=1;i<=bestHGTRed.listSource[0];i++){
		for(j=1;j<=DTSpecies[aMap.map[bestHGTRed.listSource[i]]].nbSommet;j++){
			nbSource++;
		}
	}

	bestHGT->listSource = (int*)malloc((nbSource+1)*sizeof(int));
	bestHGT->listSource[0] = nbSource;
	tmp=1;
	for(i=1;i<=bestHGTRed.listSource[0];i++){
		for(j=1;j<=DTSpecies[aMap.map[bestHGTRed.listSource[i]]].nbSommet;j++){
			bestHGT->listSource[tmp] = DTSpecies[aMap.map[bestHGTRed.listSource[i]]].Tableau[j];
			tmp++;
		}
	}
	TrierTableau(bestHGT->listSource,tmp-1);

/*	printf("\nSource :");
	for(i=1;i<=bestHGT->listSource[0];i++)
		printf(" %d",bestHGT->listSource[i]);
    */

	for(i=1;i<=bestHGTRed.listDestination[0];i++){
		for(j=1;j<=DTSpecies[aMap.map[bestHGTRed.listDestination[i]]].nbSommet;j++){
			nbDest++;
		}
	}

	bestHGT->listDestination = (int*)malloc((nbDest+1)*sizeof(int));
	bestHGT->listDestination[0] = nbDest;
	tmp=1;
	for(i=1;i<=bestHGTRed.listDestination[0];i++){
		for(j=1;j<=DTSpecies[aMap.map[bestHGTRed.listDestination[i]]].nbSommet;j++){
			bestHGT->listDestination[tmp] = DTSpecies[aMap.map[bestHGTRed.listDestination[i]]].Tableau[j];
			tmp++;
		}
	}
	TrierTableau(bestHGT->listDestination,tmp-1);

/*	printf("\nDest   :");
	for(i=1;i<=bestHGT->listDestination[0];i++)
		printf(" %d",bestHGT->listDestination[i]);
	*/

	//== find branches
	findBranch(SpeciesTree,&i,bestHGT->listSource);
	findBranch(SpeciesTree,&j,bestHGT->listDestination);
	bestHGT->source = i;
	bestHGT->destination = j;

	bestHGT->source_A = SpeciesTree.ARETE[2*i-1];
	bestHGT->source_B = SpeciesTree.ARETE[2*i-2];
	bestHGT->dest_A = SpeciesTree.ARETE[2*j-1];
	bestHGT->dest_B = SpeciesTree.ARETE[2*j-2];

	bestHGT->valide = bestHGTRed.valide;
	bestHGT->crit = bestHGTRed.crit;

	free(bestHGTRed.listSource);
	free(bestHGTRed.listDestination);
	bestHGTRed.listSource = NULL;
	bestHGTRed.listDestination = NULL;
}

//======================================================
//==
//======================================================
int isInside(struct DescTree DT1,struct DescTree DT2){

	int i;

	for(i=1;i<=DT1.nbSommet;i++){
		if(DT2.Tableau[1] == DT1.Tableau[i]) return 1;
	}
	return 0;
}

//=================================================================
//==
//=================================================================
void CreateSubStructures(struct InputTree * aTree,int inc,int binaire){

	int n = aTree->size;
	int i,j;
	int kt=0;
	inc = 10;


	//printf("n=%d",n);
	if(aTree->ARETE == NULL){
		aTree->ARETE    =(long int*)malloc(4*(2*(n+inc))*sizeof(long int));
		aTree->LONGUEUR	=(double*)malloc((4*(n+inc))*sizeof(double));
		aTree->Adjacence=(double**)malloc((2*(n+inc)+1)*sizeof(double*));
		for(i=0;i<2*(n+inc);i++)
			aTree->Adjacence[i]=(double*)malloc((2*(n+inc)+1)*sizeof(double));
	}
//	printf("\nbinaire = %d",binaire);
	kt = aTree->kt = Tree_edges (aTree->ADD,aTree->ARETE,aTree->LONGUEUR,n,binaire);

	//printf("\n\n");
	//printEdges(0,aTree->ARETE,aTree->LONGUEUR,aTree->SpeciesName,aTree->size);
	//for(i=1;i<=2*n-3-kt;i++)
	//	printf("\n%d-%d : %lf",aTree->ARETE[2*i-1],aTree->ARETE[2*i-2],aTree->LONGUEUR[i-1]);
	//scanf("%d",&i);
	loadAdjacenceMatrix(aTree->Adjacence,aTree->ARETE, aTree->LONGUEUR,n,aTree->kt);
	Floyd(aTree->Adjacence,aTree->ADD,n,aTree->kt); // 4eme fois
//  global_cpt4++;

  //===creation de degre
	aTree->degre = (int*)malloc(2*(n+inc)*sizeof(int));
	for(i=1;i<=2*n-2-kt;i++){
		aTree->degre[i]=0;
		for(j=1;j<=2*n-2-kt;j++)
			if(aTree->Adjacence[i][j] < INFINI) aTree->degre[i]++;
	}
	/*for(i=1;i<=2*n-2-kt;i++)
		printf("%d ",aTree->degre[i]);*/

}

//=================================================================
//==
//=================================================================
void ReduceTree(struct InputTree SpeciesTree,struct InputTree GeneTree,struct InputTree *SpeciesTreeRed,struct InputTree *GeneTreeRed,struct ReduceTrace *aTrace,struct DescTree *DTSpecies,struct DescTree *DTGene,int binaireSpecies,int binaireGene){


	int i,j,l=1,tmp,pos; //k;
	struct CRITERIA aCrit;

	InitCriteria(&aCrit,SpeciesTree.size);

	//printf("\nReduceTree");

	aTrace->species  = (int *)malloc((2*SpeciesTree.size+1)*sizeof(int));
	aTrace->gene     = (int *)malloc(3*(2*GeneTree.size)*sizeof(int));
	aTrace->map		 = (int *)malloc((2*SpeciesTree.size)*sizeof(int));
	//aTrace->speciesToGene = (int *)malloc((2*SpeciesTree.size-1)*sizeof(int));
	aTrace->size = 0;

	for(i=0;i<=2*SpeciesTree.size;i++){
		aTrace->species[i] = 0;
	}
/*	for(i=0;i<3*(2*GeneTree.size);i++){
		aTrace->gene[i] = 0;
	}
*/	const int gInit = 3*(GeneTree.size+1)-2;
	const int gMax = 3*(2*(GeneTree.size)-2);
	const int sMax = 2*(SpeciesTree.size)-2-SpeciesTree.kt;
	for(i=1;i<=sMax;i++){

		if(DTSpecies[i].nbSommet == 1 ){
				aTrace->species[i] = 1;
				//aTrace->gene[i] =1;
				//aTrace->speciesToGene[i] = i;
		}

		//== modifie le 19 janvier 2009 par Alix

		for(j=gInit;(j<=gMax) && (i!=SpeciesTree.size);j++){
		//for(j=1;j<=2*(GeneTree.size)-2;j++){
			//printf("[%d(%d),%d(%d)]",DTSpecies[i].nbSommet,i,DTGene[j].nbSommet,j);
			if(DTSpecies[i].nbSommet == DTGene[j].nbSommet && DTSpecies[i].nbSommet > 1 && i!= SpeciesTree.size && j!= GeneTree.size){

				//printf("\nOn a le meme nombre de sommets !!");

				if(vecteursEgaux(DTSpecies[i],DTGene[j]) == 1){
						//printf("\nLes vecteurs sont egaux !!");
						computeCriteria(DTSpecies[i].Matrice,DTGene[j].Matrice,DTSpecies[i].nbSommet+1,&aCrit,NULL,NULL,NULL,NULL);
						if(aCrit.RF == 0){
							//printf("\nLes sous arbres sont egaux !!");
							//printNoeudDT(DTSpecies,i,0,0);
							//printNoeudDT(DTGene,j,0,0);
							aTrace->species[i] = 1;
							aTrace->gene[j] = 1;
						}
				}
			}
			if(DTSpecies[i].nbSommet == 1 ){
				aTrace->species[i] = 1;
				aTrace->gene[i] =1;
				//aTrace->speciesToGene[i] = i;
			}
		}
	}



		//== ici on test les sous-ensembles 2 par 2, si un est le sous-ensemble de l'autre, il n'est pas necessaire de consid�rer
		//== ce dernier.
		for(i=1;i<=2*SpeciesTree.size-2-SpeciesTree.kt;i++){
			for(j=1;j<=2*SpeciesTree.size-2-SpeciesTree.kt;j++){
				if(aTrace->species[i] != 0 && aTrace->species[j] != 0 && i!=j && i!=SpeciesTree.Root && j!=SpeciesTree.Root)
					if(DTSpecies[i].nbSommet > DTSpecies[j].nbSommet && isInside(DTSpecies[i],DTSpecies[j])==1){
						aTrace->species[j] = 0;
					}
			}
		}



		//== creer les matrices de distances
		tmp=0;
		pos=0;
		for(i=1;i<=2*SpeciesTree.size-2-SpeciesTree.kt;i++){
			if(aTrace->species[i] != 0){
				pos++;
				aTrace->map[pos] = i; //DTSpecies[i].Tableau[1];
	//			printf("<br> %d-%d",pos,i);
			}
		}
	/*	printf("\n");
		for(i=1;i<=pos;i++){
			printf("\n%d-->%d (%d)",i,aTrace->map[i],DTSpecies[aTrace->map[i]].Tableau[1]);
		}*/

		pos++; //= racine
		if(SpeciesTreeRed->ADD == NULL){
			SpeciesTreeRed->ADD = (double**)malloc((2*pos-1)*sizeof(double*));
			GeneTreeRed->ADD = (double**)malloc((2*pos-1)*sizeof(double*));
			SpeciesTreeRed->W = (double**)malloc((2*pos-1)*sizeof(double*));
			GeneTreeRed->W = (double**)malloc((2*pos-1)*sizeof(double*));
			for(i=0;i<2*pos-1;i++){
				SpeciesTreeRed->W[i] = (double*)malloc((2*pos-1)*sizeof(double));
				GeneTreeRed->W[i] = (double*)malloc((2*pos-1)*sizeof(double));
			}
			for(i=0;i<2*pos-1;i++){
				SpeciesTreeRed->ADD[i] = (double*)malloc((2*pos-1)*sizeof(double));
				GeneTreeRed->ADD[i] = (double*)malloc((2*pos-1)*sizeof(double));
			}
		}
		for(i=1;i<pos;i++)
			if(aTrace->map[i] != SpeciesTree.Root){
				for(j=1;j<pos;j++){
					if(aTrace->map[i] != SpeciesTree.Root){
						SpeciesTreeRed->ADD[i][j] = SpeciesTree.ADD[DTSpecies[aTrace->map[i]].Tableau[1]][DTSpecies[aTrace->map[j]].Tableau[1]];
						GeneTreeRed->ADD[i][j] = GeneTree.ADD[DTSpecies[aTrace->map[i]].Tableau[1]][DTSpecies[aTrace->map[j]].Tableau[1]];
					}
				}
			}

		aTrace->species[SpeciesTree.Root] = 1;
		aTrace->map[pos] = SpeciesTree.Root;
		for(i=1;i<pos;i++){
			SpeciesTreeRed->ADD[i][pos] = SpeciesTreeRed->ADD[pos][i] = SpeciesTree.ADD[aTrace->map[pos]][DTSpecies[aTrace->map[i]].Tableau[1]];
			GeneTreeRed->ADD[i][pos] = GeneTreeRed->ADD[pos][i] = epsilon; //= modifie le 19 janvier 2010 , GeneTree.ADD[aTrace->map[pos]][DTSpecies[aTrace->map[i]].Tableau[1]];
		}

		SpeciesTreeRed->size = GeneTreeRed->size = pos;
		SpeciesTreeRed->Root = GeneTreeRed->Root = pos;

		//== ajoute le 19 janvier 2010
		GeneTreeRed->size = pos-1;

		CreateSubStructures(SpeciesTreeRed,0,binaireSpecies);
		CreateSubStructures(GeneTreeRed,0,binaireGene);

		FreeCriteria(&aCrit,SpeciesTree.size);
	//	printf("\n");
	//	printMatrix(SpeciesTree.SpeciesName,SpeciesTreeRed->ADD,SpeciesTreeRed->size);
}

//=================================================================
//==
//=================================================================
int DeleteUseLessHGT(int nbHGT,struct HGT *bestHGT, struct InputTree SpeciesTree, struct InputTree SaveTree){

	int i,j,source,destination,p,q;
	struct InputTree tmpTree;
	struct CRITERIA aCrit;
	int cpt=0,retour=0;

	InitCriteria(&aCrit,SpeciesTree.size);
	initInputTree(&tmpTree);

	for(i=1;i<=nbHGT;i++){
		if(bestHGT[i].valide == 1){
			copyInputTree(&tmpTree,SaveTree,0,0);
			for(j=1;j<=nbHGT;j++){
				if(bestHGT[j].valide == 1){
					if (i!=j){
						findBranch(tmpTree,&source,bestHGT[j].listSource);
						findBranch(tmpTree,&destination,bestHGT[j].listDestination);
						if(source != 0 && destination !=0){
							if(isAValidHGT(tmpTree,source,destination)==1){
								applyHGT(SpeciesTree.ADD,&tmpTree,source,destination);
							}
						}
					}
				}
			}
			AdjustBranchLength(&tmpTree,SpeciesTree,0,1);
			//printf("\nHGT-DETECTION : on calcule les criteres");
			computeCriteria(SpeciesTree.ADD,tmpTree.ADD,SpeciesTree.size,&aCrit,SpeciesTree.LONGUEUR,SpeciesTree.ARETE,tmpTree.LONGUEUR,tmpTree.ARETE);
			//printf("\nHGT-DETECTION : on a termine le calcul");
			if(aCrit.RF == 0){
				bestHGT[i].valide = 0;
				retour++;
			}
		}
	}
	FreeMemory_InputTree(&tmpTree,tmpTree.size);
	FreeCriteria(&aCrit,SpeciesTree.size);
	return retour;
}




//=================================================================
//==
//=================================================================
void help(){
	printf("\nHGT-DETECTION version 3.0");
	printf("\nby Alix Boc and Vladimir Makarenkov");

	printf("\n\nUsage :\nhgt -inputfile=[inputfilename] -outputfile=[outputfilename] -criterion=[rf|ls|bd]");
	printf(" -version=[web|consol] -speciesroot=[midpoint|prompt|file] -generoot=[midpoint|prompt|file]");
	printf(" -load=[no|yes] -viewTree=[no|yes] -scenario=[unique|multiple] -nbhgt=[maxhgt] -path=[path]");

	printf("\n\ncriterion          [rf] = Robinson and Foulds distance (default)");
	printf("  \n                   [ls] = Least-Square optimization");
	printf("  \n                   [bd] = Bipartition distance");
	printf("\n\nversion            [consol] (default)");
	printf("  \n                   [web] = get result file and tree files in web ");
	printf("  \n                   format for the Trex-online web site");
	printf("\n\nspeciesroot        [midpoint] = the root is selected by the midpoint ");
	printf("  \n                   method (default)");
	printf("  \n                   [prompt] = the program ask for the root branch");
	printf("  \n                   [file] = the root is in a file called speciesroot.txt");
	printf("\n\ngeneroot           [midpoint] = the root is selected by the midpoint ");
	printf("  \n                   method (default)");
	printf("  \n                   [prompt] = the program ask for the root branch");
	printf("  \n                   [file] = the root is in a file called generoot.txt");
	printf("\n\nsubtree            [yes] = use the subtree constraint (default)");
	printf("  \n                   [no]");
	printf("\n\nscenario           [unique] (default)");
	printf("  \n                   [multiple]");
	printf("\n\nnbhgt              number max of hgts for unique and multiple scenario. ");
	printf("  \n                   Default value = 50 hgts");
	printf("\n\nload               [yes] = load the tree from files speciestree.txt ");
	printf("  \n                   genetree.txt speciesroot.txt generoot.txt");
	printf("  \n                   [no] (default)");
	printf("\n\npath               path without \"/\" at the end. Default path is \".\"");
	printf("  \n                   the path will be add to the file. ex : path/input.txt");

	printf("\n\nExample. \nCompute hgt-detection with default parameters : ./hgt -inputfile=input.txt\n\n");

}


void printTransfer(FILE *out,int sortie,char ** noms,int n,int source_A,int source_B,int dest_A,int dest_B){

	if(sortie == 0)
		fprintf(out,"\nFrom branch %d--%d to branch %d--%d",source_A,source_B,dest_A,dest_B);
	else{
		if(source_A	<= n)
			fprintf(out,"\nFrom branch %s--%d",noms[source_A],source_B);
		else if(source_B <= n)
			fprintf(out,"\nFrom branch %d--%s",source_A,noms[source_B]);
		else
			fprintf(out,"\nFrom branch %d--%d",source_A,source_B);

		if(dest_A <= n)
			fprintf(out," to branch %s--%d",noms[dest_A],dest_B);
		else if(dest_B <= n)
			fprintf(out," to branch %d--%s",dest_A,noms[dest_B]);
		else
			fprintf(out," to branch %d--%d",dest_A,dest_B);
	}
}
void printListSpecies(struct HGT * bestHGT, int nbHGT, double *tabBoot, int bootmin){

		FILE * out = fopen ("result.txt","w+");
		for(int i=1; i<= nbHGT && bestHGT[i].valide == 1; i++){
			if((tabBoot != NULL && tabBoot[i] >= bootmin) || (tabBoot == NULL)){
				fprintf(out,"\n%d",bestHGT[i].listSource[0]);
				for(int j=1;j<=bestHGT[i].listSource[0];j++)
					fprintf(out," %d",bestHGT[i].listSource[j]);
				fprintf(out,"\n%d",bestHGT[i].listDestination[0]);
				for(int j=1;j<=bestHGT[i].listDestination[0];j++)
					fprintf(out," %d",bestHGT[i].listDestination[j]);
			}
		}
		fclose(out);
}


//==============================================================================
//== AFFICHAGE DES RESULTATS DANS LE FICHIER DE SORTIE
//==============================================================================
void printHGT(FILE * res,struct CRITERIA * multicheckTab,char *mode,int RFref,FILE *out,struct InputTree SpeciesTree,struct HGT *tabHGT, int nbHGT, double *tabBoot,char *subtree,int bootmin){
	int i,j,k;
	int diffHGT=0;
	FILE * res1 = res; //fopen ("result2.txt","w+");
	int nbTrivial=0;

	fprintf(res,"mode=%s\n",mode);

  if(strcmp(mode,"monocheck") == 0 || multicheckTab == NULL){
		for(i=1;i<=nbHGT;i++){
			if(i==1)
				diffHGT = RFref - tabHGT[i].crit.RF;
			else
				diffHGT = tabHGT[i-1].crit.RF - tabHGT[i].crit.RF;
		//	fprintf(out,"\n\n================ \n= HGT #%d %s \n==============",i,(diffHGT==1)?" (Trivial transfer: RF distance decreased by 1)":"");
			fprintf(res,"1");
            fprintf(res,"\nHGT %d %s",i,(diffHGT==1)?" Trivial":"");
			nbTrivial += tabHGT[i].trivial; //(diffHGT==1)?1:0;
	//		if(tabBoot != NULL)
	//			fprintf(out,"\nBootstrap value : %3.1lf%%",(double)tabBoot[i]);
			fprintf(res,"\n");
			for(int m=1;m<=tabHGT[i].listSource[0];m++)
				fprintf(res,"%s ",SpeciesTree.SpeciesName[tabHGT[i].listSource[m]]);
			//fprintf(res,"\n%d",tabHGT[i].listDestination[0]);
			fprintf(res,"\n");
			for(int m=1;m<=tabHGT[i].listDestination[0];m++)
				fprintf(res,"%s ",SpeciesTree.SpeciesName[tabHGT[i].listDestination[m]]);
	//		printTransfer(out,1,SpeciesTree.SpeciesName, SpeciesTree.size,tabHGT[i].source_A,tabHGT[i].source_B,tabHGT[i].dest_A,tabHGT[i].dest_B);
			printTransfer(res,1,SpeciesTree.SpeciesName, SpeciesTree.size,tabHGT[i].source_A,tabHGT[i].source_B,tabHGT[i].dest_A,tabHGT[i].dest_B);
	//		fprintf(out,"\nRF = %d \nLS = %lf \nBD = %lf",tabHGT[i].crit.RF,tabHGT[i].crit.LS,tabHGT[i].crit.BD,tabHGT[i].crit.QD);
			fprintf(res,"\nRF = %d , LS = %lf , BD = %lf\n",tabHGT[i].crit.RF,tabHGT[i].crit.LS,tabHGT[i].crit.BD);
		}
		if(nbHGT == 0){
	//		fprintf(out,"\n\nNo more HGT have been detected due to evolutionary constraints");
	//		fprintf(res,"\nNo more HGT have been detected due to evolutionary constraints\n");
		}
		else
			if(tabHGT[nbHGT].crit.RF != 0){
	//			fprintf(out,"\n\nNo more HGT have been detected due to evolutionary constraints ");
	//			fprintf(res,"\nNo more HGT have been detected due to evolutionary constraints\n");
      }

	}
	else{

		i=1;
		int lastRF=0;
		int trivial;

		for(j=1;j<=multicheckTab[0].m;j++){

		  if(multicheckTab[j].nbHgtFound == 0) continue;
		//	fprintf(out,"\n\n=======> %d %s been found for this iteration",multicheckTab[j].nbHgtFound,(multicheckTab[j].nbHgtFound==1)?"hgt has":"hgts have");
			fprintf(res,"%d",multicheckTab[j].nbHgtFound);

      for(k=1;k<=multicheckTab[j].nbHgtFound;k++){
			//	for(i=1;i<=nbHGT;i++){
					if(j==1)
						diffHGT = RFref - tabHGT[i].crit.RF;
					else
						diffHGT = multicheckTab[j-1].RF - tabHGT[i].crit.RF;

					//	printf("\n%d = %d - %d",multicheckTab[j-1].RF - tabHGT[i].crit.RF,multicheckTab[j-1].RF,tabHGT[i].crit.RF);
					trivial=0;
          if(
              (tabHGT[i].source_A == tabHGT[i].dest_A ) ||
              (tabHGT[i].source_A == tabHGT[i].dest_B ) ||
              (tabHGT[i].source_B == tabHGT[i].dest_A ) ||
              (tabHGT[i].source_B == tabHGT[i].dest_B )
          ){
            trivial=1;
          }
		//			fprintf(out,"\n================ \n= HGT #%d/%d %s\n================",k,multicheckTab[j].nbHgtFound,(tabHGT[i].trivial == 1 )?" (Trivial transfer: RF distance decreased by 1)":"");//((trivial==1) && (diffHGT==1) )?" (Trivial transfer: RF distance decreased by 1)":"");
					fprintf(res,"\nHGT %d / %d %s",k,multicheckTab[j].nbHgtFound,(tabHGT[i].trivial == 1 )?" Trivial":"");//((trivial==1) && (diffHGT==1) )?" (Trivial transfer: RF distance decreased by 1)":"");
					nbTrivial += tabHGT[i].trivial;
          //fprintf(out,"\n%s(%d)",(tabHGT[i].trivial == 1)?"TOTO":"TATA",tabHGT[i].trivial);
    //      if(tabBoot != NULL)
		//				fprintf(out,"\nBootstrap value : %3.1lf%%",(double)tabBoot[i]);

					if((tabBoot != NULL && tabBoot[i] >= bootmin) || (tabBoot == NULL)){
						//fprintf(res,"\n%d",tabHGT[i].listSource[0]);
						fprintf(res,"\n");
						for(int m=1;m<=tabHGT[i].listSource[0];m++)
							fprintf(res,"%s ",SpeciesTree.SpeciesName[tabHGT[i].listSource[m]]);
						//fprintf(res,"\n%d",tabHGT[i].listDestination[0]);
						fprintf(res,"\n");
						for(int m=1;m<=tabHGT[i].listDestination[0];m++)
							fprintf(res,"%s ",SpeciesTree.SpeciesName[tabHGT[i].listDestination[m]]);

					}

		//			printTransfer(out,1,SpeciesTree.SpeciesName, SpeciesTree.size,tabHGT[i].source_A,tabHGT[i].source_B,tabHGT[i].dest_A,tabHGT[i].dest_B);
					printTransfer(res1,1,SpeciesTree.SpeciesName, SpeciesTree.size,tabHGT[i].source_A,tabHGT[i].source_B,tabHGT[i].dest_A,tabHGT[i].dest_B);
		//		  fprintf(out,"\nRF = %d \nLS = %lf \nBD = %lf",tabHGT[i].crit.RF,tabHGT[i].crit.LS,tabHGT[i].crit.BD,tabHGT[i].crit.QD);
				  fprintf(res1,"\nRF = %d , LS = %lf , BD = %lf",tabHGT[i].crit.RF,tabHGT[i].crit.LS,tabHGT[i].crit.BD);
					i++;
			//	}
			}
		//	fprintf(out,"\n\n=======> After this iteration the criteria are :");
	//		fprintf(out,"\nRF = %d \nLS = %lf \nBD = %lf",multicheckTab[j].RF,multicheckTab[j].LS,multicheckTab[j].BD,multicheckTab[j].QD);
			fprintf(res1,"\nRF = %d , LS = %lf , BD = %lf\n",multicheckTab[j].RF,multicheckTab[j].LS,multicheckTab[j].BD);
			lastRF = multicheckTab[j].RF;
		}
//		fprintf(out,"\n\nTotal number of HGTs : %d",nbHGT);

		if(lastRF != 0){
		//	fprintf(out,"\n\nNo more HGT have been detected due to evolutionary constraints");
	//		fprintf(res,"\nNo more HGT have been detected due to evolutionary constraints\n");
		}
	}
	      fprintf(res1,"Stat= %d %d",nbHGT,nbTrivial);
}
//=================================================================
//==
//=================================================================
void printBestHGT_F(FILE *out,int noHGT,struct InputTree SpeciesTree,struct HGT bestHGT,int *tmp,int nbTree, int boot){

	if(bestHGT.valide == 1){

		/*printf("\n\n======== Transfer %d ==========",noHGT-(*tmp));
		printf("\n%d--%d -> %d--%d",bestHGT.source_A,bestHGT.source_B,bestHGT.dest_A,bestHGT.dest_B);
		printf("\nRF = %d \nLS = %lf \nBD = %lf",bestHGT.crit.RF,bestHGT.crit.LS,bestHGT.crit.BD);*/

		fprintf(out,"\n\n================ \n= HGT #%d \n==============",noHGT-(*tmp));
		if(boot != -1)
			fprintf(out,"\n Boottrap value : %lf%%",(double)boot/nbTree*100.0);
		printTransfer(out,1,SpeciesTree.SpeciesName, SpeciesTree.size,bestHGT.source_A,bestHGT.source_B,bestHGT.dest_A,bestHGT.dest_B);
		//printf("\n%d--%d -> %d--%d",SpeciesTree.ARETE[2*(bestHGT.source)-1],SpeciesTree.ARETE[2*(bestHGT.source)-2],SpeciesTree.ARETE[2*(bestHGT.destination)-1],SpeciesTree.ARETE[2*(bestHGT.destination)-2]);
		fprintf(out,"\nRF = %d \nLS = %lf \nBD = %lf \nQD = %d",bestHGT.crit.RF,bestHGT.crit.LS,bestHGT.crit.BD,bestHGT.crit.QD);
/*
		fprintf(out,"\n%d",bestHGT.listSource[0]);
		for(j=1;j<=bestHGT.listSource[0];j++)
			fprintf(out," %d",bestHGT.listSource[j]);
		fprintf(out,"\n%d",bestHGT.listDestination[0]);
		for(j=1;j<=bestHGT.listDestination[0];j++)
			fprintf(out," %d",bestHGT.listDestination[j]);
*/	}
	else
		(*tmp) += 1;
}




int compareLeaves(int *listLeaves,int *tab,int n){
	int i,/*val,*/cpt1=0,cpt2=0;

	for(i=1;i<=n;i++){
		if(listLeaves[i] != tab[i]) cpt1++;
		tab[i] = (tab[i]+1) % 2;
	}
	for(i=1;i<=n;i++){
		if(listLeaves[i] != tab[i]) cpt2++;
	}

	if (cpt1 < cpt2) return cpt1;
	else return cpt2;

}
//============================================================================
//==
//============================================================================
int findApproxRootBranch(int * listLeaves,struct InputTree *aTree, int maxDiff){

	int choix =-1,i,j,/*k,*/val,valmax = maxDiff;

	int * tab = (int *) malloc((aTree->size+1)*sizeof(int));

	//= recherche de la branche equivalente a listLeaves

	for(i=1;i<=2*aTree->size-3-aTree->kt;i++){
		for(j=1;j<=aTree->size;j++)
			if(aTree->ADD[j][aTree->ARETE[2*i-1]] < aTree->ADD[j][aTree->ARETE[2*i-2]])
				tab[j] = 0;
			else
				tab[j] = 1;

		val = compareLeaves(listLeaves,tab,aTree->size);
		if(val <= valmax){
			valmax = val;
			choix = i;
		/*	printf("\n");
			for(k=1;k<=aTree->size;k++){
				printf("%d",tab[k]);
			}
			printf("no branch = %d , similarite = %d",choix,val);*/
		}
	}

	free(tab);

	return choix;
}

//==========================================================================
//== selectionne la branche voisine du midpoint qui minimise la distance de
//== Robinson and Foulds
//==========================================================================
int	bestRFNeighbor(struct InputTree *geneTree,struct InputTree *speciesTree,int position){
	int choix=position;
	int v1=-1,v2=-1,v3=-1,v4=-1,i;
	struct InputTree tmpTree;	//= arbre temporaire qui va etre modifier pour calculer RF
	struct CRITERIA aCrit;		//== struture of all the criteria
	int min= (int) INFINI;

	initInputTree(&tmpTree);

	//== recherche des voisins
	for(i=1;i<=2*geneTree->size-3-geneTree->kt;i++){

		if( (geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-1] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-1] ||
		     geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-2] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-2] ) && (i!=position) && (v1==-1)){

			v1 = i; //printf("\n(v1=%d,",i);
		}
		else if( (geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-1] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-1] ||
		     geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-2] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-2] ) && (i!=position) && (v2==-1)){

			v2 = i; //printf("v2=%d,",i);
		}
		else if( (geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-1] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-1] ||
		     geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-2] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-2] ) && (i!=position) && (v3==-1)){

			v3 = i; //printf("v2=%d,",i);
		}
		else if( (geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-1] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-1] ||
		     geneTree->ARETE[2*i-1] == geneTree->ARETE[2*position-2] || geneTree->ARETE[2*i-2] == geneTree->ARETE[2*position-2] ) && (i!=position) && (v4==-1)){

			v4 = i; //printf("v3=%d)",i);
		}
	}

	//printf("\n Selection de la racine :");

	copyInputTree(&tmpTree,(*geneTree),1,1);
	addLeafAndUpdate(&tmpTree,position);
	tmpTree.Root = tmpTree.size;

	InitCriteria(&aCrit,tmpTree.size);
	computeCriteria(speciesTree->ADD,tmpTree.ADD,tmpTree.size,&aCrit,speciesTree->LONGUEUR,speciesTree->ARETE,tmpTree.LONGUEUR,tmpTree.ARETE);

	//printf("\n(%d), %d--%d = %d",position,geneTree->ARETE[2*position-1],geneTree->ARETE[2*position-2],aCrit.RF);
	if(aCrit.RF < min){
		min=aCrit.RF;
		choix = position;
	}

	copyInputTree(&tmpTree,(*geneTree),1,0);
	addLeafAndUpdate(&tmpTree,v1);
	tmpTree.Root = tmpTree.size;

	FreeCriteria(&aCrit,tmpTree.size);
	InitCriteria(&aCrit,tmpTree.size);

	computeCriteria(speciesTree->ADD,tmpTree.ADD,tmpTree.size,&aCrit,speciesTree->LONGUEUR,speciesTree->ARETE,tmpTree.LONGUEUR,tmpTree.ARETE);

	//printf("\n(%d), %d--%d = %d",v1,geneTree->ARETE[2*v1-1],geneTree->ARETE[2*v1-2],aCrit.RF);
	if(aCrit.RF < min){
		min=aCrit.RF;
		choix = v1;
	}

	copyInputTree(&tmpTree,(*geneTree),1,0);
	addLeafAndUpdate(&tmpTree,v2);
	tmpTree.Root = tmpTree.size;

	FreeCriteria(&aCrit,tmpTree.size);
	InitCriteria(&aCrit,tmpTree.size);
	computeCriteria(speciesTree->ADD,tmpTree.ADD,tmpTree.size,&aCrit,speciesTree->LONGUEUR,speciesTree->ARETE,tmpTree.LONGUEUR,tmpTree.ARETE);

	//printf("\n(%d), %d--%d = %d",v2,geneTree->ARETE[2*v2-1],geneTree->ARETE[2*v2-2],aCrit.RF);
	if(aCrit.RF < min){
		min=aCrit.RF;
		choix = v2;
	}

	copyInputTree(&tmpTree,(*geneTree),1,0);
	addLeafAndUpdate(&tmpTree,v3);
	tmpTree.Root = tmpTree.size;

	FreeCriteria(&aCrit,tmpTree.size);
	InitCriteria(&aCrit,tmpTree.size);
	computeCriteria(speciesTree->ADD,tmpTree.ADD,tmpTree.size,&aCrit,speciesTree->LONGUEUR,speciesTree->ARETE,tmpTree.LONGUEUR,tmpTree.ARETE);
	//printf("\n(%d), %d--%d = %d",v3,geneTree->ARETE[2*v3-1],geneTree->ARETE[2*v3-2],aCrit.RF);
	if(aCrit.RF < min){
		min=aCrit.RF;
		choix = v3;
	}

	copyInputTree(&tmpTree,(*geneTree),1,0);
	addLeafAndUpdate(&tmpTree,v4);
	tmpTree.Root = tmpTree.size;

	FreeCriteria(&aCrit,tmpTree.size);
	InitCriteria(&aCrit,tmpTree.size);
	computeCriteria(speciesTree->ADD,tmpTree.ADD,tmpTree.size,&aCrit,speciesTree->LONGUEUR,speciesTree->ARETE,tmpTree.LONGUEUR,tmpTree.ARETE);
	//printf("\n(%d), %d--%d = %d",v4,geneTree->ARETE[2*v4-1],geneTree->ARETE[2*v4-2],aCrit.RF);
	if(aCrit.RF < min){
		min=aCrit.RF;
		choix = v4;
	}


	FreeCriteria(&aCrit,tmpTree.size);
	return choix;
}

int compareBranches(int * tab1,int * tab2,int size){
	int score1=0,score2=0,i;

	for(int i=1;i<=size;i++){
		score1 = (tab1[i] != tab2[i])?1:0;
		score2 = (tab1[i] == tab2[i])?1:0;
	}

	return (score1<score2)?score1:score2;
}

//==========================================================================
//== selectionne la branche qui minimise la distance de Robinson and Foulds
//==========================================================================
/*int bestbipartition(struct InputTree *geneTree, struct InputTree *speciesTree){

	int indiceBranche,i,bestScore=(int)INFINI,currentScore,bestPos=-1,j;
	struct InputTree tmpTree;		//= arbre temporaire qui va etre modifier pour calculer RF
	struct CRITERIA aCrit;			//== struture of all the criteria
	int choix=-1;
	int min= (int) INFINI;

	int * PLACEk1=(int *) malloc((2*geneTree->size-2)*sizeof(int));
	int * PLACEk2=(int *) malloc((2*geneTree->size-2)*sizeof(int));

	int ** Bk1=(int **) malloc((2*geneTree->size-2)*sizeof(int*));
	int ** Bk2=(int **) malloc((2*geneTree->size-2)*sizeof(int*));

	for (i=0;i<2*geneTree->size-2;i++)
	{
		Bk1[i]=(int *) malloc((geneTree->size)*sizeof(int));
		Bk2[i]=(int *) malloc((geneTree->size)*sizeof(int));
	}

	Bipartition_Table(speciesTree->ADD,Bk1,PLACEk1,geneTree->size+1);
	Bipartition_Table(geneTree->ADD.Matrice,Bk2,PLACEk2,geneTree->size+1);

	for(i=1;i<=2*geneTree->size-3-geneTree->kt;i++){
		for(j=1;j<=2*geneTree->size-3-geneTree->kt;j++){
			currentScore = compareBranches(Bk1[i],Bk2[j],geneTree->size);
			if(currentScore < bestScore){
				bestScore = currentScore;
				bestPos = i;
				if(bestScore == 0){
					i=j=(int)INFINI;
				}
			}
		}
	}

	//============
	int pos_a=-1,pos_b=-1;
	int nbUn=0;
	int nbZero=0;
	for(i=1;i<=2*geneTree->size-3-geneTree->kt;i++){
		if()
	}


	for(i=1;i<=2*geneTree->size-3-geneTree->kt;i++){
		if(geneTree->ADD[geneTree->ARETE[2*i-1]][geneTree->ARETE[2*i-1]]){
		}
	}


	return choix;
}*/
//==========================================================================
//== selectionne la branche qui minimise la distance de Robinson and Foulds
//==========================================================================
int bestRFBranch(struct InputTree *geneTree, struct InputTree *speciesTree){

	int indiceBranche;
	struct InputTree tmpTree;	//= arbre temporaire qui va etre modifier pour calculer RF
	struct CRITERIA aCrit;		//== struture of all the criteria
	int choix=-1;
	int min= (int) INFINI;

	printf("\nTaille = %d",geneTree->size);
	initInputTree(&tmpTree);

	printf("\nTaille = %d",geneTree->size);
	for(indiceBranche=1;indiceBranche<=2*geneTree->size-3-geneTree->kt;indiceBranche++){

		copyInputTree(&tmpTree,(*geneTree),0,0);
		addLeafAndUpdate(&tmpTree,indiceBranche);
		tmpTree.Root = tmpTree.size;

		InitCriteria(&aCrit,tmpTree.size);
		computeCriteria(speciesTree->ADD,tmpTree.ADD,tmpTree.size,&aCrit,speciesTree->LONGUEUR,speciesTree->ARETE,tmpTree.LONGUEUR,tmpTree.ARETE);

		printf("\n(%d), %ld--%ld = %d",indiceBranche,geneTree->ARETE[2*indiceBranche-1],geneTree->ARETE[2*indiceBranche-2],aCrit.RF);
		if(aCrit.RF < min){
			min=aCrit.RF;
			choix = indiceBranche;
		}
		FreeCriteria(&aCrit,tmpTree.size);
	}

	return choix;
}
//==========================================================================================================
//==
//==========================================================================================================
int addRootByNoBranch(struct InputTree *aTree,int noBranch){
	addLeafAndUpdate(aTree,noBranch);
	aTree->Root = aTree->size;
	for(int i=1;i<=2*(aTree->size)-2-aTree->kt;i++){
		aTree->degre[i]=0;
		for(int j=1;j<=2*aTree->size-2-aTree->kt;j++)
			if(aTree->Adjacence[i][j] < INFINI) aTree->degre[i]++;
	}
	return 0;
}

bool file_exists(const char * filename)
{
	FILE *file;

    if ((file=fopen(filename, "r"))==NULL)
    {
		return false;
    }
	else{
        fclose(file);
        return true;
	}
}

int addRoot(struct InputTree *aTree,struct InputTree *refTree,const char * message, const char * add,char *fichier,char *fichier_leaves,int * listLeaves){

	int choix=-1,choix2=-1;
	FILE *in;
	int R1,R2,i,j;


	if(strcmp(add,"prompt") == 0 ){
		printf("\n%s",message);
		printEdges(NULL,1,aTree->ARETE,aTree->LONGUEUR,aTree->SpeciesName,aTree->size,NULL,0,aTree->kt);

		while(choix<1 || choix > 2*aTree->size-3-aTree->kt){
			printf("\n\nSelect the root branch :");
			scanf("%d",&choix2);
		}
	}
	else if(strcmp(add,"file") == 0 ){
    if(file_exists(fichier_leaves)){
			int * listLeaves2 = (int*) malloc((aTree->size+1)*sizeof(int));
			for(i=1;i<=aTree->size;i++)
				listLeaves2[i] = 2;

			if((in=fopen(fichier_leaves,"r"))==NULL){
				printf("\nHGT-DETECTION : Can't open root file (%s)",fichier_leaves);
				exit(-1);
			}

			char * tmp;
			tmp = (char*)malloc(100);

			do{
				fscanf(in,"%s",tmp);
				if(strcmp(tmp,"<>") != 0){
					for(i=1;i<=aTree->size;i++){
						if(strcmp(aTree->SpeciesName[i],tmp) == 0){
							listLeaves2[i] = 0;
						}
					}
				}
			}while(strcmp(tmp,"<>") != 0);

			while(fscanf(in,"%s",tmp) != -1){
				//printf("\n%s ",tmp);
				if(strcmp(tmp,"<>") != 0){
					for(i=1;i<=aTree->size;i++){
						if(strcmp(aTree->SpeciesName[i],tmp) == 0){
							listLeaves2[i] = 1;
							//printf(" (1) ");
						}
					}
				}
			}

			fclose(in);

			choix = choix2 = findApproxRootBranch(listLeaves2,aTree,10000);

		}else{
  		if((in=fopen(fichier,"r"))==NULL){
  			printf("\nCan't open root file (%s)",fichier);
  			exit(-1);
  		}
  		fscanf(in,"%d %d",&R1,&R2);
  		for(i=1;i<=2*aTree->size-3-aTree->kt;i++){
  			if((aTree->ARETE[2*i-1] == R1 && aTree->ARETE[2*i-2] == R2) || (aTree->ARETE[2*i-1] == R2 && aTree->ARETE[2*i-2] == R1))
  				choix2 = choix = i;
  		}
  		fclose(in);
    }
	}
	else if(strcmp(add,"bestrfbranch") == 0){
		printf("bestRfBranch = %d",choix);
		choix2 = choix = bestRFBranch(aTree,refTree);
		printf("choix = %d",choix);
		//listLeaves = NULL;
	}
	/*else if(strcmp(add,"bestbipartition") == 0){
		printf("choix = %d",choix);
		choix2 = choix = bestbipartition(aTree,refTree);
		printf("choix = %d",choix);
		//listLeaves = NULL;
	}*/
	else{
		printf("\nHGT-DETECTION : Recherche de la racine en fonction du midpoint");
		choix2 = choix = midPoint(aTree->ARETE,aTree->ADD,aTree->size,aTree->kt);
		if(refTree != NULL){
		//	choix2 = bestRFBranch(aTree,refTree);
			choix2 = choix = bestRFNeighbor(aTree,refTree,choix);
		}
	}

	//printf("\nchoix = %d",choix);
	if(listLeaves != NULL){
		//== sauvegarder la liste des feuilles de chaque cot� de l'arete.
//printf("\nje passe ici");
		if(listLeaves[0] == -1){
			for(i=1;i<=aTree->size;i++)
				if(aTree->ADD[i][aTree->ARETE[2*choix-1]] < aTree->ADD[i][aTree->ARETE[2*choix-2]])
					listLeaves[i] = 0;
				else
					listLeaves[i] = 1;
			listLeaves[0]=1;
		}
		else{

			if(strcmp(add,"bestbipartition") == 0){
				printf("\nHGT-DETECTION : Recherche de la racine en fonction de la bipartition");
				choix2 = findApproxRootBranch(listLeaves,aTree,10000);
			}
			printf("\nHGT-DETECTION : choix = %d",choix2);
		}

		if(choix2==-1)
			return -1;
		else
			choix=choix2;
	}

	printRoot(fichier,aTree->ARETE[2*choix-1],aTree->ARETE[2*choix-2]);
  printRootByLeaves(fichier_leaves,choix,aTree);

	addLeafAndUpdate(aTree,choix);

	aTree->Root = aTree->size;
//printf("\nracine = %d",aTree->Root);
	//printf("\napres addleafAndupdate");
	//printEdges(NULL,1,aTree->ARETE,aTree->LONGUEUR,aTree->SpeciesName,aTree->size,NULL,0,aTree->kt);
	//===mise a jour de degre
	for(i=1;i<=2*(aTree->size)-2-aTree->kt;i++){
		aTree->degre[i]=0;
		for(j=1;j<=2*aTree->size-2-aTree->kt;j++)
			if(aTree->Adjacence[i][j] < INFINI) aTree->degre[i]++;
	}
	/*for(i=1;i<=2*aTree->size-2-aTree->kt;i++)
		printf("%d ",aTree->degre[i]);*/


	return 0;

}


//=================================================================
//==
//=================================================================
void chargerFichier(struct InputTree *aTree,const char* branchesFile,const char* rootFile){

	FILE *branch = fopen(branchesFile,"r");
	FILE *root =  fopen(rootFile,"r");
	int aa,bb,i,noBranch=0;
	double ll;

	for(i=1;i<=2*(aTree->size)-3;i++){
		fscanf(branch,"%d %d %lf",&aa,&bb,&ll);
		aTree->ARETE[2*i-1] = aa;
		aTree->ARETE[2*i-2] = bb;
		aTree->LONGUEUR[i-1] = ll;
		//printf("\n%d %d %lf",ARETE[2*i-1],ARETE[2*i-2],LONGUEUR[i-1]);
	}
	fscanf(root,"%d %d",&aa,&bb);

	for(i=1;i<=2*(aTree->size)-3;i++)
		if((aTree->ARETE[2*i-1] == aa && aTree->ARETE[2*i-2] == bb) || (aTree->ARETE[2*i-1] == bb && aTree->ARETE[2*i-2] == aa))
			noBranch = i;
	//printf("-->nobranch=%d<--",i);
	addLeafAndUpdate(aTree,noBranch);
	aTree->Root = aTree->size;

	fclose(branch);
	fclose(root);
}


//=================================================================
//== read the species or the gene tree in the input file
//=================================================================
int readInput(int Type, const char *file,struct InputTree * aTree){

	int size,i,j;
	char name[50];
	double val;
	FILE * in;

	//= ouverture du fichier
	if((in = fopen(file,"r"))== NULL)
		return -1;

	//= lecture de la taille des matrices
	fscanf(in,"%d",&size);

	//= allocation de la m�moire
	//allocMemmory(aTree,1);
	aTree->size = size;
	size++; // more space for the root
	aTree->SpeciesName = (char **)malloc((size+1)*sizeof(char*));
	aTree->Input = (double**)malloc((2*size)*sizeof(double*));
	aTree->ADD = (double**)malloc((2*size)*sizeof(double*));
	aTree->W = (double**)malloc((size+1)*sizeof(double*));
	for(i=0;i<2*size;i++){
		aTree->ADD[i] = (double*)malloc((2*size)*sizeof(double));
		aTree->Input[i] = (double*)malloc((2*size)*sizeof(double));
		if(i<=size){
			aTree->SpeciesName[i] = (char*)malloc(SPECIES_NAME_LENGTH);
			aTree->W[i] = (double*)malloc((size+1)*sizeof(double));
		}
	}

	for(i=0;i<=size;i++)
		for(j=0;j<=size;j++)
			aTree->W[i][j] = 1.0;

	size--;
	//= reads species tree
	for(i=1;i<=size;i++)
	{
		fscanf(in,"%s",name);
		if(Type == SPECIE) strcpy(aTree->SpeciesName[i],name);
		for(j=1;j<=size;j++)
		{
			fscanf(in,"%lf",&val);
			if(Type == SPECIE) aTree->Input[i][j] = val;
		}
	}
  if(Type == SPECIE) {
    strcpy(aTree->SpeciesName[size+1],"Root");
  	fclose(in);
    return 0;
  }
	//= read gene tree
  fscanf(in,"%d",&size);

	//= allocation de la m�moire
	//allocMemmory(aTree,1);
	aTree->size = size;
	size++; // more space for the root
	aTree->SpeciesName = (char **)malloc((size+1)*sizeof(char*));
	aTree->Input = (double**)malloc((2*size)*sizeof(double*));
	aTree->ADD = (double**)malloc((2*size)*sizeof(double*));
	aTree->W = (double**)malloc((size+1)*sizeof(double*));
	for(i=0;i<2*size;i++){
		aTree->ADD[i] = (double*)malloc((2*size)*sizeof(double));
		aTree->Input[i] = (double*)malloc((2*size)*sizeof(double));
		if(i<=size){
			aTree->SpeciesName[i] = (char*)malloc(SPECIES_NAME_LENGTH);
			aTree->W[i] = (double*)malloc((size+1)*sizeof(double));
		}
	}

	for(i=0;i<=size;i++)
		for(j=0;j<=size;j++)
			aTree->W[i][j] = 1.0;

	size--;
	for(i=1;i<=size;i++)
	{
		fscanf(in,"%s",name);
    printf("->%s",name);
		if(Type == GENE) strcpy(aTree->SpeciesName[i],name);
		for(j=1;j<=size;j++)
		{
			fscanf(in,"%lf",&val);
			if(Type == GENE) aTree->Input[i][j] = val;
		}
	}

	strcpy(aTree->SpeciesName[size+1],"Root");

	fclose(in);

	return 0;
}
//===========================================================
//== print the input matrix
//===========================================================
void printMatrix(char** Name, double **Matrix,int size){

	int i,j;

	for(i=1;i<=size;i++){
//		printf("\n%s",Name[i]);
		printf("\n%d",i);
		for(j=1;j<=size;j++){
			if(Matrix[i][j] >= INFINI)
				printf("  INFINI ");
			else
				printf(" %3.5lf",Matrix[i][j]);
		}
	}
}
//===========================================================
//== print the input matrix
//===========================================================
void printMatrix2(char** Name, double **Matrix,int size){

	int i,j;

	for(i=1;i<=size;i++){
		//printf("\n%s",Name[i]);
		printf("\n%d",i);
		for(j=1;j<=size;j++){
			if(Matrix[i][j] < epsilon )
				printf("  NULL ");
			else
				printf(" %3.5lf",Matrix[i][j]);
		}
	}
}

//==============================================================
//== print branches
//==============================================================
void printBranches(FILE *out,struct InputTree aTree,const char * Message,int * BSARETE,int nbTree){

	if(out==NULL)
		printf("%s",Message);
	else
		fprintf(out,"\n\n%s",Message);
	printEdges(out,1,aTree.ARETE,aTree.LONGUEUR,aTree.SpeciesName,aTree.size,BSARETE,nbTree,aTree.kt);
}



//=====================================================================
//
//=====================================================================
int ExtraireDonnees(const char * chaine, char *champs, char * contenu){

	int cpt=0,i;
	int egale=false;
	int tailleChaine;

	if(chaine[0] != '-')
		return false;
	tailleChaine = (int)strlen(chaine);

	for(i=1;i<tailleChaine;i++){

		if (chaine[i] == '='){
			egale = true;
			champs[cpt] = '\0';
			cpt=0;
			continue;
		}
		if(egale)
			contenu[cpt++] = chaine[i];
		else
			champs[cpt++] = chaine[i];
	}
	contenu[cpt] = '\0';

	if (!egale)
		return false;

	return true;
}
//===================================================================================
//=
//===================================================================================
int readParameters(struct Parameters * param, char **argv, int nargc){

	char champs[100];
	char contenu[100];
	char input[100];
	char output[100];
	char hgtResultFile[100];
	int i;
	sprintf((*param).sort,"yes");
  sprintf((*param).printWeb,"yes");
	sprintf((*param).criterion,"bd");
	sprintf((*param).verbose,"no");
  sprintf((*param).mode,"multicheck");
	sprintf((*param).viewtree,"no");
	sprintf((*param).generoot,"midpoint");
	sprintf((*param).speciesroot,"midpoint");
	sprintf((*param).load,"no");
	sprintf((*param).version,"consol");
	sprintf((*param).multiple,"no");
	sprintf((*param).multigene,"no");
	sprintf((*param).path,".");
	sprintf(input,"_");
	sprintf(output,"output.txt");
	sprintf(hgtResultFile,"hgtresultfile.txt");
	sprintf((*param).scenario,"unique");
	sprintf((*param).subtree,"yes");
	sprintf((*param).bootstrap,"no");
  sprintf((*param).addroot,"no");
  (*param).constraints = 0;	//= 0 : pas de contraintes
	(*param).nbhgt = 50;
	(*param).bootmin = 0;
  (*param).verbose_to_screen=false;
	(*param).avgdiff = 0.35;
	(*param).avgdiffblock = 0.22;
	(*param).c1 = 100;
	(*param).c2 = 0.1;

	sprintf((*param).special,"no");

	for(i=1;i<nargc;i++){

		if(ExtraireDonnees(argv[i],champs,contenu)){

			//printf("\nchamps = %s , contenu = %s",champs,contenu);

			if(strcmp("help",champs) == 0)
				return -1;
			//======== input file ==============
			if(strcmp("printWeb",champs) == 0){
				strcpy((*param).printWeb,contenu);
			}
			//======== input file ==============
			else if(strcmp("inputfile",champs) == 0){
				strcpy(input,contenu);
			}
			//============ output file ===============
			else if(strcmp("outputfile",champs) == 0){
				strcpy(output,contenu);
			}
			//============ translations file ===============
			else if(strcmp("translationsfile",champs) == 0){
				strcpy((*param).translationsfile,contenu);
			}
			//============ output file ===============
			else if(strcmp("hgtresultfile",champs) == 0){
				strcpy(hgtResultFile,contenu);
			}
			//============ criterion ================
			else if(strcmp("criterion",champs) == 0){
				strcpy((*param).criterion,contenu);
			}
      //============= version ===============
			else if(strcmp("v",champs) == 0){
				(*param).verbose_to_screen=true;
			}
			//============= version ===============
			else if(strcmp("version",champs) == 0){
				strcpy((*param).version,contenu);
			}
			//============ species root ===============
			else if(strcmp("speciesroot",champs) == 0){
				strcpy((*param).speciesroot,contenu);
			}
			//============ gene root ===============
			else if(strcmp("generoot",champs) == 0){
				strcpy((*param).generoot,contenu);
			}
			//======== load structure ==========
			else if(strcmp("load",champs) == 0){
				strcpy((*param).load,contenu);
			}
			//============ view tree ===============
			else if(strcmp("viewtree",champs) == 0){
				strcpy((*param).viewtree,contenu);
			}
			//============ path data ===========
			else if(strcmp("path",champs) == 0){
				strcpy((*param).path,contenu);
			}
			//============ special ===========
			else if(strcmp("special",champs) == 0){
				strcpy((*param).special,contenu);
			}
			//============ subtree ===========
			else if(strcmp("subtree",champs) == 0){
				strcpy((*param).subtree,contenu);
			}
      //============ addRoot ===========
			else if(strcmp("addroot",champs) == 0){
				strcpy((*param).addroot,contenu);
			}
			//============ scenario ===========
			else if(strcmp("scenario",champs) == 0){
				strcpy((*param).scenario,contenu);
			}
			//============ nbHGT ===========
			else if(strcmp("nbhgt",champs) == 0){
				(*param).nbhgt = atoi(contenu);
			}
			else if(strcmp("constraints",champs) == 0){
				(*param).constraints = atoi(contenu);
			}
			//============ sort ===========
			else if(strcmp("sort",champs) == 0){
				sprintf((*param).sort,"%s",contenu);
				if(strcmp(contenu,"yes") != 0 && strcmp(contenu,"no") !=0){
					printf("\nincorrect value for sort [%s] ,accepted values : yes or no",(*param).sort);
					exit(-1);
				}
			}
			//============ bootstrap ===========
			else if(strcmp("bootstrap",champs) == 0){
				sprintf((*param).bootstrap,"%s",contenu);
			}
			//============ multigene ===========
			else if(strcmp("multigene",champs) == 0){
				sprintf((*param).multigene,"%s",contenu);
			}
			//============ bootmin ===========
			else if(strcmp("bootmin",champs) == 0){
				(*param).bootmin = atoi(contenu);
			}
      else if(strcmp("avg",champs) == 0){
				(*param).avgdiff = atof(contenu);
			}
			else if(strcmp("blk",champs) == 0){
				(*param).avgdiffblock = atof(contenu);
			}
      else if(strcmp("c1",champs) == 0){
				(*param).c1 = atof(contenu);
			}
			else if(strcmp("c2",champs) == 0){
				(*param).c2 = atof(contenu);
			}

			//============ mode ===========
			else if(strcmp("mode",champs) == 0){
				sprintf((*param).mode,"%s",contenu);
				if(strcmp(contenu,"multicheck") != 0 && strcmp(contenu,"monocheck") !=0){
					printf("\nincorrect value for %s , values accepted : multicheck or monocheck",(*param).mode);
					exit(-1);
				}
			}
			//============ verbose ===========
			else if(strcmp("verbose",champs) == 0){
				sprintf((*param).verbose,"%s",contenu);
				if(strcmp(contenu,"yes") != 0 && strcmp(contenu,"no") !=0){
					printf("\nincorrect value for %s , values accepted : yes or no",(*param).verbose);
					exit(-1);
				}
			}
			else{
				printf("incorrect prameter : %s",champs);
				exit(-1);
			}
		}
	}

	if(strcmp(input,"_")==0)
		return -1;

	sprintf((*param).inputfile,"%s/%s",(*param).path,input);
	sprintf((*param).input,"%s/input_.txt",(*param).path);
	sprintf((*param).outputfile,"%s/%s",(*param).path,output);
	sprintf((*param).results,"%s/results.txt",(*param).path);
	sprintf((*param).hgtResultFile,"%s/%s",(*param).path,hgtResultFile);
	sprintf((*param).speciesTree,"%s/speciesTree.txt",(*param).path);
	sprintf((*param).geneTree,"%s/geneTree.txt",(*param).path);
	sprintf((*param).speciesRootfile,"%s/speciesRoot.txt",(*param).path);
	sprintf((*param).geneRootfileLeaves,"%s/speciesRootLeaves.txt",(*param).path);
	sprintf((*param).geneRootfile,"%s/geneRoot.txt",(*param).path);
	sprintf((*param).geneRootfileLeaves,"%s/geneRootLeaves.txt",(*param).path);
	sprintf((*param).speciesTreeWeb,"%s/speciesTreeWeb.txt",(*param).path);
	sprintf((*param).geneTreeWeb,"%s/geneTreeWeb.txt",(*param).path);
	sprintf((*param).outputWeb,"%s/outputWeb.txt",(*param).path);
  sprintf((*param).filteredLanguageTree,"%s/_filteredLangueTree.new",(*param).path);

	return 0;
}

void saveTree(char * fichier,struct InputTree SpeciesTree,struct HGT *bestHGT, int addHGT,int cpt_hgt,char *subtree,char *SCENARIO,double *bootStrap){

	FILE *out;
	int i, first=0;

	if ((out=fopen(fichier,"w+"))==NULL){
		printf("can't open species web file");
		exit(-1);
	}

	fprintf(out,"maptype = VERTICAL");
	fprintf(out,"\nintVerPrint = O");
	fprintf(out,"\nproportion = N");
	fprintf(out,"\ndrawRealNames = 1");
	fprintf(out,"\ncolorobject = BLACK");
	fprintf(out,"\ncolorreticulation = MAGENTA");
	fprintf(out,"\ncoloredge = BLUE");
	fprintf(out,"\nkt = %d",SpeciesTree.kt);
	fprintf(out,"\ndrawRetic = %s",(addHGT == 1)?"2":"0");
	fprintf(out,"\nn = %d",SpeciesTree.size);
	fprintf(out,"\net = ");
	for(i=1;i<=SpeciesTree.size;i++){
		if(i!=1) fprintf(out,",");
		fprintf(out,"%s",SpeciesTree.SpeciesName[i]);
	}
	fprintf(out,"\naretes = ");
	for(i=1;i<=2*SpeciesTree.size-3-SpeciesTree.kt;i++){
		if(i!=1) fprintf(out,",");
		fprintf(out,"%ld,%ld",SpeciesTree.ARETE[2*i-1],SpeciesTree.ARETE[2*i-2]);
	}
	fprintf(out,"\nlongueur = ");
	for(i=1;i<=2*SpeciesTree.size-3-SpeciesTree.kt;i++){
		if(i!=1) fprintf(out,",");
		fprintf(out,"%lf",SpeciesTree.LONGUEUR[i-1]);
	}
	if(addHGT == 1){
		fprintf(out,"\nhgt = ");
		{
			for(i=1;i<=cpt_hgt;i++){
				if((bestHGT[i].valide ==1) && (bestHGT[i].trivial==0)){
					if(first==1) fprintf(out,",");
					first=1;
					fprintf(out,"%d,%d,%d,%d",bestHGT[i].source_A,bestHGT[i].source_B,bestHGT[i].dest_A,bestHGT[i].dest_B);
				}
			}
			if(bootStrap != NULL){
				first=0;
				fprintf(out,"\nbootHGT = ");
				for(i=1;i<=cpt_hgt;i++){
					if(first==1) fprintf(out,",");
					first=1;
					fprintf(out,"%1.0lf",bootStrap[i]);
				}
			}
		}
	}
	if(addHGT == 1)
		fprintf(out,"\nroot=%d",SpeciesTree.Root);
	else
		fprintf(out,"\nroot=0");

	fclose(out);
}



//=============================================================
//
//=============================================================
void deleteSommet(struct DescTree *tab, int t1, int t2){

	int i,j,k=0,toDel=1,tmp;

	// si les elements de tab1 se trouve dans tab2 on les supprimes dans tab2
	for(i=0;i<tab[t1].nbSommet;i++){
		tmp=0;
		for(j=0;j<tab[t2].nbSommet;j++){
			if(tab[t1].Tableau[i] == tab[t2].Tableau[j]){
				//printf("\n=>%d (%d)",tab[t1].Tableau[i],tab[t1].nbSommet);
				tmp=1;
			}
		}
		if(tmp==0) toDel=0;
	}
	for(i=0;i<tab[t1-1].nbSommet;i++){
		tmp=0;
		for(j=0;j<tab[t2].nbSommet;j++){
			if(tab[t1-1].Tableau[i] == tab[t2].Tableau[j]){
				//printf("\n=>%d (%d)",tab[t1-1].Tableau[i],tab[t1-1].nbSommet);
				tmp=1;
			}
		}
		if(tmp==0) toDel=0;
	}


	// suppression des sommets
	if(toDel == 1){
		for(i=0;i<tab[t1].nbSommet;i++){
			for(j=0;j<tab[t2].nbSommet;j++){
				if(tab[t1].Tableau[i] == tab[t2].Tableau[j]){
					tab[t2].Tableau[j]=-1;
					k=j;
				}
			}

			for(j=k;j<tab[t2].nbSommet-1;j++){
				tab[t2].Tableau[j] = tab[t2].Tableau[j+1];
			}
			tab[t2].nbSommet--;

		}
	}
}


//===================================================================
//
//===================================================================
void supprimerSousEnsemble(struct DescTree *tabDetect,int nbTableau){

	int i,j/*,k*/;

	for(i=1;i<=nbTableau;i=i+2){
		for(j=i+1;j<=nbTableau;j++){
			deleteSommet(tabDetect,i,j);
		}
	}
}

//==========================================================================================================================================
//
//==========================================================================================================================================
int formatResult(struct HGT *tabHGT,int nbHGT, struct HGT *outHGT,struct InputTree aTree){

	int nbSommet,i,j,nbTrans,sommet,tmp=1,cpt=-1,nbHGTRet=0;
	int racine;
	double max;
	int nbArete = 2*(aTree.size)-3-aTree.kt;

	struct DescTree *tabDetect  = (struct DescTree*)malloc((2*(nbHGT+1))*sizeof(struct DescTree));

	//res = fopen(result,"r");

	//==lecture des transferts

	printf("\nHGT-DETECTION : nombre de HGT avant formatRes %d",nbHGT);

	for(i=1;i<=nbHGT;i++){

		if(tabHGT[i].valide == 1){
			nbHGTRet++;
			outHGT[nbHGTRet].valide = 1;
			copyHGT(tabHGT[i],&outHGT[nbHGTRet]);
			cpt++;
			nbSommet = tabDetect[cpt].nbSommet = tabHGT[i].listSource[0];

			//printf("\ni=%d",i);
			tabDetect[cpt].Tableau = (int*)malloc(nbSommet*sizeof(int));
			//printf("\n");
			for(j=0;j<nbSommet;j++){
				tabDetect[cpt].Tableau[j] = tabHGT[i].listSource[j+1];
			//	printf("%d ",tabHGT[i].listSource[j+1]);
			}
			//printf("\ni=%d",i);
			//if(i==73) {printf("\nOn sort nbSommet=%d",nbSommet); exit(1);}
			cpt++;
			nbSommet = tabDetect[cpt].nbSommet = tabHGT[i].listDestination[0];
			tabDetect[cpt].Tableau = (int*)malloc(nbSommet*sizeof(int));
			//printf("\n");
			for(j=0;j<nbSommet;j++){
				tabDetect[cpt].Tableau[j] = tabHGT[i].listDestination[j+1];
				//printf("%d ",tabHGT[i].listDestination[j+1]);
			}
		}
		else{
			printf("\nHGT-DETECTION : Transfer non valide %d",i);
		}
	}

	outHGT[nbHGTRet+1].valide=0;
/*
	nbTrans=-1;
	do{
		val=fscanf(res,"%d",&nbSommet);
		if(val != -1){
			printf("<br>");
			nbTrans++;
			tabDetect[nbTrans].nbSommet=nbSommet;
			tabDetect[nbTrans].Tableau = (int*)malloc(nbSommet*sizeof(int));
			for(j=0;j<nbSommet;j++){
				fscanf(res,"%d",&val);
				printf("%d ",val);
				tabDetect[nbTrans].Tableau[j] = val;
			}
		}
	}while(val!=-1);*/

	nbTrans = cpt;

//	printf("\n");
	supprimerSousEnsemble(tabDetect,nbTrans);

/*	for(i=0;i<=nbTrans;i++){
		printf("\n");
		for(j=0;j<tabDetect[i].nbSommet;j++){
			printf("%d ",tabDetect[i].Tableau[j]);
		}
	}*/

//	fclose(res);

//	if(addRoot) printf("<br>Racine : %d",n);
//	else  printf("<br>Racine : %d--%d",R1,R2);

	//recherche des aretes au dessus des sous-arbres

//	if(addRoot) racine = n; else racine = R1;

	racine = aTree.Root;


	for(i=0;i<=nbTrans;i++){
		max=INFINI;
		sommet=tabDetect[i].Tableau[0];
//		printf("<br> ARETE :");
		if(tabDetect[i].nbSommet > 1){
			// calcul de la distance entre la racine et le point d'intersection de tous les sommets
			for(j=1;j<tabDetect[i].nbSommet;j++){

				if(max > ( (aTree.ADD[racine][sommet] + aTree.ADD[racine][tabDetect[i].Tableau[j]] - aTree.ADD[sommet][tabDetect[i].Tableau[j]]) / 2.0 ) )
					max = (aTree.ADD[racine][sommet] + aTree.ADD[racine][tabDetect[i].Tableau[j]] - aTree.ADD[sommet][tabDetect[i].Tableau[j]]) / 2.0 ;
			}

			// recherche du sommet
			for(j=1;j<=nbArete;j++){
				if( ((fabs(max - aTree.ADD[racine][aTree.ARETE[2*j-1]]) < 0.00001) && (aTree.ADD[racine][aTree.ARETE[2*j-2]] < aTree.ADD[racine][aTree.ARETE[2*j-1]]) && fabs(aTree.ADD[racine][sommet] - aTree.ADD[racine][aTree.ARETE[2*j-1]] - aTree.ADD[sommet][aTree.ARETE[2*j-1]])< 0.00001) ||
					((fabs(max - aTree.ADD[racine][aTree.ARETE[2*j-2]]) < 0.00001) && (aTree.ADD[racine][aTree.ARETE[2*j-1]] < aTree.ADD[racine][aTree.ARETE[2*j-2]]) && fabs(aTree.ADD[racine][sommet] - aTree.ADD[racine][aTree.ARETE[2*j-2]] - aTree.ADD[sommet][aTree.ARETE[2*j-2]])< 0.00001)  ) {

					if( fabs(aTree.ADD[racine][sommet] - aTree.ADD[racine][aTree.ARETE[2*j-1]] - aTree.ADD[aTree.ARETE[2*j-1]][sommet]) < 0.00001 ){
						if(i % 2 == 0){
							outHGT[tmp].source_A = aTree.ARETE[2*j-1];
							outHGT[tmp].source_B = aTree.ARETE[2*j-2];
						}
						else{
							outHGT[tmp].dest_A = aTree.ARETE[2*j-1];
							outHGT[tmp].dest_B = aTree.ARETE[2*j-2];
							tmp++;
						}
	//					printf(" %d--%d (%d)",ARETE[2*j-1],ARETE[2*j-2],nbTrans);
					}
				}
			}
		}
		else{
			for(j=1;j<=nbArete;j++){
				if((aTree.ARETE[2*j-1] == sommet) || (aTree.ARETE[2*j-2] == sommet)){
					if(i % 2 == 0){
						outHGT[tmp].source_A = aTree.ARETE[2*j-1];
						outHGT[tmp].source_B = aTree.ARETE[2*j-2];
					}
					else{
						outHGT[tmp].dest_A = aTree.ARETE[2*j-1];
						outHGT[tmp].dest_B = aTree.ARETE[2*j-2];
						tmp++;
					}
	//				printf(" %d--%d (%d)",ARETE[2*j-1],ARETE[2*j-2],nbTrans);
				}

			}
		}
	}

	for(i=0;i<=cpt;i++){
		free(tabDetect[i].Tableau);
	}
	free(tabDetect);

	return nbHGTRet;
}
