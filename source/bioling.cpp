bool sameGroup(char ** Names, int * source, int * destination){

  int noGroup = -1;
  int noGroup_tmp;
  char *chaine = (char*)malloc(100);
  
  for(int m=1;m<=destination[0];m++){
	  strcpy(chaine,Names[destination[m]]);
    noGroup_tmp = groupLang[atoi(strsep(&chaine,"-"))];
    if(noGroup != -1){
      if (noGroup != noGroup_tmp){
        return false;
      }
    }
    noGroup = noGroup_tmp;
	} 
  
  for(int m=1;m<=source[0];m++){
	  strcpy(chaine,Names[source[m]]);
    noGroup_tmp = groupLang[atoi(strsep(&chaine,"-"))];
    
    if (noGroup != noGroup_tmp){
      return false;
    }
    noGroup = noGroup_tmp;
	}

  return true;
}

int groupeUnique( struct DescTree *DT , struct InputTree tree, int root){
  
  int noGroup=-1,noGroup_tmp;
  char *chaine = (char*)malloc(100);
  printf("\n%d : ", root);
  for(int i=1;i<=DT[root].nbSommet;i++){
    strcpy(chaine,tree.SpeciesName[DT[root].Tableau[i]]);
    noGroup_tmp = groupLang[atoi(strsep(&chaine,"-"))];
	  printf(" %s[%d]", tree.SpeciesName[DT[root].Tableau[i]],noGroup_tmp);
    if( (noGroup_tmp != noGroup) && (noGroup != -1) ){
      return -1;
    }
    noGroup = noGroup_tmp;
  }
  return noGroup;
}

int getInternalRootNode( struct InputTree tree ){
    //== recherche de la racine interne de ROOT
    int sRoot=-1;
    for(int i=1;i<=2*tree.size-3-tree.kt;i++){
      if( tree.ARETE[2*i-1] == tree.Root) sRoot =  tree.ARETE[2*i-2];
      if( tree.ARETE[2*i-2] == tree.Root) sRoot =  tree.ARETE[2*i-1];
    }
    return sRoot;
}    

int getsubTreeRoot ( struct InputTree tree, int sRoot ){
    //== recherche de la racine du premier sous-arbre
    int str=-1;
    for(int i=1;i<=2*tree.size-3-tree.kt;i++){
      if( (tree.ARETE[2*i-1] == sRoot) &&  (tree.ARETE[2*i-2] !=  tree.Root) ) {str =  tree.ARETE[2*i-2]; break;}
      if( (tree.ARETE[2*i-2] == sRoot) &&  (tree.ARETE[2*i-1] !=  tree.Root) ) {str =  tree.ARETE[2*i-1]; break;}
    }
    return str;
}

int getsubTreeRoot ( struct InputTree tree, int sRoot, int str1 ){
    //== recherche de la racine du deuxieme sous-arbre
    int str=-1;
    for(int i=1;i<=2*tree.size-3-tree.kt;i++){
      if( (tree.ARETE[2*i-1] == sRoot) &&  (tree.ARETE[2*i-2] !=  tree.Root) &&  (tree.ARETE[2*i-2] !=  str1) ) {str =  tree.ARETE[2*i-2]; break;}
      if( (tree.ARETE[2*i-2] == sRoot) &&  (tree.ARETE[2*i-1] !=  tree.Root) &&  (tree.ARETE[2*i-1] !=  str1) ) {str =  tree.ARETE[2*i-1]; break;}
    }
    return str;
}

void saveThisHGT(struct DescTree *DT, struct InputTree tree, int saRoot1, int saRoot2){
  
  FILE *f = fopen("hgtplus.txt","w+");

  for(int i=1;i<=DT[saRoot1].nbSommet;i++){
	  fprintf(f,"%s ", tree.SpeciesName[DT[saRoot1].Tableau[i]]);
  }
	fprintf(f,"\n");
  for(int i=1;i<=DT[saRoot2].nbSommet;i++){
	  fprintf(f,"%s ", tree.SpeciesName[DT[saRoot2].Tableau[i]]);
  }
  
  fclose(f);
}

