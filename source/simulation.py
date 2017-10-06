from Bio import Phylo
from io import StringIO
import os, re, shutil

groups = {"01":'1',"02":'1',"03":'1',"04":'1',"05":'1',"06":'1',"07":'1',
		"08":'2',"09":'2',"10":'2',"11":'2',"17":'2',"18":'2',"19":'2',
		"12":'3',"13":'3',"14":'3',"15":'3',"16":'3',"20":'3',"21":'3',"22":'3',"23":'3',
		"24":'4',"25":'4',"26":'4',"27":'4',"28":'4',"29":'4',"37":'4',"38":'4',
		"30":'5',"31":'5',"32":'5',"33":'5',"34":'5',"35":'5',"36":'5',
		"39":'6',"40":'6',"41":'6',
		"42":'7',"43":'7',"44":'7',"45":'7',"46":'7',"47":'7',"48":'7',"49":'7',"50":'7',"51":'7',"52":'7',"53":'7',"54":'7',
		"55":'8',"56":'8',"57":'8',"58":'8',"59":'8',"60":'8',"61":'8',"62":'8',"63":'8',"64":'8',"65":'8',
		"73":'9',"74":'9',"75":'9',"76":'9',"77":'9',"78":'9',"79":'9',
		"80":'10',"81":'10',"82":'10',"83":'10',"84":'10',
		"66":'11',"67":'11',"68":'11',"69":'11',"70":'11',
		"71":'12',"72":'12'}

def getGroupsFromWordTree(word_tree):
    tree = Phylo.read(StringIO(word_tree), "newick")
    dict_groups = {}
    groups_to_return = ""
    for leaf in tree.get_terminals():
        tmp = re.sub(r"-[0-9]*", "", leaf.name)
        if tmp in groups:
	        key = groups[tmp]
	        if key in dict_groups:
	            dict_groups[key] += "," + leaf.name
	        else:
	            dict_groups[key] = leaf.name
    for key,value in dict_groups.items():
        groups_to_return += key + "=" + value + "\n"
    return groups_to_return

def getWordNewickString(path,word,index):
    filename = path + word + "-input-" + index + ".txt"
    with open(filename) as f:
        f.readline()
        newick = f.readline()
    newick_formated = newick.replace("-0","")
    return newick_formated

def getWordTranslations(path,word,index):
    filename = path + word + "-trans-" + index + ".txt"
    translations = ""
    with open(filename) as f:
        for line in f:
            tab = re.split(" ",line)
            translations += tab[0].replace("-0","") + "=" + tab[2]
    return translations

def getMultipleLeaves(path,word,index):
    filename = path + word + "-trans-" + index + ".txt"
    duplicate = {}
    with open(filename) as f:
        for line in f:
            tab = re.split(" ",line)
            tmp = re.sub(r"-[0-9]*", "", tab[0])
            if tmp in duplicate:
                duplicate[tmp].append(tab[0].replace("-0",""))
            else:
                duplicate[tmp] = []
    return duplicate

def getLangueNewickString(filename):
    with open(filename) as f:
        langue_newick = f.read()
    return langue_newick

def updateTree(langue_tree,duplicate):
    for key in duplicate.keys():
        if len(duplicate[key]) > 0:
            chaine = "(" + key + ":1.0," + ":1.0,".join(duplicate[key]) + ":1.0)"
            langue_tree = langue_tree.replace(key+":",chaine+":")
    return langue_tree

def createFichier(path,file):
    tab = re.split('[-.]', file)
    word = tab[0]
    index = tab[2]
    word_tree = getWordNewickString(path,word,index)
    langue_tree = getLangueNewickString("/Users/boc_a/wbe-detection/examples/langue.new")
    tranlations = getWordTranslations(path,word,index)
    duplicate = getMultipleLeaves(path,word,index)
    langue_tree = updateTree(langue_tree,duplicate)
    groupes = getGroupsFromWordTree(word_tree)
    fh = open("input.txt", "w")
    #print(word,index)
    print("language_tree:", langue_tree, sep="\n", file=fh)
    print("word_tree:", word_tree, sep="\n", file=fh)
    print("group_content:", groupes, sep="\n", file=fh)
    print("translations:", tranlations, sep="\n", file=fh)
    fh.close()

def CompileAndSave(c1,c2,blk):
    image = str(c1) + "_" + str(c2) + ".png"
    resultats = str(c1) + "_" + str(c2) + "_" + str(blk) + ".txt"
    hwt = str(c1) + "_" + str(c2) +  "_" + str(blk) + "_hwt.txt"
    os.system("perl ../statistiques/compileResultAll.pl all_hgt.txt ../statistiques/all.mots ../statistiques/db_newick_custom.txt > db_resultats.txt")
    os.system("gnuplot heatmap.gp")
    shutil.copy('heatmap_biolinguistique.png', image)
    shutil.copy('db_resultats.txt', resultats)
    shutil.copy('all_hgt.txt', hwt)

def parcourir(path, c1=1505, c2=1.5, blk=20):
    current_word = ""
    cmd = "perl run_wbe.pl -inputfile=input.txt -c1=" + str(c1) + " -c2=" + str(c2) + " -blk=" + str(blk)
    dirs = os.listdir(path)
    prog = re.compile(".*-input-[0-9]+.txt$")
    os.remove("all_hgt.txt")

    for file in sorted(dirs):
    	result = prog.match(file)
    	if result:
            if current_word != file.split("-")[0]:
                current_word = file.split("-")[0]
                os.system("echo \"=> "+ current_word +"\" >> all_hgt.txt")
            print (file)
            createFichier(path,file)
            os.system(cmd)
            if os.path.isfile("output.txt"):
                os.system("echo \"1 cognat\" >> all_hgt.txt")
                os.system("grep \"^[0-9]\" output.txt >> all_hgt.txt")

for c2 in range( 15, 16, 5):
    c2 = float(c2) / 10.0
    for c1 in range(1505, 1600, 500):
        for blk in range(20, 21, 1):
            blk = float(blk) / 100.0
            #parcourir("../data/",c1,c2,blk)
            CompileAndSave(c1,c2,blk)


print("Fin normale du script")
