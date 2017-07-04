import os, re

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

def getLangueNewickString(filename):
    with open(filename) as f:
        langue_newick = f.read()
    return langue_newick

def traitementFichier(path,file):
    tab = re.split('[-.]', file)
    word = tab[0]
    index = tab[2]
    word_tree = getWordNewickString(path,word,index)
    langue_tree = getLangueNewickString("/Users/boc_a/wbe-detection/examples/langue.new")
    tranlations = getWordTranslations(path,word,index)
    print(word,index)
    print("language_tree:",langue_tree, sep="\n")
    print("word_tree:",word_tree, sep="\n")
    print("translations:",tranlations, sep="\n")

def parcourir(path):
    dirs = os.listdir(path)
    prog = re.compile(".*-input-[0-9]+")
    for file in dirs:
        result = prog.match(file)
        if result:
            traitementFichier(path,file)

parcourir("/Users/boc_a/biolinguistique/data/")
