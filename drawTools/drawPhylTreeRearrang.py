#!/usr/bin/env python3
#  -*- coding: utf-8 -*-


""" "Dessine l'arbre phylogenetique avec des informations dessus comme des evolutions de taux"

     usage:

        ./drawPhylTreeRearrang.py data84/PhylTree.conf -lengthFile=nbDup -colorFile=nbDup +landscape
        avec nbDup = nom de gene dupliqués par branche
        format du fichier:

        Pere    Fils    Valeur(nbDup)
"""

# Librairies
import sys
import math

from LibsDyogen import myFile, myMaths, myTools, myPhylTree, myPsOutput


# Arguments
arguments = myTools.checkArgs(
	[("phylTree.conf",file)],
	[("landscape",bool,False), ("printSpecies",bool,True),
	 ("printAncestors",bool,True), ("printAges",bool,False),
	 ("lengthFile",str,""), ("colorFile",str,""), ("funcLength",str,""),
	 ("funcColor",str,""), ("root",str,""), ("min",float,None),
	 ("max",float,None)],
	__doc__
)

phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

(largeur,hauteur) = myPsOutput.printPsHeader(landscape=arguments["landscape"])

root = arguments["root"] if arguments["root"] in phylTree.items else phylTree.root
funcLength = (lambda x, a: x) if arguments["funcLength"] == "" else eval(arguments["funcLength"])
#funcColor = (lambda x, a: x) if arguments["funcColor"] == "" else eval(arguments["funcColor"])

def loadVals(filename, func):

	vals = {}
	if filename != "":
		f = myFile.openFile(filename, "r")
		for l in f:
			#print >> sys.stderr, "the line", l
			t = l.replace('\n','').split("\t")
			s1 = sys.intern(t[0])
			s2 = sys.intern(t[1])
			assert (s1 == phylTree.parent[s2][0]) or (s2 == phylTree.parent[s1][0]), (s1,s2)
			#x = func(float(t[2]) / (3. if s2 in ["Sus scrofa", "Meleagris gallopavo"] else 1.) , float(abs(phylTree.ages[s1] - phylTree.ages[s2])))
			x = func(float(t[2]), float(abs(phylTree.ages[s1] - phylTree.ages[s2])))
			vals[(s1,s2)] = x
			vals[(s2,s1)] = x
			print('%s/%s (%s)' % (s1, s2, t[2]), file=sys.stderr)
		f.close()
	return vals

lengths = loadVals(arguments["lengthFile"], funcLength)

#print >> sys.stderr, lengths

#colors = loadVals(arguments["colorFile"], funcColor)

refcolors = [(0,0,127), (0,192,192), (0,192,0), (255,255,0), (242,148,0), (255,0,0)]
inter = myMaths.myInterpolator.getMultDim(myMaths.myInterpolator.oneDimCubic, list(range(len(refcolors))), refcolors)

y = 0
dy = hauteur / (len(phylTree.species[root])+1.)

margeX1 = dy
if arguments["printSpecies"]:
	margeX1 += max([len(x) for x in phylTree.species[root]]) * 0.15

margeX2 = dy
if arguments["printAncestors"]:
	margeX2 += len(root) * 0.15

newAge = {}
#val = []
def calcNewAge(node, a):
	print(node, a, file=sys.stderr)
	newAge[node] = a
	if node in phylTree.items:
		for (e,_) in phylTree.items[node]:
			#val.append(colors[(node,e)])
			calcNewAge(e, a+lengths[(node,e)])
calcNewAge(root, 0)

m = max(newAge.values())
for a in newAge:
	newAge[a] = m-newAge[a]
print(newAge, file=sys.stderr)

#minV = float(min(val))
#maxV = float(max(val))
#print >> sys.stderr, myMaths.myStats.txtSummary(val)

#if arguments["min"] != None:
#	minV = arguments["min"]
#if arguments["max"] != None:
#	maxV = arguments["max"]

# Le degrade 
def getColor(value):
	value = min(max((value-minV) / (maxV-minV), 0), 1)
	return inter(value * (len(refcolors)-1))


#for i in xrange(101):
#	col = inter((i * (len(refcolors)-1))/100.)
#	myPsOutput.drawBox(17.5+i/10., 20.5, 0.1, 0.2, col, col)
#myPsOutput.drawText(17.5, 20, str(minV), "black")
#myPsOutput.drawText(27, 20, str(maxV), "black")

dx = (largeur-margeX1-margeX2) / m

def printSubTree(node):
	global y
	#x = margeX1 + phylTree.ages[node]*dx
	x = largeur - (margeX1 + newAge[node]*dx)
	
	if node not in phylTree.items:
		y += dy
		if arguments["printSpecies"]:
			myPsOutput.drawText(x + .3, y-.1, node, "black")
		return y
	
	mi = hauteur
	ma = 0
	todo = []
	for (e,_) in phylTree.items[node]:
		tmpY = printSubTree(e)
		if tmpY > ma:
			ma = tmpY
		if tmpY < mi:
			mi = tmpY
	
		#todo.append( (tmpY, lengths[(node,e)], getColor(colors[(node,e)])) )
		todo.append( (tmpY, lengths[(node,e)], "black" ) )

	for (tmpY,a,color) in todo:
		if color == None:
			print("0.001 cm setlinewidth")
			color = "black"
		else:
			print("0.1 cm setlinewidth")
		myPsOutput.drawLine(x, tmpY, a*dx, 0, color)
		myPsOutput.drawLine(x, tmpY, 0, (mi+ma)/2.-tmpY, color)


	if arguments["printAncestors"]:
		myPsOutput.drawText(x - 1., (mi+ma)/2. + .1, node, "black")
	if arguments["printAges"]:
		myPsOutput.drawText(x - 1., (mi+ma)/2. - .4, "(%d)" % phylTree.ages[node], "black")

	return (mi+ma)/2.

printSubTree(root)

myPsOutput.printPsFooter()

