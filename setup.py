import setuptools


try:
    with open("README.md", "r") as fh:
        long_description = fh.read()
except FileNotFoundError:
    long_description = ''


setuptools.setup(
    name="ToolsDyogen_py3",
    version="0.0.1",
    author="Matthieu MUFFATO, Alexandra LOUIS, Amélie PÉRÈS, Joseph LUCAS, Guillaume LOUVEL",
    author_email="guillaume.louvel@normalesup.org",
    description="Python 3 version of the Dyogen team scripts for comparative genomics. http://www.ibens.ens.fr/?rubrique43&lang=en",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DyogenIBENS/ToolsDyogen_py3",
    packages=setuptools.find_packages(),
    python_requires='>=3',
    install_requires=['LibsDyogen_py3'], #['enum34;python_version<"3.4"'],
    extras_require={'GFF': ['gff3'],  # Only in genomesTools/misc/convertGffToGenesST.py
                    'ODF': ['odfpy'], # Only in bin/statsOnGenesInGenomes.py,
                    'numpy': ['numpy'], # For forest_summary
                    'autoCLI': ['UItools']}, # For forest. TODO: remove
    #optional dependencies (handled in try-except)
                      #'java.lang',
                      #'psyco',  # development stopped in 2011
    scripts=[
    #    'ensemblTools/API_perl_scripts/align_cds.pl',
    #    'ensemblTools/API_perl_scripts/getChromosomeSizeForAllGenome.pl',
    #    'ensemblTools/API_perl_scripts/getChromosomeSize.pl',
    #    'ensemblTools/API_perl_scripts/getContigsSize.pl',
    #    'ensemblTools/API_perl_scripts/getDBnameFromEnsembl.pl',
    #    'ensemblTools/API_perl_scripts/downloadSpecificProteinTrees.pl',
    #    'ensemblTools/API_perl_scripts/IDmapper.pl'
    'ensemblTools/ENSEMBL.downloadProteinTreesWithSplitNodes.pl',
    'genomesTools/cmpIntervals.py',
    'genomesTools/formatTabularAncGenome.py',
    'treeTools/getInfoOnSpeciesTree.py',
    'treeTools/getSpeciesList.py',
    'statsTools/printStats.py',
    'drawTools/drawPhylTreeRearrang.py',
    'evolTools/calcDiagEvolution.py',
    'evolTools/calcGeneEvolution.py',
    'evolTools/calcGeneEvolutionOnEachBranch.py',
    'evolTools/extractGeneEvents.py',
    'evolTools/getGeneTimeline.py'
        ],
    entry_points = {
        'console_scripts': [
    #'ALL.extractStrongFamiliesFromAncGenes=ancGenesTools.ALL.extractStrongFamiliesFromAncGenes:main',
    #'ENSEMBL.biomartQueryAll:ensemblTools.ENSEMBL.biomartQueryAll:main',
    #'ENSEMBL.biomartQuery:ensemblTools.ENSEMBL.biomartQuery:main',
    #'ENSEMBL.downloadProteinTrees:ensemblTools.ENSEMBL.downloadProteinTrees:main',
    'ALL.convertNewickTree=treeTools.ALL.convertNewickTree:main',
    'ALL.cutTree=treeTools.ALL.cutTree:main',
    'ALL.extractGeneFamilies=treeTools.ALL.extractGeneFamilies:main',
    'ALL.extractNewickTrees=treeTools.ALL.extractNewickTrees:main',
    'ALL.extractOneGeneTree=treeTools.ALL.extractOneGeneTree:main',
    'ALL.extractMultipleGeneTrees=treeTools.ALL.extractMultipleGeneTrees:main',
    'ALL.getGeneHistory=treeTools.ALL.getGeneHistory:main',
    'ALL.infoProtTree=treeTools.ALL.infoProtTree:main',
    'ALL.statsKaryotype=genomesTools.ALL.statsKaryotype:main',
    'ALL.calcRearrangRates=evolTools.ALL.calcRearrangRates:main',
    'ALL.extractBranchLength=evolTools.ALL.extractBranchLength:main',
    'ENSEMBL.buildProteinTrees=treeTools.ENSEMBL.buildProteinTrees:main',
    'ENSEMBL.buildProteinTreesV1=treeTools.ENSEMBL.buildProteinTreesV1:main',
    'ENSEMBL.markLowScoreDup=treeTools.ENSEMBL.markLowScoreDup:main',
    'misc.compareGenomes=compareTools.misc.compareGenomes:main',
    'misc.convertContigsToGenome=genomesTools.misc.convertContigsToGenome:main',
    'misc.convertGffToGenesST=genomesTools.misc.convertGffToGenesST:main',
    'misc.printAllDescendants=treeTools.misc.printAllDescendants:main',
    'misc.printLastCommonAncestors=treeTools.misc.printLastCommonAncestors:main',
    'statsEvents.mkODS=statsTools.statsEvents.mkODS:main',
    'stats.getNbComparisons=statsTools.stats.getNbComparisons:main',
    'stats.mkODS=statsTools.stats.mkODS:main',
    'stats.mkODS_withSummary=statsTools.stats.mkODS_withSummary:main',
    'forest_summary=statsTools.forest_summary:main',
    'ALL.getAges=treeTools.getAges:main',
    'phyltree=treeTools.phyltree:main',
    'forest=treeTools.forest'
        ]
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        #"Operating System :: Microsoft :: Windows",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    package_data={'ensemblTools': ['*.pl']},
    include_package_data=True,
    zip_safe=False
)
