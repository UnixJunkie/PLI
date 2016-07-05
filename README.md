# Astex PLI #

PLI is a program developed at Astex Pharmaceuticals for scoring protein-ligand 
interactions. Interacting atoms are defined using a [modified][mcconkey] 
[Voronoi partitioning][voronoi] and interaction propensities (including geometric 
preferences) have been derived from over 30,000 [PDB][pdb] entries.

Much of the PLI methodology was published recently:

**Protein-Ligand-Informatics force field (PLIff): towards a fully knowledge driven “force field” for biomolecular interactions**

M.L. Verdonk, R.F. Ludlow, I. Giangreco, and P.C. Rathi

*J. Med. Chem., Just Accepted Manuscript*

DOI: [10.1021/acs.jmedchem.6b00716][doi]

We are actively developing the PLI code and will release an updated version in July 2016.

## Features ##
* Score protein-ligand interactions (for example docking/VS poses)
* Solvent Mapping
* Generate Interaction Maps for Donor/Acceptor/Lipophile etc.
* Rich atom-typing 

## Installation ##

### System Requirements ###

PLI is written in C and has been developed and tested on 64-bit Linux (CentOS) 
using the gcc compiler. It has no major external dependencies.

### Obtaining the Code ###

#### Download ####
* A zip of the current repository and tagged releases (as they become
  available) can be found at [bitbucket][downloads]. 

#### Via Git ####
* You can fork the repository on bitbucket (you will need an account)
or clone the repository from the command line
```bash
git clone https://bitbucket.org/AstexUK/pli.git
```

### Installing ###

Once you have obtained the code (either by extracting a snapshot or
cloning the repository:

* `cd` to the directory `<dir>`
* Type `make` in `<dir>`. This should compile the code and create a binary `<dir>/bin/pli`.
* Set the environment variable `PLI_DIR` to `<dir>`, e.g. for bash:
```bash
PLI_DIR=/path/to/pli
export PLI_DIR
```
* Try running one of the examples in the examples directory, e.g:
```bash
cd $PLI_DIR/examples
$PLI_DIR/bin/pli -settings $PLI_DIR/examples/1a9u_contacts.pli
```

## Modes ##
To see a list of modes and command-line options run
```sh
$PLI_DIR/bin/pli -help
```

## Get In Touch ##
If you have questions about PLI or have found it useful, let us know, on bitbucket or 
by email: "pli-at-astx-dot-com"

## Licence ##
Apache 2.0 License, see licence.txt for details. 

## Acknowledgements ##
* Brendan McConkey - for the capped and weighted Voronoi partitioning code (see the [original publication][mcconkey])
* John Burkardt - whose [SPHERE_LEBEDEV_RULE][lebedev] program was used to calculate lebedev points.
* The [PDB][pdb] - and the thousands of authors who have deposited X-ray structures of protein-ligand complexes.

[voronoi]: http://mathworld.wolfram.com/VoronoiDiagram.html
[pdb]: http://www.wwpdb.org/
[downloads]: https://bitbucket.org/AstexUK/pli/downloads
[mcconkey]: http://dx.doi.org/10.1093/bioinformatics/18.10.1365
[lebedev]: https://people.sc.fsu.edu/~jburkardt/c_src/sphere_lebedev_rule/sphere_lebedev_rule.html
[doi]: http://dx.doi.org/10.1021/acs.jmedchem.6b00716