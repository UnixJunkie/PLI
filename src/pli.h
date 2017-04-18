// Copyright 2015 Astex Therapeutics Ltd.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <ctype.h>

typedef struct Settings SETTINGS;
typedef struct Element ELEMENT;

struct Element {
  int id;
  char name[10];
  double cov_radius;
  double vdw_radius;
  unsigned int flags;
  int default_valence;
  int n_lp;
};

#define YES 1
#define NO 0
#define TRUE 1
#define FALSE 0
#define UNDEFINED -1
#define SET_BIT 1
#define UNSET_BIT 0

#define FIRST 0
#define UNIQUE 1
#define ALL 2
#define BEST -1
#define RANDOM -2
#define ORIGINAL -3
#define CURRENT -4

#define MAX_QUERY_ATOMS 20

// output bitmasks
#define OUTPUT_NONE 0
#define OUTPUT_SYSTEM 1
#define OUTPUT_ATOMS 2
#define OUTPUT_CONTACTS 4
#define OUTPUT_GEOMETRIES 8
#define OUTPUT_SCORES 16
#define OUTPUT_PDB 32
#define OUTPUT_MOL 64
#define OUTPUT_AXES 128
#define OUTPUT_VPTS 256
#define OUTPUT_COVALENT 512
#define OUTPUT_INTRAMOLECULAR 1024
#define OUTPUT_HBONDS 2048
#define OUTPUT_CLASHES 4096
#define OUTPUT_SAS_STATS 8192
#define OUTPUT_LIGAND 16384
#define OUTPUT_SIMB 32768
#define OUTPUT_STYPES 65536
#define OUTPUT_SITE_STATS 131072
#define OUTPUT_TRIPLETS 262144  

#define COVALENT_TOLERANCE 0.4
#define MAX_LINE_LEN 3000
#define sqr(x) ((x)*(x))
#define PI 3.1415926536
#define GOLDEN_RATIO 0.618034
#define MAX_RING_SIZE 500
#define DEFAULT_VDW_RADIUS 1.70
#define RT 0.59248490241
#define GAS_CONSTANT .0019872041
#define MAX_ATOM_SCORES 10
#define MAX_CONTACT_SCORES 10
#define MAX_SYSTEM_SCORES 100
#define MAX_ATOM_EDGES 12

#define MIN_COORD_SHAKE 0.0000000001
#define MAX_COORD_SHAKE 0.0000000009



// common element atomic numbers
#define HYDROGEN 1
#define CARBON 6
#define NITROGEN 7
#define OXYGEN 8
#define FLUORINE 9
#define SULFUR 16
#define CHLORINE 17
#define BROMINE 35
#define IODINE 53

// hbond geometries:
#define POINT_HBOND_GEOMETRY 0
#define VECTOR_HBOND_GEOMETRY 1
#define VECTORPAIR_HBOND_GEOMETRY 2
#define WEDGE_HBOND_GEOMETRY 3

// element bitmasks
#define METAL_ELEMENT 1

// atom flag bitmasks
#define WATER_OXYGEN 1
#define H_ATOMS_ADDED 2
#define AROMATIC_ATOM 4
#define NODE_SET_ATTEMPTED 8
#define ATOM_TYPE_ATTEMPTED 16
#define AMINO_ACID_ATOM 32
#define ATOM_VDW_RADIUS_SET 64
#define SINGLE_ATOM_PROBE 128
#define SKIP_ATOM 256
#define SELECTED_ATOM 512
#define MATCHED_ATOM 1024
#define RING_ATOM 2048

// atom status bitmasks
#define ATOM_UNCHANGED 0
#define ATOM_SCORED 0
#define ATOM_MOVED 1
#define ATOM_ROTATED 2
#define ATOM_CHANGED 4
#define ATOM_MOVED_OR_CHANGED 5
#define ATOM_NEW 7

// atom cstatus bitmasks
#define ATOM_NOTHING_CALCULATED 0
#define ATOM_CONTACTS_CALCULATED 1
#define ATOM_AREAS_CALCULATED 2
#define ATOM_GEOMETRIES_CALCULATED 4
#define ATOM_SCORES_CALCULATED 8

// atom type bitmasks
#define HBOND_DONOR_ATOM_TYPE 1
#define HBOND_ACCEPTOR_ATOM_TYPE 2
#define HBOND_DA_ATOM_TYPE 3
#define CHBOND_DONOR_ATOM_TYPE 4
#define ANY_HBOND_ATOM_TYPE 7
#define AMINO_ACID_ATOM_TYPE 8
#define METAL_ATOM_TYPE 16
#define METAL_ACCEPTOR_ATOM_TYPE 32
#define AROMATIC_ATOM_TYPE 64
#define LIPOPHILIC_ATOM_TYPE 128

// bond and bond angle bitmasks
#define SKIP_BOND 1
#define SKIP_BOND_ANGLE 1

// residue bitmasks
#define AMINO_ACID 1

// molecule bitmasks
#define LIGAND_MOLECULE 1
#define PROTEIN_MOLECULE 2
#define SYMMETRY_MOLECULE 4
#define ANY_MOLECULE 7


// atom error flag bit masks
#define ATOM_GEOMETRY_ERROR 1
#define ATOM_DICTIONARY_ERROR 2
#define ATOM_SKIP_DICTIONARY 4
#define ATOM_VORONOI_ERROR 8
#define ATOM_TYPING_ERROR 16
#define ATOM_HBOND_MISMATCH 32

#define planedef 'X'

// contact bitmasks
#define COVALENT_CONTACT 1
#define SECONDARY_CONTACT 2
#define NO_VORONOI_CONTACT 4
#define INTRAMOLECULAR_CONTACT 8
#define METAL_COORDINATING_CONTACT 16
#define HBOND_CONTACT 32

// atom sets:
#define PROTEIN_ATOM_SET 1
#define LIGAND_ATOM_SET 2
#define WATER_ATOM_SET 3
#define OTHER_ATOM_SET 4

// maximum coordination number:
#define MAX_Z 50

// maximum number of alternate atom types (for ambiguous types):
#define MAX_ALT_TYPES 10

// maximum number of virtual points on an atom:
#define MAX_VIRTUAL_POINTS 2

// extrapolation methods for histograms:
#define EXTRAPOLATE_TO_SMALLEST 0
#define EXTRAPOLATE_TO_NEAREST 1

// normalisation methods for histograms:
#define ABSOLUTE_NORMALISATION 0
#define RELATIVE_NORMALISATION 1

// geometric corrections for histograms:
#define NO_GEOMETRIC_CORRECTION 0
#define CONICAL_CORRECTION 1
#define PLIFF_R_CORRECTION 2

#define BULK_WATER_TD 2.8

// definitions for degrees of freedom:
#define DOF_LIGAND_TRANSLATION 1
#define DOF_LIGAND_ROTATION 2
#define DOF_LIGAND_TORSIONS 4
#define DOF_LIGAND_RIGID_BODY 3
#define DOF_LIGAND 7
#define DOF_LIGAND_ATOMS 8
#define DOF_LIGAND_MATCHED_ATOMS 16
#define DOF_WATER_TRANSLATION 32
#define DOF_WATER 32

// map mask types:
#define CONTACT_MASK 1
#define ATOM_MASK 2
#define VDW_MASK 3
#define CUSTOM_MASK 4
#define FIXED_MASK 4
#define PROBE_MASK 5

// atom matching
#define MAX_MATCH_ATOMS 100



enum MINIMISATION_ALGORITHM { STEEPEST_DESCENT,CONJUGATE_GRADIENT };

enum SCORE_ORDER {ASCENDING, DESCENDING};

// io definitions:

enum FILE_FORMAT { GZIPPED_FILE,ASCII_FILE };

enum OUTPUT_FORMAT { FORMATTED,JSON,CSV };

typedef unsigned char BYTE;
typedef unsigned char BOOLEAN;
typedef unsigned int FLAG;

typedef struct PLIFile {
  char filename[MAX_LINE_LEN];
  FILE *file;
  enum FILE_FORMAT format;
} PLI_FILE;

typedef struct OStyle {
  char name[20];
  char separator[5];
  char open[2];
  char close[2];
  char format[20];
} OSTYLE;

typedef void (*OBJ2TEXT)(void*,char*,char*);

typedef struct OField {
  char name[20];
  char uformat[50];
  char fformat[50];
  char format[50];
  OSTYLE *ostyle;
  OBJ2TEXT obj2text;
} OFIELD;

typedef struct Ogroup {
  char name[20];
  OFIELD *ofields;
  char fields[200];
} OGROUP;

typedef struct List {
  char name[20];
  int n_items;
  int n_alloc_items;
  int item_size;
  void *items;
} LIST;

extern PLI_FILE *PLI_STDIN;
extern PLI_FILE *PLI_STDOUT;
extern PLI_FILE *PLI_STDERR;

void obj2text(void*,char*,OFIELD*);
OSTYLE* get_ostyle(char*);
OGROUP* get_ogroup(OGROUP*,char*);
OFIELD* get_ofield(OFIELD*,char*);
void init_ofields(OFIELD*);
void set_ofields(OFIELD*,char*,char*);
void set_ofield(OFIELD*,char*);
void set_ogroups(OGROUP*,char*,char*);

// map definitions:

typedef struct Grid {
  int npoints[3]; /* number of grid points along each of the 3 axes */
  int limit[3][2]; /* limits in grid points (including sign) */
  double flimit[3][2]; /* real box limits */
  double spacing;
  double padding;
} GRID;

typedef struct Maptype {
  char name[20];
  int item_size;
  void (*init)(void*);
  void (*free)(void*);
} MAP_TYPE;

typedef struct Map {
  char name[MAX_LINE_LEN];
  MAP_TYPE *type;
  GRID *grid;
  void ***matrix;
  BYTE *mask;
} MAP;

typedef struct MapIsland {
  int id;
  int *points;
  int n_points;
  int n_frag_points;
  int n_alloc_points;
  int *edge_points;
  int n_edge_points;
  int n_alloc_edge_points;
  double maxY;
  double sumY;
  double avgY;
  double integral;
  double frag_integral;
  double volume;
  int ix;
  int iy;
  int iz;
  struct AtomList *site_atoms;
} MAP_ISLAND;

typedef struct AtomCoords {
  double position[4];
  double u[4];
  double v[4];
  double w[4];
} ATOM_COORDS;

typedef struct AtomEdge {
  struct Atom *atom;
  struct Bond *bond;
} ATOM_EDGE;

typedef struct Atom {
  int unique_id;
  int id;
  int seqid;
  ELEMENT *element;
  char name[20];
  char alt_name[10];
  char subname[10];
  int subid;
  char altloc;
  char chain;
  char icode;
  double position[4];
  double original_position[4];
  double *tether_position;
  double occupancy;
  double bfactor;
  double tdf;
  double td;
  double lnPu;
  double vdw_radius;
  double vdw_radius_H2O;
  double tether_k;
  struct AtomNode *node;
  struct AtomType *type;
  struct AtomType *ambiguous_type;
  struct AtomList *connections;
  struct AtomList *h_connections;
  struct AtomGeometry *geometry;
  struct HBondGeometry *hbond_geometry;
  struct AtomCoordination *coordination;
  struct Ring *ring;
  double type_probability;
  double u[4];
  double v[4];
  double w[4];
  double vpts[MAX_VIRTUAL_POINTS][4];
  struct Map *gridmap;
  struct AtomList *gridlist;
  struct ContactList *contactlist;
  double exposed_area;
  double contact_area;
  double intra_area;
  double da_area;
  double don_area;
  double ch_area;
  double acc_area;
  double met_area;
  double pol_area;
  double lip_area;
  double covalent_area;
  double score;
  double scores[MAX_ATOM_SCORES];
  double ref_score;
  double ref_scores[MAX_ATOM_SCORES];
  double constraint_score;
  struct Molecule *molecule;
  int set;
  int planar;
  int n_hydrogens;
  int group_status;
  int smarts_ring_label;
  unsigned int cstatus;
  unsigned int status;
  unsigned int flags;
  unsigned int error_flags;
  //
  struct TripletList *tripletlist;
  int td_converged;
  int type_id_file;
  struct AtomQuery *query;
  int n_edges;
  struct AtomEdge edges[MAX_ATOM_EDGES];
  int hybridisation;
  struct FeatureGroup *fgroup;
  MAP *depth_map;
  LIST *hydrogens;
  LIST *substructure;
} ATOM;

typedef struct AtomQuery {
  ATOM *atom;
  int ring_label;
  int aromatic;
  int hybridisation;
  int planar;
  int linear;
  int ring_atom;
  int n_edges;
  int n_hydrogens;
  int any_element;
  int n_elids;
  int elids[10];
  int mol_type;
  int n_resnames;
  char resnames[10][4];
  struct FeatureGroupType *fgtype;
} ATOM_QUERY;

struct AtomList;
struct System;
struct Torsion;

typedef struct Bond {
  int type;
  double distance;
  double Do;
  double k;
  struct Atom *atom1;
  struct Atom *atom2;
  struct Torsion *torsion;
  unsigned int flags;
} BOND;

typedef struct BondAngle {
  double angle;
  double Ao;
  double k;
  struct Atom *atom1;
  struct Atom *atom2;
  struct Atom *atom3;
  unsigned int flags;
} BOND_ANGLE;

typedef struct TorsionType {
  char name[MAX_LINE_LEN];
  int united;
  int implicit14;
  int hyb1;
  int hyb2;
  double A;
  double n;
  double Xo;
  int n_terms;
} TORSION_TYPE;

typedef struct Torsion {
  double angle;
  BOND *bond;
  ATOM *atom1;
  ATOM *atom2;
  ATOM *batom1;
  ATOM *batom2;
  struct AtomList *alist1;
  struct AtomList *alist2;
  struct AtomList *rotlist;
  struct TorsionType *type;
} TORSION;

typedef struct Molecule {
  int loaded;
  int lazy_load;
  char filename[MAX_LINE_LEN];
  PLI_FILE *file;
  long int file_pos;
  char name[MAX_LINE_LEN];
  int natoms;			/* the number of atoms in the molecule */
  int nbonds;			/* the number of bonds in the molecule */
  int n_bond_angles;		/* the number of bond angles in the molecule */
  int n_torsions;               /* the number of bond torsions in the molecule */
  ATOM* atom;			/* array of atoms (see struct ATOM) */
  BOND* bond;			/* array of bonds (see struct BOND) */
  BOND_ANGLE* bond_angles;	/* array of bond angle (see struct BOND_ANGLE) */
  TORSION *torsions;            /* array of bond torsions */
  int **nb_atom_pairs;          /* matrix describing atom pairs to be checked for nonbonded contacts */
  int n_alloc_atoms;
  int n_alloc_bonds;
  int n_alloc_bond_angles;
  MAP *covalent_map;
  MAP *contacts_map;
  MAP *conmap;
  MAP *depth_map;
  LIST *atom_ids;
  struct AtomList *selection;
  struct System *molsystem;
  int use_pdb_dict;
  int use_hydrogens;
  int use_grids;
  int rely_on_bonds;
  int prepped;
  int n_active_waters;
  int n_alloc_active_waters;
  ATOM* active_waters;
  unsigned int flags;
  double score;
  struct System *system;
  double geometric_center[4];
  double center_of_mass[4];
  double moment_of_inertia[4];
  LIST *active_atoms;
  LIST *nonh_atoms;
  LIST *hydrogens;
  LIST *hydrogen_lists;
  LIST *fgroups;
  LIST *feature_lists;
  LIST *features;
  LIST *substructures;
  LIST *subatoms;
  LIST *atom_queries;
} MOLECULE;

typedef struct Complex {
  char protein_file[MAX_LINE_LEN];
  char ligand_file[MAX_LINE_LEN];
  MOLECULE *protein;
  MOLECULE *ligand;
} COMPLEX;

typedef struct ComplexList {
  int n_complexes;
  int n_alloc_complexes;
  COMPLEX *complexes;
  PLI_FILE *file;
} COMPLEX_LIST;

typedef struct MoleculeList {
  int n_molecules;
  int n_alloc_molecules;
  MOLECULE **molecules;
  PLI_FILE *file;
} MOLECULE_LIST;

typedef struct AtomList {
  int natoms;
  int n_alloc_atoms;
  ATOM **atom;
  unsigned int *flags;
  ATOM_COORDS *coords;
} ATOMLIST;

typedef struct BondList {
  int nbonds;
  int n_alloc_bonds;
  BOND **bond;
} BONDLIST;

typedef struct Ring {
  int aromatic;
  int size;
  ATOM *atom[MAX_RING_SIZE];
  double center[4];
  double normal[4];
} RING;

typedef struct RingList {
  int n_rings;
  int n_alloc_rings;
  RING **rings;
} RING_LIST;

typedef struct AtomGeometry {
  int id;
  char name[40];
  int nbonds1;
  int nbonds2;
  int planar1;
  int planar2;
  int u_axis;
  int v_axis;
  int u_symm;
  int v_symm;
  double r_limits[2];
  double alpha_limits[2];
  double beta_limits[2];
} ATOM_GEOMETRY;

typedef struct HBondGeometry {
  int type;
  char name[40];
  int npts;
  double angle1;
  double tol11;
  double tol12;
  double angle2;
  double tol21;
  double tol22;
} HBOND_GEOMETRY;

typedef struct Residue {
  char name[10];
  int id;
  char chain;
  char icode;
  char single_letter_code;
  unsigned int flags;
} RESIDUE;

typedef struct ResidueAtom {
  int id;
  char subname[10];
  char name[10];
  char atom_type_name[20];
} RESIDUE_ATOM;

typedef struct UnitedAtom {
  int id;
  char name[10];
  char atom_node_name[10];
  struct AtomNode *atom_node;
  int hybridisation;
  int n_hydrogens;
  int formal_charge;
  double vdw_radius;
  double ideal_hbond_alpha;
  double ideal_hbond_beta;
} UNITED_ATOM;

typedef struct AtomTypeDef {
  char node_name[20];
  char element_name[2];
  struct AtomType *type;
  struct AtomNode *node;
  struct Element *element;
  int hybridisation;
  int aromatic;
  int planar;
  int n_alloc_defs;
  int n_defs;
  struct AtomTypeDef *defs;
} ATOM_TYPE_DEF;

struct AtomType;

typedef struct AltAtomType {
  struct AtomType *type;
  int group_status;
  int assigned_group_status;
  double probability;
} ALT_ATOM_TYPE;

typedef struct AtomType {
  int i;
  int id;
  char name[50];
  char united_atom_name[10];
  char atom_geometry_name[40];
  char hbond_geometry_name[40];
  UNITED_ATOM *united_atom;
  ELEMENT *element;
  ATOM_GEOMETRY *geometry;
  HBOND_GEOMETRY *hbond_geometry;
  int n_hs;
  int n_lps;
  int n_alt_types;
  ALT_ATOM_TYPE alt_types[MAX_ALT_TYPES];
  unsigned int flags;
} ATOM_TYPE;

typedef struct AtomTypingScheme {
  int n_alloc_elements;
  int n_elements;
  struct Element *elements;
  int n_alloc_atom_nodes;
  int n_atom_nodes;
  struct AtomNode *atom_nodes;
  int n_alloc_united_atoms;
  int n_united_atoms;
  struct UnitedAtom *united_atoms;
  int n_alloc_atom_geometries;
  int n_atom_geometries;
  struct AtomGeometry *atom_geometries;
  int n_alloc_atom_types;
  int n_atom_types;
  struct AtomType *atom_types;
  int n_alloc_residue_atoms;
  int n_residue_atoms;
  struct ResidueAtom *residue_atoms;
  int n_alloc_atom_type_defs;
  int n_atom_type_defs;
  struct AtomTypeDef *atom_type_defs;
} ATOM_TYPING_SCHEME;

typedef struct AtomNode {
  char name[20];
  char element_name[2];
  struct Element *element;
  int n_single;
  int n_double;
  int n_triple;
  int hybridisation;
  int formal_charge;
} ATOM_NODE;

typedef struct Contact {
  ATOM *atom1;
  ATOM *atom2;
  double distance;
  double alpha1;
  double beta1;
  double alpha2;
  double beta2;
  double area;
  double iarea;
  double contact_score;
  double r_score;
  double geometry_score;
  double score;
  double scores[MAX_CONTACT_SCORES];
  struct Polygon *poly;
  unsigned int flags;
} CONTACT;

typedef struct ContactList {
  int ncontacts;
  int n_alloc_contacts;
  CONTACT *contacts;
  struct Polyhedron *poly;
} CONTACTLIST;

typedef struct Coordination {
  CONTACT *contact;
  double s;
  double f;
} COORDINATION;

typedef struct AtomCoordination {
  int n_don;
  int n_acc;
  int n_da;
  int n_ch;
  int n_metal;
  int z;
  int max_z;
  COORDINATION list[MAX_Z];
  double intra_area;
  double inter_area;
  double exposed_polar_area;
  unsigned int error_flags;
} ATOM_COORDINATION;

typedef struct Polyhedron {
  int nvert;
  struct vertex *vert;
} POLYHEDRON;

typedef struct Polygon {
  int nvert;
  struct vertex **vert;
} POLYGON;

// Structures for Vcontacts code (McConkey et al):

typedef struct plane {
  double Ai[4];      // parameters A,B,C,D of contact plane Ax+By+Cz+D=0
  double dist;       // distance from plane to origin
  int    index;      // index to which record in PDB or ligand array.
  double area;       // contact area in square angstroms
  char   flag;       // 'X' if no contact, 'E' if an engulfing atom.
  struct Atom *atom; // MLV added - pointer to atom
} PLANE;

typedef struct vertex {
  double xi[3];      // x,y,z coordinates (x1,x2,x3)
  double dist;       // distance to origin
  int    plane[3];   // identification of intersecting planes. -1 = sphere.
} VERTEX;

typedef struct ptindex {
  int  numpts;       // number of points defining face
  int  pt[40];       // index to polyhedron points
} PTINDEX;

typedef struct edgevector {
  double V[3];       // vector for edge
  int startpt;       // initial vertex
  int endpt;         // final vertex
  int plane[2];      // planes defining edge (using single atoms for now)
  int startplane;    // third plane at start point
  int endplane;      // third plane at end point
  char arc;          // flag for arc point calculations
} EDGEVECTOR;

struct Settings;
struct System;
struct DegreeOfFreedom;

typedef struct PliMode {
  char name[MAX_LINE_LEN];
  void (*run)(struct Settings*);
  char selection[MAX_LINE_LEN];
  int selection_deviation_warning;
  int allow_pre_minimise;
} PLI_MODE;

typedef struct PocketMode {
  char name[MAX_LINE_LEN];
  struct Map* (*map)(struct System*);
} POCKET_MODE;

typedef struct PliSfunc {
  char name[MAX_LINE_LEN];
  void (*init)(void);
  double (*score_system)(struct System*);
  double (*score_gradient)(struct DegreeOfFreedom*,struct System*);
  void (*minimise_system)(struct System*);
  void (*write_system_scores)(PLI_FILE*,struct System*);
  void (*write_atom_scores)(PLI_FILE*,struct Atom*);
  void (*write_contact_scores)(PLI_FILE*,struct Contact*);
  double bad_score;
  int order;
} PLI_SFUNC;

typedef struct SysMols {
  MOLECULE *protein;
  MOLECULE *water;
  MOLECULE *symmetry;
  MOLECULE *ligand;
  MOLECULE *site_ligand;
  MOLECULE *template;
  MOLECULE_LIST *ligand_list;
  COMPLEX_LIST *complex_list;
} SYSMOLS;

struct Settings {
  char *pli_dir;
  int verbose;
  PLI_MODE *mode;
  PLI_SFUNC *sfunc;
  char jobname[MAX_LINE_LEN];
  int use_fixed_random_seeds;
  char selection[100];
  char settings_file[100];
  char types_file[1000];
  char protein_file[1000];
  char ligand_file[1000];
  char symmetry_file[1000];
  char oatoms[MAX_LINE_LEN];
  enum OUTPUT_FORMAT oformat;
  unsigned long int oflags;
  int protein_use_pdb_dict;
  int ligand_use_pdb_dict;
  int protein_use_hydrogens;
  int ligand_use_hydrogens;
  char resolve_side_chains[20];
  int resolve_metals;
  char alt_type_assignment[20];
  int report_alt_types;
  struct AtomTypingScheme *atom_typing_scheme;
  struct ForceField *force_field;
  double water_vdw_radius;
  double max_covalent_dist;
  double max_contact_dist;
  int id1,id2,eid1,eid2;
  char hist_name[20];
  SYSMOLS *sysmols;
  LIST *site_atoms;
  int minimise;
  int allow_bad_clashes;
  int allow_covalent_bonds;
  char keep_waters[200];
  unsigned int dof;
};

typedef struct Match {
  int n_atoms;
  unsigned char flags;
  int score;
  ATOM *atoms1[MAX_MATCH_ATOMS];
  ATOM *atoms2[MAX_MATCH_ATOMS];
} MATCH;

typedef struct System {
  int id;
  int prepped;
  int score_prepped;
  char name[MAX_LINE_LEN];
  struct Settings *settings;
  MOLECULE *protein;
  MOLECULE *water;
  MOLECULE *ligand;
  MOLECULE *site_ligand;
  MOLECULE *template;
  MOLECULE *symmetry;
  MOLECULE_LIST *molecule_list;
  ATOMLIST *selection;
  LIST *active_atoms;
  ATOM_COORDS *coords;
  MATCH *template_match;
  double score;
  double scores[MAX_SYSTEM_SCORES];
  double ref_score;
  double ref_scores[MAX_SYSTEM_SCORES];
  double constraint_score;
  double min_dx;
  double min_ds;
  unsigned int dof;
} SYSTEM;

typedef struct SystemList {
  int n_systems;
  SYSTEM *systems;
} SYSTEM_LIST;

typedef struct AtomFF {
  int ready;
  ATOM_TYPE *type;
  double avg_total_area;
  double std_total_area;
  double avg_inter_area;
  double std_inter_area;
  double max_polar_area;
  int Zp_histogram_id;
  int Zf_histogram_id;
  int A_histogram_ids[MAX_Z];
  struct ForceField *force_field;
} ATOM_FF;

typedef struct NonBondedFF {
  int surrogate;
  int ready;
  int usage_count;
  ATOM_TYPE *type1;
  ATOM_TYPE *type2;
  double P;
  double sP;
  double lnP;
  double slnP;
  char precision[10];
  double Dvdw;
  double Do;
  double Eo;
  double Po;
  double Dc;
  double k1;
  double k2;
  double opt_Z;
  double opt_sY;
  double LJ_n;
  double LJ_Z;
  double LJ_sY;
  double LR_A;
  double LR_B;
  double LR_D;
  double LR_n;
  double Z1;
  double Z2;
  double sY1;
  double sY2;
  char shape[10];
  int R_histogram_id;
  int ALPHA1_histogram_id;
  int BETA1_histogram_id;
  int ALPHA2_histogram_id;
  int BETA2_histogram_id;
  struct ForceField *force_field;
} NONBONDED_FF;

typedef struct ForceField {
  int AMI_HISTOGRAM_ID;
  ATOM_FF *atom;
  NONBONDED_FF **nonbonded;
  SETTINGS *settings;
} FORCE_FIELD;

typedef struct HistogramPoint {
  double X;
  double Y;
  double sY;
  double Yraw;
  double sYraw;
  int N;
} HISTOGRAM_POINT;

typedef struct Histogram {
  int id;
  int used;
  int protected;
  char name[10];
  long int sumN;
  double sumY;
  double startX;
  double endX;
  double stepX;
  double startY;
  double endY;
  double minX,minY;
  double maxX,maxY;
  int n_alloc_points;
  int n_points;
  HISTOGRAM_POINT *points;
  double *limits;
  int max_peaks;
  int extrapolation_method;
  int normalisation_method;
  int geometric_correction;
  int log_scale;
  double start_s;
  double end_s;
  double step_s;
} HISTOGRAM;

enum DOF_TYPE { RB_TRANSLATION,RB_ROTATION,BOND_ROTATION };

typedef struct DegreeOfFreedom {
  enum DOF_TYPE type;
  int axis_id;
  TORSION *torsion;
  ATOMLIST *atomlist;
  double **pos;
  double **u;
  double **v;
  double **w;
  double fd;
  double fd_shift;
  double fd_prev;
  double cg;
  double shift;
} DOF;

typedef struct DegreeOfFreedomList {
  int n_alloc_variables;
  int n_variables;
  DOF *variables;
  ATOMLIST *atomlist;
  double **pos;
  double **u;
  double **v;
  double **w;
  int cg_initialised;
} DOF_LIST;

// TODO: rename this to atom probe
typedef struct FieldProbe {
  int id;
  enum { ATOM_PROBE } type;
  char name[MAX_LINE_LEN];
  ATOM *atom;
  MOLECULE *molecule;
  double **u_axes;
} FIELD_PROBE;

union PLIParamValue {
  int i;
  double d;
  char s[MAX_LINE_LEN];
};

typedef struct PLIParam {
  char name[40];
  char type[40];
  union PLIParamValue value;
  int exposed;
  int inherit;
  char inherit_name[40];
  char default_value[MAX_LINE_LEN];
  char valid_values[MAX_LINE_LEN];
  char help_text[MAX_LINE_LEN];
  LIST *synonyms;
} PLI_PARAM;

typedef struct MoleculeProperty {
  char name[MAX_LINE_LEN];
  char type;
  union PLIParamValue value;
} MOLECULE_PROPERTY;

typedef struct ProbePoseLite {
  int orientation_index;
  double score;
} PROBE_POSE_LITE;

typedef struct PoseListLite {
  int n_poses;
  int n_alloc_poses;
  PROBE_POSE_LITE *poses;
} POSELIST_LITE;

double round(double);
double frand(void);
int irand(void);
double distance(double*,double*);
double sqr_distance(double*,double*);
int points_within_distance(double*,double*,double);
void set_vector(double*,double,double,double);
void null_vector(double*);
void copy_vector(double*,double*);
void calc_vector(double*,double*,double*);
void invert_vector(double*,double*);
double vector_length(double*);
double sqr_vector_length(double*);
double dotproduct(double*,double*);
void calc_crossproduct(double*,double*,double*);
double vector_angle(double*,double*);
void scale_vector(double*,double);
void sum_vector(double*,double*,double*);
void shift_vector(double*, double*);
void transform_vector(double*,double[4][4]);
void calc_transformed_vector(double*,double[4][4],double*);
void write_matrix(char*,double[4][4]);
void unit_matrix(double[4][4]);
void scalar_matrix(double, double[4][4]);
void copy_matrix(double[4][4], double[4][4]);
void translation_matrix(double*,double[4][4]);
void euler_matrix(double,double,double,double[4][4]);
void rotation_matrix(double*,double*,double,double[4][4]);
void calc_matrix_product(double[4][4],double[4][4],double[4][4]);
double matrix_determinant(double**,int);
double three_point_angle(double*,double*,double*);
double dihedral_angle(double*,double*,double*,double*);
double triangle_area(double*,double*,double*);
int solve_line(double,double,double,double,double*,double*);
int line_line_intersection(double,double,double,double,double*,double*);
int solve_parabola(double,double,double,double,double,double,double*,double*,double*);
int calc_parabola_vertex(double,double,double,double,double,double,double*,double*);
int calc_parabola_x_intercepts(double,double,double,double*,double*);
long int factorial(int);
double ramp_function(double,double,double,double,double);
double linear_interpolate(double, double, double);

void init_io(void);
void error_fn (char*, ...);
void warning_fn (char*, ...);
void debug_print(char*, ...);
char* get_pli_dir(void);
PLI_FILE* new_file(char*,FILE*);
PLI_FILE* open_file(char*,char*);
void close_file(PLI_FILE*);
int end_of_file(PLI_FILE*);
char* read_line(char*,int,PLI_FILE*);
void write_line(PLI_FILE*,const char*,...);
long int pli_ftell(PLI_FILE*);
int pli_fseek(PLI_FILE*,long int,int);
LIST* file2lines(char*);
LIST* filter_lines(LIST*,char*);
int read_word(char*,const char*,void*);
void substring(char*,int,int,char*);
void remove_spaces(char*,char*);
void remove_outer_spaces(char*,char*);
int read_intf(char*,int,char*,int*);
void upper_case(char*);
enum OUTPUT_FORMAT get_output_format(char*);
int is_readable_file(char *name_with_path);

void init_atom(ATOM*);
void unprep_atom(ATOM*);
void init_bond(BOND*);
void init_molecule(MOLECULE*,int);
void prep_molecule(MOLECULE*,SETTINGS*);
void unprep_molecule(MOLECULE*,SETTINGS*);
void free_molecule(MOLECULE*);
LIST* molecule_active_atoms(MOLECULE*);
void realloc_atoms(MOLECULE*);
void realloc_bonds(MOLECULE*);
ATOM* get_atom(MOLECULE*,int);
BOND* get_bond(MOLECULE*,ATOM*,ATOM*);
BOND* add_bond(MOLECULE*,ATOM*,ATOM*,int);
void delete_atom(MOLECULE*,int);
void delete_bond(MOLECULE*,int);
LIST *atoms2substructures(LIST*);
void set_molecule_connections(MOLECULE*,int);
double molecule_dict_match(MOLECULE*);
void reset_molecule_connections(MOLECULE*);
void get_connected_atoms(ATOMLIST*,ATOM*,ATOM*,int,unsigned int,unsigned int);
int same_molecule_atoms(ATOM*,ATOM*);
RESIDUE* atomlist2residues(ATOMLIST*,int*);
int same_residue_atoms(ATOM*,ATOM*);
int atom_residue_in_list(ATOM*,ATOMLIST*);
GRID* molecule2grid(MOLECULE*,double,double,GRID*);
GRID* alist2grid(LIST*,double,double,GRID*);
GRID* atomlist2grid(ATOMLIST*,double,double,GRID*);
MAP* molecule2atommap(MOLECULE*,double,MAP*,int);
MAP* molecule2contactmap(MOLECULE*);

void move_atom(ATOM*,double*);
void shift_atom(ATOM*,double*);
void transform_atom(ATOM*,double[4][4]);
void add_molecule_center(MOLECULE*, int);
void move_molecule_centre(MOLECULE*, double*);
void add_molecule_moi(MOLECULE*);
MOLECULE_LIST* alloc_molecule_list(void);
void init_molecule_list(MOLECULE_LIST*);
void add_molecule_to_list(MOLECULE_LIST*,MOLECULE*);
int remove_molecule_from_list(MOLECULE_LIST*,MOLECULE*);
void free_molecule_list(MOLECULE_LIST*);
ATOMLIST* alloc_atomlist(void);
void init_atomlist(ATOMLIST*);
void add_atom_to_list(ATOMLIST*,ATOM*,unsigned int);
int remove_atom_from_list(ATOMLIST*,ATOM*);
int atom_in_list(ATOMLIST*,ATOM*);
int atomlist_index(ATOMLIST*,ATOM*);
ATOMLIST* atom2atomlist(ATOM*,unsigned int);
ATOMLIST* molecule2atoms(MOLECULE*);
ATOM_COORDS* molecule2coords(MOLECULE*);
ATOMLIST* molecule2waters(MOLECULE*);
ATOMLIST* molecule2site(MOLECULE*);
ATOMLIST* atomlist2waters(ATOMLIST*);
void add_molecule_atoms_to_list(ATOMLIST*,MOLECULE*);
void add_atomlist_to_atomlist(ATOMLIST*,ATOMLIST*);
ATOMLIST* list2hydrogen_free_list(ATOMLIST*);
void reset_atom_status_list(ATOMLIST*);
void free_atomlist(ATOMLIST*);
void init_bondlist(BONDLIST*);
void add_bond_to_list(BONDLIST*,BOND*);
int bond_in_list(BONDLIST*,BOND*);
void free_bondlist(BONDLIST*);
BONDLIST* atomlist2bondlist(ATOMLIST*);
int linear_atom(ATOM*);
int planar_atom(ATOM*);
void set_system_vdw_radii(SYSTEM*);
void set_molecule_vdw_radii(MOLECULE*,double);
void set_atom_vdw_radii(ATOM*,double);
ATOM* atomlist2atoms(ATOMLIST*);
void atoms2atomlist(ATOMLIST*,ATOM*);
void copy_atomlist(ATOMLIST*,ATOMLIST*);
void copy_atom_coords(ATOM*,ATOM*);
double** atomlist2coordinates(ATOMLIST*);
void coordinates2atomlist(ATOMLIST*,double**);
double** molecule2coordinates(MOLECULE*);
void coordinates2molecule(MOLECULE*,double**);
int allow_covalent_bond(ATOM*,ATOM*);
void atomlist2centroid(ATOMLIST*,double*);
int heavy_atom_count(MOLECULE*);
void coords2molecule(ATOM_COORDS *coords, MOLECULE *molecule);
ATOM_COORDS** get_molecule_orientations(MOLECULE *molecule, int n_lebedev_points, double **lebdev_axes, int n_rot_lebedev_axis, int n_dihedral_steps, int *n_orientations);

void write_grid(FILE*,GRID*);
MAP *new_map(char*,char*,GRID*);
GRID* new_grid(double,double);
GRID* copy_grid(GRID*,GRID*);
MAP* copy_map(MAP*,MAP*);
int grid_points(GRID*);
void increment_map(MAP*,unsigned char*,double);
void free_map(MAP*);
void alloc_map_matrix(MAP*);
void*** alloc_3d_matrix(GRID*,int);
double*** alloc_3d_double_matrix(int,int,int);
void init_3d_matrix(MAP*);
int pos2grid(double*,int*,GRID*);
int pos2grid_round(double*,int*,GRID*);
int point_in_grid(double*,GRID*,double);
void point_to_gridpoint(double*,int*,GRID*);
void gridpoint_to_point(int*,double*,GRID*);
MAP_ISLAND* index2island(int,MAP_ISLAND*,int);
void remove_small_islands(MAP*,double,double);
void free_map_islands(MAP_ISLAND*,int);
MAP* copy_fmap(MAP*);
void inflate_map(MAP*);
void average_map(MAP*);
void smoothe_map(MAP*);
int point_in_island(int,MAP_ISLAND*);
MAP_ISLAND* map2islands(MAP*,double,int*);
void init_map_island(MAP_ISLAND*,int);
void grow_map_island(MAP_ISLAND*,MAP*,double,int);
void add_point_to_map_island(MAP_ISLAND*,int,int);
double island_mask_overlap(MAP_ISLAND*,MAP*,double);
double island_map_overlap(MAP_ISLAND*,MAP*);
void calc_map_island_properties(MAP_ISLAND*,MAP*);
int map_index2grid(int,int*,int*,int*,GRID*);
int map_grid2index(int,int,int,GRID*);
void v_init_double(void*);
void v_init_int(void*);
void v_init_atomlist(void*);
void v_init_alist(void*);
void v_free_list_items(void*);



double trilinear_interpolate(MAP*, double*, double);

int mask_bytes(GRID*);
BYTE* alloc_mask(GRID*);
void init_mask(BYTE*,GRID*);
BYTE* copy_mask(BYTE*,BYTE*,GRID*);
BYTE* mask_molecule(MOLECULE*,BYTE*,GRID*,int);
BYTE* mask_system_atoms(SYSTEM*,BYTE*,GRID*,int,double);
BYTE* mask_system_volume(SYSTEM*,BYTE*,GRID*,int,double);
BYTE* mask_atoms(LIST*,BYTE*,GRID*,int,double);
BYTE* mask_atom(ATOM*,BYTE*,GRID*,int,double);
BYTE* and_mask(BYTE*,BYTE*,BYTE*,GRID*);
BYTE* or_mask(BYTE*,BYTE*,BYTE*,GRID*);
int is_mask_bit_set(BYTE*,int);
void set_mask_bit(BYTE*,int);
void unset_mask_bit(BYTE*,int);
MAP* mask_map(BYTE*,MAP*,GRID*);

void run_field(SETTINGS*);
void init_field_settings(void);
void init_field_probes(SETTINGS*);
MAP* calculate_field(SYSTEM*,FIELD_PROBE*,MAP*,PLI_SFUNC*);
void init_field_system(SYSTEM*,FIELD_PROBE*);
void finish_field_system(SYSTEM*,FIELD_PROBE*);
FIELD_PROBE* get_field_probe(char*);
void move_probe_to_new_position(FIELD_PROBE*,double*);

void write_insight_map(char*,MAP*,int);
MAP* read_insight_map(char*); 

MOLECULE* read_pdb_molecule(char*);
int write_pdb_molecule(MOLECULE*,char*);
void write_pdb_atom_list(PLI_FILE*,ATOMLIST*,unsigned long int);
void write_pdb_atom(PLI_FILE*,ATOM*);
void write_pdb_bond(PLI_FILE*, BOND*);
void atom2pdbstr(ATOM*,char*);
MOLECULE* get_pdb_dict(char*);
MOLECULE* read_cif_molecule(char*);
ATOM* get_cif_atom(char*,MOLECULE*);
ELEMENT* pdb_atom_element(ATOM*,ATOM_TYPING_SCHEME*);

MOLECULE* read_mdl_molecule(char*);
void write_mdl_molecule(MOLECULE*,char*);
int read_mdl_molecule_fp(MOLECULE*);
MOLECULE_LIST* read_mdl_molecule_list(char*);
void write_mdl_atom_list(PLI_FILE*,ATOMLIST*, char*, int, int, ...);

MOLECULE* read_sybyl_molecule(char*);
void write_sybyl_molecule(MOLECULE*,char*);

void run_contacts(SETTINGS*);
void init_contact_settings(void);
void set_contacts_system(SYSTEM*,unsigned int);
void set_contacts_atom(ATOM*,SYSTEM*,unsigned int);
void contacts2ccg(ATOM *atom);
int metal_coordinating_contact(CONTACT*,double);
int flag_clashing_contacts(ATOMLIST*);
void remove_flagged_contacts(CONTACTLIST*,unsigned int);
int count_clashes_atom_list(ATOMLIST*);
void write_contacts_atom_list(PLI_FILE*,ATOMLIST*,enum OUTPUT_FORMAT,unsigned long int);
void write_clashes_atom_list(PLI_FILE*,ATOMLIST*,enum OUTPUT_FORMAT,unsigned long int);
void write_hbonds_atom_list(PLI_FILE*,ATOMLIST*,enum OUTPUT_FORMAT,unsigned long int);
void write_contacts_atom(PLI_FILE*,ATOM*,enum OUTPUT_FORMAT,unsigned long int);
void write_clashes_atom(PLI_FILE*,ATOM*,enum OUTPUT_FORMAT,unsigned long int);
void write_contact(PLI_FILE*,CONTACT*);
void write_clash(PLI_FILE*,CONTACT*,double,double,enum OUTPUT_FORMAT,unsigned long int);
void write_hbond(PLI_FILE*,CONTACT*,double,enum OUTPUT_FORMAT,unsigned long int);
void write_sas_stats_atom_list(PLI_FILE*,ATOMLIST*,enum OUTPUT_FORMAT);
int contact_in_list(CONTACTLIST*,CONTACT*);
CONTACT* find_contact(ATOM*,ATOM*);
int clash_contact(CONTACT*,double*,double*);
void free_contactlist(CONTACTLIST*);
double fraction_polar_point(ATOM*,double*,double,int);

void run_type(SETTINGS*);
void read_atom_typing_scheme(char*,ATOM_TYPING_SCHEME*);
ATOM_TYPE* resatoms2type(ATOM*,ATOM_TYPING_SCHEME*);
void set_system_atom_nodes(SYSTEM*);
void set_molecule_atom_nodes(MOLECULE*,ATOM_TYPING_SCHEME*);
ATOM_NODE* set_atom_node(ATOM*,ATOM_TYPING_SCHEME*);
void set_molecule_atom_types(MOLECULE*,ATOM_TYPING_SCHEME*);
void atom_type_heme_nitrogens(SYSTEM*);
void atom_type_atom(ATOM*,ATOM_TYPING_SCHEME*);
void update_group_statuses(ATOM*,ALT_ATOM_TYPE*,ATOM_TYPING_SCHEME*);
void resolve_system(SYSTEM*);
void unresolve_system(SYSTEM*);
ATOM_TYPING_SCHEME* alloc_atom_typing_scheme(void);
void set_atom_type(ATOM*,int,ATOM_TYPING_SCHEME*);
ATOM_TYPE* get_atom_type(char*,ATOM_TYPING_SCHEME*);
ATOM_TYPE* get_atom_type_by_id(int,ATOM_TYPING_SCHEME*);
ATOM_TYPE* get_atom_type_by_element_id(int,ATOM_TYPING_SCHEME*);
ALT_ATOM_TYPE *get_alt_type(ATOM_TYPE*,char*,int);
ATOM_NODE* get_atom_node(char*,ATOM_TYPING_SCHEME*);
UNITED_ATOM* get_united_atom(char*,ATOM_TYPING_SCHEME*);
double get_atom_type_vdw_radius(ATOM_TYPE*);
void change_atom_type(ATOM*,ATOM_TYPE*);
ATOM* pyridone_o_from_n(ATOM*);
ATOM* pyridone_n_from_o(ATOM*);
ATOM *carboxyl_o_from_o(ATOM*);

ATOM_TYPING_SCHEME* get_atom_typing_scheme(void);
ELEMENT* get_element(char*,ATOM_TYPING_SCHEME*);
ELEMENT* get_element_by_id(int,ATOM_TYPING_SCHEME*);
ELEMENT* get_element_from_type(ATOM_TYPE*);

RING_LIST* atomlist2ringlist(ATOMLIST*,int);
void find_ring(ATOM*, int);
void calc_ring_center(RING*);
void calc_ring_normal(RING*);
void free_ring_list(RING_LIST*);
void free_ring(RING*);
void set_molecule_rings(MOLECULE*);


void contacts2voronoi(ATOM*);
void write_system_polyhedra(FILE*,SYSTEM*);
void write_atom_polyhedron(FILE*,ATOM*,SETTINGS*);
double estimate_contact_iarea(CONTACT*);

int set_atom_axes(ATOM*);
void set_molecule_atom_geometries(MOLECULE*,ATOM_TYPING_SCHEME*);
void set_atom_geometry(ATOM*,ATOM_TYPING_SCHEME*);
HBOND_GEOMETRY *get_hbond_geometry(char*);
void calc_contact_geometries_system(SYSTEM*);
void calc_contact_geometry(CONTACT*);
void mirror_contact_geometry(CONTACT*,CONTACT*);
void write_contact_geometry(PLI_FILE*,CONTACT*,enum OUTPUT_FORMAT,unsigned long int);
ATOM_GEOMETRY* get_atom_geometry(char*,ATOM_TYPING_SCHEME*);
void read_atom_geometries(PLI_FILE*,ATOM_TYPING_SCHEME*);
double hbond_block_score(double,double,double,double);
void set_atom_virtual_points(ATOM*);

RESIDUE* get_amino_acid(char*);
int backbone_bond_type(BOND*);

void calc_atom_coordination(ATOM*);
double atom_hbond_energy(ATOM*,int);
double hbond_energy(CONTACT*);
double coord_score(ATOM_COORDINATION*,int);
double hbond_geometry_score(CONTACT*);
int hbond_flags_match(unsigned int,unsigned int,int);
void write_coordination_system(FILE*,SYSTEM*);

MOLECULE* read_molecule(char*);
void write_molecule(MOLECULE*,char*);
void read_molecule_fp(MOLECULE*);
MOLECULE_LIST* read_molecule_list(char*);
void write_atom_list(PLI_FILE*,ATOMLIST*);

// obj2text.c
void atom2text(ATOM*,char*);
void bond2text(BOND*,char*);
void residue2text(RESIDUE*,char*);
void site2text(MAP_ISLAND*,char*);
void set_atomio(char*,char*);
void set_bondio(char*,char*);
void set_residueio(char*,char*);
void set_siteio(char*,char*);

// json.c
char* atomlist2json(ATOMLIST*);
char* list2json(char**,int);

// molio.c
void write_atom(PLI_FILE*,ATOM*);
void write_bond(PLI_FILE*,BOND*);
void write_residue(PLI_FILE*,RESIDUE*);
int duplicate_atom(ATOM*,MOLECULE*);

FORCE_FIELD *get_ff(void);
void init_ff_settings(void);
void init_ff(SETTINGS*);
ATOM_FF* get_atom_ff(ATOM_FF*,ATOM_TYPE*);
NONBONDED_FF* get_nonbonded_ff(NONBONDED_FF**,ATOM_TYPE*,ATOM_TYPE*);
double r_alpha_correction(double);
double r_beta_correction(double);
double alpha_beta_correction(double);

void init_histogram_settings();
HISTOGRAM* read_histogram(PLI_FILE*,char*);
HISTOGRAM* create_histogram(void);
HISTOGRAM* add_histogram(char*,int);
void process_histogram(HISTOGRAM*);
void reset_histogram(int);
HISTOGRAM* clone_histogram(HISTOGRAM*,int);
void set_histogram_extrema(HISTOGRAM*);
double histogram_X2Y(HISTOGRAM*,double);
HISTOGRAM* get_histogram(int);
void write_histogram(HISTOGRAM*,ATOM_TYPE*,ATOM_TYPE*);
void free_histogram(HISTOGRAM*);
int calc_histogram_maximum(HISTOGRAM*,double*,double*);
void normalise_histogram(HISTOGRAM*,int);
void log_histogram(HISTOGRAM*);
int remove_local_extreme(HISTOGRAM*);
void histogram2xys(HISTOGRAM*,double,double,double,double,double*,double*,double*,int*);

void run_score(SETTINGS*);
void prep_score_system(SYSTEM*,PLI_SFUNC*);
void unprep_score_system(SYSTEM*,FLAG,PLI_SFUNC*);
void calc_system_ref_scores(SYSTEM*,PLI_SFUNC*);
void score_system(SYSTEM*,PLI_SFUNC*);
void minimise_system(SYSTEM*,PLI_SFUNC*);
void set_dof_score_gradient(DOF*,SYSTEM*,PLI_SFUNC*);
double numerical_score_gradient(DOF*,SYSTEM*,double (*score_system)(struct System*));
void init_system_scores(SYSTEM*);
void init_atom_scores(ATOM*);
void write_system_scores(PLI_FILE*,SYSTEM*);
void write_atom_scores(PLI_FILE*,ATOM*);
void write_contact_scores(PLI_FILE*,CONTACT*);
void add_system_scores(SYSTEM*,SYSTEM*);
void store_system_ref_scores(SYSTEM*);
void store_atom_ref_scores(ATOM*);
void mirror_contact_scores(CONTACT*,CONTACT*);

PLI_MODE* get_mode(char*);
PLI_SFUNC* get_sfunc(char*);

void pliff_init(void);
double pliff_score_system(SYSTEM*);
double pliff_score_gradient(DOF*,SYSTEM*);
double pliff_score_intra_covalent(SYSTEM*);
void pliff_score_molecule(MOLECULE*);
void pliff_score_contact(CONTACT*,ATOM_TYPING_SCHEME*,FORCE_FIELD*);
void pliff_write_system_scores(PLI_FILE*,SYSTEM*);
void pliff_write_atom_scores(PLI_FILE*,ATOM*);
void pliff_write_contact_scores(PLI_FILE*,CONTACT*);
void init_pliff_settings();
double pliff_score_vs_histogram(int,double,int);

void pliff_minimise_system(SYSTEM*);

double molecule_internal_energy(MOLECULE*,ATOM_TYPING_SCHEME*,FORCE_FIELD*);
double molecule_bond_energy(MOLECULE*);
double molecule_bond_angle_energy(MOLECULE*);
double molecule_torsion_energy(MOLECULE*);
double molecule_vdw_energy(MOLECULE*,ATOM_TYPING_SCHEME*,FORCE_FIELD*,int);

void minimise(SYSTEM*,PLI_SFUNC*,int);
void init_minimise_settings(void);
void read_minimise_settings(PLI_FILE*);

LIST* get_lebedev_vectors(int);
double** read_lebedev_sphere_points(int);

LIST* new_list(char*,int,int);
void init_list(LIST*,char*,int);
void reset_list(LIST*);
void* add_list_item(LIST*);
void remove_list_item(LIST*,void*);
void* get_list_item(LIST*,int);
int item_in_list(LIST*,void*);
void condense_list(LIST*);
void free_list(LIST*);
void store_list(LIST*,char*);
LIST* get_list(char*,char*);
LIST* atomlist2list(ATOMLIST*,char*);

void* mem_claim(int);
int** alloc_2d_imatrix(int,int);
char** alloc_2d_cmatrix(int,int);
double** alloc_2d_fmatrix(int,int);
void free_2d_imatrix(int**);
void free_2d_cmatrix(char**);
void free_2d_fmatrix(double**);

void normalise_bfactors(SYSTEM*);
void calc_system_simb(SYSTEM*,int);
void calc_molecule_ami(MOLECULE*,SETTINGS*);

int flag_skipped_waters(SYSTEM*);
void flag_active_waters(SYSTEM*);

SYSTEM* settings2system(SETTINGS*);
SYSTEM_LIST* settings2syslist(SETTINGS*);
void free_system(SYSTEM*);
void output_system(PLI_FILE*,SYSTEM*);
SYSTEM* molecule2system(MOLECULE*,SETTINGS*);
void reset_system(SYSTEM*,int);
void free_molsystem(SYSTEM*);
void prep_system_molecules(SYSTEM*,unsigned int);
void unprep_system_molecules(SYSTEM*,unsigned int);
GRID* system_grid(SYSTEM*,double,double);
LIST* system_active_atoms(SYSTEM*);

SETTINGS* get_settings(void);
void prep_settings(int,char**);
void unprep_settings(void);

double lennard_jones_energy(double,double,double,double);

double score_system_constraints(SYSTEM*);
double constraint_score_gradient(DOF*,SYSTEM*);
void prep_min_atom_tethers(SYSTEM*);
void unprep_min_atom_tethers(SYSTEM*);
void prep_template_atom_tethers(SYSTEM*);
void unprep_template_atom_tethers(SYSTEM*);

void set_rotatable_bonds(MOLECULE*);

SYSMOLS* load_sysmols(SETTINGS*);

unsigned int dof_get_space(char*);
unsigned int dofs_molecule_space(unsigned int,unsigned int);
void set_dof_list_shifts(DOF_LIST*,double);
void apply_dof_list_shifts(DOF_LIST*);
void apply_dof_shift(DOF*,double);
DOF_LIST* setup_dof_list(SYSTEM*);
void free_dof(DOF*);
void free_dof_list(DOF_LIST*);
void store_positions(ATOMLIST*,double**,double**,double**,double**);
void restore_positions(ATOMLIST*,double**,double**,double**,double**);
ATOM_COORDS* get_list_coords(ATOMLIST*);
void set_list_coords(ATOMLIST*,ATOM_COORDS*);

void run_help(SETTINGS*);
void params_init(int nargs,char **args);
PLI_PARAM* params_get_parameter(char*);
void params_set_parameter(PLI_PARAM*,char*);

typedef struct AtomTypeList {
  int n_types;
  int n_alloc_types;
  ATOM_TYPE **types;
} ATOM_TYPE_LIST;

typedef struct IntList {
  int *items;
  int n_items;
  int n_allocated;
} INT_LIST;

INT_LIST* parse_int_csv(char *s, char *sep);
void free_int_array (INT_LIST *array);
void skip_atom_ids(MOLECULE *molecule, INT_LIST *skip_atom_ids);
double molecule_max_size(MOLECULE *molecule);
void title_case(char*);
int word_in_text(char*,char*,char);
char* nextword(char*,char,int*);
int read_string_section(char*,char*);

// triplet
typedef struct Triplet {
  ATOM *atom1;
  ATOM *atom2;
  ATOM *atom3;
  double hbond_geom_score1;
  double hbond_geom_score2;
  double angle;
} TRIPLET;

typedef struct TripletList {
  int ntriplets;
  int n_alloc_triplets;
  TRIPLET *triplets;
} TRIPLETLIST;
//
void set_hbond_triplets_system(SYSTEM *system);
void write_hbond_triplets_system(PLI_FILE*,SYSTEM*, enum OUTPUT_FORMAT);
void write_hbond_triplets_atom(PLI_FILE*,ATOM*, enum OUTPUT_FORMAT);
void write_hbond_triplet(PLI_FILE*,TRIPLET*, enum OUTPUT_FORMAT);
int hbond_type_match(ATOM*,ATOM*);
int connected_atoms(ATOM *atom1,ATOM *atom2,ATOM *last_atom,int max_bonds);
int molecule_accessible(MAP *field, MOLECULE *molecule);
double molecule_radius(MOLECULE *molecule);
FIELD_PROBE* create_atom_type_probe(ATOM_TYPE *type, SETTINGS *settings, FIELD_PROBE *probe_in);
int position_within_dist_molecule(double*, MOLECULE*, double, double);
int uniform_rand_int(int, int);
double uniform_rand(double, double);
double normal_rand(double, double);
int map_extreme_index(MAP*,int, char*,int*, int*, int*);

BOOLEAN **alloc_2d_bool_matrix(int, int);
void free_2d_bool_matrix(BOOLEAN **);
void copy_2d_bool_matrix(BOOLEAN**, BOOLEAN**, int, int);
BOOLEAN **create_adjacency_matrix (MOLECULE *);
int is_in_iarray(int v, int *arr, int n);
void set_double_array(double*,int,double);

void get_selected_atom_indices(MOLECULE *molecule, int *ids);
BOOLEAN* alloc_bool_array(int n);
void print_grid(GRID*);
void coords_bounds(ATOM_COORDS **coords, int n_coords, MOLECULE *molecule, double bounds[3][2]);
GRID* coords2grid(ATOM_COORDS **coords, int n_coords, MOLECULE *molecule, double spacing, double padding, double *center);
ATOMLIST* get_hbond_atomlist(MOLECULE *molecule);
int position_within_dist_atomlist(double *pos, ATOMLIST *atomlist, double dist, double angle);
void mask_map_cutoff(MAP *map, double cutoff, char *keep);
ATOMLIST* molecule2hbonding_atoms(MOLECULE *molecule);
int has_hbonding_atom(MOLECULE *molecule);
MOLECULE* create_atom_type_molecules(ATOM_TYPE *type, int n_molecule, char *molecule_name, SETTINGS *settings);

void set_molecule_internal(MOLECULE*);
void init_grid(GRID*, double, double);
void update_gridlimits(GRID*, double*);
int position_within_dist_system(double*, SYSTEM*, double, double);
void prep_system(SYSTEM*, unsigned int);
void unprep_system(SYSTEM*, unsigned int);
unsigned int dofs_get_space(char*);
void shake_molecule(MOLECULE*);
void shake_atom(ATOM*,double,double);
void use_fixed_seed (void);
TORSION* get_molecule_torsion(MOLECULE*,ATOM*,ATOM*);
void set_torsion_angle(TORSION*,double);
void calc_torsion_angle(TORSION*);
void rotate_torsion(TORSION*,double);
void set_molecule_torsions(MOLECULE*);
double contact_angle(double*,double*,int);
void free_3d_matrix(void***);

double polynomial(double,double*,int);
int weighted_asymmetric_parabola_fit(double*,double*,double*,int,double,double,double*,double*);
int weighted_polynomial_fit(double*,double*,double*,int,int,double*,double*);
int weighted_power_fit(double*,double*,double*,int,double*,double*);
int weighted_linear_least_squares(double*,double*,double*,int,double*,double*,double*);
int weighted_LJ_Eo_fit(double*,double*,double*,int,int,double*,double*);
int weighted_LJ_fit(double*,double*,double*,int,int,double*,double*,double*);
double normalised_fit_error(double*,double*,double*,int);
double weighted_fit_error(double*,double*,double*,int);

MAP* ligsite_map(SYSTEM*);
MAP* ami_site_map(SYSTEM*);
MAP* digsite_map(SYSTEM*);

void pocket_run(SETTINGS*);
double pocket_score(SYSTEM*);
double atom_pocket_score(ATOM*);
POCKET_MODE* pocket_mode(char*);
void set_atom_depth_maps(SYSTEM*);
double atom_depth_weight(ATOM*);

void test_smarts(SETTINGS*);
LIST* atom_matches(ATOM*,MOLECULE*,int);
MOLECULE* smarts2mol(char*);
void atoms2pattern(LIST*,MOLECULE*);
void apply_pattern(MOLECULE*,ATOM**);

void run_resolve(SETTINGS*);



typedef struct MoleculeProbe {
  char name[MAX_LINE_LEN];
  MOLECULE *molecule;
  ATOM_COORDS** orientations;
  int n_orientations;
  MAP *allocated_maps;
  MAP **atom_maps;
  //TODO: change to atom mappings
  BOOLEAN *** auto_mappings;
  int n_mappings;
  int n_alloc_mappings;
  double score;
  int orientation_id;
  int *selected_ids;
} MOLECULE_PROBE;

typedef struct ProbePose {
  double center[4];
  int orientation_index;
  double score;
  double **position;
  double **water_positions;
  int n_added_waters;
  double water_score;
} PROBE_POSE;

typedef struct DockedWaters {
  int n_waters;
  MOLECULE *mols;
  double score;
} DOCKED_WATERS;

void run_fragmap(SETTINGS*);
void copy_field(MAP*, MAP**);
void get_min_score_point(MAP*, int*);
void copy_3d_double_matrix(double***, double***, int*);
void copy_3d_int_matrix(int***, int***, int*);
void copy_2d_int_matrix(int**, int**, int , int);
MOLECULE_PROBE* alloc_molecule_probe(void);
MOLECULE_PROBE* get_molecule_probe(char*, SETTINGS*, char*);


typedef struct PoseList {
  int max_n;
  int n_poses;
  PROBE_POSE **poses;
  PROBE_POSE *worst_score_pose;
  PROBE_POSE *best_score_pose;
} POSE_LIST;

typedef struct DivClusterSet {
  POSE_LIST **pose_lists;
  int n_alloc_clusters;
  int n_poses;
  int n_max_poses;
  POSE_LIST *worst_score_cluster;
  PROBE_POSE *allocated_poses;
} CLUSTER_SET;

enum DIST_TYPE {RMSD_DIST, CENTER_DIST};

typedef struct ClusterSettings {
  int save_poses;
  int n_saved_poses;
  int n_poses_per_orientation;
  int n_poses_per_gridpoint;
  char *saved_sdf_file;
  char *saved_data_file;
  int save_diverse_clusters;
  int n_diverse_poses;
  int diverse_cluster_size;
  double save_cutoff_atom;
  double cluster_cutoff;
  double sqr_cluster_cutoff;
  char *diverse_sdf_file_prefix;
  char *diverse_data_file;
  enum DIST_TYPE dist_type;
  int min_dist;
  int merge_clusters;
  int save;
  int save_peak_clusters;
  int n_peaks;
  int peak_cluster_size;
  int minimise;
} CLUSTER_SETTINGS;

CLUSTER_SET* alloc_div_cluster_set(int,int,int);
void free_div_cluster_set(CLUSTER_SET*);
void save_diverse_pose(CLUSTER_SET*, MOLECULE_PROBE*, DOCKED_WATERS*);
void write_cluster_set(CLUSTER_SET*, char*, MOLECULE_PROBE *probe, MOLECULE*);
PROBE_POSE* alloc_poses(int, int, int);
void free_poses(PROBE_POSE*);
POSE_LIST* alloc_pose_list(int);
void free_pose_list(POSE_LIST*);
void save_pose(POSE_LIST*, PROBE_POSE *, MOLECULE_PROBE *probe);
void write_pose_list(POSE_LIST*, MOLECULE*, ATOM_COORDS**);
CLUSTER_SETTINGS* get_cluster_settings(void);
void run_cluster(SETTINGS*);
void init_cluster_settings();
CLUSTER_SET* get_peak_clusters(MAP*, MAP*, MOLECULE_PROBE*, ATOM_COORDS**, CLUSTER_SETTINGS*, PROBE_POSE*, BOOLEAN***, int);
void init_3d_poselist_matrix(MAP*);
POSELIST_LITE*** alloc_3d_poselist_matrix(int,int,int);
void add_pose_to_poselist(POSELIST_LITE*, int, double);

typedef struct UllmanMatchList {
  int n_matches;
  int n_alloc_matches;
  int n_atoms;
  MATCH *matches;
} ULLMAN_MATCH_LIST;

typedef struct UllmanMappings {
  int n_mappings;
  int n_atoms1;
  int n_atoms2;
  int n_allocated;
  BOOLEAN ***mappings;
} ULLMAN_MAPPINGS;

enum MATCH_PROTOCOL { MATCH_ALL,MATCH_UNIQUE,MATCH_FIRST };
enum MATCH_QUALITY { MATCH_ATOM_TYPE,MATCH_NODE,MATCH_ELEMENT,MATCH_ANY };

ULLMAN_MAPPINGS* find_ullman_mappings(MOLECULE*, MOLECULE*, enum MATCH_PROTOCOL, enum MATCH_QUALITY);
ULLMAN_MATCH_LIST* find_ullman_matches(MOLECULE*,MOLECULE*,enum MATCH_PROTOCOL,enum MATCH_QUALITY);
void do_ullman (MOLECULE *mola, MOLECULE *molb, void (*func) (BOOLEAN **));



double sqr_distance_offset(double *v1,double *v2, double *offset1, double *offset2);
double rms_coords_atleat(ATOM_COORDS *coords1, ATOM_COORDS *coords2, int *selected1, int *selected2, int n_selected, BOOLEAN ***isomorphs, int n_isomprphs, double upper_bound);
double rms_coords(ATOM_COORDS *coords1, ATOM_COORDS *coords2, int *selected1, int *selected2, int n_selected, BOOLEAN ***isomorphs, int n_isomprphs, double*, double*);


void set_molecule_min_coords(MOLECULE*, double*);


int*** alloc_3d_int_matrix(int,int,int);
void free_3d_int_matrix(int***);
void move_grid_center(GRID*, double*);


void set_mask_bit(unsigned char*,int);


void mask_distant_points(MAP *map, SYSTEM *system, MOLECULE *molecule, ATOMLIST*, double max_dist, double angle);


void print_vector(double *v);


