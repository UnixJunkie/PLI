// Copyright 2015 Astex Therapautics Ltd.
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



#define YES 1
#define NO 0
#define UNDEFINED -1

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

#define COVALENT_TOLERANCE 0.4
#define MAX_LINE_LEN 3000
#define sqr(x) ((x)*(x))
#define PI 3.1415926536
#define GOLDEN_RATIO 0.618034
#define MAX_RING_SIZE 500
#define DEFAULT_VDW_RADIUS 1.70
#define RT 2.47897
#define MAX_ATOM_SCORES 10
#define MAX_CONTACT_SCORES 10
#define MAX_SYSTEM_SCORES 100

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

#define BULK_WATER_TD 2.8

// definitions for degrees of freedom:
#define DOF_LIGAND_TRANSLATION 1
#define DOF_LIGAND_ROTATION 2
#define DOF_LIGAND_TORSIONS 4
#define DOF_LIGAND_RIGID_BODY 3
#define DOF_LIGAND 7
#define DOF_WATER_TRANSLATION 8
#define DOF_WATER 8

// map mask types:
#define CONTACT_MASK 1
#define ATOM_MASK 2
#define VDW_MASK 3



enum MAP_TYPE { INTEGER_MAP,DOUBLE_MAP,ATOMLIST_MAP };

enum FILE_FORMAT { GZIPPED_FILE,ASCII_FILE };

enum OUTPUT_FORMAT { FORMATTED,JSON };

enum ATOM_STYLE { BASIC_ASTYLE,PDB_ASTYLE,PLIFF_ASTYLE,TYPE_ASTYLE };

enum MINIMISATION_ALGORITHM { STEEPEST_DESCENT,CONJUGATE_GRADIENT };



typedef struct PLIFile {
  char filename[MAX_LINE_LEN];
  FILE *file;
  enum FILE_FORMAT format;
} PLI_FILE;

extern PLI_FILE *PLI_STDIN;
extern PLI_FILE *PLI_STDOUT;
extern PLI_FILE *PLI_STDERR;

typedef struct Grid {
  int npoints[3]; /* number of grid points along each of the 3 axes */
  int limit[3][2]; /* limits in grid points (including sign) */
  double flimit[3][2]; /* real box limits */
  double spacing;
  double padding;
} GRID;

typedef struct Map {
  enum MAP_TYPE type;
  char name[MAX_LINE_LEN];
  GRID *grid;
  void ***matrix;
  unsigned char *mask;
} MAP;

typedef struct MapIsland {
  int id;
  int *points;
  int n_points;
  int n_alloc_points;
  double maxY;
  double sumY;
  double avgY;
  double integral;
  double volume;
} MAP_ISLAND;

typedef struct Element {
  int id;
  char name[10];
  double cov_radius;
  double vdw_radius;
  unsigned int flags;
} ELEMENT;

typedef struct AtomCoords {
  double position[4];
  double u[4];
  double v[4];
  double w[4];
} ATOM_COORDS;

typedef struct Atom {
  int unique_id;
  int id;
  int seqid;
  ELEMENT *element;
  char name[10];
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
  double vdw_radius;
  double vdw_radius_H2O;
  double tether_k;
  struct AtomNode *node;
  struct AtomType *type;
  struct AtomType *ambiguous_type;
  struct AtomList *connections;
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
  unsigned int cstatus;
  unsigned int status;
  unsigned int flags;
  unsigned int error_flags;
} ATOM;

struct AtomList;
struct System;
struct Torsion;

typedef struct Bond {
  int type;
  struct Atom *atom1;
  struct Atom *atom2;
  struct Torsion *torsion;
} BOND;

typedef struct TorsionType {
  char name[MAX_LINE_LEN];
  int hyb1;
  int hyb2;
  double A;
  double n;
  double Xo;
} TORSION_TYPE;

typedef struct Torsion {
  BOND *bond;
  ATOM *atom1;
  ATOM *atom2;
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
  int n_torsions;               /* the number of bond torsions in the molecule */
  ATOM* atom;			/* array of atoms (see struct ATOM) */
  BOND* bond;			/* array of bonds (see struct BOND) */
  TORSION *torsions;            /* array of bond torsions */
  int **nb_atom_pairs;          /* matrix describing atom pairs to be checked for nonbonded contacts */
  int n_alloc_atoms;
  int n_alloc_bonds;
  MAP *covalent_map;
  MAP *contacts_map;
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

typedef struct PliMode {
  char name[MAX_LINE_LEN];
  void (*run)(struct Settings*);
  char selection[MAX_LINE_LEN];
  int selection_deviation_warning;
} PLI_MODE;

typedef struct PliSfunc {
  char name[MAX_LINE_LEN];
  void (*score_system)(struct System*);
  void (*minimise_system)(struct System*);
  void (*write_system_scores)(PLI_FILE*,struct System*,enum OUTPUT_FORMAT,unsigned long int);
  void (*write_atom_scores)(PLI_FILE*,struct Atom*,enum OUTPUT_FORMAT,unsigned long int);
  void (*write_contact_scores)(PLI_FILE*,struct Contact*,enum OUTPUT_FORMAT,unsigned long int);
  double bad_score;
} PLI_SFUNC;

typedef struct SysMols {
  MOLECULE *protein;
  MOLECULE *water;
  MOLECULE *symmetry;
  MOLECULE *ligand;
  MOLECULE_LIST *ligand_list;
  COMPLEX_LIST *complex_list;
} SYSMOLS;

typedef struct Settings {
  char *pli_dir;
  int verbose;
  PLI_MODE *mode;
  PLI_SFUNC *sfunc;
  char jobname[MAX_LINE_LEN];
  char selection[100];
  char settings_file[100];
  char types_file[1000];
  char protein_file[1000];
  char ligand_file[1000];
  char symmetry_file[1000];
  char oatoms[MAX_LINE_LEN];
  enum OUTPUT_FORMAT oformat;
  unsigned long int oflags;
  enum ATOM_STYLE astyle;
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
  int id1,id2;
  char hist_name[20];
  SYSMOLS *sysmols;
  int minimise;
  int allow_bad_clashes;
  int allow_covalent_bonds;
  char keep_waters[200];
} SETTINGS;

typedef struct System {
  int id;
  char name[MAX_LINE_LEN];
  struct Settings *settings;
  MOLECULE *protein;
  MOLECULE *water;
  MOLECULE *ligand;
  MOLECULE *symmetry;
  MOLECULE_LIST *molecule_list;
  ATOMLIST *selection;
  ATOM_COORDS *coords;
  double score;
  double scores[MAX_SYSTEM_SCORES];
  double ref_score;
  double ref_scores[MAX_SYSTEM_SCORES];
  double constraint_score;
  double min_dx;
  double min_ds;
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
  int ready;
  int usage_count;
  ATOM_TYPE *type1;
  ATOM_TYPE *type2;
  double P;
  double sP;
  double lnP;
  double slnP;
  char precision[10];
  double Do;
  double Dc;
  int use_clash_potential;
  double LJ_A;
  double LJ_B;
  int R_histogram_id;
  int ALPHA1_histogram_id;
  int BETA1_histogram_id;
  int ALPHA2_histogram_id;
  int BETA2_histogram_id;
  struct ForceField *force_field;
} NONBONDED_FF;

typedef struct ForceField {
  ATOM_FF *atom;
  NONBONDED_FF **nonbonded;
  SETTINGS *settings;
} FORCE_FIELD;

typedef struct HistogramPoint {
  double X;
  double Y;
  double sY;
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
  int sharpen;
  double sharpen_sX;
  int log_scale;
  double start_s;
  double end_s;
  double step_s;
} HISTOGRAM;

enum DOF_TYPE { RB_TRANSLATION,RB_ROTATION,BOND_ROTATION };

typedef struct DegreeOfFreedom {
  enum DOF_TYPE type;
  int axis_id;
  ATOM *atom1;
  ATOM *atom2;
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

typedef struct FieldProbe {
  int id;
  enum { ATOM_PROBE,MOLECULE_PROBE } type;
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
  char default_value[MAX_LINE_LEN];
  char valid_values[MAX_LINE_LEN];
  char help_text[MAX_LINE_LEN];
} PLI_PARAM;



double round(double);
double frand(void);
int irand(void);
double distance(double*,double*);
double sqr_distance(double*,double*);
int points_within_distance(double*,double*,double);
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
void transform_vector(double*,double[4][4]);
void calc_transformed_vector(double*,double[4][4],double*);
void write_matrix(char*,double[4][4]);
void unit_matrix(double[4][4]);
void translation_matrix(double*,double[4][4]);
void euler_matrix(double,double,double,double[4][4]);
void rotation_matrix(double*,double*,double,double[4][4]);
void calc_matrix_product(double[4][4],double[4][4],double[4][4]);
double torsion_angle(double*,double*,double*,double*);
double triangle_area(double*,double*,double*);
int solve_line(double,double,double,double,double*,double*);
int solve_parabola(double,double,double,double,double,double,double*,double*,double*);
int calc_parabola_vertex(double,double,double,double,double,double,double*,double*);
long int factorial(int);
double ramp_function(double,double,double,double,double);

void init_io(void);
void error_fn (char*, ...);
void warning_fn (char*, ...);
char* get_pli_dir(void);
PLI_FILE* new_file(char*,FILE*);
PLI_FILE* open_file(char*,char*);
void close_file(PLI_FILE*);
int end_of_file(PLI_FILE*);
char* read_line(char*,int,PLI_FILE*);
void write_line(PLI_FILE*,const char*,...);
long int pli_ftell(PLI_FILE*);
int pli_fseek(PLI_FILE*,long int,int);
void read_word(char*,const char*,void*);
void substring(char*,int,int,char*);
void remove_spaces(char*,char*);
void remove_outer_spaces(char*,char*);
void upper_case(char*);
enum OUTPUT_FORMAT get_output_format(char*);

void atom2name(ATOM*,char*);
void init_atom(ATOM*);
void unprep_atom(ATOM*);
void init_bond(BOND*);
void init_molecule(MOLECULE*,int);
void prep_molecule(MOLECULE*,SETTINGS*);
void unprep_molecule(MOLECULE*,SETTINGS*);
void realloc_atoms(MOLECULE*);
void realloc_bonds(MOLECULE*);
ATOM* get_atom(MOLECULE*,int);
BOND* get_bond(MOLECULE*,ATOM*,ATOM*);
BOND* add_bond(MOLECULE*,ATOM*,ATOM*,int);
void delete_atom(MOLECULE*,int);
void delete_bond(MOLECULE*,int);
void set_system_connections(SYSTEM*);
void set_molecule_connections(MOLECULE*,int);
double molecule_dict_match(MOLECULE*);
void reset_molecule_connections(MOLECULE*);
void get_connected_atoms(ATOMLIST*,ATOM*,ATOM*,int,unsigned int,unsigned int);
int same_molecule_atoms(ATOM*,ATOM*);
int same_residue_atoms(ATOM*,ATOM*);
int atom_residue_in_list(ATOM*,ATOMLIST*);
GRID* system2grid(SYSTEM*,double,double);
MAP* system2atommap(SYSTEM*,double);
GRID* molecule2grid(MOLECULE*,double,double,GRID*);
GRID* atomlist2grid(ATOMLIST*,double,double,GRID*);
MAP* molecule2atommap(MOLECULE*,double,MAP*,int);
void move_atom(ATOM*,double*);
void shift_atom(ATOM*,double*);
void transform_atom(ATOM*,double[4][4]);
void molecule_center(MOLECULE*,double*);
ATOMLIST* molecule_sphere_to_atom_list(MOLECULE*,MOLECULE*,double);
ATOMLIST* molecule_atom_ids_to_atom_list(MOLECULE*,int*,int);
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
ATOMLIST* molecule2waters(MOLECULE*);
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
int planar_atom(ATOM*);
void set_system_vdw_radii(SYSTEM*);
void set_molecule_vdw_radii(MOLECULE*,double);
void set_atom_vdw_radii(ATOM*,double);
ATOM* atomlist2atoms(ATOMLIST*);
void atoms2atomlist(ATOMLIST*,ATOM*);
void copy_atomlist(ATOMLIST*,ATOMLIST*);
double** atomlist2coordinates(ATOMLIST*);
void coordinates2atomlist(ATOMLIST*,double**);
double** molecule2coordinates(MOLECULE*);
void coordinates2molecule(MOLECULE*,double**);
void select_molecule_atoms(MOLECULE*,char*);
int allow_covalent_bond(ATOM*,ATOM*);
void atomlist2centroid(ATOMLIST*,double*);

void write_grid(FILE*,GRID*);
MAP *new_map(char*);
MAP *new_system_map(char*,SYSTEM*,GRID*,enum MAP_TYPE);
void free_map(MAP*);
void alloc_map_matrix(MAP*);
void*** alloc_3d_matrix(GRID*,enum MAP_TYPE);
void init_3d_matrix(MAP*);
GRID *new_grid(void);
int pos2grid(double*,int*,GRID*);
int pos2grid_round(double*,int*,GRID*);
int point_in_grid(double*,GRID*,double);
void point_to_gridpoint(double*,int*,GRID*);
void gridpoint_to_point(int*,double*,GRID*);
int is_mask_bit_set(unsigned char*,int);
void mask_map_molecule(MAP*,MOLECULE*,int);
void set_mask_bit(unsigned char*,int);
void init_3d_mask(MAP*);
MAP_ISLAND* index2island(int,MAP_ISLAND*,int);
int point_in_island(int,MAP_ISLAND*);
MAP_ISLAND* map2islands(MAP*,double,int*);
double island_mask_overlap(MAP_ISLAND*,MAP*,double);
double island_map_overlap(MAP_ISLAND*,MAP*);
void calc_map_island_properties(MAP_ISLAND*,MAP*);
int map_index2grid(int,int*,int*,int*,GRID*);
int map_grid2index(int,int,int,GRID*);

void run_field(SETTINGS*);
void init_field_settings(void);
void init_field_probes(SETTINGS*);
MAP* calculate_field(SYSTEM*,FIELD_PROBE*,PLI_SFUNC*);
MAP* init_field(char*,SYSTEM*);
FIELD_PROBE* get_field_probe(char*);

void init_site(SETTINGS*);
ATOMLIST* molecule2site(MOLECULE*);
void init_site_settings(void);

void write_insight_map(char*,MAP*,int);
MAP* read_insight_map(char*); 

MOLECULE* read_pdb_molecule(char*);
int write_pdb_molecule(MOLECULE*,char*);
void write_pdb_atom_list(PLI_FILE*,ATOMLIST*,unsigned long int);
void write_pdb_atom(PLI_FILE*,ATOM*);
void atom2pdbstr(ATOM*,char*);
MOLECULE* get_pdb_dict(char*);
MOLECULE* read_cif_molecule(char*);
ATOM* get_cif_atom(char*,MOLECULE*);
ELEMENT* pdb_atom_element(ATOM*,ATOM_TYPING_SCHEME*);

MOLECULE* read_mdl_molecule(char*);
int read_mdl_molecule_fp(MOLECULE*);
MOLECULE_LIST* read_mdl_molecule_list(char*);
void write_mdl_atom_list(PLI_FILE*,ATOMLIST*);

MOLECULE* read_sybyl_molecule(char*);

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
void write_contact(PLI_FILE*,CONTACT*,enum OUTPUT_FORMAT,unsigned long int);
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

RING_LIST* atomlist2ringlist(ATOMLIST*,int);
void find_ring(ATOM*,int);
void calc_ring_center(RING*);
void calc_ring_normal(RING*);
void free_ring_list(RING_LIST*);

void contacts2voronoi(ATOM*);
void write_system_polyhedra(FILE*,SYSTEM*);
void write_atom_polyhedron(FILE*,ATOM*,SETTINGS*);
double estimate_contact_iarea(CONTACT*);

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
int hbond_flags_match(unsigned int,unsigned int);
void write_coordination_system(FILE*,SYSTEM*);

MOLECULE* read_molecule(char*);
void read_molecule_fp(MOLECULE*);
MOLECULE_LIST* read_molecule_list(char*);
void write_atom_list(PLI_FILE*,ATOMLIST*,enum ATOM_STYLE,enum OUTPUT_FORMAT,unsigned long int);
enum ATOM_STYLE get_atom_style(char*);
void write_atom(PLI_FILE*,ATOM*,enum ATOM_STYLE,enum OUTPUT_FORMAT,unsigned long int);
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
HISTOGRAM* new_histogram(void);
HISTOGRAM* add_histogram(char*,int);
void process_histogram(HISTOGRAM*);
void reset_histogram(int);
HISTOGRAM* clone_histogram(HISTOGRAM*);
void set_histogram_extrema(HISTOGRAM*);
double histogram_X2Y(HISTOGRAM*,double);
HISTOGRAM* get_histogram(int);
void write_histogram(HISTOGRAM*,ATOM_TYPE*,ATOM_TYPE*);

void run_score(SETTINGS*);
void prep_score(SYSTEM*,PLI_SFUNC*);
void unprep_score(SYSTEM*);
void calc_system_ref_scores(SYSTEM*,PLI_SFUNC*);
void score_system(SYSTEM*,PLI_SFUNC*);
void minimise_system(SYSTEM*,PLI_SFUNC*);
void init_score(SETTINGS*);
void init_system_scores(SYSTEM*);
void init_atom_scores(ATOM*);
void write_system_scores(PLI_FILE*,SYSTEM*,enum OUTPUT_FORMAT,unsigned long int);
void write_atom_scores(PLI_FILE*,ATOM*,enum OUTPUT_FORMAT,unsigned long int);
void write_contact_scores(PLI_FILE*,CONTACT*,enum OUTPUT_FORMAT,unsigned long int);
void add_system_scores(SYSTEM*,SYSTEM*);
void store_system_ref_scores(SYSTEM*);
void store_atom_ref_scores(ATOM*);
void mirror_contact_scores(CONTACT*,CONTACT*);

PLI_MODE* get_mode(char*);
PLI_SFUNC* get_sfunc(char*);
void pliff_score_system(SYSTEM*);
void pliff_score_molecule(MOLECULE*);
void pliff_score_contact(CONTACT*,ATOM_TYPING_SCHEME*,FORCE_FIELD*);
void pliff_write_system_scores(PLI_FILE*,SYSTEM*,enum OUTPUT_FORMAT,unsigned long int);
void pliff_write_atom_scores(PLI_FILE*,ATOM*,enum OUTPUT_FORMAT,unsigned long int);
void pliff_write_contact_scores(PLI_FILE*,CONTACT*,enum OUTPUT_FORMAT,unsigned long int);
void init_pliff_settings();

void pliff_minimise_system(SYSTEM*);

double molecule_internal_energy(MOLECULE*,ATOM_TYPING_SCHEME*,FORCE_FIELD*);

void minimise(SYSTEM*,PLI_SFUNC*,int);
void init_minimise_settings(void);
void read_minimise_settings(PLI_FILE*);

double** read_lebedev_sphere_points(int);

int** alloc_2d_imatrix(int,int);
double** alloc_2d_fmatrix(int,int);
void free_2d_imatrix(int**);
void free_2d_fmatrix(double**);

void normalise_bfactors(SYSTEM*);
void calc_system_simb(SYSTEM*,int);

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

SETTINGS* get_settings(int,char**);

double lennard_jones_energy(double,double,double);
double lennard_jones_energy_DE(double,double,double);

void score_system_constraints(SYSTEM*);
void setup_list_atom_tethers(ATOMLIST*,double);

void set_rotatable_bonds(MOLECULE*);

SYSMOLS* load_sysmols(SETTINGS*);

unsigned int dof_get_space(char*);
void set_dof_list_shifts(DOF_LIST*,double);
void apply_dof_list_shifts(DOF_LIST*);
void apply_dof_shift(DOF*,double);
DOF_LIST* setup_dof_list(SYSTEM*,unsigned int);
void free_dof(DOF*);
void free_dof_list(DOF_LIST*);
void store_positions(ATOMLIST*,double**,double**,double**,double**);
void restore_positions(ATOMLIST*,double**,double**,double**,double**);
ATOM_COORDS* get_list_coords(ATOMLIST*);
void set_list_coords(ATOMLIST*,ATOM_COORDS*);

void run_help(SETTINGS*);
void params_init(int nargs,char **args);
PLI_PARAM* params_get_parameter(char*);

