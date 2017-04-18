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



#include "pli.h"



#define IDEAL_HBOND_DIST 2.80
#define IDEAL_CHO_DIST 3.25
#define HBOND_DIST_TOL1 0.35
#define HBOND_DIST_TOL2 0.60
#define HBOND_ALPHA_TOL1 30.0
#define HBOND_ALPHA_TOL2 60.0
#define HBOND_BETA_TOL1 30.0
#define HBOND_BETA_TOL2 60.0
#define IDEAL_METAL_DIST 2.15
#define METAL_DIST_TOL1 0.35
#define METAL_DIST_TOL2 0.65
#define IDEAL_HBOND_ENERGY -10.0
#define HBOND_MISMATCH_ENERGY 10.0
#define IDEAL_METAL_COORD_ENERGY -20
#define HBOND_POLAR_AREA_ENERGY -0.5
//#define HBOND_CONTACT 1
//#define METAL_CONTACT 2



static double **sphere_points = NULL;



static void alloc_atom_coordination(ATOM*);
static void init_atom_coordination(ATOM_COORDINATION*);
static void write_coordination_atom(FILE*,ATOM*);
//static int hbond_type_match(ATOM*,ATOM*);
static double calc_exposed_polar_area(ATOM*);



void calc_atom_coordination(ATOM *atom1) {

  int i,j;
  double score,cscore[MAX_Z],sumf,area;
  unsigned int flags1,flags2,cflags;
  ATOM *atom2;
  ATOM_TYPE *type1,*type2;
  CONTACT *contact;
  CONTACTLIST *list;
  COORDINATION tmp_item,*item,*item1,*item2;
  ATOM_COORDINATION *coord;

  // test this is an atom we can/want to do this on:

  type1 = atom1->type;

  if (!type1) {
    
    return;
  }

  flags1 = type1->flags;

  if ((!(flags1 & ANY_HBOND_ATOM_TYPE)) && (!(flags1 & METAL_ATOM_TYPE))) {

    return;
  }

  list = atom1->contactlist;

  if (list == NULL) {

    return;
  }

  // allocate and initialise:

  alloc_atom_coordination(atom1);

  coord = atom1->coordination;

  init_atom_coordination(coord);

  coord->exposed_polar_area = calc_exposed_polar_area(atom1);

  // find the coordinating contacts:

  for (i=0,contact=list->contacts;i<list->ncontacts;i++,contact++) {

    cflags = contact->flags;

    if (!(cflags & COVALENT_CONTACT)) {

      if (cflags & INTRAMOLECULAR_CONTACT) {

	coord->intra_area += contact->area;

      } else {

	coord->inter_area += contact->area;
      }

      atom2 = contact->atom2;
      type2 = atom2->type;


      if (type2) {

	flags2 = type2->flags;

	if (hbond_flags_match(flags1,flags2,1)) { 

	  score = hbond_geometry_score(contact);

	  if ((score > 1.0E-6) && (coord->max_z < MAX_Z)) {

	    if (score > 0.2) {

	      if ((flags2 & HBOND_DA_ATOM_TYPE) == HBOND_DA_ATOM_TYPE) {

		coord->n_da++;

	      } else if (flags2 & HBOND_DONOR_ATOM_TYPE) {

		coord->n_don++;

	      } else if ((flags2 & HBOND_ACCEPTOR_ATOM_TYPE) || (flags2 & METAL_ACCEPTOR_ATOM_TYPE)) {

		coord->n_acc++;

	      } else if (flags2 & CHBOND_DONOR_ATOM_TYPE) {

		coord->n_ch++;

	      } else if (flags2 & METAL_ATOM_TYPE) {
		
		coord->n_metal++;
	      }
	    }

	    if (coord->max_z >= MAX_Z) {

	      error_fn("calc_atom_coordination: coordination number exceeds %d",MAX_Z);
	    }
	  
	    item = coord->list + coord->max_z;
	    
	    item->s = score;
	    item->contact = contact;
	    
	    coord->max_z++;
	    
	    if ((atom2->error_flags) || (atom2->occupancy < 0.99)) {
	      
	      coord->error_flags = 1;
	    }
	  }
	}
      }
    }
  }

  coord->max_z--;

  coord->z = (coord->n_acc) + (coord->n_don) + (coord->n_da) + (coord->n_ch) + (coord->n_metal);

  // sort coordinating atoms:
  
  for (i=1;i<coord->max_z;i++) {

    item1 = coord->list + i;

    for (j=i+1;j<=coord->max_z;j++) {

      item2 = coord->list + j;

      if (item2->s > item1->s) {

	memcpy(&tmp_item,item1,sizeof(COORDINATION));
	memcpy(item1,item2,sizeof(COORDINATION));
	memcpy(item2,&tmp_item,sizeof(COORDINATION));
      }
    }
  }

  for (i=0,item=coord->list;i<=coord->max_z;i++,item++) {

    cscore[i] = (i == 0) ? item->s : (item->s)*cscore[i-1];
  }

  sumf = 0.0;

  for (i=coord->max_z;i>=0;i--) {

    item = coord->list + i;

    item->f = (1.0-sumf)*cscore[i];

    sumf += item->f;
  }
}



double atom_hbond_energy(ATOM *atom,int include_exposed) {

  int i,j,metal_acc;
  double E,Q,metal_score;
  ATOM_TYPE *type,*type2;
  ALT_ATOM_TYPE *alt_type;
  CONTACTLIST *contactlist;
  CONTACT *contact;

  contactlist = atom->contactlist;

  metal_acc = 0;

  if ((contactlist) && (contactlist->ncontacts)) {

    type = atom->type;

    if (type->flags & HBOND_ACCEPTOR_ATOM_TYPE) {

      // check if atom coordinates a metal

      for (i=0,contact=contactlist->contacts;i<contactlist->ncontacts;i++,contact++) {

	type2 = contact->atom2->type;

	if ((type2) && (type2->flags & METAL_ATOM_TYPE)) {

	  if (hbond_geometry_score(contact) > 0.2) {

	    metal_acc = 1;

	  }
	}
      }
    }

    Q = 0.0;

    if (type->n_alt_types) {

      for (i=0,alt_type=type->alt_types;i<type->n_alt_types;i++,alt_type++) {

	if (alt_type->group_status == atom->group_status) {

	  atom->type = alt_type->type;

	  E = 0.0;

	  for (j=0,contact=contactlist->contacts;j<contactlist->ncontacts;j++,contact++) {

	    type2 = contact->atom2->type;

	    if (((type2) && (type2->flags & METAL_ATOM_TYPE)) || (!metal_acc)) {

	      E += hbond_energy(contact);
	    }
	  }

	  if (include_exposed) {

	    E += HBOND_POLAR_AREA_ENERGY*calc_exposed_polar_area(atom);
	  }

	  Q += (alt_type->probability)*exp(-E/RT);
	}
      }

      atom->type = type;

    } else {

      E = 0.0;

      for (i=0,contact=contactlist->contacts;i<contactlist->ncontacts;i++,contact++) {

	type2 = contact->atom2->type;

	if (((type2) && (type2->flags & METAL_ATOM_TYPE)) || (!metal_acc)) {

	  E += hbond_energy(contact);
	}
      }

      if (include_exposed) {

	E += HBOND_POLAR_AREA_ENERGY*calc_exposed_polar_area(atom);
      }

      Q += exp(-E/RT);
    }

  } else {

    E = (include_exposed) ? (HBOND_POLAR_AREA_ENERGY*calc_exposed_polar_area(atom)) : 0.0;

    Q = exp(-E/RT);
  }
  
  return(-RT*log(Q));
}



double hbond_energy(CONTACT *contact) {

  int i,j,hbond_type;
  double geom_score,Q;
  ATOM *atom1,*atom2;
  ATOM_TYPE *type1,*type2;
  ALT_ATOM_TYPE *alt_type1,*alt_type2;

  atom1 = contact->atom1;
  atom2 = contact->atom2;

  type1 = atom1->type;
  type2 = atom2->type;

  // atoms must both be H-bonding types:

  if ((!(type1->flags & ANY_HBOND_ATOM_TYPE)) && (!(type1->flags & METAL_ACCEPTOR_ATOM_TYPE)) && (!(type1->flags & METAL_ATOM_TYPE))) {

    return(0.0);
  }

  if ((!(type2->flags & ANY_HBOND_ATOM_TYPE)) && (!(type2->flags & METAL_ACCEPTOR_ATOM_TYPE)) && (!(type2->flags & METAL_ATOM_TYPE))) {

    return(0.0);
  }

  // metal acceptors must match metals:

  if ((type1->flags & METAL_ACCEPTOR_ATOM_TYPE) && (!(type1->flags & HBOND_ACCEPTOR_ATOM_TYPE)) && (!(type2->flags & METAL_ATOM_TYPE))) {

    return(0.0);
  }

  if ((type2->flags & METAL_ACCEPTOR_ATOM_TYPE) && (!(type2->flags & HBOND_ACCEPTOR_ATOM_TYPE)) && (!(type1->flags & METAL_ATOM_TYPE))) {

    return(0.0);
  }

  // allow metal-metal contacts (for metal clusters):

  if ((type1->flags & METAL_ATOM_TYPE) && (type2->flags & METAL_ATOM_TYPE)) {

    return(0.0);
  }

  geom_score = hbond_geometry_score(contact);

  Q = 0.0;

  if (type1->n_alt_types == 0) {

    if (type2->n_alt_types == 0) {

      hbond_type = hbond_flags_match(type1->flags,type2->flags,1);

      if (hbond_type == METAL_COORDINATING_CONTACT) {

	Q += exp(-(IDEAL_METAL_COORD_ENERGY*geom_score)/RT);

      } else if (hbond_type == HBOND_CONTACT) {

	Q += exp(-(IDEAL_HBOND_ENERGY*geom_score)/RT);

      } else {

	Q += exp(-(HBOND_MISMATCH_ENERGY*geom_score)/RT);
      }

    } else {

      for (i=0,alt_type2=type2->alt_types;i<type2->n_alt_types;i++,alt_type2++) {

	if (alt_type2->group_status == atom2->group_status) {

	  hbond_type = hbond_flags_match(type1->flags,alt_type2->type->flags,1);

	  if (hbond_type == METAL_COORDINATING_CONTACT) {

	    Q += (alt_type2->probability)*exp(-(IDEAL_METAL_COORD_ENERGY*geom_score)/RT);

	  } else if (hbond_type == HBOND_CONTACT) {

	    Q += (alt_type2->probability)*exp(-(IDEAL_HBOND_ENERGY*geom_score)/RT);

	  } else {

	    Q += (alt_type2->probability)*exp(-(HBOND_MISMATCH_ENERGY*geom_score)/RT);
	  }
	}
      }
    }

  } else {

    if (type2->n_alt_types == 0) {

      for (i=0,alt_type1=type1->alt_types;i<type1->n_alt_types;i++,alt_type1++) {

	if (alt_type1->group_status == atom1->group_status) {

	  hbond_type = hbond_flags_match(type2->flags,alt_type1->type->flags,1);

	  if (hbond_type == METAL_COORDINATING_CONTACT) {

	    Q += (alt_type1->probability)*exp(-(IDEAL_METAL_COORD_ENERGY*geom_score)/RT);

	  } else if (hbond_type == HBOND_CONTACT) {

	    Q += (alt_type1->probability)*exp(-(IDEAL_HBOND_ENERGY*geom_score)/RT);

	  } else {

	    Q += (alt_type1->probability)*exp(-(HBOND_MISMATCH_ENERGY*geom_score)/RT);
	  }
	}
      }

    } else {

      for (i=0,alt_type1=type1->alt_types;i<type1->n_alt_types;i++,alt_type1++) {

	if (alt_type1->group_status == atom1->group_status) {

	  for (j=0,alt_type2=type2->alt_types;j<type2->n_alt_types;j++,alt_type2++) {

	    if (alt_type2->group_status == atom2->group_status) {

	      hbond_type = hbond_flags_match(alt_type1->type->flags,alt_type2->type->flags,1);

	      if (hbond_type == METAL_COORDINATING_CONTACT) {

		Q += (alt_type1->probability)*(alt_type2->probability)*exp(-(IDEAL_METAL_COORD_ENERGY*geom_score)/RT);

	      } else if (hbond_type == HBOND_CONTACT) {

		Q += (alt_type1->probability)*(alt_type2->probability)*exp(-(IDEAL_HBOND_ENERGY*geom_score)/RT);

	      } else {

		Q += (alt_type1->probability)*(alt_type2->probability)*exp(-(HBOND_MISMATCH_ENERGY*geom_score)/RT);
	      }
	    }
	  }
	}
      }
    }
  }

  return(-RT*log(Q));
}



double coord_score(ATOM_COORDINATION *coord,int scale) {

  int i;
  double score;
  COORDINATION *item;

  score = 0.0;

  if (coord->max_z) {

    for (i=1,item=coord->list+1;i<=coord->max_z;i++,item++) {

      score += (item->s);
    }
  }

  return(score);
}



void write_coordination_system(FILE *file,SYSTEM *system) {

  int i;
  ATOM **atomp,*atom;
  ATOMLIST *list;

  list = system->selection;

  if (list == NULL) {

    error_fn("write_coordination_system: no atoms");
  }

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom = *atomp;

    if (!atom->contactlist) {

      set_contacts_atom(atom,system,1);
    }

    write_coordination_atom(file,atom);
  }
}



static void alloc_atom_coordination(ATOM *atom) {

  if (atom->coordination != NULL) {

    return;
  }

  atom->coordination = (ATOM_COORDINATION*) malloc(sizeof(ATOM_COORDINATION));

  if (atom->coordination == NULL) {

    error_fn("alloc_atom_coordination: out of memory allocating atom coordination");
  }
}



static void init_atom_coordination(ATOM_COORDINATION *coord) {

  int i;
  COORDINATION *item;

  coord->n_acc = 0;
  coord->n_don = 0;
  coord->n_da = 0;
  coord->n_ch = 0;
  coord->n_metal = 0;
  coord->z = 0;

  coord->max_z = 1;

  item = coord->list;

  item->s = 1.0;
  item->f = 0.0;

  item++;

  for (i=1;i<MAX_Z;i++,item++) {

    item->s = 0.0;
    item->f = 0.0;
  }

  coord->intra_area = 0.0;
  coord->inter_area = 0.0;
  coord->exposed_polar_area = 0.0;

  coord->error_flags = 0;
}



static void write_coordination_atom(FILE *file,ATOM *atom) {

  int i,n_don,n_acc,n_da,n_ch,n_metal,sym_atom;
  double inter_area,intra_area,score;
  unsigned int flags;
  ATOM_TYPE *type;
  CONTACT *contact;
  CONTACTLIST *list;

  if ((atom->type == NULL) || (atom->element == NULL)) {

    return;
  }

  if ((!(atom->type->flags & ANY_HBOND_ATOM_TYPE)) && (!(atom->element->flags & METAL_ELEMENT))) {

    return;
  }

  list = atom->contactlist;

  if (list == NULL) {

    return;
  }

  n_don = n_acc = n_da = n_ch = n_metal = 0;

  inter_area = intra_area = 0.0;

  for (i=0,contact=list->contacts;i<list->ncontacts;i++,contact++) {

    if (contact->flags & INTRAMOLECULAR_CONTACT) {

      if (!(contact->flags & COVALENT_CONTACT)) {

	intra_area += contact->area;
      }

    } else {

      inter_area += contact->area;
    }

    calc_contact_geometry(contact);

    if (hbond_type_match(contact->atom1,contact->atom2)) { 

      if (hbond_geometry_score(contact) > 0.2) {

	type = contact->atom2->type;

	if (type != NULL) {

	  flags = type->flags;

	  if ((flags & HBOND_DA_ATOM_TYPE) == HBOND_DA_ATOM_TYPE) {

	    n_da++;

	  } else if (flags & HBOND_DONOR_ATOM_TYPE) {

	    n_don++;

	  } else if (flags & HBOND_ACCEPTOR_ATOM_TYPE) {

	    n_acc++;

	  } else if (flags & CHBOND_DONOR_ATOM_TYPE) {

	    n_ch++;
	  }
	}
      }

    } else if (metal_coordinating_contact(contact,0.8)) {

      if (atom->element->flags & METAL_ELEMENT) {

	n_acc++;

      } else {

	n_metal++;
      }
    }
  }

  sym_atom = (atom->molecule->flags & SYMMETRY_MOLECULE) ? 1 : 0;

  fprintf(file,"COORD %6d %3d    %3d %3d %3d %3d %3d %10.4lf %10.4lf\n",atom->id,sym_atom,n_don,n_acc,n_da,n_ch,n_metal,intra_area,inter_area);
}



int hbond_type_match(ATOM *atom1,ATOM *atom2) {

  unsigned int flags1,flags2;

  if ((atom1->type == NULL) || (atom2->type == NULL)) {

    return(0);
  }

  flags1 = atom1->type->flags;
  flags2 = atom2->type->flags;

  return(hbond_flags_match(flags1,flags2,1));
}



int hbond_flags_match(unsigned int flags1,unsigned int flags2,int allow_chos) {

  if ((flags1 & METAL_ATOM_TYPE) && ((flags2 & HBOND_ACCEPTOR_ATOM_TYPE) || (flags2 & METAL_ACCEPTOR_ATOM_TYPE))) {

    return(METAL_COORDINATING_CONTACT);
  }

  if ((flags2 & METAL_ATOM_TYPE) && ((flags1 & HBOND_ACCEPTOR_ATOM_TYPE) || (flags1 & METAL_ACCEPTOR_ATOM_TYPE))) {

    return(METAL_COORDINATING_CONTACT);
  }

  if ((!(flags1 & ANY_HBOND_ATOM_TYPE)) || (!(flags2 & ANY_HBOND_ATOM_TYPE))) {

    return(0);  
  }

  if ((flags1 & HBOND_DONOR_ATOM_TYPE) && (flags2 & HBOND_ACCEPTOR_ATOM_TYPE)) {

    return(HBOND_CONTACT);
  }

  if ((flags1 & HBOND_ACCEPTOR_ATOM_TYPE) && (flags2 & HBOND_DONOR_ATOM_TYPE)) {

    return(HBOND_CONTACT);
  }

  if (!allow_chos) {

    return(0);
  }

  if ((flags1 & CHBOND_DONOR_ATOM_TYPE) && (flags2 & HBOND_ACCEPTOR_ATOM_TYPE)) {

    return(HBOND_CONTACT);
  }

  if ((flags1 & HBOND_ACCEPTOR_ATOM_TYPE) && (flags2 & CHBOND_DONOR_ATOM_TYPE)) {

    return(HBOND_CONTACT);
  }

  return(0);
}



double hbond_geometry_score(CONTACT *contact) {

  unsigned int flags1,flags2;
  double fp1,fp2,score;
  ATOM *atom1,*atom2;

  atom1 = contact->atom1;
  atom2 = contact->atom2;

  flags1 = atom1->type->flags;
  flags2 = atom2->type->flags;

  fp1 = fraction_polar_point(atom1,atom2->position,contact->alpha1,0);

  fp2 = fraction_polar_point(atom2,atom1->position,contact->alpha2,0);

  if ((flags1 & METAL_ATOM_TYPE) || (flags2 & METAL_ATOM_TYPE)) {

    score = fp1*fp2*hbond_block_score(contact->distance,IDEAL_METAL_DIST,METAL_DIST_TOL1,METAL_DIST_TOL2);

  } else if ((flags1 & CHBOND_DONOR_ATOM_TYPE) || (flags2 & CHBOND_DONOR_ATOM_TYPE)) {

    score = fp1*fp2*hbond_block_score(contact->distance,IDEAL_CHO_DIST,HBOND_DIST_TOL1,HBOND_DIST_TOL2);

  } else {

    score = fp1*fp2*hbond_block_score(contact->distance,IDEAL_HBOND_DIST,HBOND_DIST_TOL1,HBOND_DIST_TOL2);
  }

  return(score);
}



double hbond_block_score(double x,double ideal_x,double dx1,double dx2) {

  double dx;

  dx = fabs(x-ideal_x);

  if (dx < dx1) {

    return(1.0);
  }

  if (dx < dx2) {

    return(1.0-((dx-dx1)/(dx2-dx1)));
  }

  return(0.0);
}



double fraction_polar_point(ATOM *atom,double *point,double alpha,int local) {

  int geomtype;
  double x[4],y[4],v[4],l,alpha1,alpha2,f1,f2;
  HBOND_GEOMETRY *hbond_geometry;

  hbond_geometry = atom->hbond_geometry;

  if (hbond_geometry == NULL) {

    return(0.0);
  }

  geomtype = hbond_geometry->type;

  if (geomtype == POINT_HBOND_GEOMETRY) {

    // all of the atom is polar:

    return(1.0);

  } else if (geomtype == VECTOR_HBOND_GEOMETRY) {

    // use alpha to determine polar fraction:

    if (alpha < -0.001) {

      alpha = vector_angle(atom->u,point);
    }

    return(hbond_block_score(alpha,hbond_geometry->angle1,hbond_geometry->tol11,hbond_geometry->tol12));

  } else if (geomtype == WEDGE_HBOND_GEOMETRY) {

    // calculate projected angles wrt u and v:

    if (local) {

      copy_vector(point,x);

    } else {

      calc_vector(atom->position,point,x);
    }

    copy_vector(atom->v,v);

    l = dotproduct(x,v);

    scale_vector(v,l);
    
    calc_vector(v,x,y);

    alpha1 = vector_angle(y,atom->u);
    alpha2 = vector_angle(y,x);

    f1 = hbond_block_score(alpha1,hbond_geometry->angle1,hbond_geometry->tol11,hbond_geometry->tol12);
    f2 = hbond_block_score(alpha2,hbond_geometry->angle2,hbond_geometry->tol21,hbond_geometry->tol22);

    return(f1*f2);

  } else if (geomtype == VECTORPAIR_HBOND_GEOMETRY) {

    // calculate angle wrt both virtual points:

    if (local) {

      copy_vector(point,x);

    } else {

      calc_vector(atom->position,point,x);
    }

    alpha1 = vector_angle(x,atom->vpts[0]);
    alpha2 = vector_angle(x,atom->vpts[1]);

    if (alpha1 < alpha2) {

      return(hbond_block_score(alpha1,0.0,hbond_geometry->tol11,hbond_geometry->tol12));

    } else {

      return(hbond_block_score(alpha2,0.0,hbond_geometry->tol11,hbond_geometry->tol12));
    }
  }

  error_fn("fraction_polar: no such hbond geometry '%s'",hbond_geometry->name);
}



static double calc_exposed_polar_area(ATOM *atom1) {

  int i,j,buried;
  double *pos,**pointp,*point,fp,R1,R2,dist,planedist,pointdist,alpha,v[4];
  double sum_polar_exposed,area,exposed_polar_area;
  ATOM *atom2;
  CONTACTLIST *contactlist;
  CONTACT *contact;
  HBOND_GEOMETRY *hbond_geometry;

  hbond_geometry = atom1->hbond_geometry;

  if (!hbond_geometry) {

    return(0.0);
  }

  if (atom1->exposed_area < 0.001) {

    return(0.0);
  }

  // load lebedev sphere points if necessary:

  if (!sphere_points) {

    sphere_points = read_lebedev_sphere_points(110);
  }

  // loop over sphere points:

  pos = atom1->position;
  R1 = atom1->vdw_radius_H2O;

  sum_polar_exposed = 0.0;

  for (i=0,pointp=sphere_points;i<110;i++,pointp++) {

    point = *pointp;

    // test if sphere point is polar:

    fp = fraction_polar_point(atom1,point,-1,1);

    if (fp > 0.0001) {

      // test if it is buried by one of the contacts:

      contactlist = atom1->contactlist;

      if (contactlist) {

	buried = 0;

	for (j=0,contact=contactlist->contacts;(j<contactlist->ncontacts) && (!buried);j++,contact++) {

	  atom2 = contact->atom2;
	  R2 = atom2->vdw_radius_H2O;
	  dist = contact->distance;

	  calc_vector(pos,atom2->position,v);

	  alpha = vector_angle(v,point);

	  if (alpha < 89.99) {

	    planedist = (sqr(dist) + sqr(R1) - sqr(R2))/(2*dist);

	    pointdist = planedist/cos(alpha*(PI/180.0));

	    if (pointdist < R1) {

	      buried = 1;
	    }
	  }
	}

	if (!buried) {

	  sum_polar_exposed += fp;
	}
      }
    }
  }

  area = 4*PI*R1*R1;

  exposed_polar_area = area*(sum_polar_exposed/110);

  return(exposed_polar_area);
}
