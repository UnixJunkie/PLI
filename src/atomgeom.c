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



#include "pli.h"



static HBOND_GEOMETRY hbond_geometries[] = { { POINT_HBOND_GEOMETRY,      "POINT",       0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 }, // e.g. Metals
					     { VECTOR_HBOND_GEOMETRY,     "VECTOR1",     0,  0.0, 55.0, 65.0,  0.0,  0.0,  0.0 }, // e.g. Amide NH
					     { VECTOR_HBOND_GEOMETRY,     "CONE1",       0, 70.5, 55.0, 65.0,  0.0,  0.0,  0.0 }, // e.g. -NH3
					     { VECTOR_HBOND_GEOMETRY,     "CAP1",        0,  0.0, 90.0,110.0,  0.0,  0.0,  0.0 }, // e.g. =O
					     { VECTORPAIR_HBOND_GEOMETRY, "VECTORPAIR1", 2, 60.0, 40.0, 50.0,  0.0,  0.0,  0.0 }, // e.g. -NH2p
					     { VECTORPAIR_HBOND_GEOMETRY, "VECTORPAIR2", 2, 54.2, 55.0, 65.0,  0.0,  0.0,  0.0 }, // e.g. >NH2
					     { WEDGE_HBOND_GEOMETRY,      "WEDGE1",      0,  0.0,115.0,125.0,  0.0,  40.0,50.0 },
                                             { -1,                        "LAST",        0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 } };



static double contact_angle(double*,double*,int);
static double projected_contact_angle(double*,double*,double*,int);
static void alloc_atom_geometries(ATOM_TYPING_SCHEME*);
static void init_atom_geometry(ATOM_GEOMETRY*);
static int check_atom_geometry(ATOM*,ATOM_GEOMETRY*);
static int check_atom_hbond_geometry(ATOM*,HBOND_GEOMETRY*);
static int set_atom_axes(ATOM*);



void set_molecule_atom_geometries(MOLECULE *molecule,ATOM_TYPING_SCHEME *scheme) {

  int i;
  ATOM *atom;

  for (i=0,atom=molecule->atom;i<molecule->natoms;i++,atom++) {

    set_atom_geometry(atom,scheme);
  }
}



void set_atom_geometry(ATOM *atom,ATOM_TYPING_SCHEME *scheme) {

  ATOM_TYPE *type;
  ATOM_GEOMETRY *geometry;
  HBOND_GEOMETRY *hbond_geometry;

  type = atom->type;

  if ((type != NULL) && (atom->connections != NULL)) {

    geometry = type->geometry;

    if (check_atom_geometry(atom,geometry)) {

      atom->geometry = geometry;

      if (set_atom_axes(atom)) {

        hbond_geometry = type->hbond_geometry;

        if (hbond_geometry) {

    	  if (check_atom_hbond_geometry(atom,hbond_geometry)) {

	    atom->hbond_geometry = hbond_geometry;

	    set_atom_virtual_points(atom);

	  } else {

	    warning_fn("set_atom_geometry: couldn't derive hbond geometry for atom %d (%s,%d,%s,%s,%s)",
	  	       atom->id,atom->subname,atom->subid,atom->molecule->name,atom->geometry->name,hbond_geometry->name);
          }
	}

        return;
      }
    }
  }

  atom->error_flags |= ATOM_GEOMETRY_ERROR;

  atom->geometry = get_atom_geometry("X",scheme);

  set_atom_axes(atom);

  warning_fn("set_atom_geometry: unable to assign atom geometry for atom %d %s %s %d %s %s",
	     atom->id,atom->name,atom->subname,atom->subid,atom->molecule->name,atom->type->name);
}



HBOND_GEOMETRY* get_hbond_geometry(char *name) {

  int i;
  HBOND_GEOMETRY *hbond_geometry;

  if (!strcmp(name,"NONE")) {

    return(NULL);
  }

  hbond_geometry = hbond_geometries;

  while (hbond_geometry->type != -1) {

    if (!strcmp(name,hbond_geometry->name)) {

      return(hbond_geometry);
    }

    hbond_geometry++;
  }

  error_fn("get_hbond_geometry: no such hbond geometry '%s'",name);
}



void calc_contact_geometries_system(SYSTEM *system) {

  int i,j;
  ATOM **atomp,*atom;
  CONTACT *contact;
  CONTACTLIST *clist;

  for (i=0,atomp=system->selection->atom;i<system->selection->natoms;i++,atomp++) {

    atom = *atomp;

    clist = atom->contactlist;

    for (j=0,contact=clist->contacts;j<clist->ncontacts;j++,contact++) {

      calc_contact_geometry(contact);
    }
  }
}



void calc_contact_geometry(CONTACT *contact) {

  int symmetrical;
  double v12[4],v21[4];
  ATOM *atom1,*atom2;
  ATOM_GEOMETRY *geometry1,*geometry2;

  atom1 = contact->atom1;
  atom2 = contact->atom2;

  geometry1 = atom1->geometry;
  geometry2 = atom2->geometry;

  if ((geometry1 == NULL) || (geometry2 == NULL)) {

    error_fn("calc_contact_geometry: atom geometry undefined");
  }

  calc_vector(atom1->position,atom2->position,v12);

  contact->distance = vector_length(v12);

  invert_vector(v12,v21);

  if (geometry1->u_axis) {

    contact->alpha1 = contact_angle(v12,atom1->u,geometry1->u_symm);

    if (geometry1->v_axis) {

      contact->beta1 = projected_contact_angle(v12,atom1->v,atom1->w,geometry1->v_symm);
    }
  }

  if (geometry2->u_axis) {

    contact->alpha2 = contact_angle(v21,atom2->u,geometry2->u_symm);

    if (geometry2->v_axis) {

      contact->beta2 = projected_contact_angle(v21,atom2->v,atom2->w,geometry2->v_symm);
    }
  }
}



void mirror_contact_geometry(CONTACT *contact1,CONTACT *contact2) {

  contact2->distance = contact1->distance;

  contact2->alpha1 = contact1->alpha2;
  contact2->beta1 = contact1->beta2;

  contact2->alpha2 = contact1->alpha1;
  contact2->beta2 = contact1->beta1;

  contact2->area = contact1->iarea;
  contact2->iarea = contact1->area;
}



void write_contact_geometry(PLI_FILE *file,CONTACT *contact,enum OUTPUT_FORMAT oformat,unsigned long int oflags) {

  char alpha1[11],beta1[11],alpha2[11],beta2[11];
  ATOM_GEOMETRY *geometry1,*geometry2;

  geometry1 = contact->atom1->geometry;
  geometry2 = contact->atom2->geometry;

  if ((geometry1 == NULL) || (geometry2 == NULL)) {

    error_fn("write_contact_geometry: atom geometry undefined");
  }

  strcpy(alpha1,"null");
  strcpy(beta1,"null");
  strcpy(alpha2,"null");
  strcpy(beta2,"null");

  if (geometry1->u_axis) {

    sprintf(alpha1,"%.4lf",contact->alpha1);
  }
    
  if (geometry1->v_axis) {
      
    sprintf(beta1,"%.4lf",contact->beta1);
  }
      
  if (geometry2->u_axis) {
      
    sprintf(alpha2,"%.4lf",contact->alpha2);
  }
  
  if (geometry2->v_axis) {
    
    sprintf(beta2,"%.4lf",contact->beta2);
  }
  
  write_line(file,(oformat == JSON) ? 
	     "\"alpha1\":%s,\"beta1\":%s,\"alpha2\":%s,\"beta2\":%s" : 
	     " %10s %10s %10s %10s",
	     alpha1,beta1,alpha2,beta2);
}



static double contact_angle(double *v1,double *v2,int symmetry) {

  double angle;

  angle = vector_angle(v1,v2);

  if ((symmetry == 2) && (angle > 90.0)) {

    return(180.0 - angle);
  }

  return(angle);
}



static double projected_contact_angle(double *v,double *v_axis,double *w_axis,int symmetry) {

  double y,z,angle;
  
  y = dotproduct(v,v_axis);
  z = dotproduct(v,w_axis);

  if (fabs(y) < 1.0E-30) {

    angle = (z < 0.0) ? -90.0 : 90.0;

  } else {

    angle = (180.0/PI)*atan(z/y);

    if (y < 0.0) {

      angle = (z < 0.0) ? angle - 180.0 : 180.0 + angle;
    }
  }

  if (symmetry == 2) {

    angle = fabs(angle);

  } else if (symmetry == 4) {

    angle = fabs(angle);

    if (angle > 90) {

      angle = 180.0 - angle;
    }
  }

  return(angle);
}



static int check_atom_geometry(ATOM *atom1,ATOM_GEOMETRY *geometry) {

  ATOM *atom2;
  ATOMLIST *conns1,*conns2;

  conns1 = atom1->connections;

  if ((geometry->nbonds1 == UNDEFINED) || (geometry->nbonds1 == conns1->natoms)) {

    if ((geometry->planar1 == UNDEFINED) || (geometry->planar1 == planar_atom(atom1))) {

      if (conns1->natoms == 0) {

	if ((geometry->nbonds2 == UNDEFINED) && (geometry->planar2 == UNDEFINED)) {
	    
	  return(1);
	}

	return(0);
      }

      atom2 = *(conns1->atom);

      conns2 = atom2->connections;

      if (conns2 != NULL) {

	// note that here nbonds needs to be equal to or greater than the nbonds in the geometry definition

	if ((geometry->nbonds2 == UNDEFINED) || (geometry->nbonds2 <= conns2->natoms)) {

	  if ((geometry->planar2 == UNDEFINED) || (geometry->planar2 == planar_atom(atom2))) {

	    return(1);
	  }
	}
      }
    }
  }

  return(0);
}



static int check_atom_hbond_geometry(ATOM *atom,HBOND_GEOMETRY *hbond_geometry) {

  int id,u_axis,v_axis;
  ATOM_GEOMETRY *geometry;

  if (hbond_geometry == NULL) {

    return(0);
  }

  id = hbond_geometry->type;

  // point hbond geometry:

  if (id == POINT_HBOND_GEOMETRY) {

    return(1);
  }

  geometry = atom->geometry;

  if (geometry == NULL) {

    return(0);
  }

  u_axis = geometry->u_axis;
  v_axis = geometry->v_axis;

  // vector or cone hbond geometry (need u axis):

  if (id == VECTOR_HBOND_GEOMETRY) {

    if (!u_axis) {

      return(0);
    }
  }

  // vector pair hbond geometry (need u and v axis):

  if (id == VECTORPAIR_HBOND_GEOMETRY) {

    if ((!u_axis) || (!v_axis)) {

      return(0);
    }
  }

  return(1);
}



static int set_atom_axes(ATOM *atom) {

  double v1[4],v2[4],v3[4],v4[4],v5[4],v6[4],v7[4];
  ATOM *atom1,*atom2;
  ATOMLIST *conns;
  ATOM_GEOMETRY *geometry;

  geometry = atom->geometry;

  if (geometry->nbonds1 == 1) {

    if (geometry->u_axis) {

      atom1 = *(atom->connections->atom);

      calc_vector(atom1->position,atom->position,atom->u);

      scale_vector(atom->u,1.0);

      if (geometry->v_axis) {

	conns = atom1->connections;

	if (conns->natoms == 2) {

	  atom2 = (*(conns->atom) == atom) ? *(conns->atom+1) : *(conns->atom);

	  calc_vector(atom1->position,atom2->position,v1);

	  calc_crossproduct(v1,atom->u,v2);
	  calc_crossproduct(atom->u,v2,atom->v);

	  scale_vector(atom->v,1.0);

	} else if (conns->natoms == 3) {

	  if (geometry->planar2 == 1) {

	    calc_vector(atom1->position,(*(conns->atom))->position,v1);
	    calc_vector(atom1->position,(*(conns->atom+1))->position,v2);
	    calc_vector(atom1->position,(*(conns->atom+2))->position,v3);

	    calc_crossproduct(v1,v2,v4);
	    calc_crossproduct(v2,v3,v5);
	    calc_crossproduct(v3,v1,v6);

	    scale_vector(v4,1.0);
	    scale_vector(v5,1.0);
	    scale_vector(v6,1.0);

	    sum_vector(v4,v5,v7);
	    sum_vector(v7,v6,atom->v);

	    scale_vector(atom->v,1.0);

	  } else {

	    calc_vector((*(conns->atom))->position,atom1->position,v1);
	    calc_vector((*(conns->atom+1))->position,atom1->position,v2);
	    calc_vector((*(conns->atom+2))->position,atom1->position,v3);

	    scale_vector(v1,1.0);
	    scale_vector(v2,1.0);
	    scale_vector(v3,1.0);

	    sum_vector(v1,v2,v4);
	    sum_vector(v4,v3,v5);

	    calc_crossproduct(v5,atom->u,v6);
	    calc_crossproduct(atom->u,v6,atom->v);

	    scale_vector(atom->v,1.0);
	  }

	} else {

          return(0);
	}
      }
    }

  } else if (geometry->nbonds1 == 2) {

    if (geometry->u_axis) {

      conns = atom->connections;

      calc_vector((*(conns->atom))->position,atom->position,v1);
      calc_vector((*(conns->atom+1))->position,atom->position,v2);

      scale_vector(v1,1.0);
      scale_vector(v2,1.0);

      sum_vector(v1,v2,atom->u);

      scale_vector(atom->u,1.0);

      if (geometry->v_axis) {

	calc_crossproduct(v1,v2,atom->v);

	scale_vector(atom->v,1.0);
      }
    }

  } else if (geometry->nbonds1 == 3) {

    if (geometry->u_axis) {

      if (geometry->planar1 == 1) {

	calc_vector(atom->position,(*(atom->connections->atom))->position,v1);
	calc_vector(atom->position,(*(atom->connections->atom+1))->position,v2);
	calc_vector(atom->position,(*(atom->connections->atom+2))->position,v3);

	calc_crossproduct(v1,v2,v4);
	calc_crossproduct(v2,v3,v5);
	calc_crossproduct(v3,v1,v6);

	scale_vector(v4,1.0);
	scale_vector(v5,1.0);
	scale_vector(v6,1.0);

	sum_vector(v4,v5,v7);
	sum_vector(v7,v6,atom->u);

	scale_vector(atom->u,1.0);

      } else {

	calc_vector((*(atom->connections->atom))->position,atom->position,v1);
	calc_vector((*(atom->connections->atom+1))->position,atom->position,v2);
	calc_vector((*(atom->connections->atom+2))->position,atom->position,v3);
	
	scale_vector(v1,1.0);
	scale_vector(v2,1.0);
	scale_vector(v3,1.0);

	sum_vector(v1,v2,v4);
	sum_vector(v4,v3,atom->u);

	scale_vector(atom->u,1.0);
      }
    }
  }

  if (geometry->v_axis) {

    // make sure u and v are at a 90degr angle:

    calc_crossproduct(atom->u,atom->v,v1);
    calc_crossproduct(v1,atom->u,atom->v);

    scale_vector(atom->v,1.0);

    // create w:

    calc_crossproduct(atom->u,atom->v,atom->w);

    scale_vector(atom->w,1.0);
  }

  return(1);
}



void set_atom_virtual_points(ATOM *atom) {

  double angle,v0[4],m[4][4];
  ATOM_GEOMETRY *geometry;
  HBOND_GEOMETRY *hbond_geometry;

  geometry = atom->geometry;
  hbond_geometry = atom->hbond_geometry;

  if ((geometry == NULL) || (hbond_geometry == NULL)) {

    return;
  }

  if (hbond_geometry->type == VECTORPAIR_HBOND_GEOMETRY) {

    // need to generate pair of virtual points:

    angle = hbond_geometry->angle1;

    null_vector(v0);

    copy_vector(atom->u,atom->vpts[0]);
    copy_vector(atom->u,atom->vpts[1]);

    if (geometry->nbonds1 == 1) {

      // e.g. -NH2p hydrogens or =O lone pairs - coplanar with u-w plane:

      rotation_matrix(v0,atom->v,-angle,m);
      transform_vector(atom->vpts[0],m);

      rotation_matrix(v0,atom->v,angle,m);
      transform_vector(atom->vpts[1],m);

    } else if (geometry->nbonds1 == 2) {

      // e.g. >NH2 hydrogens or >O lone pairs - coplanar with u-v plane:

      rotation_matrix(v0,atom->w,-angle,m);
      transform_vector(atom->vpts[0],m);

      rotation_matrix(v0,atom->w,angle,m);
      transform_vector(atom->vpts[1],m);
    }
  }
}



ATOM_GEOMETRY* get_atom_geometry(char *name,ATOM_TYPING_SCHEME *scheme) {

  int i;
  ATOM_GEOMETRY *geometry;

  for (i=0,geometry=scheme->atom_geometries;i<scheme->n_atom_geometries;i++,geometry++) {

    if (!strcmp(geometry->name,name)) {

      return(geometry);
    }
  }

  error_fn("get_atom_geometry: no such atom geometry '%s'",name);
}




void read_atom_geometries(PLI_FILE *paramfile,ATOM_TYPING_SCHEME *scheme) {

  int flag,n_read;
  char line[MAX_LINE_LEN],word[MAX_LINE_LEN];
  ATOM_GEOMETRY *geom;

  strcpy(word,"");

  while ((!end_of_file(paramfile)) && (strcmp(word,"end_atom_geometries"))) {
  
    if (read_line(line,MAX_LINE_LEN,paramfile) == NULL)
      break;

    flag = sscanf(line,"%s",word);

    if (flag != EOF) {

      if (strcmp(word,"end_atom_geometries")) {

	alloc_atom_geometries(scheme);

	geom = scheme->atom_geometries + scheme->n_atom_geometries;

	init_atom_geometry(geom);

	n_read = sscanf(line,"%d \"%[^\"]\" %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf",
			&(geom->id),geom->name,
			&(geom->nbonds1),&(geom->planar1),
			&(geom->nbonds2),&(geom->planar2),
			&(geom->u_axis),&(geom->u_symm),
			&(geom->v_axis),&(geom->v_symm),
			&(geom->r_limits[0]),&(geom->r_limits[1]),
			&(geom->alpha_limits[0]),&(geom->alpha_limits[1]),
			&(geom->beta_limits[0]),&(geom->beta_limits[1]));

	if (n_read == 16) {

	  scheme->n_atom_geometries++;

	} else {

	  warning_fn("read_atom_geometries: skipping corrupted atom geometry line:\n%s",line);
	}
      }

    } else {

      strcpy(word,"");
    }
  }
}



static void alloc_atom_geometries(ATOM_TYPING_SCHEME *scheme) {

  if (scheme->n_alloc_atom_geometries == 0) {

    scheme->n_alloc_atom_geometries = 20;

    scheme->atom_geometries = (ATOM_GEOMETRY*) calloc(scheme->n_alloc_atom_geometries,sizeof(ATOM_GEOMETRY));

    if (scheme->atom_geometries == NULL) {

      error_fn("alloc_atom_geometries: out of memory allocating atom_geometries");
    }

  } else if (scheme->n_atom_geometries == scheme->n_alloc_atom_geometries) {

    scheme->n_alloc_atom_geometries += 20;

    scheme->atom_geometries = (ATOM_GEOMETRY*) realloc(scheme->atom_geometries,scheme->n_alloc_atom_geometries*sizeof(ATOM_GEOMETRY));

    if (scheme->atom_geometries == NULL) {

      error_fn("alloc_atom_geometries: out of memory re-allocating atom_geometries");
    }
  }
}



static void init_atom_geometry(ATOM_GEOMETRY *geom) {

  strcpy(geom->name,"");

  geom->id = -1;

  geom->nbonds1 = -1;
  geom->planar1 = -1;

  geom->nbonds2 = -1;
  geom->planar2 = -1;

  geom->u_axis = 0;
  geom->v_axis = 0;

  geom->u_symm = 0;
  geom->v_symm = 0;

  geom->r_limits[0] = 0.0;
  geom->r_limits[1] = 0.0;

  geom->alpha_limits[0] = 0.0;
  geom->alpha_limits[1] = 0.0;

  geom->beta_limits[0] = 0.0;
  geom->beta_limits[1] = 0.0;
}
