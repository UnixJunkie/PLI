// Copyright 2015 Brendan McConkey and Astex Therapeutics Ltd.
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



/**********************************************************************
 * The code in this is file was adapted directly from the Vcontacts.c
 * program, version 2.5, Dec 2014, by Brendan McConkey, which
 * calculates the atom-atom contacts within a protein, along with the
 * solvent accessible surface, and was described in the literature:
 *
 * Quantification of protein surfaces, volumes and atom-atom contacts
 * using a constrained Voronoi procedure.
 * McConkey, B.J., Sobolev, V., and Edelman, M.
 * Bioinformatics 18(10), 2002, 1365-1373
 *
 * ====================================================================

 * NOTES ON METHODOLOGY:
 *
 * The method calculates the planes of contact the center atom makes
 * will each atom in the contact list, then determines the points
 * of intersection of the planes and intersections with a sphere. The
 * sphere has a radius equal to the sum of the van der Waals radius of the
 * center atom plus the radius of a solvent atom (water). The set of
 * intersections of the planes defines a polygon surrounding the center
 * atom. This polygon is projected onto the surface of the sphere. The
 * area of each projection is calculated as the sum of spherical triangles
 * and arc segments, where an arc segment is the difference between a
 * angluar segment of a spherical cap and the corresponding spherical triangle.
 *
 * The points of intersection are calculated using a convex hull algorithm.
 * The algorithm operates with efficiency O(Nk), where N is the number of
 * input atoms and k is the number of edges on the contact polyhedron.
 *
 * The convex hull algorithm effectively describes the atom contacts only if
 * no 'engulfing' atoms are present, ie., atoms which cover more than 50% of
 * the atom contact surface. This may occur between closely bonded atoms of
 * different radii, such as some C=O bonds. A correction factor accounts for
 * this - points on the polyhedron with an engulfing contact are projected
 * to the surface of the sphere from the center of the engulfing plane,
 * preserving solvent accessible surface.
 **********************************************************************/



#include "pli.h"



static int voronoi_poly2(ATOM*,VERTEX*,VERTEX*,PLANE*,EDGEVECTOR*,unsigned int*);
static void copy_voronoi_details(ATOM*,CONTACTLIST*,PLANE*);
static char order_faces(ATOM*,int,VERTEX*,VERTEX*,PLANE*,PTINDEX*);
static void calc_areas(ATOM*,int,VERTEX*,VERTEX*,PLANE*,PTINDEX*);

static void get_firstvert(struct plane cont[], int *planeA, int *planeB, int *planeC, int NC);

static int solve_3x3(double eq0[], double eq1[], double eq2[], double pt[]);
static int solve_2xS(struct plane eq0, struct plane eq1, float rado, double pt0[], double pt1[]);
static int add_vertex(struct vertex poly[], int vn, double coor[], int planeA, int planeB, int planeC);
static void add_vedge(struct edgevector vedge[], int edgenum, struct plane cont[], int plane0, int plane1, 
		      int testplane, struct vertex poly[], int startpt);
static char test_point(double ptX[], struct plane cont[], int NC, float rado, int planeA, int planeB, int planeC);
static double cosPQR( double ptP[], double ptQ[], double ptR[]);
static void project_points(struct vertex poly[], struct vertex centerpt[], float rado, int NC, 
			   int NV, struct plane cont[]);
static double spherical_arc(struct vertex ptAo, struct vertex ptB, struct vertex ptC, float rado);



void contacts2voronoi(ATOM *atom) {

  int i,NV;
  unsigned int error_flag;
  char surfatom;         // atom type, 'I' internal, 'S' surface
  double rado;
  PTINDEX ptorder[100];  // for ordering vertices around each face
  VERTEX centerpt[100];  // center points for each contact face
  VERTEX poly[200];      // polyhedron vertices
  PLANE cont[100];       // atom and contact plane information
  EDGEVECTOR vedge[200];
  CONTACTLIST *contactlist;
  CONTACT *contacts,*contact;

  error_flag = 0;

  rado = atom->vdw_radius_H2O;

  // set solvent exposed and buried area for atom:

  atom->contact_area = 0.0;
  atom->intra_area = 0.0;
  atom->covalent_area = 0.0;

  atom->exposed_area = 4.0*PI*rado*rado;

  contactlist = atom->contactlist;

  if (contactlist->ncontacts == 0) {

    return;
  }

  if (contactlist->ncontacts > 100) {

    atom->error_flags |= ATOM_VORONOI_ERROR;

    warning_fn("contacts2voronoi: cannot calculate voronoi polyhedra for atom %d; not enough space allocated",atom->id);

    return;
  }

  NV = voronoi_poly2(atom,poly,centerpt,cont,vedge,&error_flag);

  if (error_flag) {

    atom->error_flags = ATOM_VORONOI_ERROR;

    warning_fn("contacts2voronoi: error %d produced in calculating polyhedron for atom %d - all contacts ignored",error_flag,atom->id);

    return;
  }

  surfatom = order_faces(atom,NV,poly,centerpt,cont,ptorder);

  if (surfatom == 'E') {

    warning_fn("contacts2voronoi: failed to order faces for atom %d",atom->id);

    atom->error_flags |= ATOM_VORONOI_ERROR;
   
    return;
  }

  calc_areas(atom,NV,poly,centerpt,cont,ptorder);

  copy_voronoi_details(atom,contactlist,cont);

  remove_flagged_contacts(contactlist,NO_VORONOI_CONTACT);

  if (contactlist->ncontacts) {

    contacts = contactlist->contacts;

    for (i=0,contact=contacts;i<contactlist->ncontacts;i++,contact++) {

      if (contact->flags & COVALENT_CONTACT) {

	atom->covalent_area += contact->area;

      } else if (contact->flags & INTRAMOLECULAR_CONTACT) {

	atom->intra_area += contact->area;

      } else {

	atom->contact_area += contact->area;
      }
    }
  }

  atom->exposed_area = 4.0*PI*rado*rado - atom->contact_area - atom->intra_area - atom->covalent_area;

  if (atom->exposed_area < 0.0) {

    if (atom->exposed_area < -0.01) {

      atom->error_flags |= ATOM_VORONOI_ERROR;

      warning_fn("contacts2voronoi: negative exposed area for atom %d (was %.4lf, setting to zero)",atom->id,atom->exposed_area);
    }

    atom->exposed_area = 0.0;
  }
}



double estimate_contact_iarea(CONTACT *contact) {

  double R1,R2;

  R1 = contact->atom1->vdw_radius;
  R2 = contact->atom2->vdw_radius;

  return((contact->area)*(R1/R2));
}



static void copy_voronoi_details(ATOM *atom,CONTACTLIST *list,PLANE *cont) {

  int i;
  CONTACT *contacts,*contact;
  
  contacts = list->contacts;

  for (i=0,contact=contacts;i<list->ncontacts;i++,contact++) {

    contact->area = cont[i].area;

    if ((cont[i].area < 1.0E-10) || (cont[i].flag == 'X')) {

      contact->flags |= NO_VORONOI_CONTACT;
    }
  }
}



/**************************
 * subroutine calc_areas
 * updated dec2014 BJM
 **************************/
 
 // given a polyhedron surrounding an atom, calculate the area of each
 // face projected onto the surface of a sphere radius = rad(atom) + rad(water)

static void calc_areas(ATOM *atom,int NV,VERTEX *poly,VERTEX *centerpt,PLANE *cont,PTINDEX *ptorder) {

  char   engflag;         // ='Y' if an engulfing plane is present
  int    planeX;          // current plane, 1 to NC
  int    NP;              // number of points on face
  int    vi;              // vertices counter
  double area;            // area of polygon
  int    pa1, pa2;        // possible common planes
  int    commplane;       // plane shared by adjacent points (not planeX)
  int    epi;             // engulfing plane counter
  int    engplane[4];     // index to engulfed planes
  double maxSAS;          // maximum solvent exposed surface, no atoms other than engulfing.
  struct vertex ptB, ptC; // arc point intersections
  int    currpt, nextpt;
  double cosNN1[40];      // angle between vertex N and vertex N+1
  double cosNzero[40];    // angle between vertex N and vertex N+1
  double tanprod;         // product of tangents
  int    v0, va, vb, vc;  // vertices for arc calculation
  double U,V,W,X;	        // calculations for spherical triangles
  double tansqrS, tansqrSA, tansqrSB, tansqrSC;  // calculations for spherical triangles
  double CMid[3], MidV[3]; // vectors for determining major or minor arc
  CONTACTLIST *contactlist;
  CONTACT *contact;
  double rado,crado;
  int NC;
  ATOM *catom;
  int i;

  rado = atom->vdw_radius_H2O;
  contactlist = atom->contactlist;
  NC = contactlist->ncontacts;

  engflag = 'N';
  epi = 0;

  //RESET AREAS TO ZERO, see if any planes are engulfing
  for(planeX=0; planeX < NC; ++planeX) {
    cont[planeX].area = 0.0;
    if(cont[planeX].flag == 'E') {   	// plane is engulfing, center point outside polyhedron
      engflag = 'Y';
      engplane[epi] = planeX;
      ++epi;
    }
  }

  if(engflag == 'Y') { // engulfing plane correction - project points onto sphere surface.
    // TODO: check!
    if (atom->id == 0) {

      //printf("contacts atom = %5d (%s,%d) n_contacts = %d\n",atom->id,atom->subname,atom->subid,NC);

      for(i=0,contact=contactlist->contacts; i < NC; ++i,contact++) {

	//printf("contact to atom %5d %-5s %5s %5d %10.4lf\n",contact->atom2->id,contact->atom2->name,contact->atom2->subname,contact->atom2->subid,contact->distance);

	if (contact->distance < 1.0) {

	  error_fn("%s: too short %s %s %s - %s %s %s : %.4lf\n",__func__,
		   contact->atom1->molecule->name,contact->atom1->name,contact->atom1->subname,
		   contact->atom2->molecule->name,contact->atom2->name,contact->atom2->subname,contact->distance);
	}
      }
    }

    project_points(poly, centerpt, rado, NC, NV, cont);
  }

  /* ---------------------------- */
  /* calculate area for each face */
  /* ---------------------------- */

  for(planeX=0,contact=contactlist->contacts; planeX < NC; ++planeX,contact++) {

    catom = contact->atom2;
    crado = catom->vdw_radius_H2O;

    NP = ptorder[planeX].numpts;
    area = 0.0;
    if(cont[planeX].flag == 'X') {
      continue;
    }

    // if there are no points on a valid contact, area is spherical cap
    if(NP == 0) {
      // some contacts are shielded by another atom, need to check if behind another contact plane
      if(test_point(centerpt[planeX].xi, cont, NC, rado, planeX, -1, -1) == 'Y') {
	cont[planeX].area = 2.0*PI*rado*(rado-centerpt[planeX].dist);
      }
    } else if(NP == 2) {  // only two contact points, check which part of arc
      if(test_point(centerpt[planeX].xi, cont, NC, rado, planeX, -1, -1) == 'Y') { // area is (cap - arc)
	cont[planeX].area = 2.0*PI*rado*(rado-centerpt[planeX].dist)
	  - spherical_arc(centerpt[planeX], poly[ptorder[planeX].pt[0]],
			  poly[ptorder[planeX].pt[1]], rado);
      } else {            // area is arc.
	cont[planeX].area = spherical_arc(centerpt[planeX], poly[ptorder[planeX].pt[0]],
					  poly[ptorder[planeX].pt[1]], rado);
      }
    } else {  // three or more points define face
      
      // ------ calculate cosines and angles ------
      for(vi=0; vi<NP; ++vi) {
	v0 = ptorder[planeX].pt[0];
	va = ptorder[planeX].pt[vi];
	vb = ptorder[planeX].pt[(vi+1)%NP];
	
	// calculate cosines between adjacent vertices
	cosNN1[vi] = (( poly[va].xi[0]*poly[vb].xi[0] + poly[va].xi[1]*poly[vb].xi[1]
			+ poly[va].xi[2]*poly[vb].xi[2] ) / (poly[va].dist*poly[vb].dist));
	
	// calculate cosines between vertex zero and vertex 'vi'
	if(vi != 0) {
	  cosNzero[vi] = ((poly[v0].xi[0]*poly[va].xi[0] + poly[v0].xi[1]*poly[va].xi[1]
			   + poly[v0].xi[2]*poly[va].xi[2] ) / (poly[v0].dist*poly[va].dist));
	}
      }
      
      // ----- calculate area of triangles in face -----
      for(vi=1; vi<(NP-1); ++vi) {
	U = sqrt((1+cosNzero[vi])*(1+cosNN1[vi])*(1+cosNzero[vi+1])/8.0);
	V = sqrt((1-cosNzero[vi])*(1-cosNN1[vi])*(1+cosNzero[vi+1])/8.0);
	W = sqrt((1-cosNzero[vi])*(1+cosNN1[vi])*(1-cosNzero[vi+1])/8.0);
	X = sqrt((1+cosNzero[vi])*(1-cosNN1[vi])*(1-cosNzero[vi+1])/8.0);
	tansqrS  = (1-U+V+W+X)/(1+U-V-W-X);
	tansqrSA = (1-U-V-W+X)/(1+U+V+W-X);
	tansqrSB = (1-U-V+W-X)/(1+U+V-W+X);
	tansqrSC = (1-U+V-W-X)/(1+U-V+W+X);
	tanprod = sqrt(tansqrS*tansqrSA*tansqrSB*tansqrSC);
	if(tanprod > 0.0) {
	  area += 4.0*rado*rado*atan(sqrt(tanprod));
	} 
      }
      
      // ----- add area of arc segments  -----
      for(vi=0; vi<NP; ++vi) {
	va = ptorder[planeX].pt[vi];
	vb = ptorder[planeX].pt[(vi+1)%NP];
	vc = ptorder[planeX].pt[(vi+2)%NP];
	
	//check if adjacent points are arc segments
	if((poly[va].plane[2] == -1) && (poly[vb].plane[2] == -1)) {
	  // if on different planes, area includes an arc segment
	  if((poly[va].plane[0]+poly[va].plane[1]) != (poly[vb].plane[0]+poly[vb].plane[1])) {
	    // BJM dec2014 
	    // instead of adding extra point to arc, test for major or minor arc
	    // 1. find midpoint between two arc points
	    // 2. calc vector CMid[3] from plane centerpt to line midpoint
	    // 3. calc vector MidV[3] from midpoint to any other vertex on face
	    // 4. calc dot product of CMid and MidV, if negative it's a minor arc 
	    CMid[0] = (poly[va].xi[0]+poly[vb].xi[0])/2.0 - centerpt[planeX].xi[0];
	    CMid[1] = (poly[va].xi[1]+poly[vb].xi[1])/2.0 - centerpt[planeX].xi[1];
	    CMid[2] = (poly[va].xi[2]+poly[vb].xi[2])/2.0 - centerpt[planeX].xi[2];
	    MidV[0] = poly[vc].xi[0] - (poly[va].xi[0]+poly[vb].xi[0])/2.0;
	    MidV[1] = poly[vc].xi[1] - (poly[va].xi[1]+poly[vb].xi[1])/2.0;
	    MidV[2] = poly[vc].xi[2] - (poly[va].xi[2]+poly[vb].xi[2])/2.0;
	    if ((CMid[0]*MidV[0] + CMid[1]*MidV[1] +CMid[2]*MidV[2]) < 0.0) {
	      // negative dot product is minor arc
	      area += spherical_arc(centerpt[planeX], poly[va], poly[vb], rado);
	    } else {
	      // positive dot product is major arc							
	      area += 2.0*PI*rado*(rado-centerpt[planeX].dist)
		- spherical_arc(centerpt[planeX], poly[va], poly[vb], rado);
	    }					
	  }
	}
      }
      cont[planeX].area = area;
    }
  }
  
  // --------------------------------------------------------
  //  add correction terms for engulfing planes, if required
  // --------------------------------------------------------
  
  if(engflag == 'Y') {
    for(planeX=0; planeX < NC; ++planeX) {
      if(cont[planeX].flag != 'E') {
	continue;
      }
      
      NP = ptorder[planeX].numpts;
      for(vi=0; vi<NP; ++vi) {
	currpt = ptorder[planeX].pt[vi];
	nextpt = ptorder[planeX].pt[(vi+1)%NP];
	
	// find common second plane, if any.
	if(poly[currpt].plane[0] == planeX) {
	  pa1 = poly[currpt].plane[1];
	  pa2 = poly[currpt].plane[2];
	} else {
	  pa1 = poly[currpt].plane[0];
	  if(poly[currpt].plane[1] == planeX) {
	    pa2 = poly[currpt].plane[2];
	  } else {
	    pa2 = poly[currpt].plane[1];
	  }
	}
	
	if((pa1 == poly[nextpt].plane[0]) || (pa1 == poly[nextpt].plane[1])
	   || (pa1 == poly[nextpt].plane[2])) {
	  commplane = pa1;
	} else if((pa2 == poly[nextpt].plane[0]) || (pa2 == poly[nextpt].plane[1])
		  || (pa2 == poly[nextpt].plane[2])) {
	  commplane = pa2;
	} else {
	  continue;
	}
	if((commplane != -1) && (cont[commplane].flag != 'E')) {
	  // add correction to commplane area. here centerpt is from engulfing plane.
	  cont[commplane].area += spherical_arc(centerpt[planeX], poly[currpt], poly[nextpt], rado);
	  if(NP == 2) break;  // otherwise would repeat adding area
	}
      }
    }
    
    // -----------------------------------------------
    // ------ calculate engulfed contact areas -------
    // -----------------------------------------------
    
    if(epi == 1) {
      cont[engplane[0]].area = 2.0*PI*rado*(rado+cont[engplane[0]].dist);
    } else if(epi == 2) {
      if(solve_2xS(cont[engplane[0]],cont[engplane[1]], rado, ptB.xi, ptC.xi)== -1) {
	cont[engplane[0]].area = 2.0*PI*rado*rado;
	cont[engplane[1]].area = 2.0*PI*rado*rado;
      } else {
	ptB.dist = rado;
	ptC.dist = rado;
	maxSAS = spherical_arc(centerpt[engplane[0]], ptB, ptC, rado);
	maxSAS += spherical_arc(centerpt[engplane[1]], ptB, ptC, rado);
	cont[engplane[0]].area = 2.0*PI*rado*rado - 0.5*maxSAS;
	cont[engplane[1]].area = 2.0*PI*rado*rado - 0.5*maxSAS;
      }
    } else if(epi>=3) {
      // no exposed surface if there are three or more engulfing contacts
      for(planeX=0; planeX<NC; ++planeX) {
	if(cont[planeX].flag == 'E') {
	  cont[planeX].area = 4.0*PI*rado*rado/epi;
	} else {
	  cont[planeX].area = 0.0;
	}
      }
    }
  }
  return;
}

/*********************************
 * function order_faces
 *********************************/

// output is array ptorder, the order of points around each face.
// return values are 'I' for internal atom, 'S' for surface atom.

static char order_faces(ATOM *atom,int NV,VERTEX *poly,VERTEX *centerpt,PLANE *cont,PTINDEX *ptorder) {

  int planeX;             // current plane, 1 to NC
  int planeY;             // second plane, to get adjacent points on polygon
  int surfcount;          // number of points defining a given plane of contact
  int vi, vi2;            // vertices counter
  int tempsi;             // temporary storage for exchanging indices
  double tempcos;         // for exchanging cosines
  double cos10X[50];      // cos of angle between points poly[1],poly[0],and poly[X]
  //double temppt[3];       // temp coordinates of new arc point
  char surfatom;          // return value: 'I' internal atom, 'S' surface atom
  CONTACTLIST *contactlist;
  double rado;
  int NC;

  rado = atom->vdw_radius_H2O;
  contactlist = atom->contactlist;
  NC = contactlist->ncontacts;

  surfatom = 'I'; // internal atom
  for(vi=0; vi<NV; ++vi) {
    if(poly[vi].plane[2] == -1) {
      surfatom = 'S'; // surface atom
      break;
    }
  }

  // for surface calculation only
  // if(surfatom == 'I') return('I');

  for(planeX=0; planeX < NC; ++planeX) {
    if(cont[planeX].flag == 'X') { // hidden
      ptorder[planeX].numpts = 0;
      continue;
    }
    
    surfcount=0;
    for(vi=0; vi<NV; ++vi) {
      // index all points comprising surface for planeX
      if((poly[vi].plane[0]==planeX) || (poly[vi].plane[1]==planeX) || (poly[vi].plane[2]==planeX)) {
	ptorder[planeX].pt[surfcount] = vi;
	++surfcount;
      }
    }
    
    ptorder[planeX].numpts = surfcount;
    
    if(surfcount > 3) {
      // get two points on same line (two common planes).
      // all points already share one plane (planeX), find another.
      if(poly[ptorder[planeX].pt[0]].plane[0] == planeX) {
	planeY = poly[ptorder[planeX].pt[0]].plane[1];
      } else {
	planeY = poly[ptorder[planeX].pt[0]].plane[0];
      }
      
      //find another point on the same line (2 common planes)
      for(vi=1; vi<surfcount; ++vi) {
	if((poly[ptorder[planeX].pt[vi]].plane[0]==planeY) || (poly[ptorder[planeX].pt[vi]].plane[1]==planeY)
	   || (poly[ptorder[planeX].pt[vi]].plane[2]==planeY)) {

	  break;
	}
      }
      
      // MLV 02-03-2015 return here if failed to find another point on the same line:

      if (vi == surfcount) {

	return('E');
      }

      //swap index for pt[1] and pt[vi], so points 0 and 1 are on same line
      tempsi = ptorder[planeX].pt[vi];
      ptorder[planeX].pt[vi] = ptorder[planeX].pt[1];
      ptorder[planeX].pt[1] = tempsi;

      // calculate cosine between points indexed 1,0,X
      for(vi=2; vi<surfcount; ++vi) {
	cos10X[vi] = cosPQR(poly[ptorder[planeX].pt[1]].xi, poly[ptorder[planeX].pt[0]].xi,
			    poly[ptorder[planeX].pt[vi]].xi);
      }
      
      // order by cosines, decreasing order
      for(vi=2; vi<surfcount-1; ++vi) {
	for(vi2=vi+1; vi2<surfcount; ++vi2) {
	  if(cos10X[vi] < cos10X[vi2]) {
	    // swap indices if points in wrong order
	    tempsi = ptorder[planeX].pt[vi];
	    ptorder[planeX].pt[vi] = ptorder[planeX].pt[vi2];
	    ptorder[planeX].pt[vi2] = tempsi;
	    tempcos = cos10X[vi];
	    cos10X[vi] = cos10X[vi2];
	    cos10X[vi2] = tempcos;
	  }
	}
      }
    }
  }
  return(surfatom);
}




/******************************
 * subroutine voronoi_poly2
 * created 08/07/2001  BJM
 ******************************/

static int voronoi_poly2(ATOM *atom,VERTEX *poly,VERTEX *centerpt,PLANE *cont,EDGEVECTOR *vedge,unsigned int *error_flag) {

  int    cai;          // contact atom counter
  ATOM *ca_ptr;        // pointer to pdb atom
  double atomdist;     // distance to atom
  double planedist;    // distance to plane
  double mindist;      // distance to closest plane
  int    planeA;       // closest plane to origin
  int    planeB;       // second plane, with planeA defines closest edge
  int    planeC;       // new intersection plane for edge (endpt)
  int    oldplaneC;    // old intersection plane for edge (startpt)
  double vt;           // vector parameter 't' for line x=x'+lt, y=y'+mt, z=z'+nt
  double vtmin;        // minimum value for vt
  double vtdiv;        // check for division by zero in vt calculation
  double temppt[3];
  struct vertex *stp; // pointer to start vertex, coordinates
  double *V;          // pointer to edge vector
  int    startedge;
  int    edgenum;
  int    vn = 0;
  char   edgeflag;
  int    edgei;      // edge counter
  int    vi, vj;     // vertices counters
  double arcpt0[3], arcpt1[3];
  int    testpA, testpB;
  double testvalA, testvalB;
  char   arcflag = 'N';
  
  // failsafe variables:
  char   recalc;       // flag if hull is being recalculated (orig. unbounded)
  float  origcoor[3];  // original pdb coordinates for atom.
  CONTACT *contact;
  CONTACTLIST *contactlist;
  double rad,rado;
  int NC;
  char pdbstr[MAX_LINE_LEN];
  int i;

  rad = atom->vdw_radius;
  rado = atom->vdw_radius_H2O;
  contactlist = atom->contactlist;
  NC = contactlist->ncontacts;

  recalc = 'N';
RESTART:
  planeA = -1;
  planeB = -1;
  planeC = -1;

  /* generate planes of contact with D = planedist */
  mindist = 9.9e+9;
  for(cai=0,contact=contactlist->contacts; cai<contactlist->ncontacts; ++cai,contact++) {
    ca_ptr = contact->atom2;
    atomdist = contact->distance;

    if(planedef == 'B') {  // bisection - original Voronoi procedure
      planedist = atomdist/2.0;
    } else if(planedef == 'R') { // radical plane (Gellatly and Finney) - default.
      planedist = (atomdist*atomdist + (rad*rad) - (ca_ptr->vdw_radius)*(ca_ptr->vdw_radius))/(2*atomdist);
    } else { // extended radical plane (McConkey et al), matches solvent accessible surface
      planedist = (atomdist*atomdist + (rado*rado) - (ca_ptr->vdw_radius_H2O)*(ca_ptr->vdw_radius_H2O))/(2*atomdist);
    }

    cont[cai].Ai[0] = (ca_ptr->position[0] - atom->position[0])/atomdist;
    cont[cai].Ai[1] = (ca_ptr->position[1] - atom->position[1])/atomdist;
    cont[cai].Ai[2] = (ca_ptr->position[2] - atom->position[2])/atomdist;
    cont[cai].Ai[3] = -planedist;
    cont[cai].dist  = fabs(planedist);
    cont[cai].atom = ca_ptr;
    cont[cai].flag = 'X'; // initialize contact flags to 'no contact'

    // set planeA as closest plane
    if(cont[cai].dist < mindist) {
      mindist = cont[cai].dist;
      planeA = cai;
    }
  }

  // add four planes surrounding atom, outer limit for voronoi polyhedron
  // changed to unit vectors, added distances - dec2014 BJM 
  cont[NC].Ai[0] = 0.7071;
  cont[NC].Ai[1] = 0.7071;
  cont[NC].Ai[2] = 0.0;
  cont[NC].Ai[3] = -10.0;
  cont[NC].dist = 10.0;
  cont[NC+1].Ai[0] = 0.7071;
  cont[NC+1].Ai[1] = -0.7071;
  cont[NC+1].Ai[2] = 0.0;
  cont[NC+1].Ai[3] = -10.0;
  cont[NC+1].dist = 10.0;
  cont[NC+2].Ai[0] = -0.7071;
  cont[NC+2].Ai[1] = 0.0;
  cont[NC+2].Ai[2] = 0.7071;
  cont[NC+2].Ai[3] = -10.0;
  cont[NC+2].dist = 10.0;
  cont[NC+3].Ai[0] = -0.7071;
  cont[NC+3].Ai[1] = 0.0;
  cont[NC+3].Ai[2] = -0.7071;
  cont[NC+3].Ai[3] = -10.0;
  cont[NC+3].dist = 10.0;
  // get starting vertex set of planes
  get_firstvert(cont, &planeA, &planeB, &planeC, NC);
  // calculate intersection point of three planes
  solve_3x3(cont[planeA].Ai, cont[planeB].Ai, cont[planeC].Ai, temppt);
  
  // add first vertex to vertex list
  add_vertex(poly, 0, temppt, planeA, planeB, planeC);
  
  // flag contacts as present
  cont[planeA].flag = 'Y';
  cont[planeB].flag = 'Y';
  cont[planeC].flag = 'Y';
  
  // calculate edge vectors
  add_vedge(vedge, 0, cont, planeA, planeB, planeC, poly, 0);
  add_vedge(vedge, 1, cont, planeB, planeC, planeA, poly, 0);
  add_vedge(vedge, 2, cont, planeC, planeA, planeB, poly, 0);

  startedge = 0;
  edgenum = 3;
  vn = 1;

  /* --------------------------------------------------- */
  /* Generate new polyhedron points from edge vectors    */
  /* --------------------------------------------------- */

  while(1) {
    // get next unfinished vector = startedge

    while(((edgenum-startedge) > 0) && (vedge[startedge].endpt >= 0)) {
      ++startedge;
    }
    if((edgenum-startedge) <= 0) {
      // all edges are done, polyhedron complete.
      break;
    }
    
    vtmin = 9.9e+9; // dummy value
    stp = &poly[vedge[startedge].startpt];
    V = vedge[startedge].V;
    planeA = vedge[startedge].plane[0];
    planeB = vedge[startedge].plane[1];
    oldplaneC = vedge[startedge].startplane;
    planeC = -1;
    
    // get closest positive intersection point
    for(cai=0; cai<NC+4; ++cai) {
      if((cai != planeA) && (cai != planeB) && (cai != oldplaneC)) {
	vtdiv = (cont[cai].Ai[0]*V[0] +cont[cai].Ai[1]*V[1] +cont[cai].Ai[2]*V[2]);
	if(vtdiv != 0.0) {
	  vt = -(cont[cai].Ai[0]*stp->xi[0] +cont[cai].Ai[1]*stp->xi[1]
		 +cont[cai].Ai[2]*stp->xi[2] +cont[cai].Ai[3])/vtdiv;
	  if((vt < vtmin) && (vt > 0)) {
	    vtmin = vt;
	    planeC = cai;
	  }
	}
      }
    }
    poly[vn].xi[0] = stp->xi[0] + vtmin*V[0];
    poly[vn].xi[1] = stp->xi[1] + vtmin*V[1];
    poly[vn].xi[2] = stp->xi[2] + vtmin*V[2];
    
    add_vertex(poly, vn, poly[vn].xi, planeA, planeB, planeC);
    vedge[startedge].endpt = vn;
    vedge[startedge].endplane = planeC;
    
    //flag contact as present
    cont[planeC].flag = 'Y';
    
    // ========  ADD EDGES  ========
    
    // check edge (planeA, planeC)
    edgeflag = 'Y';
    edgei = startedge+1;
    while(edgei < edgenum) {
      if(((vedge[edgei].plane[0] == planeA)&&(vedge[edgei].plane[1] == planeC)) ||
	 ((vedge[edgei].plane[0] == planeC)&&(vedge[edgei].plane[1] == planeA))) {
	// already on list, add current vertex as endpt
	vedge[edgei].endpt = vn;
	vedge[edgei].endplane = planeB;
	edgeflag = 'N';
	break;
      }
      ++edgei;
    }
    if(edgeflag == 'Y') { // add edge
      add_vedge(vedge, edgenum, cont, planeA, planeC, planeB, poly, vn);
      ++edgenum;
    }
    
    // check edge (planeB, planeC)
    edgeflag = 'Y';
    edgei = startedge+1;
    while(edgei < edgenum) {
      if(((vedge[edgei].plane[0] == planeB)&&(vedge[edgei].plane[1] == planeC)) ||
	 ((vedge[edgei].plane[0] == planeC)&&(vedge[edgei].plane[1] == planeB))) {
	// already on list, add current vertex as endpt
	vedge[edgei].endpt = vn;
	vedge[edgei].endplane = planeA;
	edgeflag = 'N';
	break;
      }
      ++edgei;
    }
    if(edgeflag == 'Y') { // add edge
      add_vedge(vedge, edgenum, cont, planeB, planeC, planeA, poly, vn);
      ++edgenum;
      
      // ===== failsafe - if solution is not converging, perturb atom  =====
      // ===== coordinates and recalculate.                            =====
      if(edgenum >= 200) {
	warning_fn("solution did not converge for atom %d, recalculating\n", atom->id);
	origcoor[0] = atom->position[0];
	origcoor[1] = atom->position[1];
	origcoor[2] = atom->position[2];
	// perturb atom coordinates
	atom->position[0] += 0.005*(float)(2*rand()-RAND_MAX)/(float)RAND_MAX;
	atom->position[1] += 0.005*(float)(2*rand()-RAND_MAX)/(float)RAND_MAX;
	atom->position[2] += 0.005*(float)(2*rand()-RAND_MAX)/(float)RAND_MAX;

	recalc = 'Y';
	
	goto RESTART;
      }
    }
    ++vn;
  }
  
  /*--------------------------------------------------*/
  /*  now have voronoi polyhedron around given atom.  */
  /*  remove vertices outside of sphere, and          */
  /*  calculate intersection points with sphere.      */
  /*--------------------------------------------------*/
  
  // flag edges that may cross sphere boundary
  for(edgei=0; edgei<edgenum; ++edgei) {
    if((rado < poly[vedge[edgei].startpt].dist) || (rado < poly[vedge[edgei].endpt].dist)) {
      // one or both vertices fall outside of sphere
      arcflag = 'Y';
      vedge[edgei].arc = '?';
    } else {
      vedge[edgei].arc = 'X';
    }
  }
  
  // calculate new arc points
  for(edgei=0; edgei<edgenum; ++edgei) {
    if(vedge[edgei].arc != 'X') {
      if(solve_2xS(cont[vedge[edgei].plane[0]], cont[vedge[edgei].plane[1]], rado, arcpt0, arcpt1) == -1) {
	vedge[edgei].arc = 'X';  // mark edge as no associated arc point
	continue;
      }
      // test new arc points vs. adjacent planes, add if ok.
      testpA = vedge[edgei].startplane;
      testpB = vedge[edgei].endplane;
      // test arc point 0
      testvalA = cont[testpA].Ai[0]*arcpt0[0] + cont[testpA].Ai[1]*arcpt0[1]
	+ cont[testpA].Ai[2]*arcpt0[2] + cont[testpA].Ai[3];
      testvalB = cont[testpB].Ai[0]*arcpt0[0] + cont[testpB].Ai[1]*arcpt0[1]
	+ cont[testpB].Ai[2]*arcpt0[2] + cont[testpB].Ai[3];
      if((testvalA < 0.0) && (testvalB < 0.0)) { // point is good
	add_vertex(poly, vn, arcpt0, vedge[edgei].plane[0], vedge[edgei].plane[1], -1);
	poly[vn].dist = rado;
	++vn;
      }
      // test arc point 1
      testvalA = cont[testpA].Ai[0]*arcpt1[0] + cont[testpA].Ai[1]*arcpt1[1]
	+ cont[testpA].Ai[2]*arcpt1[2] + cont[testpA].Ai[3];
      testvalB = cont[testpB].Ai[0]*arcpt1[0] + cont[testpB].Ai[1]*arcpt1[1]
	+ cont[testpB].Ai[2]*arcpt1[2] + cont[testpB].Ai[3];
      if((testvalA < 0.0) && (testvalB < 0.0)) { // point is good
	add_vertex(poly, vn, arcpt1, vedge[edgei].plane[0], vedge[edgei].plane[1], -1);
	poly[vn].dist = rado;
	++vn;
      }
    }
  }
  
  // reduce poly vertex list
  vj=0;
  for(vi=0; vi<vn; ++vi) {
    if(poly[vi].dist <= rado) {
      poly[vj] = poly[vi];
      ++vj;
    }
  }
  vn = vj;
  
  // ----- calculate center points and mark engulfing contacts -----
  for(cai=0; cai<NC; ++cai) {
    centerpt[cai].xi[0] = -cont[cai].Ai[0]*cont[cai].Ai[3];
    centerpt[cai].xi[1] = -cont[cai].Ai[1]*cont[cai].Ai[3];
    centerpt[cai].xi[2] = -cont[cai].Ai[2]*cont[cai].Ai[3];
    centerpt[cai].dist = fabs(cont[cai].Ai[3]);
    if(cont[cai].Ai[3] > 0.0) {
      cont[cai].flag = 'E'; // engulfing contact
    }
  }
  
  if(recalc == 'Y') {
    // reset atom coordinates to original values
    atom->position[0] = origcoor[0];
    atom->position[1] = origcoor[1];
    atom->position[2] = origcoor[2];
  }
  return(vn);
}



/****************************
 * subroutine get_firstvert
 ****************************/

void get_firstvert(struct plane cont[], int *planeA, int *planeB, int *planeC, int NC) {

  int cai;
  double mindist;
  double ptA[3];
  double vectA[4];     // 4 values so it can be used as a plane equation as well
  double vt;           // vector parameter 't' for line x=x'+lt, y=y'+mt, z=z'+nt
  double vtmin;        // minimum value for vt
  double vtdiv;        // check for division by zero in vt calculation
  double temppt[3];
  double ptdist;

  // -------------  find initial edge, on plane closest to origin  -------------
  mindist = 9.9e+9;  // dummy value
  for(cai=0; cai<NC+4; ++cai) { // added +4 BJM dec2014 
    if(cont[cai].dist < mindist) {
      mindist = cont[cai].dist;
      *planeA = cai;
    }
  }

  mindist = 9.9e+9;  // dummy value
  for(cai=0; cai<NC+4; ++cai) {
    if(cai != *planeA) {
      // define third plane through origin, perpendicular to planeA and candidate plane
      vectA[0] = cont[*planeA].Ai[1]*cont[cai].Ai[2] - cont[*planeA].Ai[2]*cont[cai].Ai[1];
      vectA[1] = cont[*planeA].Ai[2]*cont[cai].Ai[0] - cont[*planeA].Ai[0]*cont[cai].Ai[2];
      vectA[2] = cont[*planeA].Ai[0]*cont[cai].Ai[1] - cont[*planeA].Ai[1]*cont[cai].Ai[0];
      vectA[3] = 0.0;
      // get intersection of three planes, check if it is closest to origin
      if(solve_3x3(cont[*planeA].Ai, cont[cai].Ai, vectA, temppt) != -1) {
	ptdist = sqrt(temppt[0]*temppt[0] + temppt[1]*temppt[1] + temppt[2]*temppt[2]);
	if(ptdist < mindist) {
	  *planeB = cai;
	  mindist = ptdist;
	  ptA[0] = temppt[0];
	  ptA[1] = temppt[1];
	  ptA[2] = temppt[2];
	}
      }
    }
  }
  // recalc vector normal to planes A and B
  vectA[0] = cont[*planeA].Ai[1]*cont[*planeB].Ai[2] - cont[*planeA].Ai[2]*cont[*planeB].Ai[1];
  vectA[1] = cont[*planeA].Ai[2]*cont[*planeB].Ai[0] - cont[*planeA].Ai[0]*cont[*planeB].Ai[2];
  vectA[2] = cont[*planeA].Ai[0]*cont[*planeB].Ai[1] - cont[*planeA].Ai[1]*cont[*planeB].Ai[0];
  vectA[3] = 0.0;
  
  // get starting vertex on polyhedron
  vtmin = 9.9e+9;   // dummy value
  for(cai=0; cai<NC+4; ++cai) {
    if((cai != *planeA) && (cai != *planeB)) {
      vtdiv = (cont[cai].Ai[0]*vectA[0] +cont[cai].Ai[1]*vectA[1] +cont[cai].Ai[2]*vectA[2]);
      if(vtdiv != 0.0) {
	vt = -(cont[cai].Ai[0]*ptA[0] +cont[cai].Ai[1]*ptA[1] +cont[cai].Ai[2]*ptA[2] +cont[cai].Ai[3])/vtdiv;
	if(fabs(vt) < vtmin) {
	  vtmin = fabs(vt);
	  *planeC = cai;
	}
      }
    }
  }
  // three planes defining first vertex identified. return and call solve_3x3
  return;
}



/****************************
 * function solve_3x3
 ****************************/

// determines the intersection of three planes
// (solves a system of three linear equations and 3 unknowns)
// input is 3 four element arrays, representing Ax+By+Cz+D=0
// output is a three element array, (xi,xj,xk).

int solve_3x3(double eq0[], double eq1[], double eq2[], double pt[])
{
	double cof00, cof01, cof02;   // matrix cofactors
	double cof10, cof11, cof12;
	double cof20, cof21, cof22;
	double det;                   // determinant of matrix

	cof00 =  eq1[1]*eq2[2] - eq2[1]*eq1[2];
	cof01 = -eq1[0]*eq2[2] + eq2[0]*eq1[2];
	cof02 =  eq1[0]*eq2[1] - eq2[0]*eq1[1];

	cof10 = -eq0[1]*eq2[2] + eq2[1]*eq0[2];
	cof11 =  eq0[0]*eq2[2] - eq2[0]*eq0[2];
	cof12 = -eq0[0]*eq2[1] + eq2[0]*eq0[1];

	cof20 =  eq0[1]*eq1[2] - eq1[1]*eq0[2];
	cof21 = -eq0[0]*eq1[2] + eq1[0]*eq0[2];
	cof22 =  eq0[0]*eq1[1] - eq1[0]*eq0[1];

	det = eq0[0]*cof00 + eq0[1]*cof01 + eq0[2]*cof02;
	if(det == 0.0) {
		//printf("no solution for equation set, determinant is zero\n");
		return(-1);
	} else {
		pt[0] = -(eq0[3]*cof00 + eq1[3]*cof10 + eq2[3]*cof20)/det;
		pt[1] = -(eq0[3]*cof01 + eq1[3]*cof11 + eq2[3]*cof21)/det;
		pt[2] = -(eq0[3]*cof02 + eq1[3]*cof12 + eq2[3]*cof22)/det;
		return(0);
	}
}



/**********************************
 * Function solve_2xS
 * revised 19/02/2001  BJM
 ***********************************/

// determines the intersection of two planes and a sphere radius 'rad'
// input is 2 four element arrays, representing Ax+By+Cz+D=0
// plus the radius of a sphere centered on the origin.
// output is two three element arrays pt0 and pt1, (xi,xj,xk).
// return value is -1 if no real solution exists.

int solve_2xS(struct plane eq0, struct plane eq1, float rado, double pt0[], double pt1[])
{
	double eq2[3];              // eqn of plane through (0,0,0) and perp. to other two
	double cof00, cof01, cof02; // matrix cofactors
	double cof10, cof11, cof12; // (don't need cof20, cof21, cof22)
	double det;                 // determinant of matrix
	double avgpt[3];            // average of two solution points
	double t;                   // parameter in eqn of line: x=xo+At, y=yo+Bt, z=zo+Ct.
	int xi;                     // coordinate counter

	eq2[0] = eq0.Ai[1]*eq1.Ai[2] - eq0.Ai[2]*eq1.Ai[1];
	eq2[1] = eq0.Ai[2]*eq1.Ai[0] - eq0.Ai[0]*eq1.Ai[2];
	eq2[2] = eq0.Ai[0]*eq1.Ai[1] - eq0.Ai[1]*eq1.Ai[0];

	cof00 =  eq1.Ai[1]*eq2[2] - eq2[1]*eq1.Ai[2];
	cof01 = -eq1.Ai[0]*eq2[2] + eq2[0]*eq1.Ai[2];
	cof02 =  eq1.Ai[0]*eq2[1] - eq2[0]*eq1.Ai[1];

	cof10 = -eq0.Ai[1]*eq2[2] + eq2[1]*eq0.Ai[2];
	cof11 =  eq0.Ai[0]*eq2[2] - eq2[0]*eq0.Ai[2];
	cof12 = -eq0.Ai[0]*eq2[1] + eq2[0]*eq0.Ai[1];

	det = eq2[0]*eq2[0] + eq2[1]*eq2[1] + eq2[2]*eq2[2];
	if(det == 0) {
		//printf("no solution in solve_2xS\n");
		return(-1);
	}

	avgpt[0] = -(eq0.Ai[3]*cof00 + eq1.Ai[3]*cof10)/det;
	avgpt[1] = -(eq0.Ai[3]*cof01 + eq1.Ai[3]*cof11)/det;
	avgpt[2] = -(eq0.Ai[3]*cof02 + eq1.Ai[3]*cof12)/det;

	t = (rado*rado-avgpt[0]*avgpt[0]-avgpt[1]*avgpt[1]-avgpt[2]*avgpt[2])/det;
	if(t<0.0) {
		return(-1);
	} else {
		t = sqrt(t);
	}
	for(xi=0; xi<3; ++xi) {
		pt0[xi] = avgpt[xi] + t*eq2[xi];
		pt1[xi] = avgpt[xi] - t*eq2[xi];
	}
	return(0);
}



/*******************************
 * subroutine add_vertex
 *******************************/

// add polyhedron vertex to local list for current atom contacts

int add_vertex(struct vertex poly[], int vn, double coor[], int planeA, int planeB, int planeC)
{
	poly[vn].xi[0] = coor[0];
	poly[vn].xi[1] = coor[1];
	poly[vn].xi[2] = coor[2];
	poly[vn].plane[0] = planeA;
	poly[vn].plane[1] = planeB;
	poly[vn].plane[2] = planeC;
	poly[vn].dist = sqrt(coor[0]*coor[0] +coor[1]*coor[1] +coor[2]*coor[2]);
	return(0);
}



/*******************************
 * subroutine add_vedge
 *******************************/

// adds a new edge to the edge list. Direction of vector is tested so that
// it points away from the startpt, along the body of the polyhedron.

// stores information in variable vedge[edgenum].

void add_vedge(struct edgevector vedge[], int edgenum, struct plane cont[], int plane0, int plane1,
               int testplane, struct vertex poly[], int startpt)
{
	double testpt[3];
	double testval;

	vedge[edgenum].V[0] = cont[plane0].Ai[1]*cont[plane1].Ai[2] - cont[plane0].Ai[2]*cont[plane1].Ai[1];
	vedge[edgenum].V[1] = cont[plane0].Ai[2]*cont[plane1].Ai[0] - cont[plane0].Ai[0]*cont[plane1].Ai[2];
	vedge[edgenum].V[2] = cont[plane0].Ai[0]*cont[plane1].Ai[1] - cont[plane0].Ai[1]*cont[plane1].Ai[0];
	vedge[edgenum].startpt = startpt;
	vedge[edgenum].endpt = -1;  // flag, edge not completed.
	vedge[edgenum].plane[0] = plane0;
	vedge[edgenum].plane[1] = plane1;
	vedge[edgenum].startplane = testplane;
	vedge[edgenum].arc = '.'; // dummy value.

	// test direction of vector
	testpt[0] = poly[startpt].xi[0] + vedge[edgenum].V[0];
	testpt[1] = poly[startpt].xi[1] + vedge[edgenum].V[1];
	testpt[2] = poly[startpt].xi[2] + vedge[edgenum].V[2];

	testval = cont[testplane].Ai[0]*testpt[0] +cont[testplane].Ai[1]*testpt[1]
	          + cont[testplane].Ai[2]*testpt[2] +cont[testplane].Ai[3];
	if(testval > 0.0) { // vector is in wrong direction
		vedge[edgenum].V[0] = -vedge[edgenum].V[0];
		vedge[edgenum].V[1] = -vedge[edgenum].V[1];
		vedge[edgenum].V[2] = -vedge[edgenum].V[2];
	}
	return;
}



/************************
 * function test_point
 ************************/

// this function tests a given point versus a set of plane constraints
// it returns a value of 'Y' if the point is on or within the polyhedron, 
// and 'N' if the point is outside the polyhedron (violation of the plane inequality).
// ptX =(xo, yo, zo); if Axo+Byo+Czo+D > 0, point is behind plane (hidden)
// planes A,B,C are the planes that the point lies on, don't test.

char test_point(double ptX[], struct plane cont[], int NC, float rado, int planeA, int planeB, int planeC)
{
	int cp;      // counter, number of planes

	// if point is not behind any plane, keep point.
	for(cp=0; cp<NC; ++cp) {
		//if pt is on current plane, get next plane
		if((cp != planeA) && (cp != planeB) && (cp != planeC) && (cont[cp].flag != 'X')) {
			if((cont[cp].Ai[0]*ptX[0] + cont[cp].Ai[1]*ptX[1] + cont[cp].Ai[2]*ptX[2] + cont[cp].Ai[3]) > 0.0) {
				//point is behind plane, not on polyhedron.
				return('N');
			}
		}
	}
	return('Y');
}



/***************************
 * function cosPQR         *
 ***************************/

// this function returns the cosine of the angle between three points P,Q,R
// with Q being the center point.

double cosPQR( double ptP[], double ptQ[], double ptR[])
{
	double QP[3];    // vector from Q to P
	double QR[3];    // vector from Q to R
	double cosine;   // cosine of angle PQR at Q.

	// calculate vectors
	QP[0] = ptP[0] - ptQ[0];
	QP[1] = ptP[1] - ptQ[1];
	QP[2] = ptP[2] - ptQ[2];
	QR[0] = ptR[0] - ptQ[0];
	QR[1] = ptR[1] - ptQ[1];
	QR[2] = ptR[2] - ptQ[2];

	//calculate cosine
	cosine = (QP[0]*QR[0]+ QP[1]*QR[1]+ QP[2]*QR[2])
	         /sqrt((QP[0]*QP[0]+ QP[1]*QP[1]+ QP[2]*QP[2]) * (QR[0]*QR[0]+ QR[1]*QR[1]+ QR[2]*QR[2]));

	return(cosine);
}



/*****************************
 * subroutine project_points
 *****************************/

// this subroutine corrects for engulfing atoms. The center of the engulfing plane
// contact is used as the center of projection instead of the center of the atom.
// points already on the surface are not moved, preserving the SAS.

void project_points(struct vertex poly[], struct vertex centerpt[], float rado, int NC,
                    int NV, struct plane cont[])
{
	int epi;               // engulfing plane counter
	int cai;               // contact atom counter
	int engplane[4];       // index to engulfing planes
	double projpt[3];      // center point for projecting internal pts onto surface of sphere
	double pt0[3], pt1[3]; // temporary points for intersection solutions
	double V[3];           // vector from projection point to vertex
	double *P;             // pointer for vertex coordinates to be projected
	double a,b,c,k;        // for solving quadratic eqn for point projection
	int vi;                // vertex counter

	// count and mark engulfing planes
	epi = 0;
	for(cai=0; cai<NC; ++cai) {
		if(cont[cai].flag == 'E') {
			engplane[epi] = cai;
			++epi;
		}
	}

	// get projpt[] for projecting points to surface
	if(epi == 1) {
		projpt[0] = centerpt[engplane[0]].xi[0];
		projpt[1] = centerpt[engplane[0]].xi[1];
		projpt[2] = centerpt[engplane[0]].xi[2];
	} else if(epi == 2) {
		solve_2xS(cont[engplane[0]],cont[engplane[1]],rado, pt0, pt1);
		projpt[0] = (pt0[0]+pt1[0])/2;
		projpt[1] = (pt0[1]+pt1[1])/2;
		projpt[2] = (pt0[2]+pt1[2])/2;
	} else {
		solve_3x3(cont[engplane[0]].Ai, cont[engplane[1]].Ai, cont[engplane[2]].Ai, pt0);
		projpt[0] = pt0[0];
		projpt[1] = pt0[1];
		projpt[2] = pt0[2];
	}

	for(vi=0; vi<NV; ++vi) {
		if(poly[vi].plane[2] != -1) {
			// project point to intersection of engulfing plane(s) and surface of sphere
			P = poly[vi].xi;
			V[0] = P[0] - projpt[0];
			V[1] = P[1] - projpt[1];
			V[2] = P[2] - projpt[2];
			a = V[0]*V[0] + V[1]*V[1] + V[2]*V[2];
			b = 2*(P[0]*V[0] +P[1]*V[1] +P[2]*V[2]);
			c = P[0]*P[0] + P[1]*P[1] + P[2]*P[2] - rado*rado;  // c is < 0
			k = (sqrt(b*b - 4.0*a*c) - b)/(2*a);                // k is > 0
			P[0] += k*V[0];
			P[1] += k*V[1];
			P[2] += k*V[2];
			poly[vi].dist = rado;
		}
	}
	return;
}



/********************************
 * function spherical_arc       *
 ********************************/

// edited 16-Dec-2014 BJM
// Note: this calculates the minor arc, angle < 180 degrees.
// Major arc >180 degrees is area of (spherical cap - spherical arc) 

// given two points in cartesian space, the center of the spherical
// cap between atom A and the origin, and the radius of the sphere,
// this function returns the area of an arc between points B and C.
// the sides of the arc are the great circle distance (shortest distance)
// and the arc of the spherical cap centered on line AO.

double spherical_arc(struct vertex ptAo, struct vertex ptB, struct vertex ptC, float rado)
{
	// here, cosAOC=cosAOB, AOB = AOC.
	// BAC is the angle at the point on the origin, on line AO

	double cosAOB, cosBOC; // cosines of angles centered on (0,0,0)
	double cosBAC, angBAC;    // angle and cosine at vertex of cap
	double U, V, W;
	double tansqrS, tansqrSA, tansqrSB;
	double tanprod;        // product of tangents for spherical triangle
	double area;           // the value to be returned

	cosAOB = (ptAo.xi[0]*ptB.xi[0] + ptAo.xi[1]*ptB.xi[1] + ptAo.xi[2]*ptB.xi[2])/(ptAo.dist*ptB.dist);
	cosBOC = (ptB.xi[0]*ptC.xi[0] + ptB.xi[1]*ptC.xi[1] + ptB.xi[2]*ptC.xi[2])/(ptB.dist*ptC.dist);

	U = (1+cosAOB)*sqrt((1+cosBOC)/8.0);
	V = (1-cosAOB)*sqrt((1+cosBOC)/8.0);
	W = sqrt((1-cosAOB)*(1+cosAOB)*(1-cosBOC)/8.0);  // W == X
	tansqrS  = (1-U+V+W+W)/(1+U-V-W-W);
	tansqrSB = (1-U-V)/(1+U+V);
	tansqrSA = (1-U+V-W-W)/(1+U-V+W+W);
	if((tansqrS*tansqrSA) > 0.0) {
		tanprod = sqrt(tansqrSB*tansqrSB*tansqrS*tansqrSA);
	} else {
		tanprod = 0.0;
	}

	cosBAC = cosPQR(ptB.xi, ptAo.xi, ptC.xi);
	if(cosBAC>1.0) {
		angBAC = 0.0;
	} else if(cosBAC<-1.0) {
		angBAC = PI;
	} else {
		angBAC = acos(cosBAC);
	}

	area = rado*(angBAC*(rado-ptAo.dist) - 4.0*rado*atan(sqrt(tanprod)));
	return(area);
}
