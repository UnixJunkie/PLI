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

#define TFLAGS (((ATOM*) atom)->type ? ((ATOM*) atom)->type->flags : 00)

static void atom2pdb(void*,char*,char*);
static void atom2id(void*,char*,char*);
static void atom2name(void*,char*,char*);
static void atom2altloc(void*,char*,char*);
static void atom2subname(void*,char*,char*);
static void atom2chain(void*,char*,char*);
static void atom2subid(void*,char*,char*);
static void atom2icode(void*,char*,char*);
static void atom2set(void*,char*,char*);
static void atom2sym_mol(void*,char*,char*);
static void atom2x(void*,char*,char*);
static void atom2y(void*,char*,char*);
static void atom2z(void*,char*,char*);
static void atom2bfactor(void*,char*,char*);
static void atom2occupancy(void*,char*,char*);
static void atom2ami(void*,char*,char*);
static void atom2radius(void*,char*,char*);
static void atom2score(void*,char*,char*);
static void atom2carea(void*,char*,char*);
static void atom2iarea(void*,char*,char*);
static void atom2earea(void*,char*,char*);
static void atom2tid(void*,char*,char*);
static void atom2uid(void*,char*,char*);
static void atom2eid(void*,char*,char*);
static void atom2tname(void*,char*,char*);
static void atom2uname(void*,char*,char*);
static void atom2ename(void*,char*,char*);
static void atom2don(void*,char*,char*);
static void atom2acc(void*,char*,char*);
static void atom2met(void*,char*,char*);
static void atom2lip(void*,char*,char*);
static void atom2aro(void*,char*,char*);
static void atom2error_flags(void*,char*,char*);


static OFIELD atom_ofields[] = {
  { "pdb",         "",       "",         "",    NULL,   atom2pdb },
  { "id",          "%d",     "%5d",      "",    NULL,   atom2id },
  { "name",        "\"%s\"", "%4s",      "",    NULL,   atom2name },
  { "altloc",      "\"%c\"", "%c",       "",    NULL,   atom2altloc },
  { "subname",     "%s",     "%4s",      "",    NULL,   atom2subname },
  { "chain",       "\"%c\"", "%c",       "",    NULL,   atom2chain },
  { "subid",       "%d",     "%5d",      "",    NULL,   atom2subid },
  { "icode",       "\"%c\"", "%c",       "",    NULL,   atom2icode },
  { "set",         "%d",     "%5d",      "",    NULL,   atom2set },
  { "sym_mol",     "%d",     "%5d",      "",    NULL,   atom2sym_mol },
  { "bfactor",     "%.2lf",  "%6.2lf",   "",    NULL,   atom2bfactor },
  { "occupancy",   "%.2lf",  "%6.2lf",   "",    NULL,   atom2occupancy },
  { "ami",         "%.4lf",  "%10.4lf",  "",    NULL,   atom2ami },
  { "x",           "%.4lf",  "%10.4lf",  "",    NULL,   atom2x },
  { "y",           "%.4lf",  "%10.4lf",  "",    NULL,   atom2y },
  { "z",           "%.4lf",  "%10.4lf",  "",    NULL,   atom2z },
  { "radius",      "%.2lf",  "%10.2lf",  "",    NULL,   atom2radius },
  { "score",       "%.4lf",  "%10.4lf",  "",    NULL,   atom2score },
  { "carea",       "%.4lf",  "%10.4lf",  "",    NULL,   atom2carea },
  { "iarea",       "%.4lf",  "%10.4lf",  "",    NULL,   atom2iarea },
  { "earea",       "%.4lf",  "%10.4lf",  "",    NULL,   atom2earea },
  { "tid",         "%d",     "%5d",      "",    NULL,   atom2tid },
  { "uid",         "%d",     "%5d",      "",    NULL,   atom2uid },
  { "eid",         "%d",     "%5d",      "",    NULL,   atom2eid },
  { "tname",       "\"%s\"", "%-20s",    "",    NULL,   atom2tname },
  { "uname",       "\"%s\"", "%-20s",    "",    NULL,   atom2uname },
  { "ename",       "\"%s\"", "%-5s",     "",    NULL,   atom2ename },
  { "error_flags", "%d",     "%5d",      "",    NULL,   atom2error_flags },
  { "don",         "%d",     "%5d",      "",    NULL,   atom2don },
  { "acc",         "%d",     "%5d",      "",    NULL,   atom2acc },
  { "met",         "%d",     "%5d",      "",    NULL,   atom2met },
  { "lip",         "%d",     "%5d",      "",    NULL,   atom2lip },
  { "aro",         "%d",     "%5d",      "",    NULL,   atom2aro },
  { "last",        "",       "",         "",    NULL,   NULL }
};



static OGROUP atom_ogroups[] = {
  { "pliff",   atom_ofields,"id,name,altloc,subname,chain,subid,icode,set,sym_mol,xyz,radius,bfactor,occupancy,areas,eid,uid,tid,error_flags" },
  { "basic",   atom_ofields,"id,name,altloc,subname,chain,subid,icode,set,sym_mol" },
  { "types",   atom_ofields,"id,name,tid,uid,eid,tname,uname,ename" },
  { "stypes",  atom_ofields,"id,name,don,acc,met,lip,aro" },
  { "xyz",     atom_ofields,"x,y,z" },
  { "areas",   atom_ofields,"carea,iarea,earea" },
  { "last",    NULL,        "", }
};



void atom2text(ATOM *atom,char *text) {

  obj2text((void*) atom,text,atom_ofields);
}



void set_atomio(char *style,char *flags) {

  char *ostyle,*oflags;

  ostyle = (style) ? style : (params_get_parameter("atom_oformat"))->value.s;
  oflags = (flags) ? flags : (params_get_parameter("atom_oflags"))->value.s;

  init_ofields(atom_ofields);

  // recognise 'pdb' as a special output format for atoms:

  if (!strcmp(ostyle,"pdb")) {

    (get_ofield(atom_ofields,"pdb"))->ostyle= get_ostyle("formatted");

  } else {

    set_ogroups(atom_ogroups,ostyle,oflags);

    set_ofields(atom_ofields,ostyle,oflags);
  }
}



static void atom2pdb(void *vatom,char *text,char *format) {

  char subname[4],formal_charge[3];
  ATOM *atom;
  UNITED_ATOM *uatom;

  atom = (ATOM*) vatom;

  substring(atom->subname,0,3,subname);

  upper_case(subname);

  strcpy(formal_charge,"");

  if ((atom->type) && (atom->type->united_atom)) {

    uatom = atom->type->united_atom;

    if (uatom->formal_charge < 0) {

      sprintf(formal_charge,"%1d-",uatom->formal_charge);

    } else if (uatom->formal_charge > 0) {

      sprintf(formal_charge,"%1d+",uatom->formal_charge);
    }
  }

  sprintf(text,"%-6s%5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s","ATOM ",
	  atom->id,atom->name,atom->altloc,subname,atom->chain,atom->subid,atom->icode,
	  atom->position[0],atom->position[1],atom->position[2],
	  atom->occupancy,atom->bfactor,atom->element->name,formal_charge);
}



static void atom2id(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->id); }
static void atom2name(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->name); }
static void atom2altloc(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->altloc); }
static void atom2subname(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->subname); }
static void atom2chain(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->chain); }
static void atom2subid(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->subid); }
static void atom2icode(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->icode); }
static void atom2set(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->set); }
static void atom2bfactor(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->bfactor); }
static void atom2occupancy(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->occupancy); }
static void atom2ami(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->td); }
static void atom2x(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->position[0]); }
static void atom2y(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->position[1]); }
static void atom2z(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->position[2]); }
static void atom2radius(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->vdw_radius); }
static void atom2score(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->score); }
static void atom2carea(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->contact_area); }
static void atom2earea(void *atom,char *text,char *format) { sprintf(text,format,((ATOM*) atom)->exposed_area); }
static void atom2error_flags(void *atom,char *text,char *format) { sprintf(text,format,((int) ((ATOM*) atom)->error_flags)); }
static void atom2don(void *atom,char *text,char *format) { sprintf(text,format,((TFLAGS & HBOND_DONOR_ATOM_TYPE) ? 1 : 0)); }
static void atom2acc(void *atom,char *text,char *format) { sprintf(text,format,((TFLAGS & HBOND_ACCEPTOR_ATOM_TYPE) ? 1 : 0)); }
static void atom2met(void *atom,char *text,char *format) { sprintf(text,format,((TFLAGS & METAL_ATOM_TYPE) ? 1 : 0)); }
static void atom2lip(void *atom,char *text,char *format) { sprintf(text,format,((TFLAGS & LIPOPHILIC_ATOM_TYPE) ? 1 : 0)); }
static void atom2aro(void *atom,char *text,char *format) { sprintf(text,format,((TFLAGS & AROMATIC_ATOM_TYPE) ? 1 : 0)); }


static void atom2tid(void *vatom,char *text,char *format) {

  ATOM_TYPE *type;

  type = ((ATOM*) vatom)->type;

  sprintf(text,format,(type) ? type->id : -1);
}



static void atom2tname(void *vatom,char *text,char *format) {

  ATOM_TYPE *type;

  type = ((ATOM*) vatom)->type;

  sprintf(text,format,(type) ? type->name : "unknown");
}



static void atom2uid(void *atom,char *text,char *format) {

  ATOM_TYPE *type;
  UNITED_ATOM *uatom;

  type = ((ATOM*) atom)->type;

  uatom = (type) ? type->united_atom : NULL;

  sprintf(text,format,(uatom) ? uatom->id : -1);
}



static void atom2uname(void *atom,char *text,char *format) {

  ATOM_TYPE *type;
  UNITED_ATOM *uatom;

  type = ((ATOM*) atom)->type;

  uatom = (type) ? type->united_atom : NULL;

  sprintf(text,format,(uatom) ? uatom->name : "unknown");
}



static void atom2eid(void *atom,char *text,char *format) {

  ELEMENT *element;

  element = ((ATOM*) atom)->element;

  sprintf(text,format,(element) ? element->id : -1);
}



static void atom2ename(void *atom,char *text,char *format) {

  ELEMENT *element;

  element = ((ATOM*) atom)->element;

  sprintf(text,format,(element) ? element->name : "");
}



static void atom2sym_mol(void *atom,char *text,char *format) { 

  MOLECULE *molecule;

  molecule = ((ATOM*) atom)->molecule;

  sprintf(text,format,((molecule) && (molecule->flags & SYMMETRY_MOLECULE)) ? 1 : 0);
}



static void atom2iarea(void *vatom,char *text,char *format) { 

  ATOM *atom;

  atom = (ATOM*) vatom;

  sprintf(text,format,atom->intra_area + atom->covalent_area);
}
