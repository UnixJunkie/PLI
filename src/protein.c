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



#define BACKBONE_C 1
#define BACKBONE_O 2
#define BACKBONE_N 3
#define BACKBONE_CA 4



static RESIDUE AMINO_ACIDS[] = {
  { "ALA",0,' ',' ','A',1 },
  { "ARG",0,' ',' ','R',1 },
  { "ASN",0,' ',' ','N',1 },
  { "ASP",0,' ',' ','D',1 },
  { "CYS",0,' ',' ','C',1 },
  { "GLN",0,' ',' ','Q',1 },
  { "GLU",0,' ',' ','E',1 },
  { "GLY",0,' ',' ','G',1 },
  { "HIS",0,' ',' ','H',1 },
  { "ILE",0,' ',' ','I',1 },
  { "LEU",0,' ',' ','L',1 },
  { "LYS",0,' ',' ','K',1 },
  { "MET",0,' ',' ','M',1 },
  { "PHE",0,' ',' ','F',1 },
  { "PRO",0,' ',' ','P',1 },
  { "SER",0,' ',' ','S',1 },
  { "THR",0,' ',' ','T',1 },
  { "TRP",0,' ',' ','W',1 },
  { "TYR",0,' ',' ','Y',1 },
  { "VAL",0,' ',' ','V',1 },
  { "TPO",0,' ',' ','T',1 },
  { "SEP",0,' ',' ','S',1 },
  { "PTR",0,' ',' ','Y',1 },
  { "LAST",0,' ',' ',' ',0 }
};



static int backbone_atom_type(ATOM*);




RESIDUE* get_amino_acid(char *name) {

  RESIDUE *residue;

  residue = AMINO_ACIDS;

  while (strcmp(residue->name,"LAST")) {

    if ((residue->flags & AMINO_ACID) && (!strcmp(name,residue->name))) {

      return(residue);
    }

    residue++;
  }

  return(NULL);  
}



int backbone_bond_type(BOND *bond) {

  int bb_type1,bb_type2;

  bb_type1 = backbone_atom_type(bond->atom1);

  if (!bb_type1) {

    return(0);
  }

  bb_type2 = backbone_atom_type(bond->atom2);

  if (!bb_type2) {

    return(0);
  }

  if (((bb_type1 == BACKBONE_C) && (bb_type2 == BACKBONE_O)) ||
      ((bb_type2 == BACKBONE_C) && (bb_type1 == BACKBONE_O))) {

    return(2);
  }

  return(1);
}



static int backbone_atom_type(ATOM *atom) {

  if (atom->flags & AMINO_ACID_ATOM) {

    if (!strcmp(atom->name," C  ")) {

      return(BACKBONE_C);

    } else if (!strcmp(atom->name," O  ")) {

      return(BACKBONE_O);

    } else if (!strcmp(atom->name," N  ")) {

      return(BACKBONE_N);

    } else if (!strcmp(atom->name," CA ")) {

      return(BACKBONE_CA);
    }
  }

  return(0);
}

