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






char* atomlist2json(ATOMLIST *list) {

  int i;
  char **atom_text,*text;
  ATOM **atomp;

  set_atomio("json",NULL);

  atom_text = (char**) alloc_2d_cmatrix(list->natoms,MAX_LINE_LEN);

  for (i=0,atomp=list->atom;i<list->natoms;i++,atomp++) {

    atom2text((void*) *atomp,atom_text[i]);
  }

  text = list2json(atom_text,list->natoms);

  free_2d_cmatrix(atom_text);

  set_atomio(NULL,NULL);

  return(text);
}



static char* list2json(char **list,int n) {

  int i,len,pos;
  char *text;

  text = NULL;
  len = 2;

  for (i=0;i<n;i++) {

    len += strlen(list[i]) + 1;

    text = (text == NULL) ? (char*) malloc(len*sizeof(char)) : (char*) realloc(text,len*sizeof(char));
    
    if (text == NULL) {
      
      error_fn("%s: out of memory (re)allocating text",__func__);
    }
 
    (i == 0) ? sprintf(text,"[%s",list[i]) : sprintf(text+pos,",%s",list[i]);

    pos = len - 1;
  }

  text[pos] = ']';
  text[pos+1] = '\0';

  return(text);
}
