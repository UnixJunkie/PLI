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

static void init_int_array(INT_LIST *array);


INT_LIST* parse_int_csv(char *s, char *sep) {

  INT_LIST * array = malloc(sizeof(INT_LIST));
  init_int_array(array);
  char *s1,*s2,*token;
  int item,len;
  int success;
  int alloc_block_size = 1;

  len = strlen(s);

  s1 = (char*) malloc(len*sizeof(char));

  strcpy(s1,s);

  s2 = s1;

  token = strtok(s2, sep);

  while (token) {
    success = sscanf(token,"%d", &item);
    if (success) {
      if (array->n_items == array->n_allocated) {
        array->items = realloc(array->items, sizeof(int) * (array->n_allocated + alloc_block_size));
        array->n_allocated += alloc_block_size;
      }
      array->items[array->n_items] = item;
      array->n_items++;
    }
    token = strtok(NULL, ",");
  }

  free(s1);

  return array;
}



int word_in_text(char *word,char *text,char sep) {

  int pos = 0;
  char ctext[MAX_LINE_LEN],*cword;

  strcpy(ctext,text);

  while (cword = nextword(ctext,sep,&pos)) {

    if (!strcmp(cword,word)) {

      return(1);
    }
  }

  return(0);
}



char* nextword(char *text,char sep,int *pos) {

  int i;
  char *word;

  i = *pos;

  if (text[i] == '\0') {

    return(NULL);
  }

  word = text + i;

  while ((text[i] != sep) && (text[i] != '\0')) {

    i++;
  }

  if (text[i] == sep) {

    text[i] = '\0';

    *pos = i+1;

  } else {

    *pos = i;
  }

  return(word);
}



int read_string_section(char *string,char *section) {

  int i,n_open,n_close;
  char c,c_open,c_close;

  n_open = 1;
  n_close = 0;

  c_open = string[0];

  if (c_open == '(') {

    c_close = ')';

  } else if (c_open == '[') {

    c_close = ']';
    
  } else {

    error_fn("%s: unexpected error occurred reading string '%s'",__func__,string);
  }

  i = 0;

  do {

    i++;

    c = string[i];

    if (c == '\0') {

      error_fn("%s: unmatched bracket for string '%s'",__func__,string);
    }

    if (c == c_open) {

      n_open++;

    } else if (c == c_close) {

      n_close++;
    }

  } while (n_open != n_close);

  substring(string,1,i-1,section);

  return(i);
}


void init_int_array(INT_LIST *array) {
  array->items = NULL;
  array->n_items = 0;
  array->n_allocated = 0;
}


void free_int_array (INT_LIST *array) {
  free(array->items);
  free(array);
}



void title_case(char* s) {
  int i;

  for (i=0; s[i]; i++) {
    if ((i==0) || (s[i-1] == ' ')) {
      s[i] = toupper(s[i]);
    } else {
      s[i] = tolower(s[i]);
    }
  }
}

int is_in_iarray(int v, int *arr, int n) {
  int i;
  for (i=0; i<n; i++) {
    if (arr[i] == v) {
      return(1);
    }
  }
  return(0);
}



void set_double_array(double *x,int n,double value) {

  int i;

  for (i=0;i<n;i++) {

    x[i] = value;
  }
}

