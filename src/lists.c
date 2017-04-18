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



static LIST *stored_lists = NULL;




LIST* new_list(char *name,int item_size,int n) {

  LIST *list;

  list = (LIST*) malloc(sizeof(LIST));

  if (list == NULL) {

    error_fn("%s: out of memory allocating list",__func__);
  }

  init_list(list,name,item_size);

  if (n) {

    list->n_alloc_items = n;

    list->items = (void*) calloc(list->n_alloc_items,list->item_size);

    if (list->items == NULL) {

      error_fn("%s: out of memory allocating items for list %s",__func__,list->name);
    }

  } else {

    list->n_alloc_items = 0;
  }

  return(list);
}



void init_list(LIST *list,char *name,int item_size) {

  strcpy(list->name,name);
  
  list->item_size = item_size;

  list->n_items = 0;
  list->n_alloc_items = 0;
  list->items = NULL;
}



void reset_list(LIST *list) {

  if (list->items) {

    free(list->items);
  }

  list->n_items = 0;
  list->n_alloc_items = 0;
  list->items = NULL;
}



void* add_list_item(LIST *list) {

  void *item;

  if (list->items == NULL) {

    list->n_alloc_items = 1;

    list->items = (void*) calloc(list->n_alloc_items,list->item_size);

  } else {

    if (list->n_items == list->n_alloc_items) {

      list->n_alloc_items *= 2;

      list->items = (void*) realloc(list->items,(list->n_alloc_items)*(list->item_size));
    }
  }

  if (list->items == NULL) {

    error_fn("%s: out of memory adding item to list %s (n_alloc_items=%d,item_size=%d)",__func__,list->name,list->n_alloc_items,list->item_size);
  }

  item = (void*) (((char*) list->items) + ((list->n_items)*(list->item_size)));

  (list->n_items)++;

  return(item);
}



void remove_list_item(LIST *list,void *item) {

  int i,n_items,size;
  void *litem;

  i = 0;

  n_items = list->n_items;
  size = list->item_size;
  litem = list->items;

  while ((litem != item) && (i < n_items)) {

    litem = (void*) ((char*) litem + size);

    i++;
  }

  if (i < n_items-1) {

    memcpy(litem,(void*)(((char*)litem)+size),(n_items-i-1)*size);

  } else if (i == n_items) {

    error_fn("%s: item not in list",__func__);
  }

  list->n_items--;
}



void* get_list_item(LIST *list,int id) {

  void *item;

  item = (void*) (((char*) list->items) + ((id)*(list->item_size))); 

  return(item);
}



int item_in_list(LIST *list,void *item) {

  int i,size;
  char *litem;

  size = list->item_size;

  litem = (char*) list->items;

  for (i=0;i<list->n_items;i++) {

    if (!memcmp(item,(void*)litem,size)) {

      return(1);
    }

    litem += size;
  }

  return(0);
}



void condense_list(LIST *list) {

  if (list->n_alloc_items > list->n_items) {

    list->n_alloc_items = list->n_items;

    list->items = (void*) realloc(list->items,(list->n_alloc_items)*(list->item_size));

    if (list->items == NULL) {
      
      error_fn("%s: unexpected error occurred condensing list %s",__func__,list->name);
    }
  }
}



void free_list(LIST *list) {

  if (list) { 

    if (list->items) {
      
      free(list->items);
    }
    
    free(list);
  }
}



LIST* get_list(char *name,char *gname) {

  int i;
  LIST *group,**listp,*list;

  group = (gname) ? get_list(gname,NULL) : stored_lists;

  if (group) {

    for (i=0,listp=(LIST**) group->items;i<group->n_items;i++,listp++) {
      
      list = *listp;
      
      if (!strcmp(list->name,name)) {
	
	return(list);
      }
    }
  }

  return(NULL);
}



void store_list(LIST *list,char *gname) {

  LIST *group,**groupp,**listp;

  if (stored_lists == NULL) {

    stored_lists = new_list("stored lists",sizeof(LIST*),0);
  }

  group = get_list(gname,NULL);

  if (!group) {

    group = new_list(gname,sizeof(LIST*),0);

    groupp = (LIST**) add_list_item(stored_lists);

    *groupp = group;
  }

  listp = (LIST**) add_list_item(group);

  *listp = list;
}



LIST* atomlist2list(ATOMLIST *alist,char *name) {

  int i;
  ATOM **atom,**item;
  LIST *list;

  if (!alist) {

    return(NULL);
  }

  list = new_list(name,sizeof(ATOM*),alist->natoms);

  for (i=0,atom=alist->atom;i<alist->natoms;i++,atom++) {

    item = (ATOM**) add_list_item(list);

    *item = *atom;
  }

  return(list);
}
