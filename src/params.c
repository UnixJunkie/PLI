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


static LIST *pli_params = NULL;



static void params_read(char*,int);
static void params_set_defaults(void);
static void params_process_file(char*);
static void params_process_command_line(int,char**,int);
static void params_inherit(void);



void run_help(SETTINGS *settings) {

  int i;
  PLI_PARAM *param;

  printf("PLI COMMAND-LINE HELP\n\n");

  for (i=0,param=(PLI_PARAM*) pli_params->items;i<pli_params->n_items;i++,param++) {

    if (param->exposed) {

      printf("-%s\n",param->name);
      printf("  purpose     : %s\n",param->help_text);
      printf("  valid values: %s (def=%s)\n",param->valid_values,param->default_value);
      printf("\n");
    }
  }
}



void params_init(int nargs,char **args) {

  char *pli_dir,params_file_standard[MAX_LINE_LEN],params_file_custom[MAX_LINE_LEN];

  pli_dir = get_pli_dir();

  sprintf(params_file_standard,"%s/params/parameters.pli",pli_dir);
  sprintf(params_file_custom,"%s/params/parameters_custom.pli",pli_dir);

  pli_params = new_list("pli parameters",sizeof(PLI_PARAM),0);

  params_read(params_file_standard,1);
  params_read(params_file_custom,0);

  condense_list(pli_params);

  params_set_defaults();

  params_process_command_line(nargs,args,0);

  params_inherit();
}



static void params_read(char *filename,int check) {

  int n_read,pos,i;
  char line[MAX_LINE_LEN],type_line[MAX_LINE_LEN],purpose_line[MAX_LINE_LEN],valid_values_line[MAX_LINE_LEN];
  char default_value_line[MAX_LINE_LEN],exposed_line[MAX_LINE_LEN],name[40],synonyms[MAX_LINE_LEN],*word,*synonym;
  PLI_FILE *file;
  PLI_PARAM *param;

  file = open_file(filename,"r");

  if (file == NULL) {

    if (check) {

      error_fn("%s: couldn't open settings file '%s'",__func__,filename);

    } else {

      return;
    }
  }

  while (!end_of_file(file)) {

    if (read_line(line,MAX_LINE_LEN,file) == NULL)
      break;

    if (line[0] == '-') {

      // read parameter lines:

      if (read_line(type_line,MAX_LINE_LEN,file) == NULL)
	break;

      if (read_line(purpose_line,MAX_LINE_LEN,file) == NULL)
	break;
      
      if (read_line(valid_values_line,MAX_LINE_LEN,file) == NULL)
	break;
      
      if (read_line(default_value_line,MAX_LINE_LEN,file) == NULL)
	break;
      
      if (read_line(exposed_line,MAX_LINE_LEN,file) == NULL)
	break;

      param = (PLI_PARAM*) add_list_item(pli_params);

      // set param values:

      n_read = sscanf(line+1,"%s",param->name);      
      n_read += sscanf(type_line,"%*[^:]: %s",param->type);
      n_read += sscanf(purpose_line,"%*[^:]: %[^\n]",param->help_text);
      n_read += sscanf(valid_values_line,"%*[^:]: %[^\n]",param->valid_values);
      n_read += sscanf(default_value_line,"%*[^:]: %[^\n]",param->default_value);
      n_read += sscanf(exposed_line,"%*[^:]: %d",&(param->exposed));

      if (n_read != 6) {

	error_fn("%s: corrupt parameter definition for '%s'\n",__func__,param->name);
      }

      if (sscanf(line+1,"%s [%[^]]]",param->name,synonyms) == 2) {

	pos = 0;

	param->synonyms = new_list("synonyms",50*sizeof(char),0);

	while (word = nextword(synonyms,',',&pos)) {

	  synonym = (char*) add_list_item(param->synonyms);

	  strcpy(synonym,word);
	}

	condense_list(param->synonyms);

      } else {

	param->synonyms = NULL;
      }
    }
  }

  close_file(file);
}



static void params_set_defaults(void) {

  int i;
  PLI_PARAM *param;

  for (i=0,param=(PLI_PARAM*) pli_params->items;i<pli_params->n_items;i++,param++) {

    params_set_parameter(param,param->default_value);
  }
}



static void params_process_file(char *filename) {

  int n_read;  char line[MAX_LINE_LEN],name[MAX_LINE_LEN],value[MAX_LINE_LEN];
  PLI_FILE *file;
  PLI_PARAM *param;

  file = open_file(filename,"r");

  if (file == NULL) {

    error_fn("couldn't open settings file '%s'",filename);
  }

  while (!end_of_file(file)) {

    if (read_line(line,MAX_LINE_LEN,file) == NULL)
      break;

    if (line[0] != '#') {

      n_read = sscanf(line,"%s %[^\n]",name,value);

      if (n_read == 2) {
 
	param = params_get_parameter(name);

	if (param) {

	  params_set_parameter(param,value);

	} else {

	  warning_fn("params_process_file: skipping unknown parameter '%s' in file '%s'",name,filename); 
	}     
      }
    }
  }

  close_file(file);
}



static void params_process_command_line(int nargs,char **args,int i) {

  char **arg,*name;
  PLI_PARAM *param;

  if (i >= nargs) {

    return;
  }

  arg = args + i;

  if (((*arg)[0] == '-') && (!isdigit((*arg)[1]))) {

    name = (*arg) + 1;

    if (!strcmp(name,"settings")) {

      params_process_file(*(arg+1));

    } else if ((!strcmp(name,"h")) || (!strcmp(name,"help"))) {
      
      param = params_get_parameter("mode");
      
      params_set_parameter(param,"help");
      
    } else {
      
      param = params_get_parameter(name);
      
      if (param) {
	
	params_set_parameter(param,*(arg+1));
	
      } else {
	
	warning_fn("params_process_command_line: skipping unknown parameter '%s'",name); 
      }
    }

    params_process_command_line(nargs,args,i+2);

  } else {

    params_process_command_line(nargs,args,i+1);
  }
}



static void params_inherit(void) {

  int i;
  char iname[MAX_LINE_LEN];
  PLI_PARAM *param,*iparam;

  for (i=0,param=(PLI_PARAM*) pli_params->items;i<pli_params->n_items;i++,param++) {

    if ((!strncmp(param->default_value,"inherit",7)) && (!strncmp(param->value.s,"inherit",7))) {

      if (sscanf(param->default_value,"inherit %s",iname)) {

	iparam = params_get_parameter(iname);

	if ((!strcmp(param->type,"integer")) && (!strcmp(iparam->type,"integer"))){

	  param->value.i = iparam->value.i;

	} else if ((!strcmp(param->type,"double")) && (!strcmp(iparam->type,"double"))){
	  
	  param->value.d = iparam->value.d;
	  
	} else if ((!strcmp(param->type,"string")) && (!strcmp(iparam->type,"string"))) {
	  
	  strcpy(param->value.s,iparam->value.s);
	  
	} else {
	  
	  error_fn("%s: corrupt parameter definition for '%s' (1)",__func__,param->name);
	}

      } else {

	error_fn("%s: corrupt parameter definition for '%s' (2)\n",__func__,param->name);
      }
    }
  }
}



PLI_PARAM* params_get_parameter(char *name) {

  int i,j,n_matched_params;
  char *synonym;
  PLI_PARAM *param, *matched_param;

  n_matched_params = 0;

  for (i=0,param=(PLI_PARAM*) pli_params->items;i<pli_params->n_items;i++,param++) {

    if (!strcmp(param->name,name)) {

      return(param);
    }

    if (!strncmp(param->name, name, strlen(name))) {

      matched_param = param;

      n_matched_params++;
    }

    if (param->synonyms) {

      for (j=0;j<param->synonyms->n_items;j++) {

	synonym = get_list_item(param->synonyms,j);

	if (!strcmp(synonym,name)) {

	  return(param);
	}
      }
    }
  }

  if (n_matched_params == 1) {

    return(matched_param);
  }

  error_fn("params_get_parameter: no such parameter '%s'",name);
}




void params_set_parameter(PLI_PARAM *param,char *value) {

  int flag;

  // always store text value:

  strcpy(param->value.s,value);

  if (!strcmp(param->type,"integer")) {

    flag = sscanf(value,"%d",&(param->value.i));

  } else if (!strcmp(param->type,"double")) {

    flag = sscanf(value,"%lf",&(param->value.d));

  } else if (strcmp(param->type,"string")) {

    error_fn("params_set_parameter: unknown parameter type '%s'",param->type);
  }
}

