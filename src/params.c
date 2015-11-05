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



static PLI_PARAM *pli_params = NULL;
static n_pli_params = 0;
static n_alloc_pli_params = 0;



static void params_read(char*,int);
static void params_set_defaults(void);
static void params_process_file(char*);
static void params_process_command_line(int,char**);
static void params_set_parameter(PLI_PARAM*,char*);



void run_help(SETTINGS *settings) {

  int i;
  PLI_PARAM *param;

  printf("PLI COMMAND-LINE HELP\n\n");

  for (i=0,param=pli_params;i<n_pli_params;i++,param++) {

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

  params_read(params_file_standard,1);
  params_read(params_file_custom,0);

  params_set_defaults();

  params_process_command_line(nargs,args);
}



static void params_read(char *filename,int check) {

  int n_read;
  char line[MAX_LINE_LEN],type_line[MAX_LINE_LEN],purpose_line[MAX_LINE_LEN],valid_values_line[MAX_LINE_LEN];
  char default_value_line[MAX_LINE_LEN],exposed_line[MAX_LINE_LEN],name[40];
  PLI_FILE *file;
  PLI_PARAM *param;

  file = open_file(filename,"r");

  if (file == NULL) {

    if (check) {

      error_fn("couldn't open settings file '%s'",filename);

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

      // (re)alloc memory for pli_params:

      if (n_alloc_pli_params == 0) {

	n_alloc_pli_params += 100;
	
	pli_params = (PLI_PARAM*) calloc(n_alloc_pli_params,sizeof(PLI_PARAM));
	
      } else if (n_pli_params == n_alloc_pli_params) {
	
	n_alloc_pli_params += 100;
	
	pli_params = (PLI_PARAM*) realloc(pli_params,n_alloc_pli_params*sizeof(PLI_PARAM));
      }
      
      if (pli_params == NULL) {
	
	error_fn("params_read: out of memory (re)allocating pli_params");
      }

      param = pli_params + n_pli_params;

      // set param values:

      n_read = sscanf(line+1,"%s",param->name);      
      n_read += sscanf(type_line,"%*[^:]: %s",param->type);
      n_read += sscanf(purpose_line,"%*[^:]: %[^\n]",param->help_text);
      n_read += sscanf(valid_values_line,"%*[^:]: %[^\n]",param->valid_values);
      n_read += sscanf(default_value_line,"%*[^:]: %[^\n]",param->default_value);
      n_read += sscanf(exposed_line,"%*[^:]: %d",&(param->exposed));

      if (n_read == 6) {

	n_pli_params++;

      } else {

	error_fn("params_read: corrupt parameter definition for '%s'\n",param->name);
      }
    }
  }

  close_file(file);
}



static void params_set_defaults(void) {

  int i;
  PLI_PARAM *param;

  for (i=0,param=pli_params;i<n_pli_params;i++,param++) {

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



static void params_process_command_line(int nargs,char **args) {

  int i;
  char **arg,*name;
  PLI_PARAM *param;

  for (i=0,arg=args;i<nargs;i++,arg++) {

    if ((*arg)[0] == '-') {

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
    }
  }
}



PLI_PARAM* params_get_parameter(char *name) {

  int i;
  PLI_PARAM *param;

  for (i=0,param=pli_params;i<n_pli_params;i++,param++) {

    if (!strcmp(param->name,name)) {

      return(param);
    }
  }

  return(NULL);
}



static void params_set_parameter(PLI_PARAM *param,char *value) {

  int flag;

  if (!strcmp(param->type,"integer")) {

    flag = sscanf(value,"%d",&(param->value.i));

  } else if (!strcmp(param->type,"double")) {

    flag = sscanf(value,"%lf",&(param->value.d));

  } else if (!strcmp(param->type,"string")) {

    strcpy(param->value.s,value);

  } else {

    error_fn("params_set_parameter: unknown parameter type '%s'",param->type);
  }
}

