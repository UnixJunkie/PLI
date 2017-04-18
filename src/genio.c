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
#include <zlib.h>
#include <ctype.h>
#include <sys/stat.h>



PLI_FILE *PLI_STDIN;
PLI_FILE *PLI_STDOUT;
PLI_FILE *PLI_STDERR;



static OSTYLE ostyles[] = {
  { "formatted",   " ",   "",   "",   "%.0s%.0s%s" },
  { "csv",         ",",   "",   "",   "%.0s%s%.0s" },
  { "json",        ",",   "{",  "}",  "\"%s\":%s%.0s" },
  { "last",        "",    "",   "",   "" }
};



static size_t gzvprintf(gzFile,const char*,va_list);



void init_io(void) {

  PLI_STDIN = new_file("stdin",stdin);
  PLI_STDOUT = new_file("stdout",stdout);
  PLI_STDERR = new_file("stderr",stderr);
}



void error_fn (char *format, ...) {

  va_list args;
  char message[MAX_LINE_LEN];

  va_start(args,format);

  vsprintf(message,format,args);

  write_line(PLI_STDERR,"Fatal error:\n%s\n",message);

  va_end (args);

  exit (1);
}



void warning_fn (char *format, ...) {

  va_list args;
  char message[MAX_LINE_LEN];

  if (!(params_get_parameter("warnings"))->value.i) {

    return;
  }

  va_start(args,format);

  vsprintf(message,format,args);

  write_line(PLI_STDERR,"WARNING:%s\n",message);

  va_end(args);
}


void debug_print(char *format, ...) {

  va_list args;
  char message[MAX_LINE_LEN];

  if (!(params_get_parameter("debug"))->value.i) {

    return;
  }

  va_start(args,format);

  vsprintf(message,format,args);

  write_line(PLI_STDERR,"Debug print:%s\n",message);

  va_end(args);
}


char* get_pli_dir(void) {

  char *pli_dir;

  pli_dir = getenv("PLI_DIR");

  if (pli_dir == NULL) {

    error_fn("get_pli_dir: PLI_DIR undefined");
  }

  return(pli_dir);
}


int is_readable_file(char *name_with_path) {
  struct stat s;
  if( stat(name_with_path, &s) == 0)
  {
    if(s.st_mode & S_IFREG)
    {
      if(s.st_mode & S_IRUSR)
	return 1;
    }
  }
  return 0;
}


PLI_FILE* new_file(char *filename,FILE *file) {

  int len;
  PLI_FILE *pli_file;

  pli_file = (PLI_FILE*) malloc(sizeof(PLI_FILE));

  if (pli_file == NULL) {

    error_fn("new_file: out of memory allocating pli_file");
  }

  strcpy(pli_file->filename,filename);

  len = strlen(filename);

  if (!strcmp(filename+len-3,".gz")) {

    pli_file->format = GZIPPED_FILE;

  } else {

    pli_file->format = ASCII_FILE;
  }

  pli_file->file = file;

  return(pli_file);
}



PLI_FILE* open_file(char *filename,char *mode) {

  FILE* file;
  PLI_FILE *pli_file;

  pli_file = new_file(filename,NULL);

  if (pli_file->format == GZIPPED_FILE) {

    pli_file->file = (!strcmp(mode,"r")) ? (FILE*) gzopen(filename,"rb") : (FILE*) gzopen(filename,"wb");

  } else {

    pli_file->file = fopen(filename,mode);
  }

  if (pli_file->file == NULL) {

    free(pli_file);

    return(NULL);
  }

  return(pli_file);
}



void close_file(PLI_FILE *pli_file) {

  if (pli_file == NULL) {

    warning_fn("close_file: attempting to close an undefined file");

    return;
  }

  if (pli_file->format == GZIPPED_FILE) {

    gzclose(pli_file->file);

  } else {

    fclose(pli_file->file);
  }

  free(pli_file);
}



int end_of_file(PLI_FILE *pli_file) {

  if (pli_file->format == GZIPPED_FILE) {

    if (gzeof(pli_file->file)) {

      return(1);
    }

  } else {

    if (feof(pli_file->file)) {

      return(1);
    }
  }

  return(0);
}



char* read_line(char *line,int bufsize,PLI_FILE *pli_file) {

  if (pli_file->format == GZIPPED_FILE) {

    return(gzgets(pli_file->file,line,bufsize));

  } else {
    
    return(fgets(line,bufsize,pli_file->file));
  }
}



void write_line(PLI_FILE *pli_file,const char *format,...) {

  va_list arguments;

  va_start(arguments,format);

  if (pli_file->format == GZIPPED_FILE) {
    
    gzvprintf(pli_file->file,format,arguments);
    
  } else {
    
    vfprintf(pli_file->file,format,arguments);
  }
  
  va_end(arguments);
}



long int pli_ftell(PLI_FILE *pli_file) {

  if (pli_file->format == GZIPPED_FILE) {

    return(gztell(pli_file->file));
  }

  return(ftell(pli_file->file));
}



int pli_fseek(PLI_FILE *pli_file,long int offset,int origin) {

  if (pli_file->format == GZIPPED_FILE) {

    return(gzseek(pli_file->file,offset,origin));
  }

  return(fseek(pli_file->file,offset,origin));
}



LIST* file2lines(char *filename) {

  char fline[MAX_LINE_LEN],*line;
  PLI_FILE *file;
  LIST *lines;

  file = open_file(filename,"r");

  if (file == NULL) {

    error_fn("%s: failed to open ff file '%s'",__func__,filename);
  }

  lines = new_list("lines",MAX_LINE_LEN*sizeof(char),0);

  while (!end_of_file(file)) {

    if (read_line(fline,MAX_LINE_LEN,file) == NULL) {

      break;
    }

    line = (char*) add_list_item(lines);

    strcpy(line,fline);
  }

  close_file(file);

  return(lines);
}



LIST* filter_lines(LIST *lines_in,char *str) {

  int i;
  char name[MAX_LINE_LEN],word[MAX_LINE_LEN],*line_in,*line_out;
  LIST *lines_out;

  sprintf(name,"%s list",str);

  lines_out = new_list(name,lines_in->item_size,0);

  for (i=0;i<lines_in->n_items;i++) {

    line_in = (char*) get_list_item(lines_in,i);

    if ((sscanf(line_in,"%s",word) == 1) && (!strcmp(word,str))) {

      line_out = (char*) add_list_item(lines_out);

      strcpy(line_out,line_in);
    }
  }

  return(lines_out);
}



int read_word(char *s,const char *format,void *word) {

  int len,n_read;
  char subs[MAX_LINE_LEN];

  n_read = sscanf(format,"%%%d[a-z]",&len);

  if (n_read == 1) {

    substring(s,0,len,subs);

    return(sscanf(subs,format,word));

  } else {

    return(sscanf(s,format,word));
  }
}



void substring(char *s1,int start,int length,char *s2) {

  strncpy(s2,s1+start,length);
  s2[length] = '\0';
}



void remove_spaces(char *s1,char *s2) {

  int i,l1,l2;

  l1 = strlen(s1);
  l2 = 0;

  for (i=0;i<l1;i++) {

    if (s1[i] != ' ') {

      s2[l2++] = s1[i];
    }
  }

  s2[l2] = '\0';
}



void remove_outer_spaces(char *s1,char *s2) {

  int start,end,l1;

  l1 = strlen(s1);

  start = 0;

  if (start < l1) {

    end = l1 - 1;

    while ((s1[start] == ' ') && (start < l1)) {

      start++;
    }

    while ((s1[end] == ' ') && (end > 0)) {

      end--;
    }
  }

  if (start < l1) {

    substring(s1,start,end-start+1,s2);

  } else {

    strcpy(s2,"");
  }
}



void upper_case(char *s) {

  int i;

  for (i=0;s[i];i++) {

    s[i] = toupper((unsigned char) s[i]);
  }
}



enum OUTPUT_FORMAT get_output_format(char *format) {

  if ((!strcmp(format,"text")) || (!strcmp(format,"formatted"))) {

    return(FORMATTED);

  } else  if (!strcmp(format,"json")) {

    return(JSON);

  } else  if (!strcmp(format,"csv")) {

    return(CSV);
  }

  error_fn("get_output_format: no such output format '%s'",format);
}



void obj2text(void *obj,char *text,OFIELD *ofields) {

  int len,pos;
  char field_text[MAX_LINE_LEN];
  OFIELD *field;
  OSTYLE *ostyle;

  strcpy(text,"");

  ostyle = NULL;

  field = ofields;

  len = pos = 0;

  while (strcmp(field->name,"last")) {

    if (field->ostyle) {

      ostyle = field->ostyle;

      field->obj2text((void*) obj,field_text,field->format); 

      len += strlen(field_text);

      if (len > MAX_LINE_LEN-5) {

	error_fn("%s: exceeded MAX_LINE_LEN",__func__);
      }

      (pos == 0) ? sprintf(text,"%s%s",field->ostyle->open,field_text) : sprintf(text+pos,"%s%s",field->ostyle->separator,field_text);

      len = strlen(text);

      pos = len;
    }

    field++;
  }

  if (ostyle) {

    sprintf(text+pos,"%s",ostyle->close);
  }
}



OSTYLE* get_ostyle(char *name) {

  OSTYLE *style;

  style = ostyles;

  while (strcmp(style->name,"last")) {

    if (!strcmp(style->name,name)) {

      return(style);
    }

    style++;
  }

  return(NULL);
}



OGROUP* get_ogroup(OGROUP *groups,char *name) {

  OGROUP *group;

  group = groups;

  while (strcmp(group->name,"last")) {

    if (!strcmp(group->name,name)) {

      return(group);
    }

    group++;
  }

  return(NULL);
}



OFIELD* get_ofield(OFIELD *ofields,char *name) {

  OFIELD *field;

  field = ofields;

  while (strcmp(field->name,"last")) {

    if (!strcmp(field->name,name)) {

      return(field);
    }

    field++;
  }

  return(NULL);
}



void init_ofields(OFIELD *ofields) {

  OFIELD *field;

  field = ofields;

  while (strcmp(field->name,"last")) {

    field->ostyle = NULL;

    field++;
  }
}



void set_ofields(OFIELD *ofields,char *style,char *flags) {

  int pos = 0;
  char *flag,cflags[MAX_LINE_LEN];
  OFIELD *field;

  strcpy(cflags,flags);

  while (flag = nextword(cflags,',',&pos)) {

    field = get_ofield(ofields,flag);

    if (field) {

      set_ofield(field,style);
    }
  }
}



void set_ofield(OFIELD *field,char *style) {

  OSTYLE* ostyle;

  ostyle = get_ostyle(style);

  if (ostyle) {

    field->ostyle = ostyle;

    sprintf(field->format,ostyle->format,field->name,field->uformat,field->fformat);

  } else {

    error_fn("%s: no such output format '%s'",__func__,style);
  }
}



void set_ogroups(OGROUP *ogroups,char *style,char *flags) {

  int pos = 0;
  char *flag,cflags[MAX_LINE_LEN];
  OGROUP *group;

  strcpy(cflags,flags);

  while (flag = nextword(cflags,',',&pos)) {

    group = get_ogroup(ogroups,flag);

    if (group) {

      set_ofields(group->ofields,style,group->fields);
    }
  }
}



static size_t gzvprintf(gzFile gz,const char *fmt, va_list args) {

  char *buf;
  size_t len;

  if (vasprintf(&buf, fmt, args) < 0)
    return 0;

  len = strlen(buf);

  len = len == (unsigned)len ? (size_t)gzwrite(gz, buf,(unsigned)len) : 0;

  free(buf);

  return len;
}
