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



static PLI_MODE pli_modes[] = { { "contacts",     run_contacts,      "ligand",  0 },
				{ "type",         run_type,          "ligand",  0 },
				{ "score",        run_score,         "complex", 1 },
				{ "field",        run_field,         "site",    0 },
				{ "help",         run_help,          "",        0 },
				{ "last",         NULL,              "",        0 } };



PLI_MODE* get_mode(char *name) {

  PLI_MODE *mode;

  mode = pli_modes;

  while (strcmp(mode->name,"last")) {

    if (!strcmp(mode->name,name)) {

      return(mode);
    }

    mode++;
  }

  error_fn("get_mode: unknown mode '%s'",name);
}
