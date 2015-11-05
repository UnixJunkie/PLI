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



static PLI_SFUNC pli_sfuncs[] = { { "pliff",
				    pliff_score_system,pliff_minimise_system,
				    pliff_write_system_scores,pliff_write_atom_scores,pliff_write_contact_scores },
				  { "last",
				    NULL,NULL,
				    NULL,NULL,NULL } };



PLI_SFUNC* get_sfunc(char *name) {

  PLI_SFUNC *sfunc;

  sfunc = pli_sfuncs;

  while (strcmp(sfunc->name,"last")) {

    if (!strcmp(sfunc->name,name)) {

      return(sfunc);
    }

    sfunc++;
  }

  error_fn("get_sfunc: unknown scoring function '%s'",name);
}
