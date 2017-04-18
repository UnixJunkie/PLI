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



double lennard_jones_energy(double d,double Do,double Eo,double order) {

  if (d < 1.0E-10) {

    error_fn("lennard_jones_energy: d is too short (%.6lfA)",d);
  }

  return(Eo*(pow(Do/d,2.0*((double) order)) - 2.0*pow(Do/d,(double) order)));
}
