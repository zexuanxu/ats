/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  License: BSD
  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "column_subsidence_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,ColumnSubsidenceEvaluator> ColumnSubsidenceEvaluator::reg_("column subsidence");

}
}
