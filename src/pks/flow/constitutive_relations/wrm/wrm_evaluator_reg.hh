#include "wrm_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,WRMEvaluator> WRMEvaluator::factory_("WRM");

} //namespace
} //namespace
} //namespace
