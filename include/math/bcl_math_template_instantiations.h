// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_MATH_TEMPLATE_INSTANTIATIONS_H_
#define BCL_MATH_TEMPLATE_INSTANTIATIONS_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_binary_function_interface.h"
#include "bcl_math_binary_sum_function.h"
#include "bcl_math_function_interface.h"
#include "bcl_math_mutate_combine.h"
#include "bcl_math_mutate_decision_node.h"
#include "bcl_math_mutate_interface.h"
#include "bcl_math_mutate_move_wrapper.h"
#include "bcl_math_mutate_perturbation.h"
#include "bcl_math_sum_function.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_interface.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "biol/bcl_biol_membrane.h"
#include "coord/bcl_coord_geometry_interface.h"
#include "linal/bcl_linal_vector.h"
#include "score/bcl_score_protein_model.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    BCL_EXPIMP_TEMPLATE template class BCL_API FunctionInterfaceSerializable< bcl::linal::Vector< double>, double>;

    BCL_EXPIMP_TEMPLATE template class BCL_API FunctionInterfaceSerializable< bcl::assemble::ProteinModel, double>;

    BCL_EXPIMP_TEMPLATE template class BCL_API BinaryFunctionInterface< bcl::assemble::SSEGeometryInterface, bcl::assemble::SSE, double>;

    BCL_EXPIMP_TEMPLATE template class BCL_API FunctionInterfaceSerializable< double, double>;

    BCL_EXPIMP_TEMPLATE template class BCL_API FunctionInterfaceSerializable< bcl::storage::VectorND< 2, bcl::linal::Vector< double> >, bcl::storage::VectorND< 2, bcl::linal::Vector< double> > >;

    BCL_EXPIMP_TEMPLATE template class BCL_API MutateCombine< bcl::assemble::SSEPool>;

    // BCL_EXPIMP_TEMPLATE template class BCL_API MutateDecisionNode< bcl::assemble::ProteinModel>;

    BCL_EXPIMP_TEMPLATE template class BCL_API MutateDecisionNode< bcl::assemble::SSEPool>;

    BCL_EXPIMP_TEMPLATE template class BCL_API MutateInterface< bcl::assemble::Domain>;

    BCL_EXPIMP_TEMPLATE template class BCL_API MutateInterface< bcl::assemble::ProteinModel>;

	BCL_EXPIMP_TEMPLATE template class BCL_API MutateInterface< bcl::math::TransformationMatrix3D>;

    BCL_EXPIMP_TEMPLATE template class BCL_API SumFunctionMixin< bcl::score::ProteinModel>;

    BCL_EXPIMP_TEMPLATE template class BCL_API BinarySumFunction< bcl::assemble::SSEPool, bcl::biol::Membrane, double, double>;

    BCL_EXPIMP_TEMPLATE template class BCL_API MutateMoveWrapper< bcl::assemble::SSE>;

    BCL_EXPIMP_TEMPLATE template class BCL_API MutatePerturbation< bcl::assemble::ProteinModel>;

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_TEMPLATE_INSTANTIATIONS_H_
