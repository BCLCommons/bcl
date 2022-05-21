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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "math/bcl_math_template_instantiations.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    template class BCL_API FunctionInterfaceSerializable< linal::Vector< double>, double>;

    template class BCL_API FunctionInterfaceSerializable< assemble::ProteinModel, double>;

    template class BCL_API BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, double>;

    template class BCL_API FunctionInterfaceSerializable< double, double>;

    template class BCL_API FunctionInterfaceSerializable< storage::VectorND< 2, linal::Vector< double> >, storage::VectorND< 2, linal::Vector< double> > >;

    template class BCL_API MutateCombine< assemble::SSEPool>;

    template class BCL_API MutateDecisionNode< assemble::ProteinModel>;

    template class BCL_API MutateDecisionNode< assemble::SSEPool>;

    template class BCL_API MutateInterface< assemble::Domain>;

    template class BCL_API MutateInterface< assemble::ProteinModel>;

    template class BCL_API MutateInterface< TransformationMatrix3D>;

    template class BCL_API SumFunctionMixin< score::ProteinModel>;

    template class BCL_API BinarySumFunction< assemble::SSEPool, biol::Membrane, double, double>;

    template class BCL_API MutateMoveWrapper< assemble::SSE>;

    template class BCL_API MutatePerturbation< assemble::ProteinModel>;

  } // namespace math
} // namespace bcl
