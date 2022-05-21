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

#ifndef BCL_STORAGE_TEMPLATE_INSTANTIATIONS_H_
#define BCL_STORAGE_TEMPLATE_INSTANTIATIONS_H_

// include the namespace header
#include "bcl_storage.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_storage_list.h"
#include "bcl_storage_pair.h"
#include "bcl_storage_vector.h"
#include "bcl_storage_vector_nd.h"
#include "math/bcl_math_template_instantiations.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace storage
  {

    BCL_EXPIMP_TEMPLATE template class BCL_API Vector< VectorND< 2, linal::Vector< double> > >;

    BCL_EXPIMP_TEMPLATE template class BCL_API VectorND< 2, linal::Vector< double> >;

    BCL_EXPIMP_TEMPLATE template class BCL_API Pair< double, util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, math::MutateResult< assemble::ProteinModel> > > >;

    BCL_EXPIMP_TEMPLATE template class BCL_API Pair< double, util::ShPtr< math::FunctionInterfaceSerializable< assemble::SSEPool, math::MutateResult< assemble::SSEPool> > > >;

  } // namespace storage
} // namespace bcl

#endif // BCL_STORAGE_TEMPLATE_INSTANTIATIONS_H_
