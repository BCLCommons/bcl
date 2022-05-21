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

#ifndef BCL_FIND_TEMPLATE_INSTANTIATIONS_H_
#define BCL_FIND_TEMPLATE_INSTANTIATIONS_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_find_collector_criteria_combined.h"
#include "bcl_find_collector_criteria_wrapper.h"
#include "bcl_find_locator.h"
#include "bcl_find_locator_criteria.h"
#include "bcl_find_locator_criteria_wrapper.h"
#include "bcl_find_pick_criteria_interface.h"
#include "bcl_find_pick_criteria_wrapper.h"
#include "bcl_find_pick_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {

    BCL_EXPIMP_TEMPLATE template class BCL_API PickCriteriaInterface< util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::ProteinModel>;

    BCL_EXPIMP_TEMPLATE template class PickCriteriaWrapper< util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::DomainInterface>;

    BCL_EXPIMP_TEMPLATE template class BCL_API PickCriteriaInterface< util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, linal::Vector3D>;

    BCL_EXPIMP_TEMPLATE template class BCL_API LocatorCriteria< util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE, util::SiPtrList< const assemble::SSE> >;

    BCL_EXPIMP_TEMPLATE template class BCL_API CollectorCriteriaCombined< assemble::SSE>;

    BCL_EXPIMP_TEMPLATE template class BCL_API CollectorCriteriaWrapper< util::SiPtrList< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface>;

    BCL_EXPIMP_TEMPLATE template class BCL_API PickCriteriaWrapper< util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE>;

    BCL_EXPIMP_TEMPLATE template class BCL_API Locator< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, util::SiPtrList< const assemble::SSE> >;

    BCL_EXPIMP_TEMPLATE template class BCL_API LocatorCriteriaWrapper< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::SSE>;

    BCL_EXPIMP_TEMPLATE template class BCL_API LocatorCriteriaWrapper< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface>;

  } // namespace find
} // namespace bcl

#endif // BCL_FIND_TEMPLATE_INSTANTIATIONS_H_ 
