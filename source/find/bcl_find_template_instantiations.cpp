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
#include "find/bcl_find_template_instantiations.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {

    template class BCL_API PickCriteriaInterface< util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::ProteinModel>;

    template class PickCriteriaWrapper< util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::DomainInterface>;

    template class BCL_API PickCriteriaInterface< util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, linal::Vector3D>;

    template class BCL_API LocatorCriteria< util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE, util::SiPtrList< const assemble::SSE> >;

    template class BCL_API CollectorCriteriaCombined< assemble::SSE>;

    template class BCL_API CollectorCriteriaWrapper< util::SiPtrList< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface>;

    template class BCL_API PickCriteriaWrapper< util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE>;

    template class BCL_API Locator< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, util::SiPtrList< const assemble::SSE> >;

    template class BCL_API LocatorCriteriaWrapper< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::SSE>;

    template class BCL_API LocatorCriteriaWrapper< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface>;

  } // namespace find
} // namespace bcl
