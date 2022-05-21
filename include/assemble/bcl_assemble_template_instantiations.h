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

#ifndef BCL_ASSEMBLE_TEMPLATE_INSTANTIATIONS_H_
#define BCL_ASSEMBLE_TEMPLATE_INSTANTIATIONS_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_locator_sse_furthest.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    // Note: If desired (e.g. for purposes of dynamic selection on the command line), a list of all
    // Criterions can be created with this command:
    // grep 'i::Criterion[A-Za-z \n]*<[ \n]*[a-z][^_]' ./source/*/* ./include/*/* ./apps/* ./apps/*/* ./apps/*/*/*
    //      ./apps/*/*/*/* ./example/* ./example/*/* | sed 's/^.*\.[chp]*:[ \t]*//' | sed 's/\([^-]\)>[^>]*$/\1>/' |
    //      sed 's/new //' | sed 's/.*opti::/opti::/' | sed 's/> >.*/>/' | sort | uniq | sed 's/opti:://' |
    //      awk '{print "     BCL_EXPIMP_TEMPLATE template class BCL_API "$0";"}'

    BCL_EXPIMP_TEMPLATE template class BCL_API LocatorSSEFurthest< DomainInterface>;

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_TEMPLATE_INSTANTIATIONS_H_
