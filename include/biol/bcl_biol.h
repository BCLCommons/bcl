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

#ifndef BCL_BIOL_H_
#define BCL_BIOL_H_

// include the namespace forward header
#include "bcl_biol.fwd.hh"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl.h"

// external includes - sorted alphabetically

namespace bcl
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_biol.h
  //! @brief namespace for biological classes and functions in the biochemistry library
  //! @details amino acid sequence, amino acids, atoms, nucleic acids and membrane derived calsses and properties are
  //! gathered in that namespace
  //!
  //! @see @link example_biol.cpp @endlink
  //! @author karakam
  //! @date Jul 22, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace biol
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    BCL_API
    const std::string &GetNamespaceIdentifier();

    //! number of nucleic acids
    const size_t g_NumberNATypes( 9);

    //! custom type nucleic acid
    enum NATypes { e_Adenosine, e_Guanosine, e_Cytidine, e_Uridine, e_Desoxyadenosine, e_Desoxyguanosine, e_Desoxycytidine, e_Desoxythymidine, e_Undefined_na};

  } // namespace biol
} // namespace bcl

#endif //BCL_BIOL_H_
