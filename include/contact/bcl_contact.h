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

#ifndef BCL_CONTACT_H_
#define BCL_CONTACT_H_

// include the namespace forward header
#include "bcl_contact.fwd.hh"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_contact.h
  //! @brief namespace for residue - residue contact prediction and related tools
  //! @details This namespace provides classes necessary for residue-residue contact prediction, sucn as neural networks
  //! used in prediction, PredictionMap to store the predictions, ContactMap to store the actual contact information
  //! ContactOrder to calculate the contact order of a given model.
  //!
  //! @see @link example_contact.cpp @endlink
  //! @author karakam
  //! @date Jul 22, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace contact
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    BCL_API const std::string &GetNamespaceIdentifier();

    //! Number of different contacts definitions in the enumerator ContactDefinitions
    const size_t g_NumberContactDefinitions = 4;

    //! Enumerators for ContactDefinitions, used instead of 0|1 or -2|-1|0|1
    enum ContactDefinitions
    {
      e_TrueContact,
      e_TrueNonContact,
      e_WrongContact,
      e_OtherNonContact
    };

    //! vector that contains string versions of CotnactDefinitions
    const std::string g_ContactDefinitionNames[ g_NumberContactDefinitions] =
    {
      "TrueContact",
      "TrueNonContact",
      "WrongContact",
      "OtherNonContact"
    };

    //! distance to be considered as contact
    const double g_ContactCbDistanceCutoff( 8.0);

    //! Required minimal sequence separation for two residues to be considered in contact
    const size_t g_ContactMinSequenceSeparation = 6;

    //! @brief return default sequence separation range for contacts
    //! @return default sequence separation range contacts
    const math::Range< size_t> &GetDefaultSequenceSeparationRange();

    //! @brief return sequence separation range for short-range contacts
    //! @return sequence separation range for short-range contacts
    const math::Range< size_t> &GetDefaultSequenceSeparationShortRange();

    //! @brief return sequence separation range for mid-range contacts
    //! @return sequence separation range for mid-range contacts
    const math::Range< size_t> &GetDefaultSequenceSeparationMidRange();

    //! @brief return sequence separation range for long-range contacts
    //! @return sequence separation range for long-range contacts
    const math::Range< size_t> &GetDefaultSequenceSeparationLongRange();

  } // namespace contact
} // namespace bcl

#endif //  BCL_CONTACT_H_
