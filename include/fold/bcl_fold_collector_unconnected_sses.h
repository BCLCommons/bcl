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

#ifndef BCL_FOLD_COLLECTOR_UNCONNECTED_SSES_H_
#define BCL_FOLD_COLLECTOR_UNCONNECTED_SSES_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain_interface.h"
#include "biol/bcl_biol_aa_sequence_flexibility.h"
#include "find/bcl_find_collector_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @class CollectorUnconnectedSSE
    //! @brief Collects unconnected SSEs from a given domain
    //!
    //! @see @link example_fold_collector_unconnected_sses.cpp @endlink
    //! @author woetzen, alexanns, fischea
    //! @date Sep 10, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API CollectorUnconnectedSSE :
      public find::CollectorInterface< util::SiPtrList< const assemble::SSE>, assemble::DomainInterface>
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      //! sequence direction in which to collect the SSEs
      biol::AASequenceFlexibility::DirectionEnum m_Direction;

      //! whether to
      bool m_TestBondConnected;

      //! types of the SSEs to consider
      storage::Set< biol::SSType> m_SSTypes;

      //! whether to ignore terminal loops
      bool m_IgnoreTermini;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      CollectorUnconnectedSSE();

      //! @brief construct from sequence direction
      //! @param DIRECTION the sequence direction is
      //! @param TEST_BOND_CONNECTION
      //! @param SS_TYPES
      //! @param IGNORE_TERMINAL_SSE
      CollectorUnconnectedSSE
      (
        const biol::AASequenceFlexibility::SequenceDirection &DIRECTION,
        const bool TEST_BOND_CONNECTION,
        const storage::Set< biol::SSType> &SS_TYPES,
        const bool IGNORE_TERMINAL_SSE
      );

      //! @brief clone function
      //! @return pointer to a new CollectorUnconnectedSSE
      CollectorUnconnectedSSE *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief collects unconnected SSEs from the given domain
      //! @param DOMAIN_INTERFACE domain from which to collect unconnected SSEs
      //! @return unconnected SSEs in the given domain
      util::SiPtrList< const assemble::SSE> Collect( const assemble::DomainInterface &DOMAIN_INTERFACE) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class CollectorUnconnectedSSE

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_COLLECTOR_UNCONNECTED_SSES_H_
