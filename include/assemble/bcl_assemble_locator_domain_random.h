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

#ifndef BCL_ASSEMBLE_LOCATOR_DOMAIN_RANDOM_H_
#define BCL_ASSEMBLE_LOCATOR_DOMAIN_RANDOM_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_locator_interface.h"
#include "math/bcl_math_range.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorDomainRandom
    //! @brief locates a random subset of SSEs in the model and returns them as a domain
    //! @details Using the given SSType, it collects a random subset of SSEs in the given protein model and returns
    //! them as a domain
    //!
    //! @see @link example_assemble_locator_domain_random.cpp @endlink
    //! @author karakam
    //! @date Feb 3, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorDomainRandom :
      public find::LocatorInterface< util::ShPtr< Domain>, ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! range for number sses in domain to be located
      math::Range< size_t> m_DomainSizeRange;

      //! sstypes to be incldues in the domain
      biol::SSType m_SSType;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorDomainRandom();

      //! @brief constructor from a domain size range and SSTypes
      //! @param DOMAIN_SIZE_RANGE min and max sizes of the domain to be collected
      //! @param SS_TYPE SSType to be collected
      LocatorDomainRandom
      (
        const math::Range< size_t> &DOMAIN_SIZE_RANGE,
        const biol::SSType &SS_TYPE
      );

      //! @brief Clone function
      //! @return pointer to new LocatorDomainRandom
      LocatorDomainRandom *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the domain size range
      //! @return the domain size range
      const math::Range< size_t> &GetDomainSizeRange() const;

      //! @brief returns the sstype
      //! @return the sstype
      const biol::SSType &GetSSType() const;

    ////////////////
    // operations //
    ////////////////

       //! @brief locates a random domain and returns it
       //! @param PROTEIN_MODEL ProteinModel of interest
       //! @return randomly located domain
       util::ShPtr< Domain> Locate( const ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class LocatorDomainRandom

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_LOCATOR_DOMAIN_RANDOM_H_
