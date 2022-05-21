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

#ifndef BCL_ASSEMBLE_LOCATOR_SUB_DOMAIN_RANDOM_H_
#define BCL_ASSEMBLE_LOCATOR_SUB_DOMAIN_RANDOM_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "find/bcl_find_locator_interface.h"
#include "math/bcl_math_range.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorSubDomainRandom
    //! @brief class for locating a sub-domain from a given domain
    //! @details This class depending on the range of number SSEs given, constructs a sub-domain from the given domain.
    //! If the given domain has topology defined, then parts of graph that correspond to SSEs not selected for the
    //! sub-domain are pruned. Boolean member m_LocateConsecutive decides whether the subdomain should be composed
    //! of consecutive SSEs. m_UseTopologyOrrder determines whether the ordering should be according to the elements
    //! vector in the topology, this way sub-beta-sheets can be collected.
    //!
    //! @see @link example_assemble_locator_sub_domain_random.cpp @endlink
    //! @author karakam
    //! @date Mar 16, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorSubDomainRandom :
      public find::LocatorInterface< util::ShPtr< Domain>, util::ShPtr< Domain> >
    {

    private:

    //////////
    // data //
    //////////

      //! range for number of SSEs
      math::Range< size_t> m_SizeRange;

      //! boolean to whether to locate consecutive SSEs
      bool m_LocateConsecutive;

      //! boolean to order SSEs by the topology order
      bool m_UseTopologyOrder;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorSubDomainRandom();

      //! @brief construct from a range for number of SSEs in sub-domain
      //! @param SIZE_RANGE range for number of SSEs
      LocatorSubDomainRandom( const math::Range< size_t> SIZE_RANGE);

      //! @brief construct from a range for number of SSEs in sub-domain and whether they should be located consecutively
      //! @param SIZE_RANGE range for number of SSEs
      //! @param LOCATE_CONSECUTIVE boolean to whether to locate consecutive SSEs
      //! @param USE_TOPOLOGY_ORDER boolean to order SSEs by the topology order
      LocatorSubDomainRandom
      (
        const math::Range< size_t> &SIZE_RANGE,
        bool LOCATE_CONSECUTIVE,
        bool USE_TOPOLOGY_ORDER
      );

      //! @brief Clone function
      //! @return pointer to new LocatorSubDomainRandom
      LocatorSubDomainRandom *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get size range
      //! @return size range
      const math::Range< size_t> &GetSizeRange() const
      {
        return m_SizeRange;
      }

      //! @brief get locate consecutive boolean
      //! @return locate consecutive boolean
      bool GetLocateConsecutive() const
      {
        return m_LocateConsecutive;
      }

      //! @brief get use topology order boolean
      //! @return size use topology order boolean
      bool GetUseTopologyOrder() const
      {
        return m_UseTopologyOrder;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief locate and return a random sub-domain from the given domain
      //! @param SP_DOMAIN ShPtr to the domain of interest
      //! @return ShPtr to a random sub-domain located within the given SP_DOMAIN
      util::ShPtr< Domain> Locate( const util::ShPtr< Domain> &SP_DOMAIN) const;

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

    }; // class LocatorSubDomainRandom

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_LOCATOR_SUB_DOMAIN_RANDOM_H_ 
