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

#ifndef BCL_ASSEMBLE_LOCATOR_SSE_UNPAIRED_H_
#define BCL_ASSEMBLE_LOCATOR_SSE_UNPAIRED_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_collector_sse_unpaired.h"
#include "contact/bcl_contact_types.h"
#include "find/bcl_find_locator_interface.h"
#include "find/bcl_find_pick_criteria_interface.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorSSEUnpaired
    //! @brief locates an unpaired strand of a protein model.
    //! @details this class uses its member to collect all unpaired strands in the given model and uses the given
    //! pick criteria to pick one and returns it.
    //!
    //! @see @link example_assemble_locator_sse_unpaired.cpp @endlink
    //! @author alexanns, woetzen, karakam
    //! @date 03/18/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorSSEUnpaired :
      public find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface>
    {

    private:

    //////////
    // data //
    //////////

      //! the Collector which will gather the unpaired SSEs
      CollectorSSEUnpaired m_Collector;

      //! the ShPtr to a Pick object which will choose which unpaired SSE is located
      util::ShPtr< find::PickCriteriaInterface< util::SiPtr< const SSE>, util::SiPtrList< const SSE>, linal::Vector3D> > m_Pick;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorSSEUnpaired();

      //! @brief constructor taking a all members
      //! @param PICK is the ShPtr to a PickInterface which will determine which unpaired SSE is given
      //! @param CONTACT_TYPE is the contact type for which an SSE not having it is desired
      //! @param MAX_DISTANCE is the maximum distance between SSEs for them to be considered paired
      LocatorSSEUnpaired
      (
        const find::PickCriteriaInterface< util::SiPtr< const SSE>, util::SiPtrList< const SSE>, linal::Vector3D> &PICK,
        const contact::Type &CONTACT_TYPE,
        const double MAX_DISTANCE
      );

      //! @brief constructor taking a all members
      //! @param SP_PICK is the ShPtr to a PickInterface which will determine which unpaired SSE is given
      //! @param CONTACT_TYPE is the contact type for which an SSE not having it is desired
      //! @param MAX_DISTANCE is the maximum distance between SSEs for them to be considered paired
      LocatorSSEUnpaired
      (
        const util::ShPtr< find::PickCriteriaInterface< util::SiPtr< const SSE>, util::SiPtrList< const SSE>, linal::Vector3D> > &SP_PICK,
        const contact::Type &CONTACT_TYPE,
        const double MAX_DISTANCE
      );

      //! @brief Clone is the virtual Clone constructor
      //! @return a pointer to new LocatorSSEFurthest which is a copy of this
      LocatorSSEUnpaired *Clone() const;

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief SetPick changes the Pick to a different object
      //! @param PICK ShPtr to PickInterface< SSE> which "m_Pick" will be changed to
      void SetPick( const util::ShPtr< find::PickCriteriaInterface< util::SiPtr< const SSE>, util::SiPtrList< const SSE>, linal::Vector3D> > &PICK);

      //! @brief GetMaxDistance gives the maximum distance between SSEs for them to still be considered paired
      //! @return returns "m_MaxDistance" the maximum distance between SSEs for them to still be considered paired
      const CollectorSSEUnpaired &GetCollector() const;

      //! @brief SetColector Set the collector
      //! @param COLLECTOR the new collector
      void SetCollector( const CollectorSSEUnpaired &COLLECTOR);

    ////////////////
    // operations //
    ////////////////

      //! @brief Locate returns an SSE which is unpaired
      //! @param SSE_DOMAIN domain which the LocatorSSEUnpaired refers to
      //! @return returns SiPtr to the unpaired SSE furthest from the center of SSE_DOMAIN; empty if no SSEs in domain
      util::SiPtr< const SSE> Locate( const DomainInterface &SSE_DOMAIN) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief read from std::ostream
      //! @param OSTREAM input stream
      //! @return ostream which was read from
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class LocatorSSEUnpaired

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_LOCATOR_SSE_UNPAIRED_H_
