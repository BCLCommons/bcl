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
#include "assemble/bcl_assemble_locator_sse_unpaired.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LocatorSSEUnpaired::LocatorSSEUnpaired() :
      m_Collector(),
      m_Pick()
    {
    }

    //! @brief constructor taking a all members
    //! @param PICK is the ShPtr to a PickInterface which will determine which unpaired SSE is given
    //! @param CONTACT_TYPE is the contact type for which an SSE not having it is desired
    //! @param MAX_DISTANCE is the maximum distance between SSEs for them to be considered paired
    LocatorSSEUnpaired::LocatorSSEUnpaired
    (
      const find::PickCriteriaInterface< util::SiPtr< const SSE>, util::SiPtrList< const SSE>, linal::Vector3D> &PICK,
      const contact::Type &CONTACT_TYPE,
      const double MAX_DISTANCE
    ) :
      m_Collector( CONTACT_TYPE, MAX_DISTANCE),
      m_Pick( PICK.Clone())
    {
    }

    //! @brief constructor taking a all members
    //! @param SP_PICK is the ShPtr to a PickInterface which will determine which unpaired SSE is given
    //! @param CONTACT_TYPE is the contact type for which an SSE not having it is desired
    //! @param MAX_DISTANCE is the maximum distance between SSEs for them to be considered paired
    LocatorSSEUnpaired::LocatorSSEUnpaired
    (
      const util::ShPtr< find::PickCriteriaInterface< util::SiPtr< const SSE>, util::SiPtrList< const SSE>, linal::Vector3D> > &SP_PICK,
      const contact::Type &CONTACT_TYPE,
      const double MAX_DISTANCE
    ) :
      m_Collector( CONTACT_TYPE, MAX_DISTANCE),
      m_Pick( SP_PICK)
    {
    }

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorSSEFurthest which is a copy of this
    LocatorSSEUnpaired *LocatorSSEUnpaired::Clone() const
    {
      return new LocatorSSEUnpaired( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorSSEUnpaired::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief SetPick changes the Pick to a different object
    //! @param PICK ShPtr to PickInterface< SSE> which "m_Pick" will be changed to
    void LocatorSSEUnpaired::SetPick
    (
      const util::ShPtr< find::PickCriteriaInterface< util::SiPtr< const SSE>, util::SiPtrList< const SSE>, linal::Vector3D> > &PICK
    )
    {
      m_Pick = PICK;
    }

    //! @brief GetMaxDistance gives the maximum distance between SSEs for them to still be considered paired
    //! @return returns "m_MaxDistance" the maximum distance between SSEs for them to still be considered paired
    const CollectorSSEUnpaired &LocatorSSEUnpaired::GetCollector() const
    {
      return m_Collector;
    }

    //! @brief SetMaxDistance changes the type of SSE contact of interest
    //! @param COLLECTOR is the new maximum distance between SSEs for them to still be considered paired
    void LocatorSSEUnpaired::SetCollector( const CollectorSSEUnpaired &COLLECTOR)
    {
      m_Collector = COLLECTOR;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Locate returns an SSE which is unpaired
    //! @param SSE_DOMAIN domain which the LocatorSSEUnpaired refers to
    //! @return returns SiPtr to the unpaired SSE furthest from the center of SSE_DOMAIN; empty if no SSEs in domain
    util::SiPtr< const SSE> LocatorSSEUnpaired::Locate( const DomainInterface &SSE_DOMAIN) const
    {
      // return empty SiPtr to SSE if there is no unpaired SSE
      return m_Pick->Pick( m_Collector.Collect( SSE_DOMAIN), SSE_DOMAIN.GetCenter());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LocatorSSEUnpaired::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Collector, ISTREAM);
      io::Serialize::Read( m_Pick, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @param INDENT indentation
    //! @return ostream which was read from
    std::ostream &LocatorSSEUnpaired::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Collector, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Pick, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> LocatorSSEUnpaired::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorSSEUnpaired())
    );

  } // namespace assemble
} // namespace bcl
