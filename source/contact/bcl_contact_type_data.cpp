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
#include "contact/bcl_contact_type_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> TypeData::s_Instance( GetObjectInstances().AddInstance( new TypeData()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct undefined contact type
    TypeData::TypeData() :
      m_WindowRadii( util::GetUndefined< size_t>(), util::GetUndefined< size_t>()),
      m_WindowLengths( util::GetUndefined< size_t>(), util::GetUndefined< size_t>()),
      m_IsValid( false),
      m_ResidueDistanceCutoff( util::GetUndefinedDouble()),
      m_MinimalSSEDistance( util::GetUndefinedDouble()),
      m_DistanceRange(),
      m_PreferredDistanceRange(),
      m_TiltAngleRange(),
      m_MinimalFragmentInterfaceLength( util::GetUndefinedDouble()),
      m_SSTypes( biol::GetSSTypes().e_Undefined)
    {
    }

    //! @brief construct contact type from provided data
    //! @param WINDOW_RADII radius of window used for representation
    //! @param WINDOW_LENGTH_PAIR window length used for representation
    //! @param IS_VALID whether a valid contact type
    //! @param RESIDUE_DISTANCE_CUTOFF cut off distance
    //! @param MINIMAL_SSE_DISTANCE minimal distance of two sses to not clash
    //! @param MINIMAL_FRAGMENT_INTERFACE_LENGTH minimal fragment interface length
    //! @param DISTANCE_RANGE distance range
    //! @param PREFERRED_DISTANCE_RANGE preferred distance range
    //! @param TILT_ANGLE_RANGE tilt angle range
    //! @param SS_TYPES SSTypes associated with this contact type
    TypeData::TypeData
    (
      const storage::Pair< size_t, size_t> &WINDOW_RADII,
      const storage::Pair< size_t, size_t> &WINDOW_LENGTH_PAIR,
      const bool IS_VALID,
      const double RESIDUE_DISTANCE_CUTOFF,
      const double MINIMAL_SSE_DISTANCE,
      const math::Range< double> &DISTANCE_RANGE,
      const math::Range< double> &PREFERRED_DISTANCE_RANGE,
      const math::Range< double> &TILT_ANGLE_RANGE,
      const double MINIMAL_FRAGMENT_INTERFACE_LENGTH,
      const storage::Set< biol::SSType> &SS_TYPES
    ) :
      m_WindowRadii( WINDOW_RADII),
      m_WindowLengths( WINDOW_LENGTH_PAIR),
      m_IsValid( IS_VALID),
      m_ResidueDistanceCutoff( RESIDUE_DISTANCE_CUTOFF),
      m_MinimalSSEDistance( MINIMAL_SSE_DISTANCE),
      m_DistanceRange( DISTANCE_RANGE),
      m_PreferredDistanceRange( PREFERRED_DISTANCE_RANGE),
      m_TiltAngleRange( TILT_ANGLE_RANGE),
      m_MinimalFragmentInterfaceLength( MINIMAL_FRAGMENT_INTERFACE_LENGTH),
      m_SSTypes( SS_TYPES)
    {
    }

    //! @brief virtual copy constructor
    TypeData *TypeData::Clone() const
    {
      return new TypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &TypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_WindowRadii, ISTREAM);
      io::Serialize::Read( m_WindowLengths, ISTREAM);
      io::Serialize::Read( m_IsValid, ISTREAM);
      io::Serialize::Read( m_ResidueDistanceCutoff, ISTREAM);
      io::Serialize::Read( m_DistanceRange, ISTREAM);
      io::Serialize::Read( m_PreferredDistanceRange, ISTREAM);
      io::Serialize::Read( m_TiltAngleRange, ISTREAM);
      io::Serialize::Read( m_MinimalFragmentInterfaceLength, ISTREAM);
      io::Serialize::Read( m_SSTypes, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &TypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_WindowRadii, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WindowLengths, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IsValid, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ResidueDistanceCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DistanceRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PreferredDistanceRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TiltAngleRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinimalFragmentInterfaceLength, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SSTypes, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace contact
} // namespace bcl
