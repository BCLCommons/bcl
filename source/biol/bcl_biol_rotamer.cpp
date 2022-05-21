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
#include "biol/bcl_biol_rotamer.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Rotamer::s_Instance
    (
      GetObjectInstances().AddInstance( new Rotamer())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Rotamer::Rotamer() :
      m_ChiAngles()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Rotamer
    Rotamer *Rotamer::Clone() const
    {
      return new Rotamer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Rotamer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief determines if rotamer contains any chi angles or not
    //! @return boolean - true if the rotamer contains no chi angles - false otherwise
    bool Rotamer::IsEmpty() const
    {
      return m_ChiAngles.IsEmpty();
    }

    //! @brief gives iterator to first chi in rotamer
    //! @return iterator to first chi in rotamer
    Rotamer::const_iterator Rotamer::Begin() const
    {
      return m_ChiAngles.Begin();
    }

    //! @brief gives iterator to end
    //! @return iterator to end
    Rotamer::const_iterator Rotamer::End() const
    {
      return m_ChiAngles.End();
    }

    //! @brief gives the number of chi angles in the rotamer
    //! @return the number of chi angles in the rotamer
    size_t Rotamer::GetSize() const
    {
      return m_ChiAngles.GetSize();
    }

    //! @brief gives the set of chis contained in this rotamer
    //! @return the set of chis contained in this rotamer
    storage::Set< ChiAngle::ChiEnum> Rotamer::GetChis() const
    {
      // to hold the chis contained in this rotamer
      storage::Set< ChiAngle::ChiEnum> chis;

      // iterate through the chi angles in order to fill the set of chis
      for
      (
        Rotamer::const_iterator chi_itr( Begin()), chi_itr_end( End());
        chi_itr != chi_itr_end;
        ++chi_itr
      )
      {
        chis.Insert( chi_itr->GetChi());
      }

      return chis;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief inserts a chi angle into the rotamer
    //! @param CHI_ANGLE the chi angle to add to the rotamer
    //! @return pair iterator to position of inserted element and bool indicating success or not
    std::pair< Rotamer::const_iterator, bool> Rotamer::Insert( const ChiAngle &CHI_ANGLE)
    {
      return m_ChiAngles.Insert( CHI_ANGLE);
    }

    //! @brief gives the angle value of the desired chi
    //! @param CHI the chi whose angle value is desired
    //! @param ANGLE_UNIT the unit the angle should be given in
    //! @return double which is the angle of the desired chi in desird units - undefined if chi does not exist
    double Rotamer::GetAngle( const ChiAngle::Chi &CHI, const math::Angle::Unit &ANGLE_UNIT) const
    {
      // try to find the chi angle with CHI
      const_iterator chi_itr( m_ChiAngles.Find( ChiAngle( CHI)));

      // true if chi could not be found
      if( chi_itr == m_ChiAngles.End())
      {
        return util::GetUndefinedDouble();
      }

      return chi_itr->GetAngle( ANGLE_UNIT);
    }

    //! @brief determines which chi angles have the same value between this and a given rotamer
    //!        in order for a chi to match, all previous chi also must match
    //! @param ROTAMER the other Rotamer whose chi angles will be compared to this
    //! @param ANGLE_UNIT the unit the angle tolerance is provided in
    //! @return set with all Chi that match dependent on the previous chi also being equivalent
    storage::Set< ChiAngle::ChiEnum> Rotamer::ChiMatchDependent
    (
      const Rotamer &ROTAMER, const math::Angle::Unit &ANGLE_UNIT, const double TOLERANCE
    ) const
    {
      storage::Set< ChiAngle::ChiEnum> matching_chi;

      for
      (
        Rotamer::const_iterator this_itr( Begin()), this_itr_end( End()),
          other_itr( ROTAMER.Begin()), other_itr_end( ROTAMER.End());
        this_itr != this_itr_end && other_itr != other_itr_end;
        ++this_itr, ++other_itr
      )
      {
        const double chi_angle_diff( this_itr->CalculateAngleDifference( *other_itr, ANGLE_UNIT));

        BCL_MessageDbg( "chi_angle_diff is " + util::Format()( chi_angle_diff));
        if( chi_angle_diff < TOLERANCE && util::IsDefined( chi_angle_diff))
        {
          matching_chi.Insert( this_itr->GetChi());
        }
        else
        {
          break;
        }
      }

      return matching_chi;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Rotamer::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ChiAngles, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Rotamer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ChiAngles, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief reads stream formatted in way easier for user to create
    //!        format is space separated on each line
    //!        <class identifier>
    //!        <chi enum> <angle value> <angle unit>
    //!        <chi enum> <angle value> <angle unit>
    //!        continued for as many chi angles as needed for this rotamer.
    //!        An example is
    //!        bcl::biol::Rotamer
    //!        e_Two 90 degree
    //!        e_Three 1.1 radian
    //! @return istream the rotamer was read from
    std::istream &Rotamer::ReadSimple( std::istream &ISTREAM)
    {
      BCL_MessageDbg( "reading in identifier");
      util::ObjectInterface::ReadIdentifier( ISTREAM);
      BCL_MessageDbg( "done reading in identifier");
      std::string line;
      while
      (
        !ISTREAM.eof() && line != GetClassIdentifier()
      )
      {
        std::getline( ISTREAM, line);
        util::TrimString( line);
        std::stringstream read( line);
        ChiAngle current_angle;
        if( !line.empty() && line != GetClassIdentifier())
        {
          current_angle.ReadSimple( read);
          BCL_MessageDbg( "current line is |" + line + "| and GetClassIdentifier is " + GetClassIdentifier());
          BCL_Assert
          (
            m_ChiAngles.Insert( current_angle).second, "could not insert chi angle " + util::Format()( current_angle) +
            "\ninto rotamer\n" + util::Format()( *this)
          );
        }
      }
        BCL_MessageDbg( "out of while loop current line is |" + line);

      return ISTREAM;
    }

    //! @brief gives description of this in format as read by ReadSimple
    //! @return std::string gives description of this in ReadSimple format
    std::string Rotamer::WriteSimple( const math::Angle::Unit &ANGLE_UNIT) const
    {
      std::string identification( "\n" + GetClassIdentifier());
      for
      (
        Rotamer::const_iterator itr( Begin()), itr_end( End());
        itr != itr_end;
        ++itr
      )
      {
        identification +=
        (
          "\n" + itr->WriteSimple( ANGLE_UNIT)
        );
      }

      return identification;
    }

    //! @brief binary functor helper struct for determing if one chi angle is less than another
    //!        sorts according to the ChiEnum within the chi angles
    //! @param CHI_ANGLE_A first chi angle object
    //! @param CHI_ANGLE_B second chi angle object
    //! @return bool true if the chi of CHI_ANGLE_A is less than CHI_ANGLE_B - false otherwise
    bool Rotamer::ChiAngleLessThan::operator()( const ChiAngle &CHI_ANGLE_A, const ChiAngle &CHI_ANGLE_B) const
    {
      return CHI_ANGLE_A.GetChi() < CHI_ANGLE_B.GetChi();
    }

  } // namespace biol
} // namespace bcl
