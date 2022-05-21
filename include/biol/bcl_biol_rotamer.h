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

#ifndef BCL_BIOL_ROTAMER_H_
#define BCL_BIOL_ROTAMER_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_biol_chi_angle.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Rotamer
    //! @brief represents a rotational conformer of a residue side chain
    //! @details stores a set of chi angles which can define the conformation of a residue side chain
    //!
    //! @see @link example_biol_rotamer.cpp @endlink
    //! @author alexanns
    //! @date Aug 24, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Rotamer :
      public util::ObjectInterface
    {

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class ChiAngleLessThan
      //! @brief binary functor helper struct for determing if one chi angle is less than another
      //!        sorts according to the ChiEnum within the chi angles
      //! @param CHI_ANGLE_A first chi angle object
      //! @param CHI_ANGLE_B second chi angle object
      //! @author alexanns
      //! @return bool true if the chi of CHI_ANGLE_A is less than CHI_ANGLE_B - false otherwise
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct ChiAngleLessThan
      {
        bool operator()( const ChiAngle &CHI_ANGLE_A, const ChiAngle &CHI_ANGLE_B) const;
      };

    private:

    //////////
    // data //
    //////////

      //! the chi angles that make up this rotamer
      storage::Set< ChiAngle, ChiAngleLessThan> m_ChiAngles;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      typedef storage::Set< ChiAngle, ChiAngleLessThan>::const_iterator const_iterator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Rotamer();

      //! @brief Clone function
      //! @return pointer to new Rotamer
      Rotamer *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief determines if rotamer contains any chi angles or not
      //! @return boolean - true if the rotamer contains no chi angles - false otherwise
      bool IsEmpty() const;

      //! @brief gives iterator to first chi in rotamer
      //! @return iterator to first chi in rotamer
      const_iterator Begin() const;

      //! @brief gives iterator to end
      //! @return iterator to end
      const_iterator End() const;

      //! @brief gives the number of chi angles in the rotamer
      //! @return the number of chi angles in the rotamer
      size_t GetSize() const;

      //! @brief gives the set of chis contained in this rotamer
      //! @return the set of chis contained in this rotamer
      storage::Set< ChiAngle::ChiEnum> GetChis() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief inserts a chi angle into the rotamer
      //! @param CHI_ANGLE the chi angle to add to the rotamer
      //! @return pair iterator to position of inserted element and bool indicating success or not
      std::pair< const_iterator, bool> Insert( const ChiAngle &CHI_ANGLE);

      //! @brief gives the angle value of the desired chi
      //! @param CHI the chi whose angle value is desired
      //! @param ANGLE_UNIT the unit the angle should be given in
      //! @return double which is the angle of the desired chi in desird units - undefined if chi does not exist
      double GetAngle( const ChiAngle::Chi &CHI, const math::Angle::Unit &ANGLE_UNIT) const;

      //! @brief determines which chi angles have the same value between this and a given rotamer
      //!        in order for a chi to match, all previous chi also must match
      //! @param ROTAMER the other Rotamer whose chi angles will be compared to this
      //! @param ANGLE_UNIT the unit the angle tolerance is provided in
      //! @return set with all Chi that match dependent on the previous chi also being equivalent
      storage::Set< ChiAngle::ChiEnum> ChiMatchDependent
      (
        const Rotamer &ROTAMER, const math::Angle::Unit &ANGLE_UNIT, const double TOLERANCE
      ) const;

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

    //////////////////////
    // helper functions //
    //////////////////////

    public:

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
      std::istream &ReadSimple( std::istream &ISTREAM);

      //! @brief gives description of this in format as read by ReadSimple
      //! @return std::string gives description of this in ReadSimple format
      std::string WriteSimple( const math::Angle::Unit &ANGLE_UNIT) const;

    }; // class Rotamer

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_ROTAMER_H_ 
