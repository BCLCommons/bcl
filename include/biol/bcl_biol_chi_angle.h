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

#ifndef BCL_BIOL_CHI_ANGLE_H_
#define BCL_BIOL_CHI_ANGLE_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_angle.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ChiAngle
    //! @brief represents a chi angle which is most likely one of the dihedral angles of the side chain of a residue
    //! @details Chi angles in residues typically start counting closest to backbone and increase as the move out
    //!          along the atoms of the side chain away from the backbone.
    //!
    //! @see @link example_biol_chi_angle.cpp @endlink
    //! @author alexanns
    //! @date Aug 24, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ChiAngle :
      public util::ObjectInterface
    {

    public:

      //! enumerator for chi angles
      enum Chi
      {
        e_One,
        e_Two,
        e_Three,
        e_Four,
        e_Five,
        e_Undefined,
        s_NumberChi
      };

      //! @brief conversion to a string from a SequenceDirection
      //! @param SEQUENCE_DIRECTION the sequence direction to get a string for
      //! @return a string representing that sequence direction
      static const std::string &GetChiName( const Chi &CHI_ANGLE);

      //! @brief enum class wrapper for Type
      typedef util::WrapperEnum< Chi, &GetChiName, s_NumberChi> ChiEnum;

    private:

    //////////
    // data //
    //////////

      //! the chi this angle is for
      ChiEnum m_Chi;

      //! the angle of this chi
      double m_Angle;

      //! the unit type the angle is measured in
      math::Angle::UnitEnum m_Unit;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief gives the maximum magnitude that a chi angle could have
      //! @param UNIT the unit the max angle should be given in
      //! @return double which is the maximum magnitude that a chi angle could have
      static double GetMaxAngle( const math::Angle::Unit &UNIT);

      //! @brief gives the total number of angular units in a circle
      //! @param UNIT the unit the angular units are in
      //! @return double the total number of angular units in a circle
      static double GetCircle( const math::Angle::Unit &UNIT);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ChiAngle();

      //! @brief constructor from member variable parameters
      //! @param CHI the chi this will correspond to
      ChiAngle( const Chi &CHI);

      //! @brief constructor from member variable parameters
      //! @param CHI the chi this will correspond to
      //! @param ANGLE the angle of this chi
      //! @param UNIT the unit type the angle is measured in
      ChiAngle( const Chi &CHI, const double ANGLE, const math::Angle::Unit &UNIT);

      //! @brief Clone function
      //! @return pointer to new ChiAngle
      ChiAngle *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives the chi this corresponds to
      //! @return Chi which is the chi this corresponds to
      const Chi &GetChi() const;

      //! @brief gives the angular value of this angle
      //! @param ANGLE_UNIT the unit the angle should be given in
      //! @return double which is the angular value of this angle in the units desired according ANGLE_UNIT
      double GetAngle( const math::Angle::Unit &ANGLE_UNIT) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculates the absolute angular difference between this chi angle and a provided chi angle
      //! @param CHI_ANGLE the chi angle whose difference will be calculated from this
      //! @param ANGLE_UNIT the unit the angular difference should be given in
      //! @return double which is the absolute angular difference between this and the given ChiAngle in desired units
      double CalculateAngleDifference( const ChiAngle &CHI_ANGLE, const math::Angle::Unit &ANGLE_UNIT) const;

    ///////////////
    // operators //
    ///////////////

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
      //!        format is space separated
      //!        <chi enum> <angle value> <angle unit>
      //!        An example is
      //!        e_Two 90 degree
      //!        Another example is
      //!        e_Three 1.1 radian
      //! @return istream the ChiAngle was read from
      std::istream &ReadSimple( std::istream &ISTREAM);

      //! @brief gives description of this in format as read by ReadSimple
      //! @return std::string gives description of this in ReadSimple format
      std::string WriteSimple( const math::Angle::Unit &ANGLE_UNIT) const;

    }; // class ChiAngle

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_CHI_ANGLE_H_
