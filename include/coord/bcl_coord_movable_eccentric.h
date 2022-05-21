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

#ifndef BCL_COORD_MOVABLE_ECCENTRIC_H_
#define BCL_COORD_MOVABLE_ECCENTRIC_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_coord_movable_interface.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MovableEccentric
    //! @brief A movable class with the possibility of defining an external hinge
    //! @details This class allows storing a MovableInterface ( in m_Data) while providing an external hinge (m_Hinge)
    //! which is also derived from MovableInterface
    //!
    //! @see @link example_coord_movable_eccentric.cpp @endlink
    //! @author karakam
    //! @date May 6, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MovableEccentric :
      public MovableInterface
    {

    private:

    //////////
    // data //
    //////////

      //! Stores the actual MovableInterface data
      util::SiPtr< MovableInterface> m_Data;

      //! Stores the hinge information
      util::SiPtr< const MovableInterface> m_Hinge;

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
      MovableEccentric() :
        m_Data(),
        m_Hinge()
      {
      }

      //! @brief constructor from a MovableInterface and a Hinge MovableInterface
      //! @param MOVABLE_DATA MovableInterface Data
      //! @param HINGE Hinge MovableInterface
      MovableEccentric
      (
        MovableInterface &MOVABLE_DATA,
        const MovableInterface &HINGE
      ) :
        m_Data( MOVABLE_DATA),
        m_Hinge( HINGE)
      {
      }

      //! @brief constructor from SiPtrs to a MovableInterface and a Hinge MovableInterface
      //! @param SP_MOVABLE_DATA SiPtr to MovableInterface Data
      //! @param SP_HINGE SiPtr to Hinge MovableInterface
      MovableEccentric
      (
        const util::SiPtr< MovableInterface> &SP_MOVABLE_DATA,
        const util::SiPtr< const MovableInterface> &SP_HINGE
      ) :
        m_Data( SP_MOVABLE_DATA),
        m_Hinge( SP_HINGE)
      {
      }

      //! @brief copy constructor
      //! @param MOVABLE_ECCENTRIC MovableEccentric to be copied
      MovableEccentric
      (
        const MovableEccentric &MOVABLE_ECCENTRIC
      ) :
        m_Data( MOVABLE_ECCENTRIC.m_Data),
        m_Hinge( MOVABLE_ECCENTRIC.m_Hinge)
      {
      }

      //! @brief Clone function
      //! @return pointer to new MovableEccentric
      virtual MovableEccentric *Clone() const
      {
        return new MovableEccentric( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const;

      //! @brief returns the movable data
      //! @return the movable data
      const util::SiPtr< MovableInterface> &GetData() const
      {
        return m_Data;
      }

      //! @brief returns the hinge
      const util::SiPtr< const MovableInterface> &GetHinge() const
      {
        return m_Hinge;
      }

      //! @brief returns the geometric center of the object
      //! @return the geometric center of the object
      linal::Vector3D GetCenter() const;

      //! @brief returns the orientation of the object
      //! @param AXIS The requested Axis for which orientation will be returned
      //! @return the orientation of the object
      linal::Vector3D GetAxis( const Axis &AXIS) const;

      //! @brief returns the orientation and Position as TransformationMatrix3D
      //! @return the orientation and Position as TransformationMatrix3D
      const math::TransformationMatrix3D GetOrientation() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief translate the object along a given TRANSLATION vector
      //! @param TRANSLATION Vector along which the translation will occur
      void Translate( const linal::Vector3D &TRANSLATION);

      //! @brief transform the object by a given TRANSFORMATIONMATRIX3D
      //! @param TRANSFORMATIONMATRIX3D TransformationMatrix3D that will be applied
      void Transform( const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D);

      //! @brief rotate the object by a given ROTATIONMATRIX3D
      //! @brief ROTATIONMATRIX3D RotationMatrix3D that defines the rotation to be applied
      void Rotate( const math::RotationMatrix3D &ROTATIONMATRIX3D);

    ///////////////
    // operators //
    ///////////////

      //! @brief operator =
      //! @param MOVABLE_ECCENTRIC MovableEccentric to be copied
      //! @return MovableEccentric object with members updated to MOVABLE_ECCENTRIC
      MovableEccentric &operator =( const MovableEccentric &MOVABLE_ECCENTRIC);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // MovableEccentric

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_MOVABLE_ECCENTRIC_H_
