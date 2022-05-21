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

#ifndef BCL_BIOL_MEMBRANE_H_
#define BCL_BIOL_MEMBRANE_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_movable_interface.h"
#include "signal/bcl_signal_signal.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Membrane
    //! @brief This is a Membrane class.
    //! @details the membrane is located in the x-y plane, with a normal in z by default, but any other normal can be
    //! given. The membrane is movable, and environment types can be determined for passed coordinates
    //!
    //! @see @link example_biol_membrane.cpp @endlink
    //! @author woetzen, weinerbe
    //! @date 09/20/2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Membrane :
      public coord::MovableInterface
    {

    private:

    //////////
    // data //
    //////////

      math::TransformationMatrix3D m_Orientation; //!< orientation of the membrane
      storage::Vector< double> m_Thicknesses; //!< these are the z-thicknesses of the membrane in each region
      storage::Vector< double> m_Limits;      //!< these are the z-limits for regions in the membrane
      bool m_IsDefined; //!< True if the membrane is defined

      mutable signal::Signal1< const Membrane &> m_DestructorSignal; //!< signal handler for destructor

    public:

    //////////
    // data //
    //////////

      //! @brief return command line flag for defining the membrane thickness
      //! @return command line flag for defining the membrane thickness
      static util::ShPtr< command::FlagInterface> &GetFlagMembrane();

      //! @brief return command line parameter for defining the membrane core thickness
      //! @return command line parameter for defining the membrane core thickness
      static util::ShPtr< command::ParameterInterface> &GetParameterCoreThickness();

      //! @brief return command line parameter for defining the membrane transition thickness
      //! @return command line parameter for defining the membrane transition thickness
      static util::ShPtr< command::ParameterInterface> &GetParameterTransitionThickness();

      //! @brief return command line parameter for defining the membrane gap thickness
      //! @return command line parameter for defining the membrane gap thickness
      static util::ShPtr< command::ParameterInterface> &GetParameterGapThickness();

      //! @brief create a membrane object from commandline arguments
      //! @return membrane object created from commandline arguments
      static Membrane GetCommandLineMembrane();

      //! @brief static undefined membrane object
      //! @return undefined membrane object
      static const Membrane &GetUndefinedMembrane();

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Membrane();

      //! @brief constructor from membrane normal, all thicknesses and gap thickness
      //! @param THICKNESSES is a vector of membrane thicknesses
      //! @param CENTER center of the membrane
      //! @param NORMAL membrane normal
      explicit
      Membrane
      (
        const storage::Vector< double> &THICKNESSES,
        const linal::Vector3D &NORMAL = linal::Vector3D( 0.0, 0.0, 1.0),
        const linal::Vector3D &CENTER = linal::Vector3D()
      );

      //! @brief constructor from membrane normal, specified thicknesses and gap thickness
      //! @param THICKNESS_CORE thickness of membrane core
      //! @param THICKNESS_TRANSITION thickness of membrane transition region
      //! @param THICKNESS_GAP thickness of membrane gap
      //! @param NORMAL Membrane normal
      explicit
      Membrane
      (
        const double THICKNESS_CORE,
        const double THICKNESS_TRANSITION,
        const double THICKNESS_GAP,
        const linal::Vector3D &NORMAL = linal::Vector3D( 0.0, 0.0, 1.0)
      );

      //! @brief virtual copy constructor
      Membrane *Clone() const;

      //! @brief destructor
      virtual ~Membrane();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the membrane normal
      //! @return the membrane normal
      linal::Vector3D GetNormal() const;

      //! @brief return the thickness StorageVector
      //! @return the Thickness of the region
      const storage::Vector< double> &GetThicknesses() const
      {
        return m_Thicknesses;
      }

      //! @brief return the Thickness of the region
      //! @return the Thickness of the region
      double GetThickness( const EnvironmentType &ENVIRONMENT) const;

      //! @brief return the Limits StorageVector
      //! @return the Limits StorageVector
      const storage::Vector< double> &GetLimits() const
      {
        return m_Limits;
      }

      //! @brief return the Limit of the region
      //! @return the Limit of the region
      double GetLimit( const EnvironmentType &ENVIRONMENT) const;

      //! @brief access to the destructor signal handler
      //! @return ref to the SignalHandler that emits on destruction
      signal::Signal1< const Membrane &> &GetDestructorSignal() const
      {
        return m_DestructorSignal;
      }

      //! @brief returns true if the membrane is not set to static undefined membrane
      //! @return true if the membrane is not set to static undefined membrane
      bool IsDefined() const;

      //! @brief returns the geometric center of the object
      //! @return the geometric center of the object
      linal::Vector3D GetCenter() const;

      //! @brief return the orientation of the object
      //! @return orientation
      linal::Vector3D GetAxis( const coord::Axis &AXIS) const;

      //! @brief return the orientation and Position as TransformationMatrix3D
      //! @return TransformationMatrix3D that defines orientation and position
      const math::TransformationMatrix3D GetOrientation() const
      {
        return m_Orientation;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns EnvironmentType according to the coordinates
      //! @param COORDINATES coordinates of object
      //! @return EnvironmentType according to the coordinates
      EnvironmentType DetermineEnvironmentType( const linal::Vector3D &COORDINATES) const;

      //! @brief returns EnvironmentType according to the z-coordinate and weight, if it is in a gap region
      //! @param COORDINATES coordinates of object
      //! @return pair of EnvironmentType and a weight between 0 and one, 1 if it is closer to the innermost region,
      //! 0 if it is close to the outer most region
      storage::Pair< EnvironmentType, double> DetermineEnvironmentTypeAndWeight( const linal::Vector3D &COORDINATES) const;

      //! @brief  returns solvation energy for given z-coordinate and a vector containing three state solvation energies
      //! @param COORDINATES coordinates of object
      //! @return solvation energy for given z-coordinate and a vector containing three state solvation energies
      double CalculateSolvationEnergy( const linal::Vector3D &COORDINATES, const linal::Vector3D &TFE) const;

      //! @brief translate the object along a given TRANSLATION vector
      //! @param TRANSLATION Translation to be applied
      void Translate( const linal::Vector3D &TRANSLATION);

      //! @brief transform the object by a given TransformationMatrix3D
      //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
      void Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D);

      //! @brief rotate the object by a given RotationMatrix3D
      //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
      void Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D);

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
      //! @param OSTREAM output stream
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief create membrane object and transformation matrix from given pdbtm xml file
      //! only the core thickness can be retrieved from the xml file, so that transition region and gap are passed
      //! @param ISTREAM input stream of pdbtm xml file
      //! @param THICKNESS_TRANSITION thickness of the membrane transition region
      //! @param THICKNESS_GAP thickness of the gaps between the regions
      //! @return pair of membrane and transformation matrix
      static storage::Pair< Membrane, math::TransformationMatrix3D> MembraneAndTransformationFromPDBTMXML
      (
        std::istream &ISTREAM,
        const double THICKNESS_TRANSITION,
        const double THICKNESS_GAP
      );

      //! @brief returns pairs of chain id and transformation matrix that has to be applied to the chain to get one monomer for
      //! the biological relevant unit
      //! @param ISTREAM input stream of pdbtm xml file
      //! @param CHAIN_IDS chain ids in the original pdb
      //! @return vector of chain id the transformation is applied to, the new chain id and the matrix that needs to be applied order to get the relevant BIOMOLECULE
      storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> >
      static BioTransformationMatricesFromPDBTMXML( std::istream &ISTREAM, const std::string &CHAIN_IDS);

    private:

      //! @brief fills the Thicknesses with given values
      //! @param THICKNESS_CORE thickness of membrane core region
      //! @param THICKNESS_TRANSITION thickness of membrane transition region
      //! @param THICKNESS_GAP thickness of membrane gap region
      //! @return Thickness vector
      static storage::Vector< double> FillThicknessVector
      (
        const double THICKNESS_CORE,
        const double THICKNESS_TRANSITION,
        const double THICKNESS_GAP
      );

      //! @brief Calculates the Limits for the membrane regions
      //! @param THICKNESSES vector of thickness for membrane regions
      //! @return limits for the membrane regions
      static storage::Vector< double> FillLimitsVector( const storage::Vector< double> &THICKNESSES);

      //! @brief evaluates whether all entries in the vector are defined
      //! @param VECTOR vector to be evaulated
      //! @return whether all entries in the vector are defined
      static bool IsDefined( const storage::Vector< double> &VECTOR);

      //! @brief calculates the distance of a point from the center of the membrane plane
      //! @param COORDINATES coordinates to measure
      //! @return calculated distance
      double DistanceFromPlane( const linal::Vector3D &COORDINATES) const;

      //! @brief gets the membrane orientation from the given normal
      //! @param NORMAL membrane normal
      //! @param CENTER membrane center
      //! @return membrane orientation
      static math::TransformationMatrix3D OrientationFromNormal
      (
        const linal::Vector3D &NORMAL,
        const linal::Vector3D &CENTER = linal::Vector3D()
      );

    }; // class Membrane

  } // namespace biol
} // namespace bcl

#endif //BCL_BIOL_MEMBRANE_H_
