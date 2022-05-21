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

#ifndef BCL_COORD_MOMENT_OF_INERTIA_H_
#define BCL_COORD_MOMENT_OF_INERTIA_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MomentOfInertia
    //! @brief class that calculates the moment of inertia tensor for coordinates with a given weight
    //! @see @link http://en.wikipedia.org/wiki/Principal_axis_%28mechanics%29#Moment_of_inertia_tensor @endlink
    //!
    //! @see @link example_coord_moment_of_inertia.cpp @endlink
    //! @author woetzen
    //! @date Apr 11, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MomentOfInertia :
      public math::FunctionInterfaceSerializable< linal::MatrixConstInterface< double>, linal::Vector3D>
    {

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MomentOfInertia();

      //! @brief Clone function
      //! @return pointer to new MomentOfInertia
      MomentOfInertia *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate transformation that translates into the center of weights and that sorts principal axes of inertia according to principal moments of inertia x - smallest, z - largest
      //! @param COORDINATES_WEIGHT_MATRIX matrix with coordinate and weight (4 cols) in number points rows
      //! @return transformation matrix and principal moments of inertias
      storage::Pair< math::TransformationMatrix3D, linal::Vector3D> TransformationAndMoments( const linal::MatrixConstInterface< double> &COORDINATES_WEIGHT_MATRIX) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that returns the ordered principal moments of inertia
      //! @param COORDINATE_WEIGHT_MATRIX 4*n matrix with three coordinates and one weight as ciolumns, and n rows
      //! @return ordered principal moments of inertia
      linal::Vector3D operator()( const linal::MatrixConstInterface< double> &COORDINATE_WEIGHT_MATRIX) const;

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

    private:

      //! @brief calculate the center of mass
      //! @param COORDINATES_WEIGHT_MATRIX matrix with coordinate and weight (4 columns) in number points rows
      //! @return the center of mass
      static linal::Vector3D CenterOfMass( const linal::MatrixConstInterface< double> &COORDINATES_WEIGHT_MATRIX);

    }; // class MomentOfInertia

  } // namespace coord
} // namespace bcl

#endif // BCL_COORD_MOMENT_OF_INERTIA_H_ 
