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

#ifndef BCL_QUALITY_CONST_MEASURE_H_
#define BCL_QUALITY_CONST_MEASURE_H_

// include the namespace header
#include "bcl_quality.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_quality_superimpose_interface.h"
#include "math/bcl_math_transformation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConstMeasure
    //! @brief Measure that returns the given constant quality value and superimposition matrix
    //! @details This class always returns the same quality value and superimposition matrix as given in the constructor
    //! This allows us to defined default enumerators in SuperimposeMeasures class
    //!
    //! @see @link example_quality_const_measure.cpp @endlink
    //! @author karakam
    //! @date Nov 9, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConstMeasure :
      public SuperimposeInterface
    {

    private:

    //////////
    // data //
    //////////

      //! value to be returned
      double m_Value;

      //! transformation to be returned
      math::TransformationMatrix3D m_Transformation;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief singleton to return the default quality value
      //! @return the default quality value
      static double GetDefaultValue();

      //! @brief singleton to return the default transformation matrix
      //! @return the default transformation matrix
      static const math::TransformationMatrix3D &GetDefaultTransformation();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ConstMeasure();

      //! @brief constructor from a value and transformation matrix
      //! @param VALUE value to be returned by default
      //! @param TRANSFORMATION TransformationMatrix3D to be returned by default
      ConstMeasure
      (
        const double VALUE,
        const math::TransformationMatrix3D &TRANSFORMATION
      );

      //! @brief Clone function
      //! @return pointer to new ConstMeasure
      ConstMeasure *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the optimal value for that quality measurement
      //! @return the best value by which two sets of coordinates can agree
      double OptimalValue() const;

      //! @brief return the comparison function for better quality
      //! @return binary function to compare two quality measure values
      const util::BinaryFunctionInterface< double, double, bool> &GetComparisonFunction() const;

      //! @brief get value
      //! @return value
      double GetValue() const
      {
        return m_Value;
      }

      //! @brief get transformation
      //! @return ransformationMatrix3D that is returned
      const math::TransformationMatrix3D &GetTransformation() const
      {
        return m_Transformation;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief calculates a quality measure between COORDINATES and REFERENCE_COORDINATES
      //! @param COORDINATES coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return a double which is the agreement of COORDINATES between REFERENCE_COORDINATES
      double CalculateMeasure
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const
      {
        return m_Value;
      }

      //! @brief calculates the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
      //! @param COORDINATES Vector of coordinates of interest
      //! @param REFERENCE_COORDINATES reference coordinates
      //! @return transformation matrix that superimposes COORDINATES to REFERENCE_COORDINATES
      math::TransformationMatrix3D CalculateSuperimposition
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
      ) const
      {
        return m_Transformation;
      }

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

    }; // class ConstMeasure

  } // namespace quality
} // namespace bcl

#endif // BCL_QUALITY_CONST_MEASURE_H_
