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
#include "quality/bcl_quality_const_measure.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_comparisons.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ConstMeasure::s_Instance
    (
      GetObjectInstances().AddInstance( new ConstMeasure())
    );

    //! @brief singleton to return the default quality value
    //! @return the default quality value
    double ConstMeasure::GetDefaultValue()
    {
      // return 0
      return 0.0;
    }

    //! @brief singleton to return the default transformation matrix
    //! @return the default transformation matrix
    const math::TransformationMatrix3D &ConstMeasure::GetDefaultTransformation()
    {
      // initialize a static const identity matrix
      static const math::TransformationMatrix3D s_default_transformation;

      // end
      return s_default_transformation;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ConstMeasure::ConstMeasure() :
      m_Value( GetDefaultValue()),
      m_Transformation( GetDefaultTransformation())
    {
    }

    //! @brief constructor from a value and transformation matrix
    //! @param VALUE value to be returned by default
    //! @param TRANSFORMATION TransformationMatrix3D to be returned by default
    ConstMeasure::ConstMeasure
    (
      const double VALUE,
      const math::TransformationMatrix3D &TRANSFORMATION
    ) :
      m_Value( VALUE),
      m_Transformation( TRANSFORMATION)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ConstMeasure
    ConstMeasure *ConstMeasure::Clone() const
    {
      return new ConstMeasure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ConstMeasure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the optimal value for that quality measurement
    //! @return the best value by which two sets of coordinates can agree
    double ConstMeasure::OptimalValue() const
    {
      return util::GetUndefined< double>();
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &ConstMeasure::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Less;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ConstMeasure::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Value, ISTREAM);
      io::Serialize::Read( m_Transformation, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ConstMeasure::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Value, OSTREAM, INDENT);
      io::Serialize::Write( m_Transformation, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace quality

} // namespace bcl
