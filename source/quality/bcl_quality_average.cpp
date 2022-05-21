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
#include "quality/bcl_quality_average.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_statistics.h"
#include "math/bcl_math_transformation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Average::s_Instance
    (
      GetObjectInstances().AddInstance( new Average( util::ShPtrVector< SuperimposeInterface>()))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from a vector of measures
    //! @brief SUPERIMPOSE_MEASURES the measures to use - the first is used for non-average functions
    Average::Average( const util::ShPtrVector< SuperimposeInterface> &SUPERIMPOSE_MEASURES) :
      m_Measures( SUPERIMPOSE_MEASURES)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Average
    Average *Average::Clone() const
    {
      return new Average( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Average::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the optimal value for that quality measurement
    //! @return the best value by which two sets of coordinates can agree
    double Average::OptimalValue() const
    {
      return m_Measures.FirstElement()->OptimalValue();
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &Average::GetComparisonFunction() const
    {
      return m_Measures.FirstElement()->GetComparisonFunction();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates all measures for each measure between COORDINATES and REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return map of distance cutoff (key) and GDT value (value)
    storage::Vector< double> Average::CalculateMeasures
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      storage::Vector< double> measures;
      measures.AllocateMemory( m_Measures.GetSize());

      // calculate each individual measure
      for
      (
        util::ShPtrVector< SuperimposeInterface>::const_iterator itr( m_Measures.Begin()), itr_end( m_Measures.End());
        itr != itr_end;
        ++itr
      )
      {
        measures.PushBack( ( *itr)->CalculateMeasure( COORDINATES, REFERENCE_COORDINATES));
      }

      // end
      return measures;
    }

    //! @brief calculates GDT between COORDINATES and REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return GDT between COORDINATES and REFERENCE_COORDINATES
    double Average::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      const storage::Vector< double> measures( CalculateMeasures( COORDINATES, REFERENCE_COORDINATES));
      return math::Statistics::Sum( measures.Begin(), measures.End(), double( 0)) / measures.GetSize();
    }

    //! @brief calculates the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    math::TransformationMatrix3D Average::CalculateSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      return m_Measures.FirstElement()->CalculateSuperimposition( COORDINATES, REFERENCE_COORDINATES);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Average::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Measures, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Average::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Measures, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace quality
} // namespace bcl
