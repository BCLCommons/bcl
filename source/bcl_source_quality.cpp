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
#include "quality/bcl_quality.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace quality
} // namespace bcl
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
#include "quality/bcl_quality_dme.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math_comparisons.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> DME::s_Instance
    (
      GetObjectInstances().AddInstance( new DME())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DME::DME()
    {
    }

    //! @brief Clone function
    //! @return pointer to new DME
    DME *DME::Clone() const
    {
      return new DME( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &DME::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &DME::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Less;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates distance matrix error (DME) between COORDINATES and REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return distance matrix error (DME) between COORDINATES and REFERENCE_COORDINATES
    double DME::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      BCL_Assert
      (
        COORDINATES.GetSize() == REFERENCE_COORDINATES.GetSize(),
        "number of given coordinates are different: " +
        util::Format()( COORDINATES.GetSize()) + " != " + util::Format()( REFERENCE_COORDINATES.GetSize())
      );

      // calculate and store the difference of distance matrices
      const linal::Matrix< double> difference_matrix
      (
        coord::CalculateDifferenceDistanceMatrix( COORDINATES, REFERENCE_COORDINATES)
      );

      // initialize dme value to calculate
      double dme_squared( 0.0);

      // make sure the matrix is symmetrix
      BCL_Assert
      (
        difference_matrix.IsSquare(), "difference matrix is not square matrix"
      );

      // store number rows and columns
      const size_t number_rows( difference_matrix.GetNumberRows());

      // iterate over rows
      for( size_t i( 0); i < number_rows; ++i)
      {
        // iterate over the upper triangle
        for( size_t j( i + 1); j < number_rows; ++j)
        {
          // sum up the square of the distance
          dme_squared += math::Sqr( difference_matrix( i, j));
        }
      }
      // do the normalization
      dme_squared *= 2.0 / double( number_rows * ( number_rows - 1));

      // return square root of the calculated value
      return math::Sqrt( dme_squared);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DME::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DME::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace quality
} // namespace bcl
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
#include "quality/bcl_quality_dmf.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math_comparisons.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically
#include <numeric>

namespace bcl
{
  namespace quality
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> DMF::s_Instance
    (
      GetObjectInstances().AddInstance( new DMF( storage::Set< double>()))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a distance cutoff set
    //! @param DISTANCE_CUTOFFS set of distance cutoffs
    DMF::DMF( const storage::Set< double> &DISTANCE_CUTOFFS) :
      m_DistanceCutoffs( DISTANCE_CUTOFFS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new DMF
    DMF *DMF::Clone() const
    {
      return new DMF( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &DMF::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &DMF::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Greater;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates DMF between COORDINATES and REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return distance matrix fraction (DMF) between COORDINATES and REFERENCE_COORDINATES
    double DMF::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // make sure both coordinates are of same size
      BCL_Assert
      (
        COORDINATES.GetSize() == REFERENCE_COORDINATES.GetSize(),
        "The sizes of provided coordinates pair differ : " +
        util::Format()( COORDINATES.GetSize()) + " vs " + util::Format()( REFERENCE_COORDINATES.GetSize())
      );

      // calculate and store the difference of distance matrices
      const linal::Matrix< double> difference_matrix
      (
        coord::CalculateDifferenceDistanceMatrix( COORDINATES, REFERENCE_COORDINATES)
      );

      // dmf_sum
#define BCL_GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if defined(__APPLE__) && BCL_GCC_VERSION < BCL_GCC_VERSION <= 40204
      double dmf_sum( 0.0);
      for
      (
        storage::Set< double>::const_iterator
          cutoff_itr( m_DistanceCutoffs.Begin()), cutoff_itr_end( m_DistanceCutoffs.End());
        cutoff_itr != cutoff_itr_end; ++cutoff_itr
      )
      {
        // sum up the value
        dmf_sum += CalculateMeasureMatrixCutoff( difference_matrix, *cutoff_itr);
      }
#else
      std::list< double> results_each_cutoff;
      std::transform
      (
        m_DistanceCutoffs.Begin(), m_DistanceCutoffs.End(),
        std::inserter( results_each_cutoff, results_each_cutoff.begin()),
        std::bind1st( std::ptr_fun( &CalculateMeasureMatrixCutoff), difference_matrix)
      );
      const double dmf_sum( std::accumulate( results_each_cutoff.begin(), results_each_cutoff.end(), double( 0)));
#endif

      // end
      return dmf_sum / double( m_DistanceCutoffs.GetSize());
    }

    //! @brief calculates DMF between COORDINATES and REFERENCE_COORDINATES for a given distance cutoff
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @param DISTANCE_CUTOFF distance cutoff
    //! @return DMF between COORDINATES and REFERENCE_COORDINATES for a given distance cutoff
    double DMF::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES,
      const double DISTANCE_CUTOFF
    )
    {
      // make sure both coordinates are of same size
      BCL_Assert
      (
        COORDINATES.GetSize() == REFERENCE_COORDINATES.GetSize(),
        "The sizes of provided coordinates pair differ : " +
        util::Format()( COORDINATES.GetSize()) + " vs " + util::Format()( REFERENCE_COORDINATES.GetSize())
      );

      // calculate and store the difference of distance matrices
      const linal::Matrix< double> difference_matrix
      (
        coord::CalculateDifferenceDistanceMatrix( COORDINATES, REFERENCE_COORDINATES)
      );

      // calculate dmf and return it
      return CalculateMeasureMatrixCutoff( difference_matrix, DISTANCE_CUTOFF);
    }

    //! @brief calculates DMF for given difference matrix and a distance cutoff
    //! @param DIFFERENCE_MATRIX distance matrix that corresponds to two coordinate sets
    //! @param DISTANCE_CUTOFF  distance cutoff
    //! @return dmf for the given matrix
    double DMF::CalculateMeasureMatrixCutoff
    (
      const linal::Matrix< double> &DIFFERENCE_MATRIX,
      const double DISTANCE_CUTOFF
    )
    {
      // make sure the matrix is symmetrix
      BCL_Assert
      (
        DIFFERENCE_MATRIX.GetNumberRows() == DIFFERENCE_MATRIX.GetNumberCols(),
        "difference matrix is not symmetric: " + util::Format()( DIFFERENCE_MATRIX.GetNumberRows()) + " rows vs. " +
        util::Format()( DIFFERENCE_MATRIX.GetNumberCols()) + " cols!"
      );

      // initialize counts under the cutoff
      size_t counts( 0);

      // store number rows and columns
      const size_t number_rows( DIFFERENCE_MATRIX.GetNumberRows());

      // iterate over rows
      for( size_t i( 0); i < number_rows; ++i)
      {
        // iterate over the upper triangle
        for( size_t j( i + 1); j < number_rows; ++j)
        {
          // if the current cell has a value smaller than the cutoff
          if( DIFFERENCE_MATRIX( i, j) <= DISTANCE_CUTOFF)
          {
            // increment the count
            ++counts;
          }
        }
      }

      // do the normalization and return the value
      return double( 2.0) * DMF( storage::Set< double>()).OptimalValue() * double( counts) / double( number_rows * ( number_rows - 1));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DMF::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DistanceCutoffs, ISTREAM);
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DMF::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DistanceCutoffs, OSTREAM, INDENT);
      // end
      return OSTREAM;
    }

  } // namespace quality
} // namespace bcl
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
#include "quality/bcl_quality_gdt.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_range.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "quality/bcl_quality_average.h"
#include "quality/bcl_quality_lcs.h"
#include "quality/bcl_quality_rmsd.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
//    http://predictioncenter.org/casp/casp7/public/doc/LCS_GDT.README
//        LCS and GDT description
//
//    Longest Continuous Segments under specified CA RMSD cutoff (LCS).
//    The algorithm identifies all the longest continuous segments of residues
//    in the prediction deviating from the target by not more than specified
//    CA RMSD cutoff using many different superpositions.
//    Each residue in a prediction is assigned to the longest of such segments
//    provided if is a part of that segment. The absolutely longest continuous
//    segment in prediction under given RMSD cutoff is reported as well.
//    For different values of the CA RMSD cutoff (1.0 A, 2.0 A, and 5.0 A) the results
//    of the analysis are reported.
//    This measure can be used to evaluate ab initio 3D and comparative modeling
//    predictions.
//
//    Global Distance Test (GDT). The algorithm identifies in the prediction
//    the sets of residues deviating from the target by not more than
//    specified CA DISTANCE cutoff using many different superpositions.
//    Each residue in a prediction is assigned to the largest set of the residues
//    (not necessary continuous) deviating from the target by no more than a
//    specified distance cutoff.
//    This measure can be used to evaluate ab-initio 3D and comparative modeling
//    predictions.
//    For different values of DISTANCE cutoff (0.5 A, 1.0 A, 1.5 A, ... 10.0 A),
//    several measures are reported:
//    NUMBER_OF_CA_max  - the number of CA's from the "largest set" that
//           can fit under specified distance cutoff
//    PERCENT_OF_CA_Tg  - percent of CA's from the "largest set" comparing
//           to the total number of CA's in target
//    FRAGMENT: Beg-End - beginning and end of the segment containing the
//           "largest set" of CA's
//    RMS_LOCAL         - RMSD (root mean square deviation) calculated on the
//           "largest set" of CA's
//    RMS_ALL_CA        - RMSD calculated on all CA after superposition of
//           the prediction structure to the target structure
//           based on the "largest set" of CA's
//
//    The goal of introducing these two measures (GDT and LCS) is to provide a
//    tool that can be used for better detection of relatively good or bad parts
//    of the model.
//    - Using LCS we can localize the "best" continuous (along
//    the sequence) parts of the model that can fit under
//    RMSD thresholds: 1A, 2A, and 5A
//    Three blue lines represent the longest continuous sets
//    of residues that can fit under 1A, 2A, and 5A cutoff,
//    respectively.
//    - Using GDT we can localize the "best" sets of residues
//    (not necessary continuous) that can fit under DISTANCE
//    thresholds: 0.5A, 1.0A, 1.5A ,..., 10.0A
//    There are three blue lines on the GDT plot.
//    Each line represents the set of 5, 10, or 50 percent of
//    residues that can fit under specific distance cutoff (axis Y).
//    So, the lowest line represents residues (axis X) from the 5
//    percent sets of all target residues. Middle line identifies
//    those residues from the 10 percent sets, and highest from
//    50 percent sets.
//
//    The differences between LCS and GDT are the following:
//
//    1) LCS (Longest Continuous Segment) is based on RMSD cutoff.
//    2) The goal of LCS is to localize the longest continuous segment
//    of residues that can fit under RMSD cutoff.
//    3) Each residue in a prediction is assigned to the longest continuous
//    segment provided if is a part of that segment.
//    4) The data provided in the result files contains the LCS calculated
//    under three selected values of CA RMSD cutoff: 1A, 2A, and 5A
//
//    5) GDT (Global Distance Test) is based on the DISTANCE cutoff.
//    6) The goal of GDT is to localize the largest set of residues
//    (not necessary continuous) deviating from the target by no more than
//    a specified DISTANCE cutoff.
//    7) Each residue in a prediction is assigned to the largest set of the
//    residues provided if is a part of that set.
//    8) The data provided in the result files contains the GDT calculated under
//    several values of DISTANCE cutoff: 0.5, 1.0, 1.5, ... , 10.0 Angstroms.
//
//    Results of the analysis given by LCS algorithm show rather local features of
//    the model, while the residues considered in GDT come from the whole model
//    structure (they do not have to maintain the continuity along the sequence).
//
//    The GDT procedure is the following. Each three-residue segment and each
//    continuous segment found by LCS is used as a starting point to give an
//    initial equivalencies (model-target CA pairs) for a superposition.
//    The list of equivalencies is iteratively extended to produce the largest
//    set of residues that can fit under considered distance cutoff.
//    For collecting data about largest sets of residues the iterative
//    superposition procedure (ISP) is used.
//    The goal of the ISP method is to exclude from the calculations atoms
//    that are more than some threshold (cutoff) distance between the
//    model and the target structure after the transform is applied.
//    Starting from the initial set of atoms (C-alphas) the algorithm is the
//    following:
//    a) obtain the transform
//    b) apply the transform
//    c) identify all atom pairs for which distance is larger than the
//    threshold
//    d) re-obtain the transform, excluding those atoms
//    e) repeat b) - d) until the set of atoms used in calculations
//    is the same for two cycles running
//
//    -------

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> GDT::s_Instance
    (
      GetObjectInstances().AddInstance( new GDT( 0.0))
    );

    //! @brief returns the default seed length
    //! @return the default seed length
    const size_t GDT::GetDefaultSeedLength()
    {
      return 3;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from a single distance cutoff and seed length
    //! @param DISTANCE_CUTOFF distance cutoff to be used
    //! @param SEED_LENGTH length of seed
    GDT::GDT
    (
      const double &DISTANCE_CUTOFF,
      const size_t SEED_LENGTH
    ) :
      m_DistanceCutoff( DISTANCE_CUTOFF),
      m_SeedLength( GetDefaultSeedLength())
    {
    }

    //! @brief Clone function
    //! @return pointer to new GDT
    GDT *GDT::Clone() const
    {
      return new GDT( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &GDT::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &GDT::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Greater;
    }

    //! @brief get seed length
    //! @return seed length
    size_t GDT::GetSeedLength() const
    {
      return m_SeedLength;
    }

    //! @brief get distance cutoff
    //! @return distance cutoff
    double GDT::GetDistanceCutoff() const
    {
      return m_DistanceCutoff;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief collects the subset of coordinates specified by the given list of indices
    //! @param SUBSET indices of coordinates that define the subset of coordinates to be collected
    //! @param COORDINATES coordinates of interest
    //! @return vector of coordinates that correspond to the requested subset
    util::SiPtrVector< const linal::Vector3D> GDT::CollectCoordinatesSubset
    (
      const storage::List< size_t> &SUBSET,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES
    )
    {
      // initialize SiPtrVector to be returned
      util::SiPtrVector< const linal::Vector3D> coordinates_subset;
      coordinates_subset.AllocateMemory( SUBSET.GetSize());

      // iterate over the subset
      for
      (
        storage::List< size_t>::const_iterator index_itr( SUBSET.Begin()), index_itr_end( SUBSET.End());
        index_itr != index_itr_end; ++index_itr
      )
      {
        // insert the coordinates at the given index to the subset to be returned
        coordinates_subset.PushBack( COORDINATES( *index_itr));
      }

      // end
      return coordinates_subset;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculates GDT between COORDINATES and REFERENCE_COORDINATES for a given distance cutoff
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return GDT between COORDINATES and REFERENCE_COORDINATES for a given distance cutoff
    double GDT::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // return the GDT value
      return CalculateGDTAndSuperimposition( COORDINATES, REFERENCE_COORDINATES).First();
    }

    //! @brief calculates the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    math::TransformationMatrix3D GDT::CalculateSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // calculate GDT and return the superimposition
      return CalculateGDTAndSuperimposition( COORDINATES, REFERENCE_COORDINATES).Second();
    }

    //! @brief calculates GDT and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @param DISTANCE_CUTOFF distance cutoff for the superimposition
    //! @return pair of GDT and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    storage::Pair< double, math::TransformationMatrix3D> GDT::CalculateGDTAndSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // if the coordinates are empty are of different length
      if( COORDINATES.IsEmpty() || REFERENCE_COORDINATES.GetSize() != COORDINATES.GetSize())
      {
        return storage::Pair< double, math::TransformationMatrix3D>
        (
          util::GetUndefined< double>(), math::TransformationMatrix3D()
        );
      }

      // store number of points
      const size_t number_coordinates( COORDINATES.GetSize());

      // get the seed fragments
      storage::List< math::Range< size_t> > seed_fragments
      (
        GetSeedFragments( number_coordinates)
      );

      // get lcs as seeds also
      {
        const storage::List< math::Range< size_t> > lcs_seeds
        (
          LCS( m_DistanceCutoff).CalculateRanges( COORDINATES, REFERENCE_COORDINATES)
        );

        // if the lcs identified is larger than a normal seed add them
        if( !lcs_seeds.IsEmpty() && lcs_seeds.FirstElement().GetWidth() >= m_SeedLength)
        {
          seed_fragments.Append( lcs_seeds);
        }
      }

      // best indices and transformation
      storage::Pair< storage::List< size_t>, math::TransformationMatrix3D> best_indices_trans;

      // iterate over the seed fragments
      for
      (
        storage::List< math::Range< size_t> >::const_iterator itr( seed_fragments.Begin()), itr_end( seed_fragments.End());
           itr != itr_end
        && best_indices_trans.First().GetSize() != number_coordinates; // if the fragment includes all of the points meaning overall RMSD is less than the given cutoff
        ++itr
      )
      {
        // seed
        {
          const storage::Pair< storage::List< size_t>, math::TransformationMatrix3D> best_indices_trans_seed
          (
            IterativeSuperimpositionProcedure
            (
              COORDINATES, REFERENCE_COORDINATES,
              LCS::ConvertRangeToIndices( *itr)
            )
          );

          // if the longest fragment so far
          if( best_indices_trans_seed.First().GetSize() > best_indices_trans.First().GetSize())
          {
            // update it and best transformation
            best_indices_trans = best_indices_trans_seed;
          }
        }
      }

      // the longest fragment is as following
      BCL_MessageDbg
      (
        " the longest fragment and transformation is as following\n" + util::Format()( best_indices_trans)
      );

      // calculate the gdt value
      const double gdt( OptimalValue() * double( best_indices_trans.First().GetSize()) / double( number_coordinates));

      // return the pair of GDT value and the corresponding transformation matrix
      return storage::Pair< double, math::TransformationMatrix3D>( gdt, best_indices_trans.Second());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &GDT::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DistanceCutoff, ISTREAM);
      io::Serialize::Read( m_SeedLength    , ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &GDT::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DistanceCutoff, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_SeedLength    , OSTREAM, 0);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief create an average GDT object with given set of distance cutoffs
    //! @param DISTANCE_CUTOFFS set of distance cutoffs
    //! @param SEED_LENGTH length of seed
    //! @return Average object
    Average GDT::CreateAverageGDT
    (
      const storage::Set< double> &DISTANCE_CUTOFFS,
      const size_t SEED_LENGTH
    )
    {
      util::ShPtrVector< SuperimposeInterface> gdts;

      // for each distance cutoff given
      for( storage::Set< double>::const_iterator itr( DISTANCE_CUTOFFS.Begin()), itr_end( DISTANCE_CUTOFFS.End()); itr != itr_end; ++itr)
      {
        // construct a gdt measure
        gdts.PushBack( util::ShPtr< SuperimposeInterface>( new GDT( *itr, SEED_LENGTH)));
      }

      // end
      return Average( gdts);
    }

    //! @brief the iterative superimposition procedure (ISP)
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @param START_INDICES the selected coordinates as seed, do not have to fullfill the CUTOFF
    //! @return the newly selected coordinates and its transformation
    storage::Pair< storage::List< size_t>, math::TransformationMatrix3D> GDT::IterativeSuperimpositionProcedure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES,
      const storage::List< size_t> &START_INDICES
    ) const
    {
      // a)
      const math::TransformationMatrix3D start_transform( CalculateTransformation( COORDINATES, REFERENCE_COORDINATES, START_INDICES));
      // b), c)
      storage::Pair< storage::List< size_t>, math::TransformationMatrix3D> best_indices_transform
      (
        CollectIndicesBelowCutoff( COORDINATES, REFERENCE_COORDINATES, start_transform),
        math::TransformationMatrix3D()
      );
      storage::List< size_t> &best_indices( best_indices_transform.First());

      // d)
      math::TransformationMatrix3D &best_transform( best_indices_transform.Second());
      best_transform = CalculateTransformation( COORDINATES, REFERENCE_COORDINATES, best_indices_transform.First());

      // nr_extensions
      size_t nr_extension( 0);

      while( true)
      {
        ++nr_extension;
        // b), c)
        storage::List< size_t> current_indices
        (
          CollectIndicesBelowCutoff( COORDINATES, REFERENCE_COORDINATES, best_transform)
        );
        // the identified indices are as good or worse than the best
        if( current_indices.GetSize() <= best_indices.GetSize())
        {
          break;
        }
        // d)
        best_indices.InternalData().swap( current_indices.InternalData());
        best_transform = CalculateTransformation( COORDINATES, REFERENCE_COORDINATES, best_indices);

        BCL_MessageDbg
        (
          "\textension # " + util::Format()( nr_extension) +
          " fragment of size " + util::Format()( best_indices.GetSize()) +
          " start indices:\n" + util::Format()( START_INDICES)
        );
      }

      // end
      return best_indices_transform;
    }

    //! @brief returns all possible seed ranges (not necessarily valid)
    //! @param NUMBER_COORDIANTES number of coordinates
    //! @return continuous seed fragments
    storage::List< math::Range< size_t> > GDT::GetSeedFragments
    (
      const size_t NUMBER_COORDINATES
    ) const
    {
      BCL_Assert( m_SeedLength >= 3, "The seed length should be at least 3! not: " + util::Format()( m_SeedLength));

      // initialize seed fragments list to return
      storage::List< math::Range< size_t> > seed_fragments;

      // check that there are enough coordinates
      if( NUMBER_COORDINATES < m_SeedLength)
      {
        return seed_fragments;
      }

      // iterate over beginning indices
      for( size_t index( 0), index_end( NUMBER_COORDINATES + 1 - m_SeedLength); index != index_end; ++index)
      {
        // create this seed
        const math::Range< size_t> this_seed( index, index + m_SeedLength - 1);

        // add it to seeds list
        seed_fragments.PushBack( this_seed);
      }

      // end
      return seed_fragments;
    }

    //! @brief check if two coordinates and the transformation will bring them close below cutoff
    //! @param COORDINATE coord to transform
    //! @param REFERENCE_COORDINATE coord to compare against
    //! @param TRANSFORMATION transformation to apply
    //! @param DISTANCE_CUTOFF_SQUARED the distance below which they are considered good
    //! @return true if the distance between the transformed coord and the reference is below cutoff
    bool GDT::IsBelowCutoff
    (
      const linal::Vector3D &COORDINATE,
      const linal::Vector3D &REFERENCE_COORDINATE,
      const math::TransformationMatrix3D &TRANSFORMATION,
      const double DISTANCE_CUTOFF_SQUARED
    )
    {
      // make a copy of the coordinate and transform
      linal::Vector3D coordinate_copy( COORDINATE);
      coordinate_copy.Transform( TRANSFORMATION);

      // calculate the distance squared between ref coordinate transformed coordinate
      // compare to distance cutoff
      return linal::SquareDistance( REFERENCE_COORDINATE, coordinate_copy) <= DISTANCE_CUTOFF_SQUARED;
    }

    //! @brief collect all indices below the cutoff
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @param TRANSFORMATION transformation to apply before comparison
    //! @param list of coordinate indices that are below the cutoff for that transformation
    storage::List< size_t> GDT::CollectIndicesBelowCutoff
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES,
      const math::TransformationMatrix3D &TRANSFORMATION
    ) const
    {
      // store number of points and number of seed fragments
      const size_t number_points( COORDINATES.GetSize());

      // store allowed distance cutoff and its squared value
      const double distance_cutoff_squared( m_DistanceCutoff * m_DistanceCutoff);

      // create a extended subset to return
      storage::List< size_t> indices;

      // now iterate over all residues
      for( size_t index( 0); index < number_points; ++index)
      {
        if( IsBelowCutoff( *COORDINATES( index), *REFERENCE_COORDINATES( index), TRANSFORMATION, distance_cutoff_squared))
        {
          // insert it into the new subset
          indices.PushBack( index);
        }
      } // end index

      // end
      return indices;
    }

    //! @brief transformation for a new set of indices
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @param INDICES the selected coordinates
    //! @return the transformation matrix to superimpose the coordinates selected by the indices
    math::TransformationMatrix3D GDT::CalculateTransformation
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES,
      const storage::List< size_t> &INDICES
    )
    {
      return RMSD::SuperimposeCoordinates
      (
        CollectCoordinatesSubset( INDICES, REFERENCE_COORDINATES),
        CollectCoordinatesSubset( INDICES, COORDINATES)
      );
    }

  } // namespace quality
} // namespace bcl
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
#include "quality/bcl_quality_lcs.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_range.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "quality/bcl_quality_rmsd.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_si_ptr_vector.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> LCS::s_Instance
    (
      GetObjectInstances().AddInstance( new LCS())
    );

    //! @brief returns the default Rmsd cutoff
    //! @return the default Rmsd cutoff
    double LCS::GetDefaultRmsdCutoff()
    {
      return 5.0;
    }

    //! @brief returns the default seed length
    //! @return the default seed length
    size_t LCS::GetDefaultSeedLength()
    {
      return 3;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from a RMSD cutoff and a seed length
    //! @param RMSD_CUTOFF distance cutoff
    //! @param SEED_LENGTH length of seeds
    LCS::LCS
    (
      const double RMSD_CUTOFF,
      const size_t SEED_LENGTH
    ) :
      m_RmsdCutoff( RMSD_CUTOFF),
      m_SeedLength( SEED_LENGTH)
    {
      BCL_Assert( SEED_LENGTH >= 3, "The seed length has to be at least 3 not :" + util::Format()( SEED_LENGTH));
    }

    //! @brief Clone function
    //! @return pointer to new LCS
    LCS *LCS::Clone() const
    {
      return new LCS( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &LCS::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the optimal value for that quality measurement
    //! @return the best value by which two sets of coordinates can agree
    double LCS::OptimalValue() const
    {
      return util::GetUndefined< double>();
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &LCS::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Greater;
    }

    //! @brief return rmsd cutoff
    //! @return rmsd cutoff
    double LCS::GetCutoff() const
    {
      return m_RmsdCutoff;
    }

    //! @brief get seed length
    //! @return seed length
    size_t LCS::GetSeedLength() const
    {
      return m_SeedLength;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief find a larger range by extending the given one that has a RMSD below the cutoff
    //! @param RANGE range of coordinates that are used to be as seed and to be extended
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return extended range
    math::Range< size_t> LCS::ExtendRange
    (
      const math::Range< size_t> &RANGE,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // store number of points and number of seed fragments
      const size_t number_points( COORDINATES.GetSize());

      // initialize rmsd object
      const RMSD rmsd( true);

      // calculate range length
      const size_t initial_range_length( RANGE.GetWidth() + 1);

      // get the initial coordinates and allocate memory for expansion
      util::SiPtrVector< const linal::Vector3D>
        coords( COORDINATES.SubSiPtrVector( RANGE.GetMin(), initial_range_length));
      coords.AllocateMemory( number_points - RANGE.GetMin());
      util::SiPtrVector< const linal::Vector3D>
        reference_coords( REFERENCE_COORDINATES.SubSiPtrVector( RANGE.GetMin(), initial_range_length));
      reference_coords.AllocateMemory( number_points - RANGE.GetMin());

      // make a copy of the range
      math::Range< size_t> longest_range( RANGE);
      math::Range< size_t> current_range( RANGE);

      // while the range can still be extended
      while( current_range.GetMax() < number_points - 1)
      {
        // increment the current range length
        current_range.SetMax( current_range.GetMax() + 1);

        // add the new coordinates
        coords.PushBack( COORDINATES( current_range.GetMax()));
        reference_coords.PushBack( REFERENCE_COORDINATES( current_range.GetMax()));

        // if the current RMSD is smaller
        if( rmsd.CalculateMeasure( coords, reference_coords) < m_RmsdCutoff)
        {
          // update the longest range
          longest_range = current_range;
        }
      }

      // return the longest range found so far
      return longest_range;
    }

    //! @brief returns the ranges of longest continuous segments that can be superimposed below cutoff for given coordinates
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the range of longest continuous segment that can be superimposed below cutoff for given coordinates
    storage::List< math::Range< size_t> > LCS::CalculateRanges
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // initialize variable to store the longest fragments
      storage::List< math::Range< size_t> > longest_fragments;

      // number points
      const size_t number_points( COORDINATES.GetSize());

      size_t max_frag_length( number_points), min_frag_length( 3);
      static size_t s_last_super_imp_size( size_t( math::Sqrt( max_frag_length * min_frag_length) + 1));
      size_t test_frag_length
      (
        s_last_super_imp_size + 1 < max_frag_length
        ? s_last_super_imp_size
        : size_t( math::Sqrt( max_frag_length * min_frag_length) + 1)
      );
      while( max_frag_length > min_frag_length)
      {
        bool found_good_range( false);
        for
        (
          size_t index( 0), mx_index( number_points - test_frag_length + 1);
          index < mx_index;
          ++index
        )
        {
          // range
          const math::Range< size_t> range( index, index + test_frag_length - 1);
          // fragment should be superimposable below cutoff
          if( IsGoodRange( range, COORDINATES, REFERENCE_COORDINATES))
          {
            found_good_range = true;
            break;
          }
        }
        if( found_good_range)
        {
          min_frag_length = test_frag_length;
        }
        else
        {
          max_frag_length = test_frag_length - 1;
        }
        test_frag_length = size_t( math::Sqrt( max_frag_length * min_frag_length) + 1);
      }

      s_last_super_imp_size = max_frag_length;
      // iterate over starting ranges
      for( size_t index( 0), mx_index( number_points - max_frag_length + 1); index < mx_index; ++index)
      {
        // range
        const math::Range< size_t> range( index, index + max_frag_length - 1);
        // fragment should be superimposable below cutoff
        if( IsGoodRange( range, COORDINATES, REFERENCE_COORDINATES))
        {
          longest_fragments.PushBack( range);
        }
      }

      // if message level is verbose or less
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        size_t ctr( 1);
        // iterate over each range
        for
        (
          storage::List< math::Range< size_t> >::const_iterator
            range_itr( longest_fragments.Begin()), range_itr_end( longest_fragments.End());
          range_itr != range_itr_end;
          ++range_itr, ++ctr
        )
        {
          // get the coordinates
          const util::SiPtrVector< const linal::Vector3D>
            coords( COORDINATES.SubSiPtrVector( range_itr->GetMin(), range_itr->GetWidth() + 1));
          const util::SiPtrVector< const linal::Vector3D>
            reference_coords( REFERENCE_COORDINATES.SubSiPtrVector( range_itr->GetMin(), range_itr->GetWidth() + 1));
          // make a copy of COORDINATES
          const math::TransformationMatrix3D transformation
          (
            RMSD::SuperimposeCoordinates( reference_coords, coords)
          );
          // make a copy of the actual coordinates so that they can be transformed and create SiPtrVector
          storage::Vector< linal::Vector3D>
            copy_coords( util::ConvertToStorageVector< linal::Vector3D, const linal::Vector3D>( COORDINATES));
          util::SiPtrVector< linal::Vector3D>
            copy_sp_coords( util::ConvertToSiPtrVector< linal::Vector3D>( copy_coords));
          // transform all the coordinates with the calculated translation
          coord::TransformCoordinates( copy_sp_coords, transformation);
          const double local_rmsd( RMSD( true).CalculateMeasure( coords, reference_coords));
          const double global_rmsd( RMSD( false).CalculateMeasure( copy_sp_coords, reference_coords));

          // print out the values
          BCL_MessageDbg
          (
            "the longest fragment # " + util::Format()( ctr) + " => " + range_itr->GetString() +
            " local_rmsd:\t" + util::Format()( local_rmsd) + " global_rmsd:\t" + util::Format()( global_rmsd)
          );
        }
      }

      // end
      return longest_fragments;
    }

    //! @brief returns the indices to the coordinates of longest continuous segments that can be superimposed below cutoff for given coordinates
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the indices of longest continuous segments that can be superimposed below cutoff for given coordinates
    storage::List< storage::List< size_t> > LCS::CalculateIndices
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      return ConvertRangesToLists( CalculateRanges( COORDINATES, REFERENCE_COORDINATES));
    }

    //! @brief calculates LCS between given coordinates
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return LCS between COORDINATES and REFERENCE_COORDINATES
    double LCS::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // calculate the LCS and store the range for the first one
      const storage::List< math::Range< size_t> > longest_ranges( CalculateRanges( COORDINATES, REFERENCE_COORDINATES));
      if( longest_ranges.IsEmpty())
      {
        return 0.0;
      }

      // take first of the longest ranges, since all are equally long, they might just be at different ranges
      const math::Range< size_t> &lcs( longest_ranges.FirstElement());

      // return the length of the segment
      return lcs.GetWidth() + 1;
    }

    //! @brief calculates the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    math::TransformationMatrix3D LCS::CalculateSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // calculate the LCS and store the range for the first one
      const storage::List< math::Range< size_t> > longest_ranges( CalculateRanges( COORDINATES, REFERENCE_COORDINATES));
      if( longest_ranges.IsEmpty())
      {
        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      // take first of the longest ranges, since all are equally long, they might just be at different ranges
      const math::Range< size_t> &lcs( longest_ranges.FirstElement());

      // get the coordinates for the subset
      const util::SiPtrVector< const linal::Vector3D>
        coords( COORDINATES.SubSiPtrVector( lcs.GetMin(), lcs.GetWidth() + 1));
      const util::SiPtrVector< const linal::Vector3D>
        reference_coords( REFERENCE_COORDINATES.SubSiPtrVector( lcs.GetMin(), lcs.GetWidth() + 1));

      // calculate and return the transformation
      return RMSD::SuperimposeCoordinates( reference_coords, coords);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LCS::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_RmsdCutoff, ISTREAM);
      io::Serialize::Read( m_SeedLength, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LCS::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_RmsdCutoff, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_SeedLength, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief check if a given range superimposes below the cutoff
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return true, if coordinates within the given range are superimposable below the cutoff
    bool LCS::IsGoodRange
    (
      const math::Range< size_t> &RANGE,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // calculate range length
      const size_t initial_range_length( RANGE.GetWidth() + 1);

      // coordinates for that range
      const util::SiPtrVector< const linal::Vector3D>
        coords( COORDINATES.SubSiPtrVector( RANGE.GetMin(), initial_range_length));
      const util::SiPtrVector< const linal::Vector3D>
        reference_coords( REFERENCE_COORDINATES.SubSiPtrVector( RANGE.GetMin(), initial_range_length));

      // if the current RMSD is smaller
      return RMSD( true).CalculateMeasure( coords, reference_coords) < m_RmsdCutoff;
    }

    //! @brief converts a range to list of indices
    //! @param RANGE list of indices
    //! @return list that contains the indices in the range
    storage::List< size_t> LCS::ConvertRangeToIndices
    (
      const math::Range< size_t> &RANGE
    )
    {
      // initialize list to return
      storage::List< size_t> list;

      // iterate over the number
      for( size_t index( RANGE.GetMin()), index_end( RANGE.GetMax() + 1); index < index_end; ++index)
      {
        list.PushBack( index);
      }

      // end
      return list;
    }

    //! @brief converts a vector of ranges to list of list of indices
    //! @param RANGES vector of range of indices
    //! @return vector of lists that contains the indices in the range vector
    storage::List< storage::List< size_t> > LCS::ConvertRangesToLists
    (
      const storage::List< math::Range< size_t> > &RANGES
    )
    {
      // initialize list to return
      storage::List< storage::List< size_t> > list;

      // iterate over vector
      std::transform( RANGES.Begin(), RANGES.End(), std::inserter( list.InternalData(), list.Begin()), std::ptr_fun( &ConvertRangeToIndices));

      // end
      return list;
    }

  } // namespace quality
} // namespace bcl
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
#include "quality/bcl_quality_maxsub.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "quality/bcl_quality_rmsd.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MaxSub::s_Instance
    (
      GetObjectInstances().AddInstance( new MaxSub())
    );

    //! @brief returns default RMSD cutoff
    //! @return default RMSD cutoff
    double MaxSub::GetDefaultRMSDCutoff()
    {
      return 3.5;
    }

    //! @brief returns default length of the seed subset of coordinates
    //! @return default seed length of the seed subset of coordinates
    size_t MaxSub::GetDefaultSeedLength()
    {
      return 3;
    }

    //! @brief returns default number of iterations to extend an individual seed subset
    //! @return default number of iterations to extend an individual seed subset
    size_t MaxSub::GetDefaultNumberIterations()
    {
      return 4;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a RMSD cutoff, seed length and number of iterations
    //! @param RMSD_CUTOFF RMSD cutoff for the longest subset of superimposed coordinates
    //! @param SEED_LENGTH length of the seed subset of coordinates
    //! @param NUMBER_ITERATIONS number of iterations to extend an individual seed subset
    MaxSub::MaxSub
    (
      const double RMSD_CUTOFF,
      const size_t SEED_LENGTH,
      const size_t NUMBER_ITERATIONS
    ) :
      m_RMSDCutoff( RMSD_CUTOFF),
      m_SeedLength( SEED_LENGTH),
      m_NumberIterations( NUMBER_ITERATIONS)
    {
      // check the given seed length
      BCL_Assert
      (
        SEED_LENGTH >= 3,
        "The seed length has to be at least 3 for transformations to work; not " + util::Format()( SEED_LENGTH)
      );
    }

    //! @brief Clone function
    //! @return pointer to new MaxSub
    MaxSub *MaxSub::Clone() const
    {
      return new MaxSub( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MaxSub::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &MaxSub::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Greater;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief find a larger subset by extending the given one that has a RMSD below the cutoff
    //! @param SUBSET subset of coordinates that are used to be as seed and to be extended
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return transformation for the largest subset
    math::TransformationMatrix3D MaxSub::ExtendSubset
    (
      storage::List< size_t> &SUBSET,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // store number of points and number of seed fragments
      const size_t number_points( COORDINATES.GetSize());

      // make a copy of the seed subset
      const storage::List< size_t> seed_subset( SUBSET);

      // iterate over coordinates to extend the subset
      for( size_t nr_itr( 1); nr_itr <= m_NumberIterations; ++nr_itr)
      {
        // double check that there are at least 3 members in the subset
        if( SUBSET.GetSize() < 3)
        {
          BCL_MessageVrb
          (
            "the size of subset is smaller than 3: " + util::Format()( SUBSET.GetSize()) +
            " skipping more iterations"
          )
          return math::TransformationMatrix3D();
        }

        // store allowed distance cutoff and its squared value
        const double this_distance_cutoff( double( nr_itr) * m_RMSDCutoff / double( m_NumberIterations));
        const double distance_cutoff_squared( this_distance_cutoff * this_distance_cutoff);

        // superimpose the coordinates for the subset of coordinates
        const math::TransformationMatrix3D transformation
        (
          RMSD::SuperimposeCoordinates
          (
            CollectCoordinatesSubset( SUBSET, REFERENCE_COORDINATES),
            CollectCoordinatesSubset( SUBSET, COORDINATES)
          )
        );

        // if the transformation was not defined
        if( !transformation.IsDefined())
        {
          BCL_MessageVrb
          (
            "undefined transformation in Maxsub calculation at iteration " + util::Format()( nr_itr) +
            " skipping more iterations"
          )
          return transformation;
        }

        // initialize the extended subset to be returned
        storage::List< size_t> extended_subset;

        // now iterate over all residues
        for( size_t index( 0); index < number_points; ++index)
        {
          // if it's already
          if( index >= seed_subset.FirstElement() && index <= seed_subset.LastElement())
          {
            extended_subset.PushBack( index);
            continue;
          }

          // make a copy of the coordinates from COORDINATES_B
          linal::Vector3D coordinate_copy( *COORDINATES( index));

          // calculate the distance squared between the coordinate A and transformed coordinate B
          const double this_distance_squared
          (
            ( ( *REFERENCE_COORDINATES( index)) - coordinate_copy.Transform( transformation)).SquareNorm()
          );

          // if the distance between the transformed coordinates is smaller than the RMSD cutoff
          if( this_distance_squared < distance_cutoff_squared)
          {
            // insert it into the new subset
            extended_subset.PushBack( index);
          }
        } // end index

        // update the subset the extended subset
        SUBSET = extended_subset;

      } // end iterations

      // now recompute the transformation with the latest SUBSET
      // superimpose the coordinates for the subset of coordinates
      const math::TransformationMatrix3D final_transformation
      (
        RMSD::SuperimposeCoordinates
        (
          CollectCoordinatesSubset( SUBSET, REFERENCE_COORDINATES),
          CollectCoordinatesSubset( SUBSET, COORDINATES)
        )
      );

      // create a new subset
      storage::List< size_t> filtered_subset;

      // calculate the cutoff square
      const double rmsd_cutoff_squared( m_RMSDCutoff * m_RMSDCutoff);

      // iterate over the indices in SUBSET
      for
      (
        storage::List< size_t>::iterator index_itr( SUBSET.Begin()), index_itr_end( SUBSET.End());
        index_itr != index_itr_end;
        ++index_itr
      )
      {
        // make a copy of the coordinates from COORDINATES_B
        linal::Vector3D coordinate_copy( *COORDINATES( *index_itr));

        // calculate the distance squared between coordinate A and transformed coordinate B
        const double this_distance_squared
        (
          ( ( *REFERENCE_COORDINATES( *index_itr)) - coordinate_copy.Transform( final_transformation)).SquareNorm()
        );

        // if the distance is not greater than the cutoff
        if( this_distance_squared <= ( rmsd_cutoff_squared))
        {
          // add it to filtered subset
          filtered_subset.PushBack( *index_itr);
        }
      }

      // assign the filtered subset
      SUBSET = filtered_subset;

      // return final transformation
      return final_transformation;
    }

    //! @brief collects the subset of coordinates specified by the given list of indices
    //! @param SUBSET indices of coordinates that define the subset of coordinates to be collected
    //! @param COORDINATES util::SiPtrVector< Vector3D> which will be measured for agreement with COORDINATES_B
    //! @return vector of coordinates that correspond to the requested subset
    util::SiPtrVector< const linal::Vector3D> MaxSub::CollectCoordinatesSubset
    (
      const storage::List< size_t> &SUBSET,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES
    ) const
    {
      // initialize SiPtrVector to be returned
      util::SiPtrVector< const linal::Vector3D> coordinates_subset;
      coordinates_subset.AllocateMemory( SUBSET.GetSize());

      // iterate over the subset
      for
      (
        storage::List< size_t>::const_iterator index_itr( SUBSET.Begin()), index_itr_end( SUBSET.End());
        index_itr != index_itr_end;
        ++index_itr
      )
      {
        // insert the coordinates at the given index to the subset to be returned
        coordinates_subset.PushBack( COORDINATES( *index_itr));
      }

      // end
      return coordinates_subset;
    }

    //! @brief calculates MaxSub between given coordinates
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return MaxSub between given coordinates
    double MaxSub::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      return CalculateMaxSubAndSuperimposition( COORDINATES, REFERENCE_COORDINATES).First();
    }

    //! @brief calculates the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    math::TransformationMatrix3D MaxSub::CalculateSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      return CalculateMaxSubAndSuperimposition( COORDINATES, REFERENCE_COORDINATES).Second();
    }

    //! @brief calculates the MaxSub and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return pair of the MaxSub and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    storage::Pair< double, math::TransformationMatrix3D> MaxSub::CalculateMaxSubAndSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // make sure the sizes of the coordinates are the same
      if( COORDINATES.GetSize() != REFERENCE_COORDINATES.GetSize())
      {
        BCL_MessageStd( "The number of coordinates given to MaxSub calculation differ!");
        return
          storage::Pair< double, math::TransformationMatrix3D>
          (
            util::GetUndefinedDouble(), math::TransformationMatrix3D()
          );
      }

      // initialize the variable to hold the indices for the largest subset found so far
      storage::List< size_t> largest_subset;

      // initialize best transformation
      math::TransformationMatrix3D best_transformation;

      // store number of points and number of seed fragments
      const size_t number_points( COORDINATES.GetSize());
      const size_t number_seed_subsets( number_points - m_SeedLength + 1);

      // iterate over coordinates to create seed segments
      for( size_t index( 0); index < number_seed_subsets; ++index)
      {
        // initialize seed subset
        storage::List< size_t> seed_subset;

        // form the seed subset
        for( size_t seed_index( index); seed_index < index + m_SeedLength; ++seed_index)
        {
          seed_subset.PushBack( seed_index);
        }

        // extend the seed and store the transformation
        const math::TransformationMatrix3D transformation
        (
          ExtendSubset( seed_subset, COORDINATES, REFERENCE_COORDINATES)
        );

        // if the extended subset leads to size smaller than 3
        if( seed_subset.GetSize() < 3)
        {
          continue;
        }

        // if the length of the extended seed is larger then the longest seed so far
        if( seed_subset.GetSize() > largest_subset.GetSize())
        {
          // update this
          largest_subset = seed_subset;

          // also update the best transformation
          best_transformation = transformation;
        }
      }
      // calculate MaxSub value
      const double maxsub( OptimalValue() * double( largest_subset.GetSize()) / double( number_points));

      // return the MaxSub value and the corresponding transformation matrix
      return storage::Pair< double, math::TransformationMatrix3D>( maxsub, best_transformation);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MaxSub::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_RMSDCutoff      , ISTREAM);
      io::Serialize::Read( m_SeedLength      , ISTREAM);
      io::Serialize::Read( m_NumberIterations, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MaxSub::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_RMSDCutoff      , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SeedLength      , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumberIterations, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace quality
} // namespace bcl
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
#include "quality/bcl_quality_measures.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "quality/bcl_quality_average.h"
#include "quality/bcl_quality_dme.h"
#include "quality/bcl_quality_dmf.h"
#include "quality/bcl_quality_gdt.h"
#include "quality/bcl_quality_rmsd.h"
#include "quality/bcl_quality_superimpose_measures.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {

  //////////
  // data //
  //////////

    //! @brief returns the default distance cutoff set for HA high accuracy
    //! @return the default distance cutoff set HA
    const storage::Set< double> &Measures::GetDistanceCutoffsHA()
    {
      // initialize default static distance cutoff vector
      static const storage::Set< double> s_cutoffs( storage::Set< double>::Create( 0.5, 1.0, 2.0, 4.0));
      return s_cutoffs;
    }

    //! @brief returns the default distance cutoff set for TS (total score)
    //! @return the default distance cutoff set TS
    const storage::Set< double> &Measures::GetDistanceCutoffsTS()
    {
      // initialize default static distance cutoff vector
      static const storage::Set< double> s_cutoffs( storage::Set< double>::Create( 1.0, 2.0, 4.0, 8.0));
      return s_cutoffs;
    }

    //! @brief return command line flag for defining the quality measures to be calculated
    //! @return command line flag for defining the quality measures to be calculated
    util::ShPtr< command::FlagInterface> &Measures::GetFlagQualityMeasures()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "quality",
          "\tlist of quality measures to be calculated",
          command::Parameter
          (
            "quality_measure",
            "\tquality measure to be calculated",
            command::ParameterCheckEnumerate< Measures>()
          ),
          0,
          GetMeasures().GetEnumCount()
        )
      );

      // end
      return s_flag;
    }

    //! @brief function to return the list of quality measures in a set defined by the command line flag
    //! @return the list of quality measures in a set defined by the command line flag
    storage::Set< Measure> Measures::GetCommandLineQualityMeasures()
    {
      return GetFlagQualityMeasures()->GetObjectSet< Measure>();
    }

    //! @brief construct all Measures
    Measures::Measures() :
      e_RMSD(                   AddEnum( "RMSD"                  , *GetSuperimposeMeasures().e_RMSD)),
      e_RMSD_NoSuperimposition( AddEnum( "RMSD_NoSuperimposition", util::ShPtr< MeasureInterface>( new RMSD( false)))),
      e_RMSD_XYSuperimposition( AddEnum( "RMSD_XYSuperimposition", *GetSuperimposeMeasures().e_RMSD_XYSuperimposition)),
      e_DME(                    AddEnum( "DME"                   , util::ShPtr< MeasureInterface>( new DME()))),
      e_DMF_HA(                 AddEnum( "DMF_HA"                , util::ShPtr< MeasureInterface>( new DMF( GetDistanceCutoffsHA())))),
      e_DMF_TS(                 AddEnum( "DMF_TS"                , util::ShPtr< MeasureInterface>( new DMF( GetDistanceCutoffsTS())))),
      e_LCS(                    AddEnum( "LCS"                   , *GetSuperimposeMeasures().e_LCS)),
      e_GDT_HA(                 AddEnum( "GDT_HA"                , util::ShPtr< MeasureInterface>( GDT::CreateAverageGDT( GetDistanceCutoffsHA()).Clone()))),
      e_GDT_TS(                 AddEnum( "GDT_TS"                , util::ShPtr< MeasureInterface>( GDT::CreateAverageGDT( GetDistanceCutoffsTS()).Clone()))),
      e_GDT_1A(                 AddEnum( "GDT_1A"                , *GetSuperimposeMeasures().e_GDT_1A)),
      e_GDT_2A(                 AddEnum( "GDT_2A"                , *GetSuperimposeMeasures().e_GDT_2A)),
      e_GDT_4A(                 AddEnum( "GDT_4A"                , *GetSuperimposeMeasures().e_GDT_4A)),
      e_GDT_8A(                 AddEnum( "GDT_8A"                , *GetSuperimposeMeasures().e_GDT_8A)),
      e_MaxSub(                 AddEnum( "MaxSub"                , *GetSuperimposeMeasures().e_MaxSub)),
      e_Zero(                   AddEnum( "Zero"                  , *GetSuperimposeMeasures().e_NoSuperimpose))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Measures::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief function that returns the static instance of the Measures class
    //! @return the static instance of the Measures class
    Measures &GetMeasures()
    {
      return Measures::GetEnums();
    }

  } // namespace quality

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< quality::MeasureInterface>, quality::Measures>;

  } // namespace util
} // namespace bcl
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
#include "quality/bcl_quality_rmsd.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RMSD::s_Instance
    (
      GetObjectInstances().AddInstance( new RMSD())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a boolean to whether superimpose coordinates (set to true by default)
    //! @param SUPERIMPOSE_COORDINATES boolean to whether superimpose coordinates before calculating RMSD
    //! @param IGNORE_Z_COORDINATES boolean whether to superimpose using the Z-coordinates
    RMSD::RMSD( const bool SUPERIMPOSE_COORDINATES, const bool IGNORE_Z_COORDINATES) :
      m_SuperimposeCoordinates( SUPERIMPOSE_COORDINATES),
      m_IgnoreZCoordinates( IGNORE_Z_COORDINATES)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new RMSD
    RMSD *RMSD::Clone() const
    {
      return new RMSD( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &RMSD::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &RMSD::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Less;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief creates and returns a new coordinate set with all Z coordinates set to 0 for the given coordinate set
    //! @param COORDINATES vector of coordinates of interest
    //! @return new vector of coordinates with all Z coordinates set to 0
    storage::Vector< linal::Vector3D> RMSD::RemoveZCoordinates
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES
    ) const
    {
      // initialize vector
      storage::Vector< linal::Vector3D> xy_coordinates;
      xy_coordinates.AllocateMemory( COORDINATES.GetSize());

      // iterate over given coordinates
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          itr( COORDINATES.Begin()), itr_end( COORDINATES.End());
        itr != itr_end; ++itr
      )
      {
        // create new coordinate set without the Z coordinate and push it back
        xy_coordinates.PushBack( linal::Vector3D( ( *itr)->X(), ( *itr)->Y(), 0.0));
      }

      // end
      return xy_coordinates;
    }

    //! @brief calculates root mean square deviation between given coordinates
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return root mean square deviation between given coordinates
    double RMSD::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // if superimpose coordinates is set
      if( m_SuperimposeCoordinates)
      {
        // if Z coordinates are to be used for the superimposition
        if( !m_IgnoreZCoordinates)
        {
          return SuperimposedRMSD( COORDINATES, REFERENCE_COORDINATES);
        }
        // Z coordinates are not to be used for the superimposition
        else
        {
          // create vectors for the new coordinates
          storage::Vector< linal::Vector3D> xy_coordinates( RemoveZCoordinates( COORDINATES));
          storage::Vector< linal::Vector3D> reference_xy_coordinates( RemoveZCoordinates( REFERENCE_COORDINATES));

          // calculate the superimposition
          const math::TransformationMatrix3D transform
          (
            SuperimposeCoordinates
            (
              util::ConvertToConstSiPtrVector( reference_xy_coordinates),
              util::ConvertToConstSiPtrVector( xy_coordinates)
            )
          );

          // create a vector to hold the transformed coordinates
          storage::Vector< linal::Vector3D> new_coordinates;

          // apply the transformation to the b coordinates
          for
          (
            util::SiPtrVector< const linal::Vector3D>::const_iterator coord_itr( COORDINATES.Begin()),
              coord_itr_end( COORDINATES.End());
            coord_itr != coord_itr_end; ++coord_itr
          )
          {
            linal::Vector3D this_coordinate( **coord_itr);
            this_coordinate.Transform( transform);
            new_coordinates.PushBack( this_coordinate);
          }

          // calculate the RMSD and return
          return RMSD::RealSpaceRMSD( util::ConvertToConstSiPtrVector( new_coordinates), REFERENCE_COORDINATES);
        }
      }
      // don't superimpose coordinates
      else
      {
        return RMSD::RealSpaceRMSD( COORDINATES, REFERENCE_COORDINATES);
      }
    }

    //! @brief calculates the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    math::TransformationMatrix3D RMSD::CalculateSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // if Z coordinates are to be ignore for the superimposition
      if( m_IgnoreZCoordinates)
      {
        // make copies of the coordinates without Z coordinates
         storage::Vector< linal::Vector3D> xy_coordinates( RemoveZCoordinates( COORDINATES));
         storage::Vector< linal::Vector3D> reference_xy_coordinates( RemoveZCoordinates( REFERENCE_COORDINATES));

         // calculate the transformation matrix and return it
         return
           SuperimposeCoordinates
           (
             util::ConvertToConstSiPtrVector( reference_xy_coordinates),
             util::ConvertToConstSiPtrVector( xy_coordinates)
           );
      }
      // otherwise do regular superimposition
      else
      {
        return SuperimposeCoordinates( REFERENCE_COORDINATES, COORDINATES);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RMSD::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SuperimposeCoordinates, ISTREAM);
      io::Serialize::Read( m_IgnoreZCoordinates, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @param INDENT indentation
    //! @return ostream which was read from
    std::ostream &RMSD::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SuperimposeCoordinates, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IgnoreZCoordinates, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief determine the transformation matrix to optimally (lowest RMSD_ superimpose two set of coordinates)
    //! @param COORDINATES_A set of coordinates A
    //! @param COORDINATES_B set of coordinates B
    //! @return Transformation matrix that superimposes B onto A
    math::TransformationMatrix3D RMSD::SuperimposeCoordinates
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_A,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_B
    )
    {
      // check size of vectors
      if( COORDINATES_A.GetSize() != COORDINATES_B.GetSize())
      {
        BCL_MessageCrt
        (
          "number of points differs: " + util::Format()( COORDINATES_A.GetSize()) +
          " != " + util::Format()( COORDINATES_B.GetSize())
        );

        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      const size_t number( COORDINATES_A.GetSize());

      // check for minimal size
      if( number < 3)
      {
//        BCL_Message
//        (
//          util::Message::e_Critical,
//          "number points should be at least 3, otherwise ambiguous solution for super imposition"
//        );

        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      // calculate the center of mass
      const linal::Vector3D shift_a( coord::CenterOfMass( COORDINATES_A, false));
      const linal::Vector3D shift_b( coord::CenterOfMass( COORDINATES_B, false));

      // Calculate the covariance matrix
      linal::Matrix3x3< double> moment( BuildCovarianceMatrix( COORDINATES_A, COORDINATES_B, shift_a, shift_b));
      if( !moment.IsDefined())
      {
        BCL_MessageVrb( "covariance matrix is undefined");
        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      return CovarianceToTransformationMatrix( moment, shift_a, shift_b);
    }

    //! @brief calculate the real space rmsd of two sets of coordinates
    //! uses the coordinates as they are passed
    //! @param COORDINATES_A set of coordinates A
    //! @param COORDINATES_B set of coordinates B
    //! @return the rmsd between the passed coordinates
    double RMSD::RealSpaceRMSD
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_A,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_B
    )
    {
      // check size of vectors
      if( COORDINATES_A.GetSize() != COORDINATES_B.GetSize())
      {
        return util::GetUndefined< double>();
      }

      // compute RMSD
      double sum_square_rmsd( 0);
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          itr_a( COORDINATES_A.Begin()), itr_a_end( COORDINATES_A.End()),
          itr_b( COORDINATES_B.Begin()), itr_b_end( COORDINATES_B.End());
        itr_a != itr_a_end && itr_b != itr_b_end;
        ++itr_a, ++itr_b
      )
      {
        sum_square_rmsd += ( **itr_a - **itr_b).SquareNorm();
      }

      // end
      return math::Sqrt( sum_square_rmsd / COORDINATES_A.GetSize());
    }

    //! @brief calculate the real space rmsd for a set of coordinates - compares each coordinate with every other
    //! uses the coordinates as they are passed
    //! @param COORDINATES set of coordinates that will be compared with themselves
    //! @return first the rmsd between the passed coordinates and the standard deviation of the rmsd
    //!         the RunningAverageSD< double> has the mean distance and the standard deviation in the distances
    storage::Pair< storage::VectorND< 2, double>, math::RunningAverageSD< double> > RMSD::RealSpaceRMSDPairwise
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES
    )
    {
      // true if only one coordinate or less
      if( COORDINATES.GetSize() < 2)
      {
        return storage::Pair< storage::VectorND< 2, double>, math::RunningAverageSD< double> >
        (
          storage::VectorND< 2, double>( util::GetUndefinedDouble(), util::GetUndefinedDouble()),
          math::RunningAverageSD< double>()
        );
      }

      math::RunningAverageSD< double> mean_sd;

      // compute RMSD
      double sum_square_rmsd( 0);

      // iterate over coordinates
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          itr_a( COORDINATES.Begin()), itr_a_end( COORDINATES.End());
        itr_a != itr_a_end;
        ++itr_a
      )
      {
        // iterate over coordinates
        for
        (
          util::SiPtrVector< const linal::Vector3D>::const_iterator itr_b( itr_a + 1);
          itr_b != itr_a_end;
          ++itr_b
        )
        {
          const double distance( ( **itr_a - **itr_b).Norm());
          BCL_MessageDbg( " distance " + util::Format()( distance));
          sum_square_rmsd += math::Sqr( distance);
          mean_sd += distance;
        }
      }

      // number pairwise comparisons
      const size_t num_comparisons( COORDINATES.GetSize() * ( COORDINATES.GetSize() - 1) / 2);
      BCL_MessageDbg( " sum_square_rmsd " + util::Format()( sum_square_rmsd));
      BCL_MessageDbg( " num_comparisons " + util::Format()( num_comparisons));
      // calculate standard deviation of rmsd
      const double rmsd_std_dev
      (
        math::Sqrt( std::abs( sum_square_rmsd / num_comparisons - math::Sqr( mean_sd.GetAverage())))
      );

      // calculate the rmsd
      const double rmsd( math::Sqrt( sum_square_rmsd / num_comparisons));

      // end
      return storage::Pair< storage::VectorND< 2, double>, math::RunningAverageSD< double> >
      (
        storage::VectorND< 2, double>( rmsd, rmsd_std_dev), mean_sd
      );
    }

    //! @brief calculate the rmsd of two sets of coordinates if they are optimally
    //! they are never superimposed - the rmsd is derived from the eigenvalues of the SVD of the covariance matrix
    //! @param COORDINATES_A set of coordinates A
    //! @param COORDINATES_B set of coordinates B
    //! @return the rmsd of the coordinates
    double RMSD::SuperimposedRMSD
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_A,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_B
    )
    {
      if( COORDINATES_A.IsEmpty() || COORDINATES_B.IsEmpty())
      {
        return util::GetUndefined< double>();
      }

      double square_norm_centered_a( 0.0);
      double square_norm_centered_b( 0.0);

      // Calculate the covariance matrix
      linal::Matrix3x3< double> moment
      (
        BuildCovarianceMatrix
        (
          COORDINATES_A,
          COORDINATES_B,
          coord::CenterOfMass( COORDINATES_A, false),
          coord::CenterOfMass( COORDINATES_B, false),
          &square_norm_centered_a,
          &square_norm_centered_b
        )
      );

      if( !moment.IsDefined())
      {
        BCL_MessageCrt( "covariance matrix is undefined");

        return util::GetUndefinedDouble();
      }

      // determine sign of last element
      static const double s_chi_threshold( 1e-10);
      const int chi( moment.Determinant() < s_chi_threshold ? -1 : 1);

      moment = linal::MatrixTimesItselfTransposed( moment);
      // sort diagonal
      linal::Vector< double> eigenvalues( moment.EigenValues());
      // handle numerical issues that could cause one of the eigenvalues to become slightly less than 0; also take
      // square root
      for
      (
        linal::Vector< double>::iterator itr( eigenvalues.Begin()), itr_end( eigenvalues.End());
        itr != itr_end; ++itr
      )
      {
        *itr = math::Sqrt( std::max( *itr, double( 0.0)));
      }
      std::sort( eigenvalues.Begin(), eigenvalues.End());
      eigenvalues( 0) *= chi;

      // calculate the square deviation
      double square_deviation( 0.0);
      square_deviation += square_norm_centered_a;
      square_deviation += square_norm_centered_b;
      square_deviation -= 2 * eigenvalues.Sum();

      // root mean and return
      return math::Sqrt( std::max( square_deviation, double( 0)) / double( COORDINATES_A.GetSize()));
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief compute covariance matrix of two sets of coordinates COORDINATES_A on COORDINATES_B
    //! both coordinate sets are translated to the center of mass
    //! @param COORDINATES_A set of coordinates
    //! @param COORDINATES_B set of coordinates
    //! @param CENTER_A the center of COORDINATES_A
    //! @param CENTER_B the center of COORDINATES_B
    //! @param SQUARE_NORM_CENTERED_COORDINATES_A optional pointer to which the square norm of the centered coordinates a will be deposited
    //! @param SQUARE_NORM_CENTERED_COORDINATES_B optional pointer to which the square norm of the centered coordinates b will be deposited
    //! @return COORDINATES_A * COORDINATES_B
    linal::Matrix3x3< double> RMSD::BuildCovarianceMatrix
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_A,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_B,
      const linal::Vector3D &CENTER_A,
      const linal::Vector3D &CENTER_B,
      double *SQUARE_NORM_CENTERED_COORDINATES_A,
      double *SQUARE_NORM_CENTERED_COORDINATES_B
    )
    {
      // check the centers are defined
      if( !CENTER_A.IsDefined() || !CENTER_B.IsDefined())
      {
        BCL_MessageVrb
        (
          "given centers are undefined:\n" + util::Format()( CENTER_A) + "\n" + util::Format()( CENTER_B)
        );

        // return empty matrix
        return linal::Matrix3x3< double>( util::GetUndefined< double>());
      }
      BCL_Assert
      (
        COORDINATES_A.GetSize() == COORDINATES_B.GetSize(),
        "Non-equal sized sets of coordinates cannot be used to build a covariance matrix!"
      );

      const size_t number( COORDINATES_A.GetSize());

      // initialize 2 matrices with coordinates and 2 to keep the start coordinates
      linal::Matrix< double> coord_a( number, 3);
      linal::Matrix< double> coord_b( number, 3);

      // copy and shift coordinates

      for( size_t row( 0); row < number; ++row)
      {
        linal::VectorReference< double> row_a( coord_a.GetRow( row)), row_b( coord_b.GetRow( row));
        row_a.CopyValues( *COORDINATES_A( row));
        row_a -= CENTER_A;
        row_b.CopyValues( *COORDINATES_B( row));
        row_b -= CENTER_B;
      }

      // make cross moments matrix
      const linal::Matrix< double> covariance_matrix( linal::MatrixTransposeTimesMatrix( coord_a, coord_b));

      // set the argument centered coordinates if desired
      if( SQUARE_NORM_CENTERED_COORDINATES_A != NULL)
      {
        *SQUARE_NORM_CENTERED_COORDINATES_A = coord_a.AsVector().SquareNorm();
      }
      if( SQUARE_NORM_CENTERED_COORDINATES_B != NULL)
      {
        *SQUARE_NORM_CENTERED_COORDINATES_B = coord_b.AsVector().SquareNorm();
      }

      // end
      return covariance_matrix;
    }

    //! @brief Transformation matrix from Covariance matrix
    //! @param MOMENT covariance matrix
    //! @param CENTER_COORDINATES of coordinates
    //! @param CENTER_REFERENCE_COORDINATES center of reference coordinates
    //! @return transformation matrix
    math::TransformationMatrix3D RMSD::CovarianceToTransformationMatrix
    (
      const linal::Matrix3x3< double> &MOMENT,
      const linal::Vector3D &CENTER_COORDINATES,
      const linal::Vector3D &CENTER_REFERENCE_COORDINATES
    )
    {
      // diagonalization
      linal::Matrix3x3< double> rotate( MOMENT);
      linal::Matrix3x3< double> moment_t( MOMENT);
      moment_t.Transpose();
      rotate *= moment_t;

      // solve Eigensystem
      linal::Vector< double> eigenvalues( 3, 0.0);
      linal::Matrix3x3< double> eigenvectors;
      if( !rotate.EigenVectorsSymmetric( eigenvectors, eigenvalues))
      {
        BCL_MessageCrt( "Non-symmetric eigenmatrix! Check: " + util::Format()( eigenvectors) + "; should be A^T*A of " + util::Format()( MOMENT));
        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      if( !util::IsDefined( eigenvalues( 0)))
      {
        BCL_MessageCrt( "Undefined principle eigenvalue (degenerate solutions)!");
        return math::TransformationMatrix3D( util::UndefinedObject());
      } // error

      // sort eigenvectors by eigenvalues
      eigenvectors.Transpose();
      eigenvectors.SortRowsAndVector( eigenvalues);

      // check second eigenvalue
      if( eigenvalues( 1) <= 0.0 || eigenvalues( 0) <= 0.0)
      {
        BCL_MessageCrt( "Undefined secondary eigenvalue (degenerate solutions)!");
        return math::TransformationMatrix3D();
      } // error

      // build rotation matrix
      eigenvectors.ReplaceRow( 2, linal::CrossProduct( linal::Vector3D( eigenvectors[ 0]), linal::Vector3D( eigenvectors[ 1])));
      rotate = eigenvectors * MOMENT;

      std::transform( eigenvalues.Begin(), eigenvalues.End(), eigenvalues.Begin(), &math::Sqrt< double>);
      std::transform( eigenvalues.Begin(), eigenvalues.End(), eigenvalues.Begin(), std::bind1st( std::divides< double>(), 1.0));
      for( double *ptr( rotate.Begin()), *ptr_e( eigenvalues.Begin()), *ptr_end( rotate.End()); ptr < ptr_end; ptr += rotate.GetNumberCols(), ++ptr_e)
      {
        std::transform( ptr, ptr + rotate.GetNumberCols(), ptr, std::bind2nd( std::multiplies< double>(), *ptr_e));
      }

      rotate.ReplaceRow( 2, linal::CrossProduct( linal::Vector3D( rotate[ 0]), linal::Vector3D( rotate[ 1])));
      rotate.Transpose();
      rotate *= eigenvectors;

      // shift and rotate molecule
      math::TransformationMatrix3D transform;
      transform( -CENTER_REFERENCE_COORDINATES);
      transform( math::RotationMatrix3D( rotate));
      transform( CENTER_COORDINATES);

      // return transformation matrix calculated
      return transform;
    }

  } // namespace quality
} // namespace bcl
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
#include "quality/bcl_quality_rmsd_preprocessor.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "quality/bcl_quality_rmsd.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RMSDPreprocessor::s_Instance
    (
      GetObjectInstances().AddInstance( new RMSDPreprocessor())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a boolean to whether superimpose coordinates (set to true by default)
    RMSDPreprocessor::RMSDPreprocessor( const util::SiPtrVector< const linal::Vector3D> &COORDINATES, const bool &RECENTER) :
      m_SquareCenteredNorm( 0.0),
      m_CoordinatesTransposed( size_t( 3), COORDINATES.GetSize()),
      m_Center( size_t( 3), float( 0.0)),
      m_RecenterForSuperimposition( RECENTER)
    {
      size_t coord_id( 0);
      for
      (
        auto itr_coords( COORDINATES.Begin()), itr_coords_end( COORDINATES.End());
        itr_coords != itr_coords_end;
        ++itr_coords, ++coord_id
      )
      {
        linal::Vector< float> coords( ( *itr_coords)->Begin(), ( *itr_coords)->End());
        m_Center += coords;
        m_CoordinatesTransposed( 0, coord_id) = coords( 0);
        m_CoordinatesTransposed( 1, coord_id) = coords( 1);
        m_CoordinatesTransposed( 2, coord_id) = coords( 2);
      }
      m_Center /= float( COORDINATES.GetSize());
      if( m_RecenterForSuperimposition)
      {
        auto x_trans( m_CoordinatesTransposed.GetRow( 0));
        auto y_trans( m_CoordinatesTransposed.GetRow( 1));
        auto z_trans( m_CoordinatesTransposed.GetRow( 2));
        x_trans -= m_Center( 0);
        y_trans -= m_Center( 1);
        z_trans -= m_Center( 2);
        m_SquareCenteredNorm = m_CoordinatesTransposed.AsVector().SquareNorm();
      }
    }

    //! @brief virtual copy constructor
    //! @return pointer to new RMSD
    RMSDPreprocessor *RMSDPreprocessor::Clone() const
    {
      return new RMSDPreprocessor( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &RMSDPreprocessor::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RMSDPreprocessor::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Center, ISTREAM);
      io::Serialize::Read( m_CoordinatesTransposed, ISTREAM);
      m_SquareCenteredNorm = m_CoordinatesTransposed.AsVector().SquareNorm();

      // end
      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @param INDENT indentation
    //! @return ostream which was read from
    std::ostream &RMSDPreprocessor::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Center, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CoordinatesTransposed, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief determine the transformation matrix to optimally (lowest RMSD_ superimpose two set of coordinates)
    //! @param COORDINATES_A set of coordinates A
    //! @param COORDINATES_B set of coordinates B
    //! @return Transformation matrix that superimposes B onto A
    math::TransformationMatrix3D RMSDPreprocessor::SuperimposeCoordinates
    (
      const RMSDPreprocessor &COORDINATES_B
    ) const
    {
      BCL_Assert
      (
        m_RecenterForSuperimposition && COORDINATES_B.m_RecenterForSuperimposition,
        "Tried to superimpose coordinates without recentering; indicates constructor was called improperly"
      );
      // check size of vectors
      if( m_CoordinatesTransposed.GetNumberCols() != COORDINATES_B.m_CoordinatesTransposed.GetNumberCols())
      {
        BCL_MessageCrt
        (
          "number of points differs: " + util::Format()( m_CoordinatesTransposed.GetNumberCols()) +
          " != " + util::Format()( COORDINATES_B.m_CoordinatesTransposed.GetNumberCols())
        );

        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      // check for minimal size
      if( m_CoordinatesTransposed.GetNumberCols() < 3)
      {
        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      // Calculate the covariance matrix
      linal::Matrix3x3< float> moment( BuildCovarianceMatrix( COORDINATES_B));
      if( !moment.IsDefined())
      {
        BCL_MessageVrb( "covariance matrix is undefined");
        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      linal::Matrix3x3< double> moment_dbl;
      std::copy( moment.Begin(), moment.End(), moment_dbl.Begin());
      linal::Vector3D center_dbl( m_Center.Begin(), m_Center.End());
      linal::Vector3D b_center_dbl( COORDINATES_B.m_Center.Begin(), COORDINATES_B.m_Center.End());
      return RMSD::CovarianceToTransformationMatrix( moment_dbl, center_dbl, b_center_dbl);
    }

    //! @brief calculate the rmsd of two sets of coordinates if they are optimally
    //! they are never superimposed - the rmsd is derived from the eigenvalues of the SVD of the covariance matrix
    //! @param COORDINATES_A set of coordinates A
    //! @param COORDINATES_B set of coordinates B
    //! @return the rmsd of the coordinates
    double RMSDPreprocessor::SuperimposedRMSD
    (
      const RMSDPreprocessor &COORDINATES_B
    ) const
    {
      BCL_Assert
      (
        m_RecenterForSuperimposition && COORDINATES_B.m_RecenterForSuperimposition,
        "Tried to superimpose coordinates without recentering; indicates constructor was called improperly"
      );
      if( !GetSize() || !COORDINATES_B.GetSize())
      {
        return util::GetUndefined< double>();
      }

      // Calculate the covariance matrix
      linal::Matrix3x3< float> moment( BuildCovarianceMatrix( COORDINATES_B));

      if( !moment.IsDefined())
      {
        BCL_MessageCrt( "covariance matrix is undefined");

        return util::GetUndefinedDouble();
      }

      // determine sign of last element
      static const float s_chi_threshold( 1e-10);
      const int chi( moment.Determinant() < s_chi_threshold ? -1 : 1);

      moment = linal::MatrixTimesItselfTransposed( moment);
      // sort diagonal
      linal::Vector< float> eigenvalues( moment.EigenValues());
      // handle numerical issues that could cause one of the eigenvalues to become slightly less than 0; also take
      // square root
      for
      (
        linal::Vector< float>::iterator itr( eigenvalues.Begin()), itr_end( eigenvalues.End());
        itr != itr_end; ++itr
      )
      {
        *itr = math::Sqrt( std::max( *itr, float( 0.0)));
      }
      std::sort( eigenvalues.Begin(), eigenvalues.End());
      eigenvalues( 0) *= chi;

      // calculate the square deviation
      float square_deviation( m_SquareCenteredNorm + COORDINATES_B.m_SquareCenteredNorm - 2.0 * eigenvalues.Sum());

      // root mean and return
      return math::Sqrt( std::max( square_deviation, float( 0)) / float( m_CoordinatesTransposed.GetNumberCols()));
    }

    //! @brief calculate the rmsd of two sets of coordinates if they are optimally
    //! they are never superimposed - the rmsd is derived from the eigenvalues of the SVD of the covariance matrix
    //! @param COORDINATES_A set of coordinates A
    //! @param COORDINATES_B set of coordinates B
    //! @return the rmsd of the coordinates
    double RMSDPreprocessor::RMSD( const RMSDPreprocessor &COORDINATES_B) const
    {
      BCL_Assert
      (
        !m_RecenterForSuperimposition && !COORDINATES_B.m_RecenterForSuperimposition,
        "Tried to superimpose coordinates but did recentering; indicates constructor was called improperly"
      );
      // check size of vectors
      if( GetSize() != COORDINATES_B.GetSize())
      {
        return util::GetUndefined< double>();
      }

      // compute RMSD
      double sum_square_rmsd( 0);
      for
      (
        auto itr_a( m_CoordinatesTransposed.Begin()), itr_a_end( m_CoordinatesTransposed.End()),
             itr_b( COORDINATES_B.m_CoordinatesTransposed.Begin());
        itr_a != itr_a_end;
        ++itr_a, ++itr_b
      )
      {
        sum_square_rmsd += math::Sqr( *itr_a - *itr_b);
      }

      // end
      return math::Sqrt( sum_square_rmsd / float( 3 * GetSize()));
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief compute covariance matrix of two sets of coordinates COORDINATES_A on COORDINATES_B
    //! both coordinate sets are translated to the center of mass
    linal::Matrix3x3< float> RMSDPreprocessor::BuildCovarianceMatrix( const RMSDPreprocessor &COORDINATES_B) const
    {
      // check the centers are defined
      if( !m_Center.IsDefined() || !COORDINATES_B.m_Center.IsDefined())
      {
        BCL_MessageVrb
        (
          "given centers are undefined:\n" + util::Format()( m_Center) + "\n" + util::Format()( COORDINATES_B.m_Center)
        );

        // return empty matrix
        return linal::Matrix3x3< float>( util::GetUndefined< float>());
      }
      BCL_Assert
      (
        m_CoordinatesTransposed.GetNumberCols() == COORDINATES_B.m_CoordinatesTransposed.GetNumberCols(),
        "Non-equal sized sets of coordinates cannot be used to build a covariance matrix!"
      );

      // make cross moments matrix
      const linal::Matrix< float> covariance_matrix
      (
        linal::MatrixTimesMatrixTranspose( m_CoordinatesTransposed, COORDINATES_B.m_CoordinatesTransposed)
      );

      // end
      return covariance_matrix;
    }

  } // namespace quality
} // namespace bcl
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
#include "quality/bcl_quality_superimpose_measures.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "quality/bcl_quality_const_measure.h"
#include "quality/bcl_quality_gdt.h"
#include "quality/bcl_quality_lcs.h"
#include "quality/bcl_quality_maxsub.h"
#include "quality/bcl_quality_rmsd.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
  //////////
  // data //
  //////////

    //! @brief return command line flag for defining the measures to calculate superimposition
    //! @return command line flag for defining the measures to calculate superimposition
    util::ShPtr< command::FlagInterface> &SuperimposeMeasures::GetFlagSuperimposeMeasure()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "superimpose",
          "\tFlag for defining the quality measure to use for superimposing models onto template/native model",
          command::Parameter
          (
            "superimpose_measure",
            "\tquality measure to use for superimposition",
            command::ParameterCheckEnumerate< SuperimposeMeasures>(),
            GetSuperimposeMeasures().e_NoSuperimpose.GetName()
          )
        )
      );

      // end
      return s_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct all Measures
    SuperimposeMeasures::SuperimposeMeasures() :
      e_RMSD(                   AddEnum( "RMSD",                   util::ShPtr< SuperimposeInterface>( new RMSD( true)))),
      e_RMSD_XYSuperimposition( AddEnum( "RMSD_XYSuperimposition", util::ShPtr< SuperimposeInterface>( new RMSD( true, true)))),
      e_GDT_1A(                 AddEnum( "GDT_1A",                 util::ShPtr< SuperimposeInterface>( new GDT( 1.0)))),
      e_GDT_2A(                 AddEnum( "GDT_2A",                 util::ShPtr< SuperimposeInterface>( new GDT( 2.0)))),
      e_GDT_4A(                 AddEnum( "GDT_4A",                 util::ShPtr< SuperimposeInterface>( new GDT( 4.0)))),
      e_GDT_8A(                 AddEnum( "GDT_8A",                 util::ShPtr< SuperimposeInterface>( new GDT( 8.0)))),
      e_LCS(                    AddEnum( "LCS",                    util::ShPtr< SuperimposeInterface>( new LCS()))),
      e_MaxSub(                 AddEnum( "MaxSub",                 util::ShPtr< SuperimposeInterface>( new MaxSub()))),
      e_NoSuperimpose(          AddEnum( "NoSuperimpose",          util::ShPtr< SuperimposeInterface>( new ConstMeasure())))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SuperimposeMeasures::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief function that returns static instance of SuperimposeMeasures
    //! @return static instance of SuperimposeMeasures
    SuperimposeMeasures &GetSuperimposeMeasures()
    {
      return SuperimposeMeasures::GetEnums();
    }

  } // namespace quality

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< quality::SuperimposeInterface>, quality::SuperimposeMeasures>;

  } // namespace util
} // namespace bcl
