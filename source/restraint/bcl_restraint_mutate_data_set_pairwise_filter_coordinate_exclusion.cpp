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
#include "math/bcl_math_mutate_result.h"
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_coordinate_exclusion.h"
#include "score/bcl_score_data_set_pairwise_coordinate_exclusion.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateDataSetPairwiseFilterCoordinateExclusion::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDataSetPairwiseFilterCoordinateExclusion())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &MutateDataSetPairwiseFilterCoordinateExclusion::GetDefaultScheme()
    {
      static const std::string s_scheme( "filter_coord_exclusion");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateDataSetPairwiseFilterCoordinateExclusion::MutateDataSetPairwiseFilterCoordinateExclusion( const std::string &SCHEME) :
        m_Radius(),
        m_DistanceMap(),
        m_Scheme( SCHEME)
      {
      }

    //! @brief constructor taking members
    //! @param EXCLUSION_RADIUS data with atoms inside this radius of any exclusion coordinates will be penalized
    //! @param X_COORD_COLUMN column in istream that has the x-coordinate - columns start at 0
    //! @param Y_COORD_COLUMN column in istream that has the y-coordinate - columns start at 0
    //! @param Z_COORD_COLUMN column in istream that has the z-coordinate - columns start at 0
    //! @param ALL_POSSIBLE_DATA_POINTS the set of data points that should be subjected to this filter
    //! @param ENSEMBLE ensemble for which coordinates will be checked to make sure they aren't near exclusion coords
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterCoordinateExclusion::MutateDataSetPairwiseFilterCoordinateExclusion
    (
      const double &EXCLUSION_RADIUS,
      std::istream &READ,
      const size_t X_COORD_COLUMN,
      const size_t Y_COORD_COLUMN,
      const size_t Z_COORD_COLUMN,
      const storage::Set
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
        assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > &ALL_POSSIBLE_DATA_POINTS,
      const assemble::ProteinEnsemble &ENSEMBLE,
      const std::string &SCHEME
    ) :
      m_Radius( EXCLUSION_RADIUS),
      m_DistanceMap(),
      m_Scheme( SCHEME)
    {
      m_DistanceMap = score::DataSetPairwiseCoordinateExclusion::FillDistanceMap
        ( READ, X_COORD_COLUMN, Y_COORD_COLUMN, Z_COORD_COLUMN, ALL_POSSIBLE_DATA_POINTS, ENSEMBLE);
    }

    //! @brief Clone function
    //! @return pointer to new MutateDataSetPairwiseFilterCoordinateExclusion
    MutateDataSetPairwiseFilterCoordinateExclusion *MutateDataSetPairwiseFilterCoordinateExclusion::Clone() const
    {
      return new MutateDataSetPairwiseFilterCoordinateExclusion( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDataSetPairwiseFilterCoordinateExclusion::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &MutateDataSetPairwiseFilterCoordinateExclusion::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an DATA and returning a mutated object of DATA
    //! @param DATA DATA of interest that will be mutated
    //! @return MutateResult that results from mutating to the argument DATA
    math::MutateResult< DataSetPairwise>
    MutateDataSetPairwiseFilterCoordinateExclusion::operator()( const DataSetPairwise &DATA) const
    {
      static util::ShPtr< DataSetPairwise> s_empty;

      if( DATA.IsEmpty())
      {
        BCL_MessageStd( "MutateDataSetPairwiseFilterCoordinateExclusion data set is empty");
        return math::MutateResult< DataSetPairwise>( s_empty, *this);
      }

      // the dataset filtered
      util::ShPtr< DataSetPairwise> filtered_data( new DataSetPairwise());

      // iterate through the sorted dataset
      for
      (
        DataSetPairwise::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        const double score_first
        (
          score::DataSetPairwiseCoordinateExclusion::CalculateExclusionScore( data_itr->First(), m_DistanceMap, m_Radius)
        );
        const double score_second
        (
          score::DataSetPairwiseCoordinateExclusion::CalculateExclusionScore( data_itr->Second(), m_DistanceMap, m_Radius)
        );

        // true if both the scores are defined and zero - meaning they meet the radius cutoff around exclusion coords
        if
        (
          util::IsDefined( score_first) && util::IsDefined( score_second) &&
          score_first == double( 0) && score_first == score_second
        )
        {
          // add data pair to filtered_data
          filtered_data->Insert( *data_itr);
        }
      }

      return math::MutateResult< DataSetPairwise>( filtered_data, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDataSetPairwiseFilterCoordinateExclusion::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Radius, ISTREAM);
      io::Serialize::Read( m_DistanceMap, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDataSetPairwiseFilterCoordinateExclusion::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Radius, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DistanceMap, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
  
} // namespace bcl
