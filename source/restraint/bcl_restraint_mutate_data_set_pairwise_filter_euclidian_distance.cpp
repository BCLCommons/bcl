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
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_euclidian_distance.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_data_set_pairwise.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateDataSetPairwiseFilterEuclidianDistance::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDataSetPairwiseFilterEuclidianDistance())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &MutateDataSetPairwiseFilterEuclidianDistance::GetDefaultScheme()
    {
      static const std::string s_scheme( "distance_range");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterEuclidianDistance::MutateDataSetPairwiseFilterEuclidianDistance
    (
      const std::string &SCHEME
    ) :
      m_DistanceRange(),
      m_Ensemble(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking optional scheme
    //! @param RANGE desired range of euclidian separation between data points
    //! @param ENSEMBLE ensemble distances will be calculated from
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterEuclidianDistance::MutateDataSetPairwiseFilterEuclidianDistance
    (
      const math::Range< double> &RANGE,
      const util::ShPtr< assemble::ProteinEnsemble> &ENSEMBLE,
      const std::string &SCHEME
    ) :
      m_DistanceRange( RANGE),
      m_Ensemble( ENSEMBLE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateDataSetPairwiseFilterEuclidianDistance
    MutateDataSetPairwiseFilterEuclidianDistance *MutateDataSetPairwiseFilterEuclidianDistance::Clone() const
    {
      return new MutateDataSetPairwiseFilterEuclidianDistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDataSetPairwiseFilterEuclidianDistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &MutateDataSetPairwiseFilterEuclidianDistance::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an DataSetPairwise and returning a mutated object
    //! @param DATA DataSetPairwise of interest that will be mutated
    //! @return MutateResult that results from mutating to the argument DATA_SET
    math::MutateResult< DataSetPairwise>
    MutateDataSetPairwiseFilterEuclidianDistance::operator()( const DataSetPairwise &DATA) const
    {

      static util::ShPtr< DataSetPairwise> s_empty;

      if( DATA.IsEmpty())
      {
        BCL_MessageDbg( "MutateDataSetPairwiseFilterEuclidianDistance data set is empty");
        return math::MutateResult< DataSetPairwise>( s_empty, *this);
      }

      // the data set filtered
      util::ShPtr< DataSetPairwise> filtered_data( new DataSetPairwise());

      // iterate through the data set to filter out data pairs that don't meet exposure criteria
      for
      (
        DataSetPairwise::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // get the average distance of the current data pair in m_Ensemble
        const math::RunningAverageSD< double> mean_sd( data_itr->EuclidianDistance( *m_Ensemble).First());

        // if there are counts get the mean, otherwise mean is undefined
        const double mean_distance( mean_sd.GetWeight() ? mean_sd.GetAverage() : util::GetUndefinedDouble());

        // true if size is within the desired range and defined
        if( m_DistanceRange.IsWithin( mean_distance) && util::IsDefined( mean_distance))
        {
          // add data pair to filtered_data
          filtered_data->Insert( *data_itr);
        }

      }

      BCL_MessageDbg
      (
        "MutateDataSetPairwiseFilterEuclidianDistance data set size " +
        util::Format()( filtered_data->GetSize())
      );

      return math::MutateResult< DataSetPairwise>( filtered_data, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDataSetPairwiseFilterEuclidianDistance::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DistanceRange, ISTREAM);
      io::Serialize::Read( m_Ensemble, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDataSetPairwiseFilterEuclidianDistance::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DistanceRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Ensemble, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
  
} // namespace bcl
