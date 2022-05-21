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
#include "assemble/bcl_assemble_protein_model.h"
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "restraint/bcl_restraint_mutate_data_set_pairwise_filter_triangulation.h"
#include "score/bcl_score_data_set_pairwise_coordinate_triangulation.h"

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
    const util::SiPtr< const util::ObjectInterface> MutateDataSetPairwiseFilterTriangulation::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDataSetPairwiseFilterTriangulation())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &MutateDataSetPairwiseFilterTriangulation::GetDefaultScheme()
    {
      static const std::string s_scheme( "triangulation");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterTriangulation::MutateDataSetPairwiseFilterTriangulation( const std::string &SCHEME) :
      m_RadiusCutoff(),
      m_Ensemble(),
      m_SortedData(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking members
    //! @param RADIUS_CUTOFF coordinates are considered far enough apart when above this distance
    //! @param ENSEMBLE ensemble for which exposures will be calculated
      //! @param SORTED_DATA data sorted in a desired order in which it will be filtered out
    //! @param SCHEME the scheme for this scoring function
    MutateDataSetPairwiseFilterTriangulation::MutateDataSetPairwiseFilterTriangulation
    (
      const double &RADIUS_CUTOFF,
      const util::ShPtr< assemble::ProteinEnsemble> &ENSEMBLE,
      const storage::List< DataPairwise> &SORTED_DATA,
      const std::string &SCHEME
    ) :
      m_RadiusCutoff( RADIUS_CUTOFF),
      m_Ensemble( ENSEMBLE),
      m_SortedData( SORTED_DATA),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateDataSetPairwiseFilterTriangulation
    MutateDataSetPairwiseFilterTriangulation *MutateDataSetPairwiseFilterTriangulation::Clone() const
    {
      return new MutateDataSetPairwiseFilterTriangulation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateDataSetPairwiseFilterTriangulation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &MutateDataSetPairwiseFilterTriangulation::GetScheme() const
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
    //! @param DATA_SET DataSetPairwise of interest that will be mutated
    //! @return MutateResult that results from mutating to the argument DATA_SET
    math::MutateResult< DataSetPairwise>
    MutateDataSetPairwiseFilterTriangulation::operator()( const DataSetPairwise &DATA_SET) const
    {
      static util::ShPtr< DataSetPairwise> s_empty;

      if( DATA_SET.IsEmpty())
      {
        BCL_MessageStd( "MutateDataSetPairwiseFilterTriangulation data set is empty");
        return math::MutateResult< DataSetPairwise>( s_empty, *this);
      }

      // the dataset filtered
      util::ShPtr< DataSetPairwise> filtered_data( new DataSetPairwise());

      storage::Set
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > used_atoms;

      // iterate through the sorted dataset
      for
      (
        storage::List< DataPairwise>::const_iterator data_itr( m_SortedData.Begin()), data_itr_end( m_SortedData.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // try to find the data pair in DATA
        DataSetPairwise::const_iterator found_data_itr( DATA_SET.Find( *data_itr));

        // true if the current data is not found in DATA
        if( found_data_itr == DATA_SET.End())
        {
          // go to next data in sorted list of possible data
          continue;
        }
        BCL_MessageDbg
        (
          "checking " + data_itr->First()->GetIdentification()   + " and " + data_itr->Second()->GetIdentification()
        );

        // check if the data pair is usable considering the atoms that have been used already
        const bool usable_a_a
        (
          score::DataSetPairwiseCoordinateTriangulation::UsableLocator
          (
            used_atoms, *m_Ensemble, found_data_itr->First(), m_RadiusCutoff
          )
        );
        const bool usable_a_b
        (
          score::DataSetPairwiseCoordinateTriangulation::UsableLocator
          (
            used_atoms, *m_Ensemble, found_data_itr->Second(), m_RadiusCutoff
          )
        );

        // true if either aa or ab are not usable
        if( !usable_a_a || !usable_a_b)
        {
          // go to next iteration
          continue;
        }

        // add data pair to filtered_data
        filtered_data->Insert( *found_data_itr);

        // iterate through dataset again
        for
        (
          storage::List< DataPairwise>::const_iterator data_itr_b( ++storage::List< DataPairwise>::const_iterator( data_itr));
          data_itr_b != data_itr_end;
          ++data_itr_b
        )
        {
          // try to find the data pair in DATA
          DataSetPairwise::const_iterator found_data_itr_b( DATA_SET.Find( *data_itr_b));

          // true if the current data is not found in DATA
          if( found_data_itr_b == DATA_SET.End())
          {
            // go to next data in sorted list of possible data
            continue;
          }

          // check if the data pair is usable considering the atoms that have been used already
          const bool usable_b_a
          (
            score::DataSetPairwiseCoordinateTriangulation::UsableLocator
            (
              used_atoms, *m_Ensemble, found_data_itr_b->First(), m_RadiusCutoff
            )
          );
          const bool usable_b_b
          (
            score::DataSetPairwiseCoordinateTriangulation::UsableLocator
            (
              used_atoms, *m_Ensemble, found_data_itr_b->Second(), m_RadiusCutoff
            )
          );

          // true if either ba or bb are not usable
          if( !usable_b_a || !usable_b_b)
          {
            // go to next iteration
            continue;
          }

          // add data pair to filtered_data
          filtered_data->Insert( *found_data_itr_b);
        }
        BCL_MessageDbg
        (
          "finished checking " + data_itr->First()->GetIdentification() +
          " and " + data_itr->Second()->GetIdentification()
        );
      }

      BCL_MessageDbg
      (
        "MutateDataSetPairwiseFilterTriangulation data set size " +
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
    std::istream &MutateDataSetPairwiseFilterTriangulation::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_RadiusCutoff, ISTREAM);
      io::Serialize::Read( m_Ensemble, ISTREAM);
      io::Serialize::Read( m_SortedData, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDataSetPairwiseFilterTriangulation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_RadiusCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Ensemble, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SortedData, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
  
} // namespace bcl
