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
#include "score/bcl_score_data_set_pairwise_coordinate_triangulation.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_data_set_pairwise.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DataSetPairwiseCoordinateTriangulation::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwiseCoordinateTriangulation())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &DataSetPairwiseCoordinateTriangulation::GetDefaultScheme()
    {
      static const std::string s_scheme( "triangulation");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking optional scheme
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseCoordinateTriangulation::DataSetPairwiseCoordinateTriangulation( const std::string &SCHEME) :
      m_RadiusCutoff(),
      m_Ensemble(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking members
    //! @param RADIUS_CUTOFF coordinates are considered far enough apart when above this distance
    //! @param ENSEMBLE ensemble for which exposures will be calculated
    //! @param SCHEME the scheme for this scoring function
    DataSetPairwiseCoordinateTriangulation::DataSetPairwiseCoordinateTriangulation
    (
      const double &RADIUS_CUTOFF,
      const util::ShPtr< assemble::ProteinEnsemble> &ENSEMBLE,
      const std::string &SCHEME
    ) :
      m_RadiusCutoff( RADIUS_CUTOFF),
      m_Ensemble( ENSEMBLE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new DataSetPairwiseCoordinateTriangulation
    DataSetPairwiseCoordinateTriangulation *DataSetPairwiseCoordinateTriangulation::Clone() const
    {
      return new DataSetPairwiseCoordinateTriangulation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwiseCoordinateTriangulation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DataSetPairwiseCoordinateTriangulation::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the score of a data set
    //! @param DATA data set to be scored
    //! @return the score of the current data set
    double DataSetPairwiseCoordinateTriangulation::operator()( const restraint::DataSetPairwise &DATA) const
    {
      // initialize score to zero
      double score( 0);

      // will hold the data points that have already been used to help determine which new ones should be added
      storage::Set
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > used_atoms;

      // iterate through the data set
      for
      (
        restraint::DataSetPairwise::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // check if the data pair is usable
        const bool usable_a_a( UsableLocator( used_atoms, *m_Ensemble, data_itr->First(), m_RadiusCutoff));
        const bool usable_a_b( UsableLocator( used_atoms, *m_Ensemble, data_itr->Second(), m_RadiusCutoff));

        // true if either aa or ab are not usable
        if( !usable_a_a || !usable_a_b)
        {
          // increment the score and go to next iteration
          ++score;
          continue;
        }

        // iterate through dataset again
        for
        (
          restraint::DataSetPairwise::const_iterator data_itr_b( ++restraint::DataSetPairwise::const_iterator( data_itr));
          data_itr_b != data_itr_end;
          ++data_itr_b
        )
        {
          // check if the data pair is usable
          const bool usable_b_a( UsableLocator( used_atoms, *m_Ensemble, data_itr_b->First(), m_RadiusCutoff));
          const bool usable_b_b( UsableLocator( used_atoms, *m_Ensemble, data_itr_b->Second(), m_RadiusCutoff));

          // true if either ba or bb are not usable
          if( !usable_b_a || !usable_b_b)
          {
            // increment the score and go to next iteration
            ++score;
            continue;
          }
        }
      }

      // return score
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetPairwiseCoordinateTriangulation::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_RadiusCutoff, ISTREAM);
      io::Serialize::Read( m_Ensemble, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataSetPairwiseCoordinateTriangulation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_RadiusCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Ensemble, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief determines if a locator can be used or if one is already similar to it so don't use it
    //! @param USED_LOCATORS locators that already exist in a data set
    //! @param ENSEMBLE structures that will be used to determine average spatial proximity
    //! @param LOCATOR locator which will be checked to see if it can be used or not based on already used locators
    //! @param DISTANCE_THRESHOLD the minimum spatial distance allowed between LOCATOR and any other used locator
    bool DataSetPairwiseCoordinateTriangulation::UsableLocator
    (
      storage::Set
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      > &USED_LOCATORS,
      const assemble::ProteinEnsemble &ENSEMBLE,
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR,
      const double DISTANCE_THRESHOLD
    )
    {
      // check to see if the locator_a is in used atoms
      storage::Set
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, assemble::LocatorAtomCoordinatesInterface::PtrLessThan
      >::const_iterator used_itr_a( USED_LOCATORS.Find( LOCATOR));

      // assume at first that LOCATOR can be used
      bool usable( true);

      // true if LOCATOR is not in used atoms
      if( used_itr_a == USED_LOCATORS.End())
      {
        // check to see if it is too close to one of the other atoms already in used_atoms
        // iterate through used atoms
        for
        (
          storage::Set
          <
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>,
            assemble::LocatorAtomCoordinatesInterface::PtrLessThan
          >::const_iterator
            used_itr( USED_LOCATORS.Begin()), used_itr_end( USED_LOCATORS.End());
          used_itr != used_itr_end;
          ++used_itr
        )
        {
          // calculate the mean and standard deviation of distances between LOCATOR and used_itr over models in ENSEMBLE
          const math::RunningAverageSD< double> mean_sd( ENSEMBLE.GetDistanceStatistics( restraint::DataPairwise( LOCATOR, *used_itr)));

          // get the mean distance between the used locator and the current locator
          const double distance_mean( mean_sd.GetWeight() > 0 ? mean_sd.GetAverage() : util::GetUndefinedDouble());

          // true if the two locators are on average too close or the distance mean is undefined
          if( distance_mean < DISTANCE_THRESHOLD || !util::IsDefined( distance_mean))
          {
            usable = false;
          }

          BCL_MessageDbg
          (
            "data pair distance " +
            restraint::DataPairwise( LOCATOR, *used_itr).GetIdentification() + " is " + util::Format()( distance_mean)
          );
          BCL_Assert( util::IsDefined( distance_mean) && !usable, "unusable but is defined");
        }

        // true if current locator is not close to any other used coordinates so far
        if( usable)
        {
          // insert the locator into the set of used data points
          BCL_Assert( USED_LOCATORS.Insert( LOCATOR).second, "could not insert " + LOCATOR->GetIdentification());
        }
      }

      // return bool indicating if the data point LOCATOR is usable or not
      return usable;
    }

  } // namespace score
  
} // namespace bcl
