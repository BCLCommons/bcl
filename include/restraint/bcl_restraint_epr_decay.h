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

#ifndef BCL_RESTRAINT_EPR_DECAY_H_
#define BCL_RESTRAINT_EPR_DECAY_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EPRDecay
    //! @brief Stores observed EPR decay patterns for a spin-labeling pair
    //! @detail This class stores the the results of EPR decay measurements for a certain spin-labeling pair. The
    //! observed data consists of a point in time and the corresponding decay. This class is intended to store
    //! for only one spin-labeling pair.
    //!
    //! @see @link example_restraint_epr_decay.cpp @endlink
    //! @author fischea
    //! @date Nov 10, 2016
    //!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API EPRDecay :
      public util::SerializableInterface
    {

    //////////////
    // typedefs //
    //////////////

    public:

      //! type containing a decay measurement at a certain point in time
      typedef storage::Pair< double, double> Measurement;

      //! type describing one spin-labeling site
      typedef storage::Pair< storage::Pair< char, int>, storage::Pair< char, int> > SLPair;

    //////////
    // data //
    //////////

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      //! chain and sequence IDs of the first residue
      storage::Pair< char, int> m_FirstResidueIDs;

      //! chain and sequence IDs of the second residue
      storage::Pair< char, int> m_SecondResidueIDs;

      //! list of EPR decay measurements for the given spin-labeling pair
      storage::Vector< storage::Pair< double, double> > m_Measurements;

      //! path to the file containing the results of the EPR decay measurements
      std::string m_MeasurementsFilePath;

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

    public:

      //! @brief default constructor
      EPRDecay();

      //! @brief construct from file name
      //! @param MEASUREMENT_FILE_PATH path to the file containing the results of the EPR decay measurements
      EPRDecay( const std::string &MEASUREMENT_FILE_PATH);

      //! @brief construct from members
      //! @param CHAIN_1 chain ID of the first residue
      //! @param SEQ_1 sequence ID of the first residue
      //! @param CHAIN_2 chain ID of the second residue
      //! @param SEQ_2 sequence ID of the second residue
      EPRDecay( const char CHAIN_1, int SEQ_1, char CHAIN_2, int SEQ_2);

      //! @brief clone function
      //! @return pointer to a new EPRDecay
      EPRDecay *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief returns the spin-labeling sites these measurements are describing
      //! @return pair of spin-labeling sites described by chain ID and sequence ID
      SLPair GetSpinLabelingSites() const;

      //! @brief adds a data point to the measurements
      //! @param TIME time of the measurement
      //! @DECAY observed decay
      void AddMeasurement( double TIME, double DECAY);

      //! @brief returns the EPR decay measurements for this spin-labeling pair
      //! @return the EPR decay measurements for this spin-labeling pair
      const storage::Vector< Measurement> &GetMeasurements() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      //! @return return code indicating success or failure
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class EPRDecay

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_EPR_DECAY_H_
