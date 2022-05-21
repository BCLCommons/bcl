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
#include "restraint/bcl_restraint_epr_decay.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> EPRDecay::s_Instance( GetObjectInstances().AddInstance( new EPRDecay()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    EPRDecay::EPRDecay() :
      m_FirstResidueIDs(),
      m_SecondResidueIDs(),
      m_Measurements(),
      m_MeasurementsFilePath()
    {
    }

    //! @brief construct from file name
    //! @param MEASUREMENT_FILE_PATH path to the file containing the results of the EPR decay measurements
    EPRDecay::EPRDecay( const std::string &MEASUREMENT_FILE_PATH) :
      m_FirstResidueIDs(),
      m_SecondResidueIDs(),
      m_Measurements(),
      m_MeasurementsFilePath( MEASUREMENT_FILE_PATH)
    {
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief construct from members
    //! @param CHAIN_1 chain ID of the first residue
    //! @param SEQ_1 sequence ID of the first residue
    //! @param CHAIN_2 chain ID of the second residue
    //! @param SEQ_2 sequence ID of the second residue
    EPRDecay::EPRDecay( const char CHAIN_1, int SEQ_1, char CHAIN_2, int SEQ_2) :
      m_FirstResidueIDs( CHAIN_1, SEQ_1),
      m_SecondResidueIDs( CHAIN_2, SEQ_2),
      m_Measurements(),
      m_MeasurementsFilePath()
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new EPRDecay
    EPRDecay *EPRDecay::Clone() const
    {
      return new EPRDecay( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &EPRDecay::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &EPRDecay::GetAlias() const
    {
      static const std::string s_alias( "EPRDecay");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer EPRDecay::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Stores EPR decay measurements.");
      serializer.AddInitializer
      (
        "file path",
        "path to the file containing the results of the EPR decay measurements",
        io::Serialization::GetAgent( &m_MeasurementsFilePath)
      );

      return serializer;
    }

    //! @brief returns the spin-labeling sites these measurements are describing
    //! @return pair of spin-labeling sites described by chain ID and sequence ID
    EPRDecay::SLPair EPRDecay::GetSpinLabelingSites() const
    {
      return SLPair( m_FirstResidueIDs, m_SecondResidueIDs);
    }

    //! @brief adds a data point to the measurements
    //! @param TIME time of the measurement
    //! @DECAY observed decay
    void EPRDecay::AddMeasurement( double TIME, double DECAY)
    {
      m_Measurements.PushBack( Measurement( TIME, DECAY));
    }

    //! @brief returns the EPR decay measurements for this spin-labeling pair
    //! @return the EPR decay measurements for this spin-labeling pair
    const storage::Vector< EPRDecay::Measurement> &EPRDecay::GetMeasurements() const
    {
      return m_Measurements;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool EPRDecay::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
} // namespace bcl
