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

#ifndef BCL_SSPRED_BOCTOPUS_H_
#define BCL_SSPRED_BOCTOPUS_H_

// include the namespace header
#include "bcl_sspred.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sspred_method_interface.h"
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BOCTOPUS
    //! @brief stores prediction for BOCTOPUS
    //! @details stores prediction for BOCTOPUS, which is a TM-helix topology predictor. http://BOCTOPUS.cbr.su.se/
    //!
    //! @see @link example_sspred_boctopus.cpp @endlink
    //! @author mendenjl
    //! @date Aug 25, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BOCTOPUS :
      public MethodInterface
    {

    public:

      //! Boctopus prediction types
      enum Prediction
      {
        e_Inside    = int( 'i'),
        e_Membrane  = int( 'M'),
        e_Outside   = int( 'o')
      };

    private:

    //////////
    // data //
    //////////

      //! actual BOCTOPUS prediction
      Prediction m_Prediction;

    public:

      //! single instance of that class
      static const util::ObjectInterface *s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      BOCTOPUS();

      //! @brief constructor from the BOCTOPUS prediction
      //! @param PREDICTION prediction from BOCTOPUS file
      BOCTOPUS( const char &PREDICTION);

      //! @brief Clone function
      //! @return pointer to new Boctopus
      BOCTOPUS *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get file extension associated with this Method
      //! @return file extension associated with this Method
      const std::string &GetFileExtension() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get three state, environment independent secondary structure prediction
      //! @return three state, environment independent secondary structure prediction
      linal::Vector3D GetThreeStatePrediction() const;

      //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
      //! @return three state secondary structure prediction
      linal::Matrix< double> GetNineStatePrediction() const;

      //! @brief get the BOCTOPUS prediction
      Prediction GetPrediction() const;

      //! @brief find the TMTypes with highest prediction and returns it
      //! @return TMTYpe with highest prediction
      biol::EnvironmentType GetOneStateTMPrediction() const;

      //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
      //! @param ISTREAM input stream
      //! @param AMINO_ACID amino acid into which sspredictions will be read
      //! @return std::istream which was read from
      std::istream &ReadPredictionsForAA( std::istream &ISTREAM, biol::AABase &AMINO_ACID) const;

      //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
      //! @param ISTREAM input stream
      //! @param AA_SEQUENCE AASequence into which sspredictions will be read
      //! @return std::istream which was read from
      std::istream &ReadPredictionsForAASequence( std::istream &ISTREAM, biol::AASequence &AA_SEQUENCE) const;

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

    }; // class BOCTOPUS

  } // namespace sspred
} // namespace bcl

#endif // BCL_SSPRED_BOCTOPUS_H_
