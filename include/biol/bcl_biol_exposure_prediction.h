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

#ifndef BCL_BIOL_EXPOSURE_PREDICTION_H_
#define BCL_BIOL_EXPOSURE_PREDICTION_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ExposurePrediction
    //! @brief Predicts exposure for residues using models trained with machine learning algorithms
    //! @details generates descriptors for a sequence, runs them through a trained model, and reports predicted exposure
    //!
    //! @see @link example_biol_exposure_prediction.cpp @endlink
    //! @author weinerbe, lib14
    //! @date Aug 27, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ExposurePrediction :
      public util::ObjectInterface
    {
    public:

        //! exposure types
        enum ExposureType
        {
          e_ProtomericContactNumber, // contact number of a residue when the protein is treated as a monomer
          e_OligomericContactNumber, // contact number of a residue when the protein is treated as a multimer
          e_ProtomericRSA, // relative solvent accessible surface area of a residue when the protein is treated as a monomer
          e_OligomericRSA, // relative solvent accessible surface area of a residue when the protein is treated as a multimer
          s_NumberExposureTypes
        };

        //! @brief return model path as a string
        //! @param ExposureType one of the exposure types
        //! @return model path as a string
        static const std::string &GetModelPath( const ExposureType &EXPOSURE_TYPE);

        //! @brief get file extension associated with the exposure type
        //! @param EXPOSURE_TYPE one of the exposure types
        //! @return file extension associated with this Method
        static const std::string &GetFileExtension( const ExposureType &EXPOSURE_TYPE);

        //! @brief OutputOptionEnum enum I/O helper
        typedef util::WrapperEnum< ExposureType, &GetModelPath, s_NumberExposureTypes> ExposureTypeEnum;

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ExposurePrediction();

      //! @brief Clone function
      //! @return pointer to new ExposurePrediction
      ExposurePrediction *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief iterates over the sequences and calculates the predictions for every residue in the sequence
      //! @param SEQUENCE sequence of interest
      //! @param EXPOSURE_TYPE one of the exposure types
      static void Calculate( AASequence &SEQUENCE, const ExposureType &EXPOSURE_TYPE = e_ProtomericContactNumber);

      //! @brief iterates over the sequences and calculates the predictions for every residue in the sequence
      //! @param MODEL protein model of interest
      //! @param EXPOSURE_TYPE one of the exposure types
      static void Calculate
      (
        assemble::ProteinModel &MODEL,
        const ExposureType &EXPOSURE_TYPE = e_ProtomericContactNumber
      );

      //! @brief iterates over the sequences and calculates the predictions for every residue in the sequence
      //! @param MODEL protein model of interest
      //! @param EXPOSURE_TYPE one of the exposure types
      static void Calculate
      (
        assemble::ProteinModelWithCache &MODEL,
        const ExposureType &EXPOSURE_TYPE = e_ProtomericContactNumber
      );

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief reads the exposure predictions for a model given a path and prefix
      //! @param PROTEIN_MODEL protein model to contain predictions
      //! @param PREFIX prefix of the sequence ( usually pdb id)
      //! @param PATH path where exposure prediction files can be found
      //! @param EXPOSURE_TYPE one of the exposure types
      static void ReadPredictions
      (
        assemble::ProteinModel &PROTEIN_MODEL,
        const std::string &PREFIX,
        const std::string &PATH,
        const ExposureType &EXPOSURE_TYPE = e_ProtomericContactNumber
      );

      //! @brief reads the exposure predictions from a file
      //! @param ISTREAM stream to read from
      //! @param PROTEIN_MODEL protein model to contain predictions
      static void ReadPredictions( std::istream &ISTREAM, assemble::ProteinModel &PROTEIN_MODEL);

      //! @brief reads the exposure predictions from a file
      //! @param ISTREAM stream to read from
      //! @param SEQUENCE sequence to contain predictions
      static void ReadPredictions( std::istream &ISTREAM, AASequence &SEQUENCE);

      //! @brief writes out the exposure predictions to a file
      //! @param OSTREAM stream to write to
      //! @param SEQUENCE sequence containing predictions
      static void WritePredictions( std::ostream &OSTREAM, const AASequence &SEQUENCE);

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

    }; // class ExposurePrediction

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_EXPOSURE_PREDICTION_H_
