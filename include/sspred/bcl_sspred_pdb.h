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

#ifndef BCL_SSPRED_PDB_H_
#define BCL_SSPRED_PDB_H_

// include the namespace header
#include "bcl_sspred.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_sspred_method_interface.h"
#include "biol/bcl_biol_environment_types.h"
#include "biol/bcl_biol_ss_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PDB
    //! @brief MethodInterface derived class to store the SSType and the EnvironmentType of the native/template structure
    //! @details This MethodInterface derived class allows storing the SSType and the EnvironmentType as defined by
    //! the native/template structure at AAData level.
    //!
    //! @see @link example_sspred_pdb.cpp @endlink
    //! @author karakam
    //! @date Jan 14, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PDB :
      public MethodInterface
    {

    private:

    //////////
    // data //
    //////////

      //! native SSType defined by the PDB
      biol::SSType m_SSType;

      //! native environment type defined by the PDB
      biol::EnvironmentType m_EnvironmentType;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PDB();

      //! @brief constructor from a given SSType enum and an EnvironmentType enum
      //! @param SS_TYPE SSType enum of interest
      //! @param ENVIRONMENT_TYPE SSType EnvironmentType enum
      PDB( const biol::SSType &SS_TYPE, const biol::EnvironmentType &ENVIRONMENT_TYPE);

      //! @brief Clone function
      //! @return pointer to new PDB
      PDB *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get file extension associated with this Method
      //! @return file extension associated with this Method
      const std::string &GetFileExtension() const
      {
        // initialize static variable and return it
        static const std::string s_file_extension( ".pdb");
        return s_file_extension;
      }

      //! @brief get whether this method determined the secondary structure / membrane environment from the structure
      //! @return true if this method determined the secondary structure / membrane environment from the structure
      bool GetIsDeterminedFromSturcture() const
      {
        return true;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief get three state, environment independent secondary structure prediction
      //! @return three state, environment independent secondary structure prediction
      linal::Vector3D GetThreeStatePrediction() const;

      //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
      //! @return three state secondary structure prediction
      linal::Matrix< double> GetNineStatePrediction() const;

      //! @brief find the SSType with highest prediction and returns it
      //! @return SSType with highest prediction
      biol::SSType GetOneStateSSPrediction() const;

      //! @brief find the TMTypes with highest prediction and returns it
      //! @return TMTYpe with highest prediction
      biol::EnvironmentType GetOneStateTMPrediction() const;

      //! @brief find the SSType-TMType pair with highest prediction and returns it
      //! @return SSType-TMType pair with highest prediction
      storage::Pair< biol::SSType, biol::EnvironmentType> GetOneStateSSTMPrediction() const;

      //! @brief set the environment types for the residues in the protein model
      //! @param PROTEIN_MODEL model containing residues to be set
      //! @param USE_PDBTM_MEMBRANE_THICKNESS true to use membrane thickness from the PDBTM if it was read in
      //! @return true if the membrane environment was found; false otherwise
      static bool SetEnvironmentTypes
      (
        assemble::ProteinModel &PROTEIN_MODEL,
        const bool &USE_PDBTM_MEMBRANE_THICKNESS = false
      );

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
      //! @param ISTREAM input stream
      //! @param AMINO_ACID amino acid into which sspredictions will be read
      //! @return std::istream which was read from
      std::istream &ReadPredictionsForAA
      (
        std::istream &ISTREAM,
        biol::AABase &AMINO_ACID
      ) const;

      //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
      //! @param ISTREAM input stream
      //! @param AA_SEQUENCE AASequence into which sspredictions will be read
      //! @return std::istream which was read from
      std::istream &ReadPredictionsForAASequence
      (
        std::istream &ISTREAM,
        biol::AASequence &AA_SEQUENCE
      ) const;

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

    }; // class PDB

  } // namespace sspred
} // namespace bcl

#endif // BCL_SSPRED_PDB_H_
