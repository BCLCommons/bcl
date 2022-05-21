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

#ifndef BCL_SSPRED_MAHSSMI_H_
#define BCL_SSPRED_MAHSSMI_H_

// include the namespace header
#include "bcl_sspred.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_sspred_method_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Mahssmi
    //! @brief identifies secondary structure and membrane topology based on a hybrid method combining membrane location
    //!        from method e_PDB, and secondary structure from StrideDSSP (soluble regions) and PALSSE (membrane only)
    //! @details Membrane Aware Hybrid Secondary Structure & Membrane topology Identification
    //!
    //! @see @link example_sspred_mahssmi.cpp @endlink
    //! @author mendenjl
    //! @date Dec 05, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Mahssmi :
      public MethodInterface
    {

    private:

    //////////
    // data //
    //////////

      //! native SSType defined by the MAHSSMI algorithm
      biol::SSType m_SSType;

      //! native environment type defined by MAHSSMI
      biol::EnvironmentType m_EnvironmentType;

      //! for residues in a membrane beta barrel, true if the residues' CA-CB points towards the pore
      //! as opposed to membrane or neighboring member of an oligomer
      bool m_BetaBarrelResidueFacesPore;

      //! true if the residue is on a TM Helix/Strand that originates (has N-term) in the cytosol vs false for those w/ N-term
      //! in the extracellular space (for sequences). As per the OPM convention, cytosol is indicated by Z-coordinate < 0, so
      //! when computed from a 3D structure, TM-helices/strands that transit at least 4 angstroms and go positive will also
      //! be considered exiting cytosol
      bool m_OriginatesInCytosol;

    public:

      //! single instance of that class
      static const util::ObjectInterface *s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Mahssmi();

      //! @brief constructor from members
      //! @param SS_TYPE type of sse
      //! @param ENVIRONMENT environment of the residue
      //! @param FACES_PORE true if the residue faces the pore of the membrane
      Mahssmi
      (
        const biol::SSType &SS_TYPE,
        const biol::EnvironmentType &ENVIRONMENT,
        const bool &FACES_PORE,
        const bool &ORIGINATES_CYTOSOL = false
      );

      //! @brief Clone function
      //! @return pointer to new Mahssmi
      Mahssmi *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get file extension associated with this Method
      //! @return file extension associated with this Method
      const std::string &GetFileExtension() const;

      //! @brief get whether this method determined the secondary structure / membrane environment from the structure
      //! @return true if this method determined the secondary structure / membrane environment from the structure
      bool GetIsDeterminedFromSturcture() const
      {
        return true;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief get whether the residue faces the pore
      //! @return true if the residue is in a beta barrel and faces the pore
      bool DoesBetaBarrelResidueFacePore() const;

      //! @brief get whether the residue is in a TM-region whose N-term is towards the cytosol
      //! @return true if the residue is in a TM-region whose N-term is towards the cytosol
      bool OriginatesInCytosol() const;

      //! @brief get three state, environment independent secondary structure prediction
      //! @return three state, environment independent secondary structure prediction
      linal::Vector3D GetThreeStatePrediction() const;

      //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
      //! @return three state secondary structure prediction
      linal::Matrix< double> GetNineStatePrediction() const;

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

      //! @brief iterates over the sequences in ProteinModel and calculates MAHSSMI for every residue in the sequence
      //! @param PROTEIN_MODEL ProteinModel for which MAHSSMI will be calculated
      //! @param USE_PDBTM_MEMBRANE_THICKNESS true to use membrane thickness from the PDBTM if it was read in
      static void Calculate( assemble::ProteinModel &PROTEIN_MODEL, const bool &USE_PDBTM_MEMBRANE_THICKNESS);

      //! @brief iterates over the sequence and calculates MAHSSMI for every residue in the sequence
      //! @param SEQUENCE Sequence of interest
      //! @param MEMBRANE_PTR pointer to membrane object
      static void Calculate( biol::AASequence &SEQUENCE, const util::SiPtr< const biol::Membrane> &MEMBRANE_PTR);

      //! @brief iterates over the sequence and writes MAHSSMI for every residue in every chain
      //! @param STREAM stream to write the analysis to
      //! @param PROTEIN_MODEL Model to write the analysis for
      static void WriteAnalysis( std::ostream &STREAM, const assemble::ProteinModel &PROTEIN_MODEL);

      //! @brief iterates over the protein and calculates MAHSSMI for every residue in every chain
      //! @param STREAM stream to write the analysis to
      //! @param SEQUENCE Sequence to write the analysis for
      //! @param CHAIN_ID Id of the chain that the sequence represents
      static void WriteAnalysis( std::ostream &STREAM, const biol::AASequence &SEQUENCE, const char &CHAIN_ID);

      //! @brief Compute which TM regions originate in cytosol
      static void ComputeCytosolicTMOrientation( biol::AASequence &AA_SEQUENCE);

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

    }; // class Mahssmi

  } // namespace sspred
} // namespace bcl

#endif // BCL_SSPRED_MAHSSMI_H_
