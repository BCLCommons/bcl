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

#ifndef BCL_SSPRED_CI_PHI_PSI_H_
#define BCL_SSPRED_CI_PHI_PSI_H_

// include the namespace header
#include "bcl_sspred.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_sspred_method_interface.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CIPhiPsi
    //! @brief identifies secondary structure and membrane topology based primarily on phi-psi angles and AA distances,
    //!        and membrane location for TM strands.  Specifically, TM-strand residues are relabeled as coil whenever
    //!        the alternating pattern between pore and membrane orientation is broken.
    //!
    //! @see @link example_sspred_ci_phi_psi.cpp @endlink
    //! @author mendenjl
    //! @date Apr 09, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CIPhiPsi :
      public MethodInterface
    {

    public:

      //! Direction of TM-ss elements
      enum TMDirection
      {
        e_LeavingCytosol,
        e_EnteringCytosol,
        e_Amphipathic,       //!< always helix
        e_Pore,              //!< could be any ss
        e_Inside,            //!< Inside cytosol
        e_Outside,           //!< Outside cytosol
        e_NonMembrane,
        e_Unknown,           //!< Unknown direction, appropriate for residues with undefined positions
        s_NumberTMDirections //!< Non-membrane protein
      };

      //! @brief Data as string
      //! @param DATA the data whose name is desired
      //! @return the name as string
      static const std::string &GetTMDirectionType( const TMDirection &DATA);

      //! DataEnum simplifies the usage of the Data enum of this class
      typedef util::WrapperEnum< TMDirection, &GetTMDirectionType, s_NumberTMDirections>
        TMDirectionTypeEnum;

    private:

    //////////
    // data //
    //////////

      //! native SSType defined by the MAHSSMI algorithm
      biol::SSType m_SSType;

      //! native environment type defined by MAHSSMI
      biol::EnvironmentType m_EnvironmentType;

      //! TM direction
      TMDirectionTypeEnum m_TMDirection;

      //! for residues in a membrane beta barrel, true if the residues' CA-CB points towards the pore
      //! as opposed to membrane or neighboring member of an oligomer
      bool m_BetaBarrelResidueFacesPore;

      //! sse number, stored for convenience to calling classes
      size_t m_SSENumber;

      typedef storage::Triplet
      <
        util::SiPtr< biol::AABase>,
        storage::VectorND< 2, double>,
        double
      > t_SegmentStorage;

    public:

      //! single instance of that class
      static const util::ObjectInterface *s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CIPhiPsi();

      //! @brief constructor from members
      //! @param SS_TYPE type of sse
      //! @param ENVIRONMENT environment of the residue
      //! @param FACES_PORE true if the residue faces the pore of the membrane
      //! @param DIRECTION Direction of trans-membrane secondary structure elements. Undefined for soluble elements
      CIPhiPsi
      (
        const biol::SSType &SS_TYPE,
        const biol::EnvironmentType &ENVIRONMENT,
        const bool &FACES_PORE,
        const TMDirection &DIRECTION = e_NonMembrane
      );

      //! @brief Clone function
      //! @return pointer to new CIPhiPsi
      CIPhiPsi *Clone() const;

      ~CIPhiPsi();

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

      //! @brief get whether the residue faces the pore
      //! @return true if the residue is in a beta barrel and faces the pore
      bool DoesBetaBarrelResidueFacePore() const;

      //! @brief get whether the residue is in a TM-region whose N-term is towards the cytosol
      //! @return true if the residue is in a TM-region whose N-term is towards the cytosol
      TMDirectionTypeEnum GetTMDirection() const;

      //! @brief get three state, environment independent secondary structure prediction
      //! @return three state, environment independent secondary structure prediction
      linal::Vector3D GetThreeStatePrediction() const;

      //! @brief get the SSE Number
      //! @return the sse number
      const size_t &GetSSENumber() const;

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
      void Calculate( assemble::ProteinModel &PROTEIN_MODEL, const bool &USE_PDBTM_MEMBRANE_THICKNESS) const;

      //! @brief iterates over the sequence and calculates MAHSSMI for every residue in the sequence
      //! @param SEQUENCE Sequence of interest
      //! @param MEMBRANE_PTR pointer to membrane object
      void Calculate( biol::AASequence &SEQUENCE, const util::SiPtr< const biol::Membrane> &MEMBRANE_PTR) const;

      //! @brief iterates over the sequence and writes MAHSSMI for every residue in every chain
      //! @param STREAM stream to write the analysis to
      //! @param PROTEIN_MODEL Model to write the analysis for
      void WriteAnalysis( std::ostream &STREAM, const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief iterates over the protein and calculates MAHSSMI for every residue in every chain
      //! @param STREAM stream to write the analysis to
      //! @param SEQUENCE Sequence to write the analysis for
      //! @param CHAIN_ID Id of the chain that the sequence represents
      void WriteAnalysis( std::ostream &STREAM, const biol::AASequence &SEQUENCE, const char &CHAIN_ID) const;

      //! @brief set SSE Numbers for a given AA Sequence
      static void SetSSENumbers( biol::AASequence &AA_SEQUENCE);

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

      //! @brief iterates over the sequences in ProteinModel and calculates CIPhiPsi for every residue in the sequence
      //! @param PROTEIN_MODEL ProteinModel for which CIPhiPsi will be calculated
      void CalculateSolubleSS( assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief assign SS for all residues, ignoring the membrane
      //! @param SEQUENCE Sequence of interest
      //! @param MEMBRANE_PTR pointer to membrane object
      void CalculateSolubleSS( biol::AASequence &SEQUENCE, const double &MEM_THICKNESS = 0.0) const;

      //! @brief assign ss for the contiguous segment of residues given
      //! @param SEGMENT_AA_PHI_PSI aas, phi, psi angles in the segment
      //! @return vector of predicted SS types
      storage::Vector< biol::SSType> CalculateSolubleSSForSegment
      (
        storage::Vector< t_SegmentStorage> &SEGMENT_AA_PHI_PSI,
        const double &MEM_THICKNESS = 0.0
      ) const;

      //! @brief assign ss for the contiguous segment of residues given
      //! @param SEGMENT_AA_PHI_PSI aas, phi, psi angles in the segment
      //! @return vector of predicted SS types
      storage::Vector< TMDirectionTypeEnum> CalculateTMDirectionsForSegment
      (
        storage::Vector< t_SegmentStorage> &SEGMENT_AA_PHI_PSI,
        const double &MEM_THICKNESS
      ) const;

    }; // class CIPhiPsi

  } // namespace sspred
} // namespace bcl

#endif // BCL_SSPRED_CI_PHI_PSI_H_
