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

#ifndef BCL_BIOL_AA_SEQUENCE_FACTORY_H_
#define BCL_BIOL_AA_SEQUENCE_FACTORY_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_biol_aa_classes.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AASequenceFactory
    //! @brief Generates new sequences by reading from FASTA files or combining/fitting existing sequences
    //! @details Generates new AASequences either from FASTA files or from combining or fitting existing AA sequences.
    //!          Additionally contains some AASequence-specific functions.
    //!
    //! @see @link example_biol_aa_sequence_factory.cpp @endlink
    //! @author weinerbe, alexanns
    //! @date Jan 25, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AASequenceFactory :
      public util::ObjectInterface
    {

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
      AASequenceFactory();

      //! @brief Clone function
      //! @return pointer to new AASequenceFactory
      AASequenceFactory *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief appends a residue to the C-terminal end of a sequence
      //! @param SEQUENCE AASequence to be extended
      //! @param AMINO_ACID to be appended
      //! @param PHI phi angle for added residue
      //! @param PSI psi angle for last residue in original sequence
      static void AppendAA
      (
        AASequence &SEQUENCE,
        const AABase &AMINO_ACID,
        const double PHI,
        const double PSI
      );

      //! @brief appends a sequence to the C-terminal end of a sequence
      //! @param N_TERMINAL_SEQUENCE AASequence to be extended
      //! @param C_TERMINAL_SEQUENCE AASequence to be appended
      //! @param PHI phi angle for first residue in C-terminal sequence
      //! @param PSI psi angle for last residue in N-terminal sequence
      static void AppendSequence
      (
        AASequence &N_TERMINAL_SEQUENCE,
        const AASequence &C_TERMINAL_SEQUENCE,
        const double PHI,
        const double PSI
      );

      //! @brief transformation to append a c term sequence onto the given n-term sequence
      //! @param N_TERMINAL_SEQUENCE AASequence fixed in space
      //! @param C_TERMINAL_AA amino acid to be appended
      //! @param PHI angle for first residue in cterm
      static math::TransformationMatrix3D TransformationAppend
      (
        const AASequence &N_TERMINAL_SEQUENCE,
        const AABase &C_TERMINAL_AA,
        const double PHI
      );

      //! @brief prepends a residue to the N-terminal end of a sequence
      //! @param AMINO_ACID to be prepended
      //! @param SEQUENCE AASequence to be extended
      //! @param PHI phi angle for first residue in original sequence
      //! @param PSI psi angle for added residue
      static void PrependAA
      (
        const AABase &AMINO_ACID,
        AASequence &SEQUENCE,
        const double PHI,
        const double PSI
      );

      //! @brief prepends a sequence to the N-terminal end of a sequence
      //! @param N_TERMINAL_SEQUENCE AASequence to be prepended
      //! @param C_TERMINAL_SEQUENCE AASequence to be extended
      //! @param PHI phi angle for first residue in C-terminal sequence
      //! @param PSI psi angle for last residue in N-terminal sequence
      static void PrependSequence
      (
        const AASequence &N_TERMINAL_SEQUENCE,
        AASequence &C_TERMINAL_SEQUENCE,
        const double PHI,
        const double PSI
      );

      //! @brief transformation to append a c term sequence onto the given n-term sequence
      //! @param N_TERMINAL_AA amino acid to be appended
      //! @param C_TERMINAL_SEQUENCE AASequence fixed in space
      //! @param PHI angle for first residue in cterm
      static math::TransformationMatrix3D TransformationPrepend
      (
        const AABase &N_TERMINAL_AA,
        const AASequence &C_TERMINAL_SEQUENCE,
        const double PHI
      );

      //! @brief fits the passed sequence to the phi/psi angles
      //! @param SEQUENCE AASequence to be used
      //! @param PHI_PSI phi and psi values to be applied
      //! @param SS_TYPE sstype to be applied if the sequence is longer than the given phi/psi information
      static void FitSequence
      (
        AASequence &SEQUENCE,
        const AASequencePhiPsi &PHI_PSI,
        const SSType &SS_TYPE
      );

      //! @brief Idealize the coordinates for the given AASequence according to given SSType
      //! @param AA_SEQUENCE AASequence to be idealized
      //! @param SS_TYPE SSType to be used for idealization
      static void IdealizeSequence( AASequence &AA_SEQUENCE, const SSType &SS_TYPE);

      //! @brief read fasta from given ISTREAM and initialize the sequence with instances of AA_CLASS
      //! @param ISTREAM input stream
      //! @param AA_CLASS type of amino acid to be used
      //! @param CHAIN_ID chain id to be used for the sequence
      //! @return newly constructed AASequence
      static AASequence BuildSequenceFromFASTA
      (
        std::istream &ISTREAM,
        const AAClass &AA_CLASS = GetAAClasses().e_AA,
        const char CHAIN_ID = 'A'
      );

      //! @brief read fasta from provided SEQUENCE_STRING and initialize the sequences with amino acids of type AA_CLASS
      //! @param SEQUENCE_STRING sequence string that has the one letter code sequence
      //! @param AA_CLASS type of amino acid to be used
      //! @param CHAIN_ID chain id to be used for the sequence
      //! @return AA sequence built from the string
      static AASequence BuildSequenceFromFASTAString
      (
        const std::string &SEQUENCE_STRING,
        const AAClass &AA_CLASS = GetAAClasses().e_AA,
        const char CHAIN_ID = 'A'
      );

      //! @brief read sequence data, fasta and pdb files, from the given filenames and returns their sequences
      //! @param FILENAMES a vector of sequence filenames
      //! @return a ShPtrVector of sequences
      static util::ShPtrVector< AASequence> ReadSequenceData( const storage::Vector< std::string> &FILENAMES);

      //! @brief calculates the transformation needed for superimposing SEQUENCE_A to SEQUENCE_B
      //! @param AA_SEQUENCE_A AASequence which will be superimposed to AA_SEQUENCE_B
      //! @param AA_SEQUENCE_B AASequence to which superimposition will be calculated against
      //! @return transformation matrix that gives the best superimposition
      static math::TransformationMatrix3D CalculateSuperimposition
      (
        const AASequence &AA_SEQUENCE_A,
        const AASequence &AA_SEQUENCE_B
      );

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

    }; // class AASequenceFactory

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_AA_SEQUENCE_FACTORY_H_ 
