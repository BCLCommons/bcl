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

#ifndef BCL_APP_MSA2PSSM_H_
#define BCL_APP_MSA2PSSM_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MSA2PSSM
    //! @brief Reads a CLUSTAL alignment and generates a PSSM, optionally using constant pseudocounts and background
    //!        probabilities
    //!
    //! @author mendenjl
    //! @date Oct 31, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MSA2PSSM :
      public Interface
    {
    private:

    //////////
    // data //
    //////////

      //! input fasta file, alignment file, and output file
      util::ShPtr< command::ParameterInterface> m_AlignmentFile;
      util::ShPtr< command::ParameterInterface> m_OutputFile;

      //! Number of pseudocounts to use
      util::ShPtr< command::FlagInterface> m_PseudocountFlag;

      //! Whether to use MSA to compute background frequencies for pseudocounts
      util::ShPtr< command::FlagInterface> m_UseMSABackgroundFlag;

      //! If true, number of pseudocounts is interpreted as a fraction of the residues that aligned at a given position
      util::ShPtr< command::FlagInterface> m_FractionalPseudocountsFlag;

      //! If true, number of pseudocounts is interpreted as a fraction of the residues that aligned at a given position
      util::ShPtr< command::FlagInterface> m_EffFractionalPseudocountsFlag;

      //! If set, all pssm positions are simply log(pseudocount+real_count), without any row-based normalization
      util::ShPtr< command::FlagInterface> m_LogOnlyFlag;

      //! Pseudocount value
      mutable double m_Pseudocount;

      //! Fractional pseudocount flag value
      mutable double m_FractionalPseudocount;

      //! Effective fractional pseudocount flag value
      mutable double m_EffFractionalPseudocount;

      //! number of sequences in the alignment; used to compute alignment weight
      mutable size_t m_NumberSeqs;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! @brief default constructor
      MSA2PSSM();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      MSA2PSSM *Clone() const
      {
        return new MSA2PSSM( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      // instantiate enumerator for MSA2PSSM class
      static const ApplicationType MSA2PSSM_Instance;

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! Main
      int Main() const;

      //! @brief compute a blast profile
      //! @param ALIGNED AA-types that aligned at this position (should contain only upper case letters, no gaps).
      //!        First letter should be the sequence character itself
      //! @param BACKGROUND_AA_FREQ frequency of AAs in the MSA
      //! @return BlastProfile
      biol::BlastProfile ComputeBlastProfile
      (
        const std::string &ALIGNED,
        const linal::Vector< double> &BACKGROUND_AA_FREQ
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    private:

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class MSA2PSSM

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MSA2PSSM_H_
