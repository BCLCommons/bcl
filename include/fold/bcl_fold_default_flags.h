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

#ifndef BCL_FOLD_DEFAULT_FLAGS_H_
#define BCL_FOLD_DEFAULT_FLAGS_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "sspred/bcl_sspred.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DefaultFlags
    //! @brief This class provides access to all the flags used in default Fold protocol
    //! @details This class has static access functions to all the flags used in Fold protocol and serves as a
    //! supplement to Fold setup class.
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date Nov 11, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DefaultFlags :
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
      DefaultFlags();

      //! @brief Clone function
      //! @return pointer to new DefaultFlags
      DefaultFlags *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////
    // flags //
    ///////////

      //! @brief returns all flags that are used by default
      //! @return all flags that are used by default
      static const util::ShPtrVector< command::FlagInterface> &GetAllFlags();

      //! @brief return command line flag for writing the minimization to files - only of given step statuses
      //! @return command line flag for writing the minimization to files - only of given step statuses
      static util::ShPtr< command::FlagInterface> &GetFlagPrintMinimization();

      //! @brief return command line parameter for specifying the path for where tracker files are generated
      //! @return command line parameter for specifying the path for where tracker files are generated
      static util::ShPtr< command::ParameterInterface> &GetParameterPrintTrackerHistoryPath();

      //! @brief return command line flag for a detailed analysis of each step in a tabulated format
      //! @return command line flag for a detailed analysis of each step in a tabulated format
      static util::ShPtr< command::FlagInterface> &GetFlagPrintTrackerHistory();

      //! @brief return command line flag for a detailed analysis of success rates of each mutate
      //! @return return command line flag for a detailed analysis of success rates of each mutate
      static util::ShPtr< command::FlagInterface> &GetFlagPrintStepCounts();

      //! @brief return command line flag for monte carlo minimization, max number of rejected steps and max iterations
      //! @return command line flag for monte carlo minimization, max number of rejected steps and max iterations
      static util::ShPtr< command::FlagInterface> &GetFlagMCNumberIterations();

      //! @brief return command line flag for specifying whether or not to print the start model
      //! @return command line flag for specifying  whether or not to print the start model
      static util::ShPtr< command::FlagInterface> &GetFlagPrintStartModel();

      //! @brief return command line flag for specifying a prefix to be used for writing files
      //! @return command line flag for specifying a prefix to be used for writing files
      static util::ShPtr< command::FlagInterface> &GetFlagPrefix();

      //! @brief return command line flag for specifying the number of models to build
      //! @return command line flag for specifying the number of models to build
      static util::ShPtr< command::FlagInterface> &GetFlagNumberModels();

      //! @brief return command line flag for providing one or more fasta files for complete de novo folding
      //! @return command line flag for providing one or more fasta files for complete de novo folding
      static util::ShPtr< command::FlagInterface> &GetFlagFastaRead();

      //! @brief return command line flag for providing one or more chain ids for fasta files for complete de novo folding
      //! @return command line flag for providing one or more chain ids for fasta files for complete de novo folding
      static util::ShPtr< command::FlagInterface> &GetFlagChainIdRead();

      //! @brief return command line flag for providing a native model for comparison
      //! @return command line flag for providing a native model for comparison
      static util::ShPtr< command::FlagInterface> &GetFlagNativeModel();

      //! @brief return command line flag for providing a starting model
      //! @return command line flag for providing a starting model
      static util::ShPtr< command::FlagInterface> &GetFlagStartModel();

      //! @brief return command line flag for using native SSE definitions as the SSE pool
      //! @return command line flag for using native SSE definitions as the SSE pool
      static util::ShPtr< command::FlagInterface> &GetFlagUseNativeSSEsAsPool();

      //! @brief return command line flag for enable resizing of SSEs
      //! @return command line flag for enable resizing of SSEs
      static util::ShPtr< command::FlagInterface> &GetFlagEnableSSEResize();

      //! @brief return command line flag for separating adjoining SSE pools with specified number of loop residues
      //! @return command line flag for separating adjoining SSE pools with specified number of loop residues
      static util::ShPtr< command::FlagInterface> &GetFlagPoolSeparate();

      //! @brief return command line flag for specifying the path where the ss predictions should be read from
      //! @return command line flag for specifying the path where the ss predictions should be read from
      static util::ShPtr< command::FlagInterface> &GetFlagReadSequenceDataPath();

      //! @brief return command line flag for reading scores from a file
      //! @return command line flag for reading scores from a file
      static util::ShPtr< command::FlagInterface> &GetFlagScoreRead();

      //! @brief return command line flag for writing scores to a file
      //! @return command line flag for writing scores to a file
      static util::ShPtr< command::FlagInterface> &GetFlagScoreWrite();

      //! @brief return command line flag for reading score weight set from a file
      //! @return command line flag for reading score weight set from a file
      static util::ShPtr< command::FlagInterface> &GetFlagScoreWeightSetRead();

      //! @brief return command line flag for writing score weight set to a file
      //! @return command line flag for writing score weight set to a file
      static util::ShPtr< command::FlagInterface> &GetFlagScoreWeightSetWrite();

      //! @brief return command line flag for reading mutates from a file
      //! @return command line flag for reading mutates from a file
      static util::ShPtr< command::FlagInterface> &GetFlagMutateRead();

      //! @brief return command line flag for writing mutates to a file
      //! @return command line flag for writing mutates to a file
      static util::ShPtr< command::FlagInterface> &GetFlagMutateWrite();

      //! @brief return command line flag for reading mutate weight set from a file
      //! @return command line flag for reading mutate weight set from a file
      static util::ShPtr< command::FlagInterface> &GetFlagMutateWeightSetRead();

      //! @brief return command line flag for writing mutate weight set to a file
      //! @return command line flag for writing mutate weight set to a file
      static util::ShPtr< command::FlagInterface> &GetFlagMutateWeightSetWrite();

      //! @brief return command line flag for defining the filename to read the stages from
      //! @return command line flag for defining the filename to read the stages from
      static util::ShPtr< command::FlagInterface> &GetFlagStagesFileRead();

      //! @brief return command line flag for defining the filename to write the stages to
      //! @return command line flag for defining the filename to write the stages to
      static util::ShPtr< command::FlagInterface> &GetFlagStagesFileWrite();

      //! @brief return command line flag for defining the number of cycles stages should go through
      //! @return command line flag for defining the number of cycles stages should go through
      static util::ShPtr< command::FlagInterface> &GetFlagStagesNumberCycles();

      //! @brief return command line flag for defining whether to fit swapped sses
      //! @return command line flag for defining whether to fit swapped sses
      static util::ShPtr< command::FlagInterface> &GetFlagFitSwappedSSEs();

      //! @brief command line flag indicating the input file numbering is PDB numbering instead of sequence id numbering
      //! @return flag indicating the input file numbering is PDB numbering instead of sequence id numbering
      static util::ShPtr< command::FlagInterface> &GetFlagPDBIDNumbering();

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

    }; // class DefaultFlags

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_DEFAULT_FLAGS_H_
