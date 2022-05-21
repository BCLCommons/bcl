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

#ifndef BCL_ASSEMBLE_PRINTER_PROTEIN_MODEL_MULTIMER_H_
#define BCL_ASSEMBLE_PRINTER_PROTEIN_MODEL_MULTIMER_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "mc/bcl_mc.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "mc/bcl_mc_print_interface.h"
#include "quality/bcl_quality_measures.h"
#include "quality/bcl_quality_superimpose_measures.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrinterProteinModelMultimer
    //! @brief Prints quality measures and optional PDB file for multimeric proteins
    //! @details Appends multimer quality measures to pdb file.  Also optionally superimposes the multimer with a
    //!          multimer native and writes a separate pdb file.
    //!
    //! @see @link example_assemble_printer_protein_model_multimer.cpp @endlink
    //! @author weinerbe, alexanns
    //! @date May 6, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PrinterProteinModelMultimer :
      public mc::PrintInterface< ProteinModel, double>
    {

    private:

    //////////
    // data //
    //////////

      //! prefix prepended to all model names
      std::string m_Prefix;

      //! superimposition measure to use
      quality::SuperimposeMeasure m_Superimpose;

      //! round number
      size_t m_RoundNumber;

      //! stage number
      size_t m_StageNumber;

      //! native multimer
      util::ShPtr< ProteinModel> m_NativeMultimer;

      //! storage
      mutable util::ShPtr< ProteinStorageFile> m_Storage;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PrinterProteinModelMultimer();

      //! @brief construct with all member variables
      //! @param PREFIX prefix string
      //! @param NATIVE_MULTIMER_MODEL native multimer model
      //! @param STORAGE protein storage to use
      //! @param SUPERIMPOSE measure to use for superimposition
      PrinterProteinModelMultimer
      (
        const std::string &PREFIX,
        const util::ShPtr< ProteinModel> &NATIVE_MULTIMER_MODEL,
        const util::ShPtr< ProteinStorageFile> &STORAGE,
        const quality::SuperimposeMeasure &SUPERIMPOSE = quality::GetSuperimposeMeasures().e_NoSuperimpose
      );

      //! @brief Clone function
      //! @return pointer to new PrinterProteinModelMultimer
      PrinterProteinModelMultimer *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return prefix
      //! @return prefix
      const std::string &GetPrefix() const;

      //! @brief set prefix to given PREFIX
      //! @param PREFIX new prefix
      void SetPrefix( const std::string &PREFIX);

    ////////////////
    // operations //
    ////////////////

      //! @brief reset and initialize the printer
      //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
      //! @return true if initialization was successful
      void Initialize( const size_t &ROUND_NUMBER);

      //! @brief reset and initialize the printer
      //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
      //! @param STAGE_NUMBER for multiple stage optimizations, a different round and stage number will be passed
      //! @return true if initialization was successful
      void Initialize( const size_t &ROUND_NUMBER, const size_t &STAGE_NUMBER);

      //! @brief prints information concerning the approximation process based on the status of the given tracker
      //! @param TRACKER holds the status of the approximation process
      void Print( const opti::Tracker< ProteinModel, double> &TRACKER) const;

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief function to write given model
      //! @param ARGUMENT_RESULT_PAIR Pair of argument and corresponding result
      //! @param STORAGE_STRINGS source and key strings from the minimization information
      //! @return whether write was successful
      bool WriteToStorage
      (
        const util::ShPtr< storage::Pair< ProteinModel, double> > &ARGUMENT_RESULT_PAIR,
        const storage::VectorND< 2, std::string> &STORAGE_STRINGS
      ) const;

    public:

      //! @brief determines the best combination of chain ids to match the native model
      //! @param MODEL model whose chains will be changed in order to find the best combination for superimposing
      //! @param NATIVE the model that will be used to calculate quality against
      //! @param QUALITY the quality that will be used to calculate the best model
      //! @param MULTIPLIER multiplier used to create the model
      //! @return protein model which can be best superimposed onto NATIVE
      static ProteinModel CalculateBestMultimer
      (
        const ProteinModel &MODEL,
        const ProteinModel &NATIVE,
        const quality::Measure &QUALITY,
        const util::ShPtr< ProteinModelMultiplier> &MULTIPLIER
      );

    }; // class PrinterProteinModelMultimer

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_PRINTER_PROTEIN_MODEL_MULTIMER_H_
