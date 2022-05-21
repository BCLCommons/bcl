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

#ifndef BCL_ASSEMBLE_PRINTER_PROTEIN_MODEL_ENSEMBLE_H_
#define BCL_ASSEMBLE_PRINTER_PROTEIN_MODEL_ENSEMBLE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_protein_model.h"
#include "mc/bcl_mc_print_interface.h"
#include "opti/bcl_opti_tracker.h"
#include "quality/bcl_quality_superimpose_measures.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrinterProteinModelEnsemble
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_assemble_printer_protein_model_ensemble.cpp @endlink
    //! @author alexanns
    //! @date Apr 12, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PrinterProteinModelEnsemble :
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

      //! storage
      mutable util::ShPtr< ProteinStorageFile> m_Storage;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PrinterProteinModelEnsemble();

      //! @brief construct with prefix
      //! @param PREFIX prefix string
      PrinterProteinModelEnsemble( const std::string &PREFIX);

      //! @brief construct with all member variables
      //! @param PREFIX prefix string
      //! @param STORAGE protein storage to use
      //! @param SUPERIMPOSE measure to use for superimposition
      PrinterProteinModelEnsemble
      (
        const std::string &PREFIX,
        const util::ShPtr< ProteinStorageFile> &STORAGE,
        const quality::SuperimposeMeasure &SUPERIMPOSE = quality::GetSuperimposeMeasures().e_NoSuperimpose
      );

      //! @brief Clone function
      //! @return pointer to new PrinterProteinModelEnsemble
      PrinterProteinModelEnsemble *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return const reference to the prefix member variable
      //! @return prefix as string
      const std::string &GetPrefix() const
      {
        return m_Prefix;
      }

      //! @brief set prefix to given PREFIX
      //! @param PREFIX new prefix
      void SetPrefix( const std::string &PREFIX)
      {
        m_Prefix = PREFIX;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief reset and initialize the printer
      //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
      //! @return true if initialization was successful
      void Initialize( const size_t &ROUND_NUMBER);

      //! @brief reset and initialize the printer with the given round and stage number
      //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
      //! @param STAGE_NUMBER for multiple optimizations, a different stage number will be passed
      //! @return true if initialization was successful
      void Initialize( const size_t &ROUND_NUMBER, const size_t &STAGE_NUMBER);

      //! @brief prints information concerning the approximation process based on the status of the given tracker
      //! @param TRACKER holds the status of the approximation process
      void Print( const opti::Tracker< ProteinModel, double> &TRACKER) const;

    ///////////////
    // operators //
    ///////////////

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
      //! @param ARGUMENT protein model to print,
      //! @param STORAGE_STRINGS source and key strings from the minimization information
      //! @return whether write was successful
      bool WriteToStorage
      (
        const ProteinModel &ARGUMENT,
        const storage::VectorND< 2, std::string> &STORAGE_STRINGS
      ) const;

      //! @brief writes a single model to file
      //! @param MODEL the model to write
      //! @param CONFORMATION_NUMBER the identifying number string
      //! @param STORAGE_STRINGS the other strings for storage
      void WriteModel
      (
        const ProteinModel &MODEL, const std::string &CONFORMATION_NUMBER,
        const storage::VectorND< 2, std::string> &STORAGE_STRINGS
      ) const;

    }; // class PrinterProteinModelEnsemble

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_PRINTER_PROTEIN_MODEL_ENSEMBLE_H_
