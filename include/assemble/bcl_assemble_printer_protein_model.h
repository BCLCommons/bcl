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

#ifndef BCL_ASSEMBLE_PRINTER_PROTEIN_MODEL_H_
#define BCL_ASSEMBLE_PRINTER_PROTEIN_MODEL_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
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
    //! @class PrinterProteinModel
    //! @brief Prints out a protein model during a MC minimization
    //! @details Prints start and end models from a Monte Carlo minimization using a storage object
    //!
    //! @see @link example_assemble_printer_protein_model.cpp @endlink
    //! @author weinerbe, fischea
    //! @date Nov 15, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PrinterProteinModel :
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

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PrinterProteinModel();

      //! @brief construct with prefix
      //! @param PREFIX prefix string
      PrinterProteinModel( const std::string &PREFIX);

      //! @brief construct with all member variables
      //! @param PREFIX prefix string
      //! @param STORAGE protein storage to use
      //! @param SUPERIMPOSE measure to use for superimposition
      PrinterProteinModel
      (
        const std::string &PREFIX,
        const util::ShPtr< ProteinStorageFile> &STORAGE,
        const quality::SuperimposeMeasure &SUPERIMPOSE = quality::GetSuperimposeMeasures().e_NoSuperimpose
      );

      //! @brief Clone function
      //! @return pointer to new PrinterProteinModel
      PrinterProteinModel *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
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

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief gets the source and key strings from the minimization information
      //! @param PREFIX prefix
      //! @param TAG tag for the step status if applicable
      //! @param ROUND_NUMBER round number
      //! @param STAGE_NUMBER stage number
      //! @return source and key strings, respectively
      static storage::VectorND< 2, std::string> GetStorageStrings
      (
        const std::string &PREFIX,
        const std::string &TAG,
        const size_t ROUND_NUMBER,
        const size_t STAGE_NUMBER
      );

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

    }; // class PrinterProteinModel

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_PRINTER_PROTEIN_MODEL_H_
