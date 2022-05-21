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

#ifndef BCL_PDB_PRINTER_QUALITY_DOCKING_H_
#define BCL_PDB_PRINTER_QUALITY_DOCKING_H_

// include the namespace header
#include "bcl_pdb_line.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "quality/bcl_quality_measures.h"
#include "score/bcl_score_protein_model_score_sum.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @class PrinterQualityDocking
    //! @brief Prints PDB lines that are quality measures for docking, for example RMSDs over the ligand
    //! @details Evaluates the quality of a docked model of a protein complex based on the given set of quality
    //!          measures and chain IDs of the ligand. Writes out lines containing detailed quality information
    //!          to PDB file.
    //!
    //! @see @link example_pdb_printer_quality_docking.cpp @endlink
    //! @author lib14
    //! @date May 2, 2017
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PrinterQualityDocking :
      public util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< Line> >
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      //! quality measures to be printed out for a more detailed information about the protein model
      storage::Set< quality::Measure> m_QualityMeasures;

      //! ID of the chain over which to evaluate the quality measures
      std::string m_ChainIDs;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PrinterQualityDocking();

      //! @brief construct from given quality measures and chain IDs
      //! @param QUALITY_MEASURES quality measures to use
      //! @param CHAIN_IDS chain IDs of the ligand
      PrinterQualityDocking
      (
        const storage::Set< quality::Measure> &QUALITY_MEASURES,
        const std::string &CHAIN_IDS
      );

      //! @brief Clone function
      //! @return pointer to a copy of this PrinterQualityDocking object
      PrinterQualityDocking *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief get quality measures
      //! @return quality measures
      const storage::Set< quality::Measure> &GetQualityMeasures() const;

      //! @brief set quality measures
      void SetQualityMeasures( const storage::Set< quality::Measure> &QUALITY_MEASURES);

      //! @brief get chain IDs of the ligand
      //! @return chain IDs of the ligand
      const std::string &GetChainIDs() const;

      //! @brief set chain IDs to retrieve for the ligand
      void SetChainIDs( const std::string &CHAIN_IDS);

    ///////////////
    // operators //
    ///////////////

      //! @brief write lines containing model quality details
      //! @param MODEL the protein model of which to write quality measure lines
      //! @return lines containing model quality details
      util::ShPtrList< Line> operator()( const assemble::ProteinModel &MODEL) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    };

  } // namespace pdb
} // namespace bcl

#endif // BCL_PDB_PRINTER_QUALITY_DOCKING_H_
