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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "pdb/bcl_pdb_printer_quality_multimer.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_printer_protein_model_multimer.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "assemble/bcl_assemble_quality_batch.h"
#include "pdb/bcl_pdb_line.h"
#include "pdb/bcl_pdb_printer_score.h"
#include "storage/bcl_storage_table.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PrinterQualityMultimer::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterQualityMultimer())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterQualityMultimer::PrinterQualityMultimer() :
      m_Qualities(),
      m_NativeMultimer()
    {
    }

    //! @brief construct from members
    //! @param QUALITIES quality measures to calculate
    //! @param NATIVE native multimer model
    PrinterQualityMultimer::PrinterQualityMultimer
    (
      const storage::Set< quality::Measure> &QUALITIES,
      const util::ShPtr< assemble::ProteinModel> &NATIVE
    ) :
      m_Qualities( QUALITIES),
      m_NativeMultimer( NATIVE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterQualityMultimer
    PrinterQualityMultimer *PrinterQualityMultimer::Clone() const
    {
      return new PrinterQualityMultimer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterQualityMultimer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a protein model and returns PDB lines
    //! @param PROTEIN_MODEL protein model to print
    //! @return PDB lines
    util::ShPtrList< Line> PrinterQualityMultimer::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // cast a pointer to the multiplier data if any
      const util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Multiplier)
      );

      // if no multiplier data, return empty list
      if( !sp_multiplier.IsDefined() || !m_NativeMultimer.IsDefined())
      {
        return util::ShPtrList< Line>();
      }

      // model to hold the model multimer
      assemble::ProteinModel model_multimer( sp_multiplier->operator()( PROTEIN_MODEL));

      // set model as best multimer based on lowest RMSD
      model_multimer = assemble::PrinterProteinModelMultimer::CalculateBestMultimer
      (
        model_multimer,
        *m_NativeMultimer,
        quality::GetMeasures().e_RMSD,
        sp_multiplier
      );

      // construct the table
      const storage::Table< double> table
      (
        assemble::QualityBatch( m_Qualities, biol::GetAtomTypes().CA, "_mult").ConstructTable
        (
          model_multimer,
          *m_NativeMultimer
        )
      );

      // return pdb lines
      return PrinterScore::WriteTableToLines( table, false);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterQualityMultimer::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterQualityMultimer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
