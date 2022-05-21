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
#include "pdb/bcl_pdb_printer_quality_membrane.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality_batch.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
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
    const util::SiPtr< const util::ObjectInterface> PrinterQualityMembrane::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterQualityMembrane())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterQualityMembrane::PrinterQualityMembrane() :
      m_Qualities(),
      m_Environments(),
      m_NativeModel()
    {
    }

    //! @brief construct from members
    //! @param QUALITIES quality measures to calculate
    //! @param ENVIRONMENTS environment types
    //! @param NATIVE native model
    PrinterQualityMembrane::PrinterQualityMembrane
    (
      const storage::Set< quality::Measure> &QUALITIES,
      const storage::Set< biol::EnvironmentType> &ENVIRONMENTS,
      const util::ShPtr< assemble::ProteinModel> &NATIVE
    ) :
      m_Qualities( QUALITIES),
      m_Environments( ENVIRONMENTS),
      m_NativeModel( NATIVE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterQualityMembrane
    PrinterQualityMembrane *PrinterQualityMembrane::Clone() const
    {
      return new PrinterQualityMembrane( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterQualityMembrane::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a protein model and returns PDB lines
    //! @param PROTEIN_MODEL protein model to print
    //! @return PDB lines
    util::ShPtrList< Line> PrinterQualityMembrane::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // get the membrane
      const util::ShPtr< biol::Membrane> sp_membrane
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
      );

      // if no membrae data, return empty list
      if( !sp_membrane.IsDefined() || !m_NativeModel.IsDefined())
      {
        return util::ShPtrList< Line>();
      }

      // initialize undefined CA atom
      const biol::Atom undefined_ca( linal::Vector3D( util::GetUndefined< double>()), biol::GetAtomTypes().CA);

      // hardcopy the native
      util::ShPtr< assemble::ProteinModel> native_copy( m_NativeModel->HardCopy());

      // iterate over the SSEs
      const util::SiPtrVector< const assemble::SSE> model_sses( native_copy->GetSSEs());
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( model_sses.Begin()),
          sse_itr_end( model_sses.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // hardcopy the sse
        util::ShPtr< assemble::SSE> sp_sse( ( *sse_itr)->HardCopy());

        // iterate over the SSE
        for
        (
          util::ShPtrVector< biol::AABase>::iterator aa_itr( sp_sse->Begin()), aa_itr_end( sp_sse->End());
          aa_itr != aa_itr_end; ++aa_itr
        )
        {
          // get the CA
          biol::Atom ca( ( *aa_itr)->GetCA());

          // if it is not in the right environment
          if( !m_Environments.Contains( sp_membrane->DetermineEnvironmentType( ca.GetCoordinates())))
          {
            // set the CA to undefined
            ( *aa_itr)->SetAtom( undefined_ca);
          }
        }

        // replace the SSE in the model
        native_copy->Replace( sp_sse);
      }

      // construct the table
      const storage::Table< double> table
      (
        assemble::QualityBatch( m_Qualities, biol::GetAtomTypes().CA, "_mem").ConstructTable
        (
          PROTEIN_MODEL,
          *native_copy
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
    std::istream &PrinterQualityMembrane::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterQualityMembrane::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
