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
#include "pdb/bcl_pdb_printer_quality_docking.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_quality_batch.h"
#include "pdb/bcl_pdb_printer_score.h"

namespace bcl
{
  namespace pdb
  {
    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterQualityDocking())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterQualityDocking::PrinterQualityDocking() :
      m_QualityMeasures(),
      m_ChainIDs()
    {
    }

    //! @brief construct from given quality measures and chain IDs
    //! @param QUALITY_MEASURES quality measures to use
    //! @param CHAIN_IDS chain IDs of the ligand
    PrinterQualityDocking::PrinterQualityDocking
    (
      const storage::Set< quality::Measure> &QUALITY_MEASURES,
      const std::string &CHAIN_IDS
    ) :
      m_QualityMeasures( QUALITY_MEASURES),
      m_ChainIDs( CHAIN_IDS)
    {
    }

    //! @brief Clone function
    //! @return pointer to a copy of this PrinterQualityDocking object
    PrinterQualityDocking *PrinterQualityDocking::Clone() const
    {
      return new PrinterQualityDocking( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get the name of this class
    //! @return the name of this class
    const std::string &PrinterQualityDocking::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get quality measures
    //! @return quality measures
    const storage::Set< quality::Measure> &PrinterQualityDocking::GetQualityMeasures() const
    {
      return m_QualityMeasures;
    }

    //! @brief set quality measures
    void PrinterQualityDocking::SetQualityMeasures
    (
      const storage::Set< quality::Measure> &QUALITY_MEASURES
    )
    {
      m_QualityMeasures = QUALITY_MEASURES;
    }

    //! @brief get chain IDs of the ligand
    //! @return chain IDs of the ligand
    const std::string &PrinterQualityDocking::GetChainIDs() const
    {
      return m_ChainIDs;
    }

    //! @brief set chain IDs to retrieve for the ligand
    void PrinterQualityDocking::SetChainIDs( const std::string &CHAIN_IDS)
    {
      m_ChainIDs = CHAIN_IDS;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief write lines containing model quality details
    //! @param MODEL the protein model of which to write quality measure lines
    //! @return lines containing model quality details
    util::ShPtrList< Line> PrinterQualityDocking::operator()
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      // create a ProteinModel object from the docked ligand chains
      assemble::ProteinModel docked_ligand_model( MODEL.GetChains( m_ChainIDs));

      // create a ProteinModel object from the native ligand chains
      util::ShPtr< assemble::ProteinModel> sp_native_model
      (
        MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_NativeModel)
      );

      // if native model is not defined, return nan
      if( !sp_native_model.IsDefined())
      {
        // get table headers
        storage::Table< double> undefined_table
        (
          math::SumFunctionMixin< score::ProteinModel>::GetValueTableVerticalColumnNames()
        );
        for
        (
          storage::Set< quality::Measure>::const_iterator
            measure_itr( m_QualityMeasures.Begin()), measure_itr_end( m_QualityMeasures.End());
          measure_itr != measure_itr_end;
          ++measure_itr
        )
        {
          undefined_table.InsertRow
          (
            measure_itr->GetName(),
            storage::Vector< double>::Create( 1.0, util::GetUndefinedDouble(), util::GetUndefinedDouble())
          );
        }
        return PrinterScore::WriteTableToLines( undefined_table);
      }

      // get native ligand model
      const assemble::ProteinModel native_ligand_model( sp_native_model->GetChains( m_ChainIDs));

      // add the native ligand model into docked ligand model
      util::ShPtr< assemble::ProteinModelData> sp_data( new assemble::ProteinModelData());
      sp_data->Insert
      (
        assemble::ProteinModelData::e_NativeModel,
        util::ShPtr< assemble::ProteinModel>( native_ligand_model.Clone())
      );
      sp_data->Insert
      (
        assemble::ProteinModelData::e_NativeFilteredModel,
        util::ShPtr< assemble::ProteinModel>( native_ligand_model.Clone())
      );
      docked_ligand_model.SetProteinModelData( sp_data);

      // construct table for qualities
      const assemble::QualityBatch qualities( m_QualityMeasures, biol::GetAtomTypes().CA);
      const storage::Table< double> table( qualities.ConstructTable( docked_ligand_model));

      // write table to PDB lines
      return PrinterScore::WriteTableToLines( table, false);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterQualityDocking::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterQualityDocking::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace pdb
} // namespace bcl
