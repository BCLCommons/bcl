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
#include "assemble/bcl_assemble_printer_protein_model_multimer.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_alignment_node.h"
#include "assemble/bcl_assemble_printer_protein_model.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "assemble/bcl_assemble_protein_storage_file.h"
#include "assemble/bcl_assemble_quality.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PrinterProteinModelMultimer::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterProteinModelMultimer())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterProteinModelMultimer::PrinterProteinModelMultimer() :
      m_Prefix( ""),
      m_Superimpose(),
      m_RoundNumber( 0),
      m_StageNumber( util::GetUndefined< size_t>()),
      m_NativeMultimer(),
      m_Storage( ProteinStorageFile::GetDefaultStorage())
    {
    }

    //! @brief construct with all member variables
    //! @param PREFIX prefix string
    //! @param NATIVE_MULTIMER_MODEL native multimer model
    //! @param STORAGE protein storage to use
    //! @param SUPERIMPOSE measure to use for superimposition
    PrinterProteinModelMultimer::PrinterProteinModelMultimer
    (
      const std::string &PREFIX,
      const util::ShPtr< ProteinModel> &NATIVE_MULTIMER_MODEL,
      const util::ShPtr< ProteinStorageFile> &STORAGE,
      const quality::SuperimposeMeasure &SUPERIMPOSE
    ) :
      m_Prefix( PREFIX),
      m_Superimpose( SUPERIMPOSE),
      m_RoundNumber( 0),
      m_StageNumber( util::GetUndefined< size_t>()),
      m_NativeMultimer( NATIVE_MULTIMER_MODEL),
      m_Storage( STORAGE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterProteinModelMultimer
    PrinterProteinModelMultimer *PrinterProteinModelMultimer::Clone() const
    {
      return new PrinterProteinModelMultimer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterProteinModelMultimer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return prefix
    //! @return prefix
    const std::string &PrinterProteinModelMultimer::GetPrefix() const
    {
      return m_Prefix;
    }

    //! @brief set prefix to given PREFIX
    //! @param PREFIX new prefix
    void PrinterProteinModelMultimer::SetPrefix( const std::string &PREFIX)
    {
      m_Prefix = PREFIX;
    }

    //! @brief reset and initialize the printer
    //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
    //! @return true if initialization was successful
    void PrinterProteinModelMultimer::Initialize( const size_t &ROUND_NUMBER)
    {
      // update round number
      m_RoundNumber = ROUND_NUMBER;

      // set the stage number to undefined
      m_StageNumber = util::GetUndefined< size_t>();
    }

    //! @brief reset and initialize the printer with the given round and stage number
    //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
    //! @param STAGE_NUMBER for multiple optimizations, a different stage number will be passed
    //! @return true if initialization was successful
    void PrinterProteinModelMultimer::Initialize( const size_t &ROUND_NUMBER, const size_t &STAGE_NUMBER)
    {
      // initialize according to round number
      Initialize( ROUND_NUMBER);

      // update stage number
      m_StageNumber = STAGE_NUMBER;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief prints information concerning the approximation process based on the status of the given tracker
    //! @param TRACKER holds the status of the approximation process
    void PrinterProteinModelMultimer::Print( const opti::Tracker< ProteinModel, double> &TRACKER) const
    {
      // call WriteToStorage function
      WriteToStorage
      (
        TRACKER.GetBest(),
        PrinterProteinModel::GetStorageStrings( m_Prefix, TRACKER.GetTag() + "_mult", m_RoundNumber, m_StageNumber)
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterProteinModelMultimer::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterProteinModelMultimer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief function to write given model
    //! @param ARGUMENT_RESULT_PAIR Pair of argument and corresponding result
    //! @param STORAGE_STRINGS source and key strings from the minimization information
    //! @return whether write was successful
    bool PrinterProteinModelMultimer::WriteToStorage
    (
      const util::ShPtr< storage::Pair< ProteinModel, double> > &ARGUMENT_RESULT_PAIR,
      const storage::VectorND< 2, std::string> &STORAGE_STRINGS
    ) const
    {
      // copy the protein model data
      util::ShPtr< ProteinModelData> mult_pmd( ARGUMENT_RESULT_PAIR->First().GetProteinModelData()->HardCopy());

      // cast a pointer to the multiplier data if any
      util::ShPtr< ProteinModelMultiplier> sp_multiplier( mult_pmd->GetData( ProteinModelData::e_Multiplier));

      // model to hold the model multimer
      ProteinModel model_multimer;

      // if the pointer is defined
      if( sp_multiplier.IsDefined())
      {
        // apply the multiplier
        model_multimer = sp_multiplier->operator()( ARGUMENT_RESULT_PAIR->First());

        // set the pointer to null
        mult_pmd->Replace( ProteinModelData::e_Multiplier, util::ShPtr< ProteinModelMultiplier>());
      }

      // set the protein identification
      const util::ShPtr< util::Wrapper< std::string> > &sp_id
      (
        mult_pmd->GetData( ProteinModelData::e_Identification)
      );
      BCL_Assert( sp_id.IsDefined(), "No identification stored for protein model");
      mult_pmd->Replace
      (
        ProteinModelData::e_Identification,
        util::ShPtr< util::Wrapper< std::string> >( new util::Wrapper< std::string>( std::string( *sp_id) + "_mult"))
      );

      // if a native multimer model was supplied
      if( m_NativeMultimer.IsDefined())
      {
        // determine the best multimer
        model_multimer = CalculateBestMultimer( model_multimer, *m_NativeMultimer, quality::GetMeasures().e_RMSD, sp_multiplier);

        // create a new model using the chains from the native
        ProteinModel new_model( m_NativeMultimer->GetEmptyChains());

        // iterate over the SSEs in the multimer
        const util::SiPtrVector< const SSE> multimer_sses( model_multimer.GetSSEs());
        for
        (
          util::SiPtrVector< const SSE>::const_iterator sse_itr( multimer_sses.Begin()),
            sse_itr_end( multimer_sses.End());
          sse_itr != sse_itr_end; ++sse_itr
        )
        {
          // insert the sse
          new_model.Insert( util::ShPtr< SSE>( ( *sse_itr)->Clone()));
        }
        model_multimer = new_model;

        if( !mult_pmd->Insert( ProteinModelData::e_NativeModel, m_NativeMultimer))
        {
          mult_pmd->Replace( ProteinModelData::e_NativeModel, m_NativeMultimer);
        }

        // make a copy of the native model
        util::ShPtr< ProteinModel> sp_native_model_copy( m_NativeMultimer->Clone());

        // filter the SSEs by given pool min sse sizes
        sp_native_model_copy->FilterByMinSSESizes( SSEPool::GetCommandLineMinSSELengths());

        // add filtered model
        if( !mult_pmd->Insert( ProteinModelData::e_NativeFilteredModel, sp_native_model_copy))
        {
          mult_pmd->Replace( ProteinModelData::e_NativeFilteredModel, sp_native_model_copy);
        }
      }

      model_multimer.ConnectSSEToChainData();

      // if there is a native and the model should be superimposed
      if( m_Superimpose.IsDefined() && m_NativeMultimer.IsDefined())
      {
        // get the alignment between the current model and the native
        const storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > alignments
        (
          Quality::CreateAlignmentProteinModels( model_multimer, *m_NativeMultimer)
        );

        const storage::Set< biol::AtomType> &atom_type( biol::GetAtomTypes().CA);

        // initialize coordinates vector
        const storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > coords
        (
          Quality::CoordinatesFromAlignments( alignments, atom_type)
        );

        // transform the model to superimpose it to the native
        const math::TransformationMatrix3D transformation( ( *m_Superimpose)->CalculateSuperimposition( coords.Second(), coords.First()));
        model_multimer.Transform( transformation);
      }

      // set the protein model data
      model_multimer.SetProteinModelData( mult_pmd);

      // store the model
      const std::string filename( STORAGE_STRINGS.First() + STORAGE_STRINGS.Second() + ".pdb");
      if( m_Storage->Store( model_multimer, STORAGE_STRINGS.First(), STORAGE_STRINGS.Second()))
      {
        BCL_MessageStd( "pdb written to " + filename);
      }
      else
      {
        BCL_MessageCrt
        (
          "Unable to store file, " + filename + ", file may already exist and overwrite is not set"
        );
      }

      // end
      return true;
    }

    //! @brief determines the best combination of chain ids to match the native model
    //! @param MODEL model whose chains will be changed in order to find the best combination for superimposing
    //! @param NATIVE the model that will be used to calculate quality against
    //! @param QUALITY the quality that will be used to calculate the best model
    //! @param MULTIPLIER multiplier used to create the model
    //! @return protein model which can be best superimposed onto NATIVE
    ProteinModel PrinterProteinModelMultimer::CalculateBestMultimer
    (
      const ProteinModel &MODEL,
      const ProteinModel &NATIVE,
      const quality::Measure &QUALITY,
      const util::ShPtr< ProteinModelMultiplier> &MULTIPLIER
    )
    {
      // quality calculated on CAs
      const storage::Set< biol::AtomType> atom_type( biol::GetAtomTypes().CA);

      // get the alignment between the model and the native
      storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > alignments
      (
        Quality::CreateAlignmentProteinModels( MODEL, NATIVE)
      );

      // initialize coordinates vector
      storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > coords
      (
        Quality::CoordinatesFromAlignments( alignments, atom_type)
      );

      // calculate the quality measure for the model
      const double quality( ( *QUALITY)->CalculateMeasure( coords.First(), coords.Second()));

      // now reverse the chain ordering (i.e. counter-clockwise to clockwise) and see if the quality improves

      // initialize chains
      util::ShPtrVector< Chain> chains;

      // iterate over the chain mappings
      const storage::Map< char, std::string> chain_id_mappings( MULTIPLIER->GetTargetChains());
      for
      (
        storage::Map< char, std::string>::const_iterator map_itr( chain_id_mappings.Begin()),
          map_itr_end( chain_id_mappings.End());
        map_itr != map_itr_end; ++map_itr
      )
      {
        // iterate forwards and backwards through the target chains
        std::string::const_reverse_iterator rev_itr( map_itr->second.rbegin());
        const std::string::const_reverse_iterator rev_itr_end( map_itr->second.rend());
        for
        (
          std::string::const_iterator forward_itr( map_itr->second.begin()), forward_itr_end( map_itr->second.end());
          forward_itr != forward_itr_end && rev_itr != rev_itr_end; ++forward_itr, ++rev_itr
        )
        {
          // get the chain
          util::ShPtr< Chain> sp_chain( MODEL.GetChain( *forward_itr)->HardCopy());

          // set the chain id
          sp_chain->SetChainID( *rev_itr);

          // add to chains vector
          chains.PushBack( sp_chain);
        }
      }

      // sort
      chains.Sort( ChainLessThan());

      // create a new model
      ProteinModel model_copy( chains);

      // get the alignment between the current model and the native
      alignments = Quality::CreateAlignmentProteinModels( model_copy, NATIVE);

      // get coordinates vector
      coords = Quality::CoordinatesFromAlignments( alignments, atom_type);

      // return whichever model has the best quality
      return ( *QUALITY)->CalculateMeasure( coords.First(), coords.Second()) < quality ? model_copy : MODEL;
    }

  } // namespace assemble

} // namespace bcl
