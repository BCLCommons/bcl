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
#include "assemble/bcl_assemble_protein_ensemble.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_printer_protein_model.h"
#include "assemble/bcl_assemble_printer_protein_model_ensemble.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "assemble/bcl_assemble_protein_storage_file.h"
#include "assemble/bcl_assemble_quality.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PrinterProteinModelEnsemble::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterProteinModelEnsemble())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterProteinModelEnsemble::PrinterProteinModelEnsemble() :
      m_Prefix( ""),
      m_Superimpose( quality::GetSuperimposeMeasures().e_NoSuperimpose),
      m_RoundNumber( 0),
      m_StageNumber( util::GetUndefined< size_t>()),
      m_Storage( ProteinStorageFile::GetDefaultStorage())
    {
    }

    //! @brief construct with prefix
    //! @param PREFIX prefix string
    PrinterProteinModelEnsemble::PrinterProteinModelEnsemble( const std::string &PREFIX) :
      m_Prefix( PREFIX),
      m_Superimpose( quality::GetSuperimposeMeasures().e_NoSuperimpose),
      m_RoundNumber( 0),
      m_StageNumber( util::GetUndefined< size_t>()),
      m_Storage( ProteinStorageFile::GetDefaultStorage())
    {
    }

    //! @brief construct with all member variables
    //! @param PREFIX prefix string
    //! @param STORAGE protein storage to use
    //! @param SUPERIMPOSE measure to use for superimposition
    PrinterProteinModelEnsemble::PrinterProteinModelEnsemble
    (
      const std::string &PREFIX,
      const util::ShPtr< ProteinStorageFile> &STORAGE,
      const quality::SuperimposeMeasure &SUPERIMPOSE
    ) :
      m_Prefix( PREFIX),
      m_Superimpose( SUPERIMPOSE),
      m_RoundNumber( 0),
      m_StageNumber( util::GetUndefined< size_t>()),
      m_Storage( STORAGE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterProteinModelEnsemble
    PrinterProteinModelEnsemble *PrinterProteinModelEnsemble::Clone() const
    {
      return new PrinterProteinModelEnsemble( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterProteinModelEnsemble::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reset and initialize the printer
    //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
    //! @return true if initialization was successful
    void PrinterProteinModelEnsemble::Initialize( const size_t &ROUND_NUMBER)
    {
      // call initialize with undefined stage number
      Initialize( m_RoundNumber, util::GetUndefined< size_t>());
    }

    //! @brief reset and initialize the printer with the given round and stage number
    //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
    //! @param STAGE_NUMBER for multiple optimizations, a different stage number will be passed
    //! @return true if initialization was successful
    void PrinterProteinModelEnsemble::Initialize( const size_t &ROUND_NUMBER, const size_t &STAGE_NUMBER)
    {
      // update round number
      m_RoundNumber = ROUND_NUMBER;

      // set the stage number
      m_StageNumber = STAGE_NUMBER;
    }

    //! @brief print function taking a tracker
    //! @param TRACKER the tracker to print
    void PrinterProteinModelEnsemble::Print( const opti::Tracker< ProteinModel, double> &TRACKER) const
    {
      WriteToStorage
      (
        TRACKER.GetCurrent()->First(),
        PrinterProteinModel::GetStorageStrings( m_Prefix, TRACKER.GetTag(), m_RoundNumber, m_StageNumber)
      );
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterProteinModelEnsemble::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterProteinModelEnsemble::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief function to write given model
    //! @param ARGUMENT protein model to print,
    //! @param STORAGE_STRINGS source and key strings from the minimization information
    //! @return whether write was successful
    bool PrinterProteinModelEnsemble::WriteToStorage
    (
      const ProteinModel &ARGUMENT,
      const storage::VectorND< 2, std::string> &STORAGE_STRINGS
    ) const
    {
      // get the conformational ensemble of the protein
      const ProteinEnsemble &ensemble( ARGUMENT.GetConformationalEnsemble());
      BCL_MessageDbg
      (
        "ensemble size is " + util::Format()( ensemble.GetSize())
      );

      for
      (
        ProteinEnsemble::const_iterator ensemble_itr( ensemble.Begin()), ensemble_itr_end( ensemble.End());
        ensemble_itr != ensemble_itr_end;
        ++ensemble_itr
      )
      {
        // copy the model
        util::ShPtr< ProteinModel> copy( ( *ensemble_itr)->HardCopy());

        // get the number identifier for this model
        const size_t conformation_number( ensemble_itr - ensemble.Begin());

        // model id string based on its number in the ensemble
        const std::string model_num_id( GetRoundNumberFormat()( conformation_number));

        WriteModel( *copy, model_num_id, STORAGE_STRINGS);
      }

      {
        // copy the model
        util::ShPtr< ProteinModel> copy( ARGUMENT.HardCopy());
        WriteModel( *copy, GetRoundNumberFormat()( ensemble.GetSize()), STORAGE_STRINGS);
      }

      // end
      return true;
    }

    //! @brief writes a single model to file
    //! @param MODEL the model to write
    //! @param CONFORMATION_NUMBER the identifying number string
    //! @param STORAGE_STRINGS the other strings for storage
    void PrinterProteinModelEnsemble::WriteModel
    (
      const ProteinModel &MODEL, const std::string &CONFORMATION_NUMBER,
      const storage::VectorND< 2, std::string> &STORAGE_STRINGS
    ) const
    {
      util::ShPtr< ProteinModel> model_copy( MODEL.HardCopy());

      // initialize transformation matrix for returning multiplier to original state
      math::TransformationMatrix3D multiplier_transform;

      // cast a pointer to the multiplier data if any
      util::ShPtr< ProteinModelMultiplier> sp_multiplier
      (
        model_copy->GetProteinModelData()->GetData( ProteinModelData::e_Multiplier)
      );

      // superimpose if specified
      if( m_Superimpose.IsDefined() && m_Superimpose != quality::GetSuperimposeMeasures().e_NoSuperimpose)
      {
        // superimpose the model and store the transformation to be used for the multiplier
        const math::TransformationMatrix3D transform
        (
          Quality::SuperimposeModel( m_Superimpose, *model_copy, biol::GetAtomTypes().CA).Second()
        );

        // if the pointer is defined
        if( sp_multiplier.IsDefined())
        {
          // set the transformation
          sp_multiplier->Transform( transform);
          multiplier_transform = math::Inverse( transform);
        }
      }

      // store the model
      const std::string filename // can get from storage object
      (
        m_Storage->GetInitializer() + "/" + STORAGE_STRINGS.First() + "_" + CONFORMATION_NUMBER + STORAGE_STRINGS.Second()
         + ".pdb"
      );
      if( m_Storage->Store( *model_copy, STORAGE_STRINGS.First() + "_" + CONFORMATION_NUMBER, STORAGE_STRINGS.Second()))
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

      // if the multiplier pointer is defined
      if( sp_multiplier.IsDefined())
      {
        // move it back to where it was originally
        sp_multiplier->Transform( multiplier_transform);
      }
    }

  } // namespace assemble
} // namespace bcl
