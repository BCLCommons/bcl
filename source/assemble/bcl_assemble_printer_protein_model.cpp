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
#include "assemble/bcl_assemble_printer_protein_model.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "assemble/bcl_assemble_protein_storage_file.h"
#include "assemble/bcl_assemble_quality.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PrinterProteinModel::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterProteinModel())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterProteinModel::PrinterProteinModel() :
      m_Prefix( ""),
      m_Superimpose( quality::GetSuperimposeMeasures().e_NoSuperimpose),
      m_RoundNumber( 0),
      m_StageNumber( util::GetUndefined< size_t>()),
      m_Storage( ProteinStorageFile::GetDefaultStorage())
    {
    }

    //! @brief construct with prefix
    //! @param PREFIX prefix string
    PrinterProteinModel::PrinterProteinModel( const std::string &PREFIX) :
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
    PrinterProteinModel::PrinterProteinModel
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
    //! @return pointer to new PrinterProteinModel
    PrinterProteinModel *PrinterProteinModel::Clone() const
    {
      return new PrinterProteinModel( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &PrinterProteinModel::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reset and initialize the printer
    //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
    //! @return true if initialization was successful
    void PrinterProteinModel::Initialize( const size_t &ROUND_NUMBER)
    {
      // call initialize with undefined stage number
      Initialize( ROUND_NUMBER, util::GetUndefined< size_t>());
    }

    //! @brief reset and initialize the printer with the given round and stage number
    //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
    //! @param STAGE_NUMBER for multiple optimizations, a different stage number will be passed
    //! @return true if initialization was successful
    void PrinterProteinModel::Initialize( const size_t &ROUND_NUMBER, const size_t &STAGE_NUMBER)
    {
      // update round number
      m_RoundNumber = ROUND_NUMBER;

      // set the stage number
      m_StageNumber = STAGE_NUMBER;
    }

    //! @brief prints information concerning the approximation process based on the status of the given tracker
    //! @param TRACKER holds the status of the approximation process
    void PrinterProteinModel::Print( const opti::Tracker< ProteinModel, double> &TRACKER) const
    {
      // call WriteToStorage function
      WriteToStorage
      (
        TRACKER.GetPhase() == opti::e_End ? TRACKER.GetBest() : TRACKER.GetCurrent(),
        GetStorageStrings( m_Prefix, TRACKER.GetTag(), m_RoundNumber, m_StageNumber)
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterProteinModel::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &PrinterProteinModel::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief gets the source and key strings from the minimization information
    //! @param PREFIX prefix
    //! @param TAG tag for the step status if applicable
    //! @param ROUND_NUMBER round number
    //! @param STAGE_NUMBER stage number
    //! @return source and key strings, respectively
    storage::VectorND< 2, std::string> PrinterProteinModel::GetStorageStrings
    (
      const std::string &PREFIX,
      const std::string &TAG,
      const size_t ROUND_NUMBER,
      const size_t STAGE_NUMBER
    )
    {
      // initialize strings
      std::string tag( TAG);
      if( tag.compare( "") != 0)
      {
        tag = "_" + tag;
      }
      storage::VectorND< 2, std::string> strings( PREFIX + TAG, "_" + GetRoundNumberFormat()( ROUND_NUMBER));

      // if the stage is defined
      if( util::IsDefined( STAGE_NUMBER))
      {
        strings.Second() += "_" + GetStageNumberFormat()( STAGE_NUMBER);
      }

      // end
      return strings;
    }

    //! @brief function to write given model
    //! @param ARGUMENT_RESULT_PAIR Pair of argument and corresponding result
    //! @param STORAGE_STRINGS source and key strings from the minimization information
    //! @return whether write was successful
    bool PrinterProteinModel::WriteToStorage
    (
      const util::ShPtr< storage::Pair< ProteinModel, double> > &ARGUMENT_RESULT_PAIR,
      const storage::VectorND< 2, std::string> &STORAGE_STRINGS
    ) const
    {
      // copy the model
      util::ShPtr< ProteinModel> copy_ptr( ARGUMENT_RESULT_PAIR->First().HardCopy());
      ProteinModel &copy( *copy_ptr);

      // initialize transformation matrix for returning multiplier to original state
      math::TransformationMatrix3D orig_transform;

      // cast a pointer to the multiplier data if any
      util::ShPtr< ProteinModelMultiplier> sp_multiplier
      (
        copy.GetProteinModelData()->GetData( ProteinModelData::e_Multiplier)
      );

      // cast a pointer to the membrane data if any
      util::ShPtr< biol::Membrane> sp_membrane
      (
        copy.GetProteinModelData()->GetData( ProteinModelData::e_Membrane)
      );

      // superimpose if specified
      if( m_Superimpose.IsDefined() && m_Superimpose != quality::GetSuperimposeMeasures().e_NoSuperimpose)
      {
        // superimpose the model and store the transformation to be used for the multiplier
        const math::TransformationMatrix3D transform
        (
          Quality::SuperimposeModel( m_Superimpose, copy, biol::GetAtomTypes().CA).Second()
        );

        // if the pointer is defined
        if( sp_multiplier.IsDefined())
        {
          // set the transformation
          sp_multiplier->Transform( transform);
          orig_transform = math::Inverse( transform);
        }

        // if the pointer is defined
        if( sp_membrane.IsDefined())
        {
          // set the transformation
          sp_membrane->Transform( transform);
          orig_transform = math::Inverse( transform);
        }
      }

      // store the model
      const std::string filename // can get from storage object
      (
        m_Storage->GetInitializer() + "/" + STORAGE_STRINGS.First() + STORAGE_STRINGS.Second() + ".pdb"
      );

      if( m_Storage->Store( copy, STORAGE_STRINGS.First(), STORAGE_STRINGS.Second()))
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
        sp_multiplier->Transform( orig_transform);
      }

      // if the membrane pointer is defined
      if( sp_membrane.IsDefined())
      {
        // move it back to where it was originally
        sp_membrane->Transform( orig_transform);
      }

      // end
      return true;
    }

  } // namespace assemble
} // namespace bcl
