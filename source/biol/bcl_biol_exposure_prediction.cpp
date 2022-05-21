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
#include "biol/bcl_biol_exposure_prediction.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "descriptor/bcl_descriptor_base.h"
#include "descriptor/bcl_descriptor_iterator.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_running_average.h"
#include "model/bcl_model_interface_retrieve_from_file.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ExposurePrediction::s_Instance
    (
      GetObjectInstances().AddInstance( new ExposurePrediction())
    );

    //! @brief return model path as a string
    //! @param ExposureType one of the exposure types
    //! @return model path as a string
    const std::string &ExposurePrediction::GetModelPath( const ExposureType &EXPOSURE_TYPE)
    {
      static std::string s_paths[] =
      {
          "exposure/contact_number/protomeric",
          "exposure/contact_number/oligomeric",
          "exposure/rsa/protomeric",
          "exposure/rsa/oligomeric",
          GetStaticClassName< ExposurePrediction::ExposureType>()
      };
      return s_paths[ size_t( EXPOSURE_TYPE)];
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &ExposurePrediction::GetFileExtension( const ExposureType &EXPOSURE_TYPE)
    {
      static const std::string s_extensions[] =
      {
          ".proto_cn",
          ".oligo_cn",
          ".proto_rsa",
          ".oligo_rsa",
          GetStaticClassName< ExposurePrediction::ExposureType>()
      };
      return s_extensions[ size_t( EXPOSURE_TYPE)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ExposurePrediction::ExposurePrediction()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ExposurePrediction
    ExposurePrediction *ExposurePrediction::Clone() const
    {
      return new ExposurePrediction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ExposurePrediction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief iterates over the sequences and calculates the jufo predictions for every residue in the sequence
    //! @param SEQUENCE sequence of interest
    //! @param EXPOSURE_TYPE one of the exposure types
    void ExposurePrediction::Calculate( AASequence &SEQUENCE, const ExposureType &EXPOSURE_TYPE)
    {
      // create a dummy model for the prediction
      assemble::ProteinModel model
      (
        util::ShPtr< assemble::Chain>
        (
          new assemble::Chain
          (
            util::ShPtr< AASequence>( new AASequence( SEQUENCE))
          )
        )
      );
      ExposurePrediction::Calculate( model, EXPOSURE_TYPE);
    }

    //! @brief iterates over the sequences in ProteinModel and calculates the jufo predictions for every residue in the sequence
    //! @param PROTEIN_MODEL ProteinModel for which JUFO will be calculated
    //! @param EXPOSURE_TYPE one of the exposure types
    void ExposurePrediction::Calculate( assemble::ProteinModel &PROTEIN_MODEL, const ExposureType &EXPOSURE_TYPE)
    {
      // create a protein-model-with-cache
      assemble::ProteinModelWithCache pmwc( PROTEIN_MODEL, false);
      Calculate( pmwc, EXPOSURE_TYPE);
    }

    //! @brief iterates over the sequences in ProteinModel and calculates the jufo predictions for every residue in the sequence
    //! @param PROTEIN_MODEL ProteinModel for which JUFO will be calculated
    //! @param EXPOSURE_TYPE one of the exposure types
    void ExposurePrediction::Calculate( assemble::ProteinModelWithCache &PROTEIN_MODEL, const ExposureType &EXPOSURE_TYPE)
    {
      // make sure the blast profile exists for this sequence
      BCL_Assert
      (
        PROTEIN_MODEL.GetIterator()->GetBlastProfilePtr().IsDefined(),
        "Blast Profile is not available"
      );

      // path to the model
      const std::string path_to_model( GetModelPath( EXPOSURE_TYPE));

      // create the descriptor to generate the dataset
      util::Implementation< descriptor::Base< AABase, float> > aa_descriptor
      (
        "PredictionMean(storage=File(directory=" + model::Model::AddModelPath( path_to_model) + ", prefix=model))"
      );

      // set the object up
      aa_descriptor->SetObject( PROTEIN_MODEL);

      // set the dimension (1 because we operate on elements of the sequence)
      aa_descriptor->SetDimension( 1);

      // create a descriptor iterator
      descriptor::Iterator< AABase> itr( PROTEIN_MODEL.GetIterator());

      // create a non-const iterator over the same amino acids
      iterate::Generic< AABase> itr_non_const( PROTEIN_MODEL.GetIteratorNonConst());

      // iterate over the amino acids
      for( ; itr.NotAtEnd(); ++itr, ++itr_non_const)
      {
        // set the prediction
        itr_non_const->SetExposurePrediction( aa_descriptor->operator ()( itr)( 0));
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads the exposure predictions for a model given a path and prefix
    //! @param PROTEIN_MODEL protein model to contain predictions
    //! @param PREFIX prefix of the sequence ( usually pdb id)
    //! @param PATH path where exposure prediction files can be found
    void ExposurePrediction::ReadPredictions
    (
      assemble::ProteinModel &PROTEIN_MODEL,
      const std::string &PREFIX,
      const std::string &PATH,
      const ExposureType &EXPOSURE_TYPE
    )
    {
      //iterate over all chains in the sequence
      for
      (
        util::ShPtrVector< assemble::Chain>::iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // construct the filename
        std::string filename
        (
          PATH + PATH_SEPARATOR + PREFIX + ( *chain_itr)->GetChainID() + GetFileExtension( EXPOSURE_TYPE)
        );

        if( !io::DirectoryEntry( filename).DoesExist())
        {
          filename = PATH + PATH_SEPARATOR + PREFIX + GetFileExtension( EXPOSURE_TYPE);
        }
        // open the file
        io::IFStream read;
        io::File::MustOpenIFStream( read, filename);

        // read predictions for the sequence
        ReadPredictions( read, *( *chain_itr)->GetSequence());
        io::File::CloseClearFStream( read);
      }
    }

    //! @brief reads the exposure predictions from a file
    //! @param ISTREAM stream to read from
    //! @param PROTEIN_MODEL protein model to contain predictions
    void ExposurePrediction::ReadPredictions( std::istream &ISTREAM, assemble::ProteinModel &PROTEIN_MODEL)
    {
      //iterate over all chains in the sequence
      for
      (
        util::ShPtrVector< assemble::Chain>::iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // read predictions for the sequence
        ReadPredictions( ISTREAM, *( *chain_itr)->GetSequence());
      }
    }

    //! @brief reads the exposure predictions from a file
    //! @param ISTREAM stream to read from
    //! @param SEQUENCE sequence to contain predictions
    void ExposurePrediction::ReadPredictions( std::istream &ISTREAM, AASequence &SEQUENCE)
    {
      // iterate over the sequence
      for
      (
        AASequence::iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        // read in the seq id and exposure
        int seq_id;
        std::string exposure;
        ISTREAM >> seq_id >> exposure;

        BCL_Assert
        (
          ( *aa_itr)->GetSeqID() == seq_id,
          "Mismatch in seq id for exposure predictions. " + ( *aa_itr)->GetIdentification() + " does not match " +
            util::Format()( seq_id)
        );

        // set the prediction
        ( *aa_itr)->SetExposurePrediction( util::ConvertStringToNumericalValue< double>( exposure));
      }
    }

    //! @brief writes out the exposure predictions to a file
    //! @param OSTREAM stream to write to
    //! @param SEQUENCE sequence containing predictions
    void ExposurePrediction::WritePredictions( std::ostream &OSTREAM, const AASequence &SEQUENCE)
    {
      // iterate over the sequence
      for
      (
        AASequence::const_iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        OSTREAM << ( *aa_itr)->GetSeqID() << '\t' << ( *aa_itr)->GetExposurePrediction() << '\n';
      }
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ExposurePrediction::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ExposurePrediction::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
