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
#include "sspred/bcl_sspred_jufo9d.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "descriptor/bcl_descriptor_base.h"
#include "descriptor/bcl_descriptor_iterator.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "math/bcl_math_running_average.h"
#include "model/bcl_model_interface_retrieve_from_file.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {

  //////////
  // data //
  //////////

    //! initialize base probability
    const double JUFO9D::s_BaseProbability( 1.0 / 9.0);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    JUFO9D::JUFO9D() :
      m_Prediction( GetDefaultPredictionMatrix())
    {
    }

    //! @brief constructor from a vector of 9 values
    //! @param VECTOR vector of values
    JUFO9D::JUFO9D( const linal::Vector< double> &VECTOR) :
      m_Prediction( 3, 3, VECTOR.Begin())
    {
    }

    //! @brief constructor from a linal::Matrix
    //! @param MATRIX linal::Matrix of probabilities
    JUFO9D::JUFO9D( const linal::Matrix< double> &MATRIX) :
      m_Prediction( MATRIX)
    {
      // check the size
      BCL_Assert
      (
        MATRIX.GetNumberRows() == 3 &&
        MATRIX.GetNumberCols() == 3,
        "The matrix should be 3 by 3 instead it is  " +
        util::Format()( MATRIX.GetNumberRows()) + " by " + util::Format()( MATRIX.GetNumberCols())
      )
    }

    //! @brief Clone function
    //! @return pointer to new JUFO
    JUFO9D *JUFO9D::Clone() const
    {
      return new JUFO9D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &JUFO9D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get file extension associated with this Method
    //! @return file extension associated with this Method
    const std::string &JUFO9D::GetFileExtension() const
    {
      // static extension string
      static const std::string s_file_extension( ".jufo9d");

      // end
      return s_file_extension;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief probability that this is a TM-span of type given
    //! @return probability that this is a TM-span of type given
    double JUFO9D::TMTypeProbability( const biol::SSType &SS_TYPE) const
    {
      BCL_Assert( SS_TYPE <= biol::GetSSTypes().COIL, "No prediction data stored for type " + SS_TYPE.GetName());

      // get the tm-type prediction
      const double prob_tm_type
      (
        m_Prediction( biol::GetEnvironmentTypes().e_MembraneCore->GetReducedIndex(), SS_TYPE)
      );

      // end
      return prob_tm_type + ( 7.0 * s_BaseProbability / 2.0);
    }

    //! @brief get three state, environment independent secondary structure prediction
    //! @return three state, environment independent secondary structure prediction
    linal::Vector3D JUFO9D::GetThreeStatePrediction() const
    {
      return ConvertNineStateToThreeState( m_Prediction);
    }

    //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
    //! @return three state secondary structure prediction
    linal::Matrix< double> JUFO9D::GetNineStatePrediction() const
    {
      return m_Prediction;
    }

    //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AMINO_ACID amino acid into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &JUFO9D::ReadPredictionsForAA
    (
      std::istream &ISTREAM,
      biol::AABase &AMINO_ACID
    ) const
    {
      // initialize matrix of predictions
      linal::Matrix< double> prediction( 3, 3, double( 0.0));

      // temporary variable to store one state
      std::string one_state;

      // read one state ss prediction and one state tm prediction
      ISTREAM >> one_state >> one_state;

      // read COIL, HELIX and STRAND
      ISTREAM >> prediction( biol::GetEnvironmentTypes().e_MembraneCore->GetReducedIndex(), biol::GetSSTypes().COIL);
      ISTREAM >> prediction( biol::GetEnvironmentTypes().e_Transition->GetReducedIndex()  , biol::GetSSTypes().COIL);
      ISTREAM >> prediction( biol::GetEnvironmentTypes().e_Solution->GetReducedIndex()    , biol::GetSSTypes().COIL);
      ISTREAM >> prediction( biol::GetEnvironmentTypes().e_MembraneCore->GetReducedIndex(), biol::GetSSTypes().HELIX);
      ISTREAM >> prediction( biol::GetEnvironmentTypes().e_Transition->GetReducedIndex()  , biol::GetSSTypes().HELIX);
      ISTREAM >> prediction( biol::GetEnvironmentTypes().e_Solution->GetReducedIndex()    , biol::GetSSTypes().HELIX);
      ISTREAM >> prediction( biol::GetEnvironmentTypes().e_MembraneCore->GetReducedIndex(), biol::GetSSTypes().STRAND);
      ISTREAM >> prediction( biol::GetEnvironmentTypes().e_Transition->GetReducedIndex()  , biol::GetSSTypes().STRAND);
      ISTREAM >> prediction( biol::GetEnvironmentTypes().e_Solution->GetReducedIndex()    , biol::GetSSTypes().STRAND);
      // read characters until the end of the line
      std::getline( ISTREAM, one_state);

      // normalize the matrix
      prediction.AsVector().SetToSum( 1.0);

      // set the predictions for this amino acid
      AMINO_ACID.SetSSPrediction( GetMethods().e_JUFO9D, JUFO9D( prediction));

      // end
      return ISTREAM;
    }

    //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
    //! @param ISTREAM input stream
    //! @param AA_SEQUENCE AASequence into which sspredictions will be read
    //! @return std::istream which was read from
    std::istream &JUFO9D::ReadPredictionsForAASequence
    (
      std::istream &ISTREAM,
      biol::AASequence &AA_SEQUENCE
    ) const
    {
      // call standard read function and return it
      return ReadStandardPredictionsForAASequence( ISTREAM, AA_SEQUENCE, GetMethods().e_JUFO9D);
    }

    //! @brief helper function to retrieve the jufo9d descriptors
    //! @param MULTIMER true if the sequence is multimeric
    util::Implementation< descriptor::Base< biol::AABase, float> > JUFO9D::GetJufo9DANNDescriptors( const bool &MULTIMER)
    {
      return
        util::Implementation< descriptor::Base< biol::AABase, float> >
        (
          "Combine("
          "  DefineNaN("
          "    Window("
          "    Combine("
          "      AA_StericalParameter, AA_Polarizability, AA_Volume, AA_IsoelectricPoint, "
          "      AA_SASA, AA_FreeEnergyHelix, AA_FreeEnergyStrand, AA_FreeEnergyCoil, "
          "      AA_FreeEnergyCore, AA_FreeEnergyTransition, AA_FreeEnergySolution, "
          "      AA_FreeEnergyCoreHelix, AA_FreeEnergyTransitionHelix, AA_FreeEnergySolutionHelix, "
          "      AA_FreeEnergyCoreStrand,  AA_FreeEnergyTransitionStrand, AA_FreeEnergySolutionStrand,"
          "      AA_FreeEnergyCoreCoil, AA_FreeEnergyTransitionCoil, AA_FreeEnergySolutionCoil"
          "    ),"
          "    size=15,"
          "    alignment=JufoCenter"
          "    ),"
          "    replacement=Constant(2.165,0.147,3.33,6.141,263.191,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)"
          "  ),"
          "  DefineNaN(Window(AABlastProfile,size=15,alignment=JufoCenter),replacement=Constant(0)),"
          "  Constant(" + std::string( MULTIMER ? "1" : "0" ) +"),                          "
          "  NElements,                                                                     "
          "  SequenceMean(DefineNaN(AA_StericalParameter,replacement=Constant(2.165))),     "
          "  SequenceMean(DefineNaN(AA_Polarizability,replacement=Constant(0.147))),        "
          "  SequenceMean(DefineNaN(AA_Volume,replacement=Constant(3.33))),                 "
          "  SequenceMean(DefineNaN(AA_IsoelectricPoint,replacement=Constant(6.141))),      "
          "  SequenceMean(DefineNaN(AA_SASA,replacement=Constant(263.191))),                "
          "  SequenceMean(DefineNaN(AA_FreeEnergyHelix,replacement=Constant(0))),           "
          "  SequenceMean(DefineNaN(AA_FreeEnergyStrand,replacement=Constant(0))),          "
          "  SequenceMean(DefineNaN(AA_FreeEnergyCoil,replacement=Constant(0))),            "
          "  SequenceMean(DefineNaN(AA_FreeEnergyCore,replacement=Constant(0))),            "
          "  SequenceMean(DefineNaN(AA_FreeEnergyTransition,replacement=Constant(0))),      "
          "  SequenceMean(DefineNaN(AA_FreeEnergySolution,replacement=Constant(0))),        "
          "  SequenceMean(DefineNaN(AA_FreeEnergyCoreHelix,replacement=Constant(0))),       "
          "  SequenceMean(DefineNaN(AA_FreeEnergyTransitionHelix,replacement=Constant(0))), "
          "  SequenceMean(DefineNaN(AA_FreeEnergySolutionHelix,replacement=Constant(0))),   "
          "  SequenceMean(DefineNaN(AA_FreeEnergyCoreStrand,replacement=Constant(0))),      "
          "  SequenceMean(DefineNaN(AA_FreeEnergyTransitionStrand,replacement=Constant(0))),"
          "  SequenceMean(DefineNaN(AA_FreeEnergySolutionStrand,replacement=Constant(0))),  "
          "  SequenceMean(DefineNaN(AA_FreeEnergyCoreCoil,replacement=Constant(0))),        "
          "  SequenceMean(DefineNaN(AA_FreeEnergyTransitionCoil,replacement=Constant(0))),  "
          "  SequenceMean(DefineNaN(AA_FreeEnergySolutionCoil,replacement=Constant(0))),    "
          "  SequenceMean(DefineNaN(AABlastProfile,replacement=Constant(0)))                "
          ")"
        );
    }

    //! @brief iterates over the sequences in ProteinModel and calculates the jufo predictions for every residue in the sequence
    //! @param PROTEIN_MODEL ProteinModel for which JUFO will be calculated
    //! @param MULTIMER true if the sequence is multimeric
    void JUFO9D::Calculate( assemble::ProteinModel &PROTEIN_MODEL, const bool MULTIMER)
    {
      // create a protein-model-with-cache
      assemble::ProteinModelWithCache pmwc( PROTEIN_MODEL, false);
      Calculate( pmwc, MULTIMER);
    }

    //! @brief iterates over the sequences in ProteinModel and calculates the jufo predictions for every residue in the sequence
    //! @param PROTEIN_MODEL ProteinModel for which JUFO will be calculated
    //! @param MULTIMER true if the sequence is multimeric
    void JUFO9D::Calculate( assemble::ProteinModelWithCache &PROTEIN_MODEL, const bool MULTIMER)
    {
      // make sure the blast profile exists for this sequence
      BCL_Assert
      (
        PROTEIN_MODEL.GetIterator()->GetBlastProfilePtr().IsDefined(),
        "Blast Profile is not available"
      );

      // create the descriptor to generate the jufo dataset
      util::Implementation< descriptor::Base< biol::AABase, float> > aa_descriptor( GetJufo9DANNDescriptors( MULTIMER));

      // set the object up
      aa_descriptor->SetObject( PROTEIN_MODEL);

      // set the dimension (1 because we operate on elements of the sequence)
      aa_descriptor->SetDimension( 1);

      // create a descriptor iterator
      descriptor::Iterator< biol::AABase> itr( PROTEIN_MODEL.GetIterator());

      // when jufo 9d was trained, columns (descriptors) that had no range were given a value of 0.5.  Now such columns
      // are given a value of 0; but the trained networks rely on the old behavior
      model::RescaleFeatureDataSet::GetDefaultValueForEmptyRangedColumns() = 0.5;

      // retrieve the ANNs
      const model::InterfaceRetrieveFromFile storage( model::Model::AddModelPath( "jufo"), "jufo9d_");
      const model::RetrieveInterface::t_Container first_layer_anns( storage.RetrieveEnsemble());
      BCL_Assert
      (
        !first_layer_anns.IsEmpty(),
        "No JUFO9D models found! Files should be in " + model::Model::AddModelPath( "jufo") + " with a prefix, jufo9d_"
      );

      // create a non-const iterator over the same amino acids
      iterate::Generic< biol::AABase> itr_non_const( PROTEIN_MODEL.GetIteratorNonConst());

      // iterate over the amino acids
      for( ; itr.NotAtEnd(); ++itr, ++itr_non_const)
      {
        math::RunningAverage< linal::Vector< float> > average_prediction;
        linal::Vector< float> descriptor( aa_descriptor->GetSizeOfFeatures(), ( *aa_descriptor)( itr).Begin());
        const model::FeatureDataSet< float> descriptor_data
        (
          linal::Matrix< float>( 1, descriptor.GetSize(), descriptor.Begin())
        );
        // create predictions from all networks
        for
        (
          model::RetrieveInterface::t_Container::const_iterator
            itr_model( first_layer_anns.Begin()), itr_model_end( first_layer_anns.End());
          itr_model != itr_model_end;
          ++itr_model
        )
        {
          average_prediction += ( **itr_model)( descriptor_data).GetMatrix().GetRow( 0);
        }

        // get the average
        linal::Vector< float> prediction_ave( average_prediction.GetAverage());
        prediction_ave += float( 0.15);

        // normalize, then store it in ssprediction
        prediction_ave.SetToSum( 1);

        // get a pointer to the amino acid

        // set the prediction
        itr_non_const->SetSSPrediction( GetMethods().e_JUFO9D, JUFO9D( prediction_ave));
      }
      model::RescaleFeatureDataSet::GetDefaultValueForEmptyRangedColumns() = util::GetUndefined< float>();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator = for assigning given JUFO_RHS to this object
    //! @param JUFO_RHS JUFO to be assigned
    //! @return this object after assignment
    JUFO9D &JUFO9D::operator =( const JUFO9D &JUFO_RHS)
    {
      // update members
      m_Prediction = JUFO_RHS.m_Prediction;

      // end
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &JUFO9D::Read( std::istream &ISTREAM)
    {
      // read members
      ISTREAM >> m_Prediction;

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &JUFO9D::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_Prediction;

      // end
      return OSTREAM;
    };

    //! @brief write two state TM-helix prediction from sequence to stream
    //! @param OSTREAM stream to write to
    //! @param AA_SEQUENCE sequence containing predictions
    //! @param SS_TYPE sstype for the 2 state prediction, must be helix or strand
    //! @return output stream
    std::ostream &JUFO9D::WriteTwoStateTMPredictions
    (
      std::ostream &OSTREAM,
      const biol::AASequence &AA_SEQUENCE,
      const biol::SSType &SS_TYPE
    )
    {
      // check the ss type
      if( SS_TYPE != biol::GetSSTypes().HELIX && SS_TYPE != biol::GetSSTypes().STRAND)
      {
        // don't write anything
        return OSTREAM;
      }

      // iterate over amino acids in the sequence
      for
      (
        biol::AASequence::const_iterator aa_itr( AA_SEQUENCE.Begin()), aa_itr_end( AA_SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // check existence of requested method
        const util::SiPtr< const JUFO9D> sp_prediction( ( *aa_itr)->GetSSPrediction( GetMethods().e_JUFO9D));
        if( !sp_prediction.IsDefined())
        {
          BCL_MessageCrt
          (
            "secondary structure prediction method " +
            util::Format()( GetMethods().e_JUFO9D) + " is not stored for this residue"
          );
        }

        // output the residue seqid and the residue one letter code
        OSTREAM << util::Format().W( 4)( ( *aa_itr)->GetSeqID()) << ' '
               << ( *aa_itr)->GetType()->GetOneLetterCode() << ' ';

        // write out 2 state prediction
        const double tm_prob( sp_prediction->TMTypeProbability( SS_TYPE));
        const double other_prob( 1.0 - tm_prob);
        util::Format prob_format;
        prob_format.W( 7).FFP( 3);
        if( tm_prob >= 0.5)
        {
          OSTREAM << SS_TYPE->GetOneLetterCode() << " MC ";
        }
        else
        {
          OSTREAM << "- -- ";
        }
        OSTREAM << prob_format( tm_prob) << prob_format( other_prob) << '\n';
      }

      // end
      return OSTREAM;
    }

    //! @brief write three state TM prediction from sequence to stream
    //! @param OSTREAM stream to write to
    //! @param AA_SEQUENCE sequence containing predictions
    //! @return output stream
    std::ostream &JUFO9D::WriteThreeStateTMPredictions( std::ostream &OSTREAM, const biol::AASequence &AA_SEQUENCE)
    {
      // iterate over amino acids in the sequence
      for
      (
        biol::AASequence::const_iterator aa_itr( AA_SEQUENCE.Begin()), aa_itr_end( AA_SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // check existence of requested method
        const util::SiPtr< const JUFO9D> sp_prediction( ( *aa_itr)->GetSSPrediction( GetMethods().e_JUFO9D));
        if( !sp_prediction.IsDefined())
        {
          BCL_MessageCrt
          (
            "secondary structure prediction method " +
            util::Format()( GetMethods().e_JUFO9D) + " is not stored for this residue"
          );
        }

        // output the residue seqid and the residue one letter code
        OSTREAM << util::Format().W( 4)( ( *aa_itr)->GetSeqID()) << ' '
               << ( *aa_itr)->GetType()->GetOneLetterCode() << ' ';

        // get three state predictions
        const linal::Vector3D tm_prediction( sp_prediction->GetThreeStateTMPrediction());

        // write data
        OSTREAM << sp_prediction->GetOneStateTMPrediction()->GetTwoLetterCode() << ' '
               << util::Format().W( 7).FFP( 3)( tm_prediction( biol::GetEnvironmentTypes().e_MembraneCore->GetReducedIndex()))
               << util::Format().W( 7).FFP( 3)( tm_prediction( biol::GetEnvironmentTypes().e_Transition->GetReducedIndex()))
               << util::Format().W( 7).FFP( 3)( tm_prediction( biol::GetEnvironmentTypes().e_Solution->GetReducedIndex()))
               << '\n';
      }

      // end
      return OSTREAM;
    }

  } // namespace sspred
} // namespace bcl
