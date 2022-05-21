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
#include "score/bcl_score_phi_psi_with_sspred.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "math/bcl_math_bicubic_spline.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PhiPsiWithSSPred::s_Instance
    (
      GetObjectInstances().AddInstance( new PhiPsiWithSSPred())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &PhiPsiWithSSPred::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "phipsi_sspred");

      // end
      return s_default_scheme;
    }

    //! @brief get the name of the object
    //! @return the name of the object
    const std::string &PhiPsiWithSSPred::GetAlias() const
    {
      static const std::string s_name( "PhiPsiWithSSPred");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PhiPsiWithSSPred::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "scores the phi and psi angles of an SSE taking into account SS predictions");
      serializer.AddInitializer
      (
        "histogram filename",
        "path to where the statistics and energy potentials are read from",
        io::Serialization::GetAgent( &m_HistogramFileName)
       );
      serializer.AddInitializer
      (
        "ss methods",
        "set of SS methods to use in evaluation",
        io::Serialization::GetAgent( &m_SSMethods)
       );

      return serializer;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PhiPsiWithSSPred::PhiPsiWithSSPred() :
      m_Scheme( GetDefaultScheme()),
      m_HistogramFileName(),
      m_EnergyMap(),
      m_SSMethods()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param HISTOGRAM_FILENAME path to file where statistics and in consequence the energy potentials are read from
    //! @param SS_METHODS set of SSMethods to use in evaluation
    //! @param SCHEME scheme to be used in outputting
    PhiPsiWithSSPred::PhiPsiWithSSPred
    (
      const storage::Set< sspred::Method> &SS_METHODS,
      const std::string &SCHEME,
      const std::string &HISTOGRAM_FILENAME
    ) :
      m_Scheme( SCHEME),
      m_HistogramFileName( HISTOGRAM_FILENAME),
      m_EnergyMap(),
      m_SSMethods( SS_METHODS)
    {
      m_EnergyMap = PhiPsi( "", m_HistogramFileName).GetEnergyFunctions();
    }

    //! @brief Clone function
    //! @return pointer to new PhiPsiWithSSPred
    PhiPsiWithSSPred *PhiPsiWithSSPred::Clone() const
    {
      return new PhiPsiWithSSPred( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PhiPsiWithSSPred::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read from this class
    //! @param ERROR_STREAM stream with which to write errors
    bool PhiPsiWithSSPred::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      m_EnergyMap = PhiPsi( "", m_HistogramFileName).GetEnergyFunctions();
      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that calculates the score for a given SSE
    //! @param SSE SSE of interest
    //! @param MEMBRANE membrane object
    //! @return score calculated for the given SSE
    storage::Pair< double, size_t> PhiPsiWithSSPred::operator()
    (
      const assemble::SSE &SSE,
      const biol::Membrane &MEMBRANE
    ) const
    {
      // make sure the ssmethods set is not empty
      if( m_SSMethods.IsEmpty())
      {
        static bool s_have_warned( false);
        if( !s_have_warned)
        {
          BCL_MessageStd( "No SS prediction methods set for PhiPsiWithSSPred");
          s_have_warned = true;
        }
        return storage::Pair< double, size_t>( 0.0, 0);
      }

      // ensure that the amino acid chain length is at least three
      if( SSE.GetSize() < 3)
      {
        return storage::Pair< double, size_t>( double( 0.0), 0);
      }

      // varables for total score of sse and number residues scored
      double scoresum( 0.0);
      size_t scored_entities( 0);

      // iterate over amino acids
      for
      (
        biol::AASequence::const_iterator
          aa_itr( SSE.Begin() + 1), aa_itr_end( SSE.End() - 1);
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        BCL_MessageDbg( "aa " + util::Format()( aa_itr - SSE.Begin()));

        // get references to the previous and next residues
        const biol::AABase &previous_aa( **( aa_itr - 1));
        const biol::AABase &next_aa( **( aa_itr + 1));

        // will hold the best score for the current amino acid out of all the ss types
        double best_score( util::GetUndefinedDouble());

        // iterate over the SStypes
        for
        (
          storage::Map
          <
            biol::SSType, storage::Map< biol::AAType, math::BicubicSpline>
          >::const_iterator ss_type_itr( m_EnergyMap.Begin()), ss_type_itr_end( m_EnergyMap.End());
          ss_type_itr != ss_type_itr_end;
          ++ss_type_itr
        )
        {
          // checks whether amino acid type is contained within the energy map.
          const storage::Map< biol::AAType, math::BicubicSpline>::const_iterator spline_itr
          (
            ss_type_itr->second.Find( ( *aa_itr)->GetType())
          );

          // make sure a bicubic spline exists for the selected amino acid type
          if( spline_itr == ss_type_itr->second.End())
          {
            BCL_MessageStd
            (
              "There is no spline stored for aa type " + ( *aa_itr)->GetType().GetName()
            )
            continue;
          }

          // get the phi psi score for the current residue for the current ss type
          double sstype_phi_psi_score( PhiPsi::ScoreAAPhiPsi( spline_itr->second, **aa_itr, previous_aa, next_aa));

          // if the phi psi score is not defined go to next ss type
          if( !util::IsDefined( sstype_phi_psi_score))
          {
            BCL_MessageDbg( "sstype_phi_psi_score is undefined");
            continue;
          }

          // will hold the sum of all phi psi score weighted with each ss prediction for th current aa and ss type
          double weighted_scoresum( 0.0);

          // iterate over the ss prediction methods
          for
          (
            storage::Set< sspred::Method>::const_iterator
              method_itr( m_SSMethods.Begin()), method_itr_end( m_SSMethods.End());
            method_itr != method_itr_end; ++method_itr
          )
          {
            // store the predicted one state
            const util::SiPtr< const sspred::MethodInterface> this_prediction( ( *aa_itr)->GetSSPrediction( *method_itr));

            // TODO: the SS-predictions are not ordinarily stored on the AAs that are in the SSE since these change during
            //       folding; instead, they're stored on the main sequence. This code presupposes that they're stored on the
            //       SSE too, which is generally not true.
            // check that ss prediction is defined
            if( !this_prediction.IsDefined())
            {
              BCL_MessageVrb
              (
                "SS Prediction is not defined " + method_itr->GetName()
                + " for residue " + ( *aa_itr)->GetIdentification()
                + " PhiPsiWithSSPred will be 0!"
              );
            }
            else
            {

              // square root of the prediction for the current method and current ss type
              const double prediction_value( this_prediction->GetThreeStatePrediction()( ss_type_itr->first));

              // add the current weighted score to the weighted score sum
              weighted_scoresum += prediction_value * sstype_phi_psi_score;
            }
          }

          // true if the weighted score sum for the current ss type is better than for the other ss types
          if( !util::IsDefined( best_score))
          {
            best_score = weighted_scoresum;
          }
          else if( weighted_scoresum < best_score)
          {
            // best score is now the score for obtained for the current ss type
            best_score = weighted_scoresum;
          }
        }

        // add the score for the current residue to the total scoresum of the amino acid sequence
        if( util::IsDefined( best_score))
        {
          scoresum += best_score;
        }

        // increase the number of scored residues
        ++scored_entities;
      }

      // normalize by the number of methods
      scoresum /= double( m_SSMethods.GetSize());

      // return the result of score and number scored residues
      return storage::Pair< double, size_t>( scoresum, scored_entities);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PhiPsiWithSSPred::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Scheme, ISTREAM);
      io::Serialize::Read( m_HistogramFileName, ISTREAM);
      io::Serialize::Read( m_SSMethods, ISTREAM);

      m_EnergyMap = PhiPsi( "", m_HistogramFileName).GetEnergyFunctions();

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PhiPsiWithSSPred::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HistogramFileName, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SSMethods, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
