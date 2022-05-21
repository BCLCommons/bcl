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
#include "score/bcl_score_sse_pair_contact.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "math/bcl_math_binary_sum_function.h"
#include "score/bcl_score_aa_pair_contact.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SSEPairContact::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEPairContact())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEPairContact::SSEPairContact()
    {
    }

    //! @brief constructor with a prediction map
    //! @param SP_PREDICTION_MAP ShPtr to PredictionMap to be used
    SSEPairContact::SSEPairContact( const util::ShPtr< contact::PredictionMap> &SP_PREDICTION_MAP) :
      m_ScoringFunctions( contact::Types::s_NumberValidTypes)
    {

      // constructor the scoring function for HELIX_HELIX
      m_ScoringFunctions( contact::GetTypes().HELIX_HELIX) =
        util::ShPtr< math::BinaryFunctionInterface< biol::AABase, biol::AABase, double> >
        (
          new AAPairContact( contact::GetTypes().HELIX_HELIX,   SP_PREDICTION_MAP)
        );

      // constructor the scoring function for HELIX_SHEET
      m_ScoringFunctions( contact::GetTypes().HELIX_SHEET) =
        util::ShPtr< math::BinaryFunctionInterface< biol::AABase, biol::AABase, double> >
        (
          new AAPairContact( contact::GetTypes().HELIX_SHEET, SP_PREDICTION_MAP)
        );

      // constructor the scoring function for SHEET_HELIX
      m_ScoringFunctions( contact::GetTypes().SHEET_HELIX) =
        util::ShPtr< math::BinaryFunctionInterface< biol::AABase, biol::AABase, double> >
        (
          new AAPairContact( contact::GetTypes().SHEET_HELIX, SP_PREDICTION_MAP)
        );

      // constructor the scoring function for STRAND_STRAND
      m_ScoringFunctions( contact::GetTypes().STRAND_STRAND) =
        util::ShPtr< math::BinaryFunctionInterface< biol::AABase, biol::AABase, double> >
        (
          new AAPairContact( contact::GetTypes().STRAND_STRAND, SP_PREDICTION_MAP)
        );

      // constructor the scoring function for SHEET_SHEET
      m_ScoringFunctions( contact::GetTypes().SHEET_SHEET) =
        util::ShPtr< math::BinaryFunctionInterface< biol::AABase, biol::AABase, double> >
        (
          new AAPairContact( contact::GetTypes().SHEET_SHEET,   SP_PREDICTION_MAP)
        );
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new SSEPairContact copied from this one
    SSEPairContact *SSEPairContact::Clone() const
    {
      return new SSEPairContact( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPairContact::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &SSEPairContact::GetScheme() const
    {
      // initialize static string to store scheme
      static const std::string s_scheme( "ssecont");
      return s_scheme;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the score for the given SSE pair
    //! @param SSE_A first SSE of interest
    //! @param SSE_B second SSE of interest
    //! @return the score for the given SSE pair
    double SSEPairContact::operator()( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
    {
      BCL_MessageDbg( "score sse pair contact");

      // Retrieve Contacttype and distance_angleX_angleZ of these two SSEs
      const assemble::SSEGeometryPacking ssepack( SSE_A, SSE_B);

      // util::ShPtr to function which will be assigned according to the contact type and the relative position of
      // the two sses to each other
      util::ShPtr
      <
        math::BinaryFunctionInterface< biol::AABase, biol::AABase, double>
      > merged_scoring;

      // store contact type
      contact::Type contact_type( ssepack.GetContactType());

      // switch over contact type
      // if UNKNOWN
      if( contact_type == contact::GetTypes().e_Undefined)
      {
        return 0.0;
      }
      // if HELIX_HELIX, HELIX_SHEET, SHEET_HELIX or STRAND_STRAND
      else if
      (
        contact_type == contact::GetTypes().HELIX_HELIX  ||
        contact_type == contact::GetTypes().HELIX_SHEET ||
        contact_type == contact::GetTypes().SHEET_HELIX ||
        contact_type == contact::GetTypes().STRAND_STRAND
      )
      {
        merged_scoring = m_ScoringFunctions( ssepack.GetContactType());
      }
      // if UNDEFINED_HELIX_STRAND
      else if( contact_type == contact::GetTypes().UNDEFINED_HELIX_STRAND)
      {
        merged_scoring =
          util::ShPtr< math::BinaryFunctionInterface< biol::AABase, biol::AABase, double> >
          (
            (
              ( *m_ScoringFunctions( contact::GetTypes().HELIX_SHEET)) *
              math::WeightBetweenZeroAndPi_ThreeSections( ssepack.GetRelativePosition() * 2)
            ).Clone()
          );
      }
      // if UNDEFINED_STRAND_HELIX
      else if( contact_type == contact::GetTypes().UNDEFINED_STRAND_HELIX)
      {
        merged_scoring =
          util::ShPtr< math::BinaryFunctionInterface< biol::AABase, biol::AABase, double> >
          (
            (
              ( *m_ScoringFunctions( contact::GetTypes().SHEET_HELIX)) *
              math::WeightBetweenZeroAndPi_ThreeSections( ssepack.GetRelativePosition() * 2)
            ).Clone()
          );
      }
      // if SHEET_SHEET
      else if( contact_type == contact::GetTypes().SHEET_SHEET)
      {
        merged_scoring = m_ScoringFunctions( contact::GetTypes().SHEET_SHEET);
      }
      // if UNDEFINED_STRAND_STRAND
      else if( contact_type == contact::GetTypes().UNDEFINED_STRAND_STRAND)
      {
        // calculate the weight
        const double weight( math::WeightBetweenZeroAndPi_ThreeSections( ssepack.GetRelativePosition()));

        merged_scoring =
          util::ShPtr< math::BinaryFunctionInterface< biol::AABase, biol::AABase, double> >
          (
            (
              weight * ( *m_ScoringFunctions( contact::GetTypes().SHEET_SHEET)) +
              ( 1 - weight) * ( *m_ScoringFunctions( contact::GetTypes().STRAND_STRAND))
            ).Clone()
          );
      }
      // if none above
      else
      {
        BCL_MessageCrt( "undefined contacttype has been found");
        return util::GetUndefined< double>();
      }

      double score( 0);

      // iterate over every residue in first SSE against every residue in the other SSE
      for
      (
        biol::AASequence::const_iterator
          aa_itr_a( SSE_A.GetData().Begin()),
          aa_itr_a_end( SSE_A.GetData().End());
        aa_itr_a != aa_itr_a_end; ++aa_itr_a
      )
      {
        for
        (
          biol::AASequence::const_iterator
            aa_itr_b( SSE_B.GetData().Begin()),
            aa_itr_b_end( SSE_B.GetData().End());
          aa_itr_b != aa_itr_b_end; ++aa_itr_b
        )
        {
          // and sum over the scores
          score += merged_scoring->operator()( **aa_itr_a, **aa_itr_b);
        }
      }

      // return the total score times the interaction weight if the packing is not optimal (orthogonal)
      return ssepack.GetInteractionWeight() * score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &SSEPairContact::Read( std::istream &ISTREAM)
    {
      // read members
      ISTREAM >> m_ScoringFunctions;

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &SSEPairContact::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      OSTREAM << m_ScoringFunctions;

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param SSE_A first SSE of interest
    //! @param SSE_B second SSE of interest
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &SSEPairContact::WriteDetailedSchemeAndValues
    (
      const assemble::SSE &SSE_A,
      const assemble::SSE &SSE_B,
      std::ostream &OSTREAM
    ) const
    {
      // Retrieve Contacttype and distance_angleX_angleZ of these two SSEs
      const assemble::SSEGeometryPacking ssepack( SSE_A, SSE_B);

      //util::ShPtr to function which will be assigned according to the contacttype and the relative position of the two sses to each other
      util::ShPtr
      <
        math::BinaryFunctionInterface< biol::AABase, biol::AABase, double>
      > merged_scoring;

      // store contact type
      contact::Type contact_type( ssepack.GetContactType());

      // switch over contact type
      // if UNKNOWN
      if( contact_type == contact::GetTypes().e_Undefined)
      {
        return OSTREAM;
      }
      // if HELIX_HELIX, HELIX_SHEET, SHEET_HELIX or STRAND_STRAND
      else if
      (
        contact_type == contact::GetTypes().HELIX_HELIX  ||
        contact_type == contact::GetTypes().HELIX_SHEET ||
        contact_type == contact::GetTypes().SHEET_HELIX ||
        contact_type == contact::GetTypes().STRAND_STRAND
      )
      {
        merged_scoring = m_ScoringFunctions( ssepack.GetContactType());
      }
      // if UNDEFINED_HELIX_STRAND
      else if( contact_type == contact::GetTypes().UNDEFINED_HELIX_STRAND)
      {
        merged_scoring
          = util::ShPtr< math::BinaryFunctionInterface< biol::AABase, biol::AABase, double> >
            (
              (
                ( *m_ScoringFunctions( contact::GetTypes().HELIX_SHEET))
                * math::WeightBetweenZeroAndPi_ThreeSections( ssepack.GetRelativePosition() * 2)
              ).Clone()
            );
      }
      // if UNDEFINED_STRAND_HELIX
      else if( contact_type == contact::GetTypes().UNDEFINED_STRAND_HELIX)
      {
        merged_scoring
          = util::ShPtr< math::BinaryFunctionInterface< biol::AABase, biol::AABase, double> >
            (
              (
                ( *m_ScoringFunctions( contact::GetTypes().SHEET_HELIX))
                * math::WeightBetweenZeroAndPi_ThreeSections( ssepack.GetRelativePosition() * 2)
              ).Clone()
            );
      }
      // if UNDEFINED_STRAND_HELIX
      else if( contact_type == contact::GetTypes().SHEET_SHEET)
      {
        merged_scoring = m_ScoringFunctions( contact::GetTypes().SHEET_SHEET);
      }
      // if UNDEFINED_STRAND_HELIX
      else if( contact_type == contact::GetTypes().UNDEFINED_STRAND_STRAND)
      {
        const double weight( math::WeightBetweenZeroAndPi_ThreeSections( ssepack.GetRelativePosition()));
        merged_scoring
          = util::ShPtr< math::BinaryFunctionInterface< biol::AABase, biol::AABase, double> >
            (
              (
                weight * ( *m_ScoringFunctions( contact::GetTypes().SHEET_SHEET)) + ( 1 - weight)
                * ( *m_ScoringFunctions( contact::GetTypes().STRAND_STRAND))
              ).Clone()
            );
      }
      // if none above
      else
      {
        BCL_MessageCrt( "undefined contacttype has been found");
        return OSTREAM;
      }

      //iterate over all pairs of SSElements
      for
      (
        biol::AASequence::const_iterator
          aa_itr_a( SSE_A.GetData().Begin()),
          aa_itr_a_end( SSE_A.GetData().End());
        aa_itr_a != aa_itr_a_end; ++aa_itr_a
      )
      {
        for
        (
          biol::AASequence::const_iterator
            aa_itr_b( SSE_B.GetData().Begin()),
            aa_itr_b_end( SSE_B.GetData().End());
          aa_itr_b != aa_itr_b_end;
          ++aa_itr_b
        )
        {
          merged_scoring->WriteDetailedSchemeAndValues( **aa_itr_a, **aa_itr_b, OSTREAM);
        }
      }

      //end
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl

