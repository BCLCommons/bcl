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
#include "score/bcl_score_aa_assignment_ss_prediction.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAAssignmentSSPrediction::s_Instance
    (
      GetObjectInstances().AddInstance( new AAAssignmentSSPrediction())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AAAssignmentSSPrediction::AAAssignmentSSPrediction() :
      m_SSMethod( sspred::GetMethods().e_Undefined)
    {
    }

    //! @brief construct from SSMethod
    //! @param SS_METHOD method to be scored
    AAAssignmentSSPrediction::AAAssignmentSSPrediction( const sspred::Method SS_METHOD) :
      m_SSMethod( SS_METHOD)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AAAssignmentSSPrediction copied from this one
    AAAssignmentSSPrediction *AAAssignmentSSPrediction::Clone() const
    {
      return new AAAssignmentSSPrediction( *this);
    }

    //! @brief destructor
    AAAssignmentSSPrediction::~AAAssignmentSSPrediction()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAAssignmentSSPrediction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief operator that calculates the score between two assigned members
    //! @param MEMBER_A amino acid A that is compared
    //! @param MEMBER_B amino acid A that is compared
    //! @return logarithm of the scalar product of the consensus SSPrediction for the given m_SSMethod normalized by the NumberSSTypes - at least 1
    double AAAssignmentSSPrediction::operator()( const biol::AABase &MEMBER_A, const biol::AABase &MEMBER_B) const
    {
      util::SiPtr< const sspred::MethodInterface> prediction_a( MEMBER_A.GetSSPrediction( m_SSMethod));
      util::SiPtr< const sspred::MethodInterface> prediction_b( MEMBER_B.GetSSPrediction( m_SSMethod));

      // check definition
      BCL_Assert
      (
        prediction_a.IsDefined() && prediction_b.IsDefined(),
        "secondary structure prediction (" + m_SSMethod.GetName() + ") is not stored for amino acids!"
      );

      // compute scalar product of vectors
      double scalar
      (
        ( biol::GetSSTypes().COIL.GetIndex() + 1) *
        ( prediction_a->GetThreeStatePrediction() * prediction_b->GetThreeStatePrediction())
      );

      // ensure scalar > 0 and scale scalar between GetSize() and 1 / GetSize()
      scalar = std::max< double>( scalar, double( 1) / ( biol::GetSSTypes().COIL.GetIndex() + 1));

      // logarithmize and return
      return std::log10( scalar);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &AAAssignmentSSPrediction::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSMethod, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AAAssignmentSSPrediction::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSMethod, ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace score
} // namespace bcl
