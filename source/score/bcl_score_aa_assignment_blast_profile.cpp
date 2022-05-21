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
#include "score/bcl_score_aa_assignment_blast_profile.h"

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
    const util::SiPtr< const util::ObjectInterface> AAAssignmentBlastProfile::s_Instance
    (
      GetObjectInstances().AddInstance( new AAAssignmentBlastProfile())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    AAAssignmentBlastProfile::AAAssignmentBlastProfile()
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AAAssignmentBlastProfile copied from this one
    AAAssignmentBlastProfile *AAAssignmentBlastProfile::Clone() const
    {
      return new AAAssignmentBlastProfile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAAssignmentBlastProfile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief operator that calculates the score between two assigned members
    //! @param MEMBER_A amino acid A that is compared
    //! @param MEMBER_B amino acid A that is compared
    //! @return logarithm of the scalar product of the blast profiles multiplied by the blastprofile length
    double AAAssignmentBlastProfile::operator()( const biol::AABase &MEMBER_A, const biol::AABase &MEMBER_B) const
    {
      // check definition
      BCL_Assert
      (
        MEMBER_A.GetBlastProfilePtr().IsDefined() && MEMBER_B.GetBlastProfilePtr().IsDefined(),
        "blast profile is not stored for amino acid!"
      );

      const biol::BlastProfile &profile_a( MEMBER_A.GetBlastProfile());
      const biol::BlastProfile &profile_b( MEMBER_B.GetBlastProfile());

      // compute scalar product of vectors
      double scalar
      (
        profile_a.GetProfile().GetSize() * ( profile_a.GetProbabilities() * profile_b.GetProbabilities())
      );

      // ensure scalar > 0 and scale scalar between GetSize() and 1 / GetSize()
      scalar = std::max< double>( scalar, double( 1) / profile_a.GetProfile().GetSize());
      return std::log10( scalar);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &AAAssignmentBlastProfile::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AAAssignmentBlastProfile::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

  } // namespace score
} // namespace bcl
