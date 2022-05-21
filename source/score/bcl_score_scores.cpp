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
#include "score/bcl_score_scores.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_default_scores.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> Scores::s_Instance
    (
      GetObjectInstances().AddInstance( new Scores())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Scores::Scores()
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new Scores
    Scores *Scores::Clone() const
    {
      return new Scores( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @returns single instance of this class
    //! @return single instance of this class
    Scores &Scores::GetInstance()
    {
      static Scores s_single_instance;
      return s_single_instance;
    }

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &Scores::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initializes the mutates and adds them to the enumerator
    void Scores::Initialize()
    {
      // initialize scoring functions
      fold::DefaultScores::GetInstance().InitializeScores();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads members from a given input stream
    //! @param ISTREAM input stream to read members from
    //! @return input stream which members were read from
    std::istream &Scores::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief writes members into a given output stream
    //! @param OSTREAM output stream to write members into
    //! @param INDENT number of indentations
    //! @return output stream into which members were written
    std::ostream &Scores::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace score
} // namespace bcl
