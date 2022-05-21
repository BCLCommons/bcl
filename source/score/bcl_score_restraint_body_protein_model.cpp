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
#include "score/bcl_score_restraint_body_protein_model.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"
#include "restraint/bcl_restraint_assignment.h"
#include "restraint/bcl_restraint_body.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RestraintBodyProteinModel::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new RestraintBodyProteinModel())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RestraintBodyProteinModel::RestraintBodyProteinModel() :
      m_Restraint(),
      m_Score()
    {
    }

    //! @brief construct from a ShPtrList of Restraints and a score object
    //! @param RESTRAINT is the ShPtrList of Restraints which will be "m_Restraint"
    //! @param SCORE is the FunctionInterface object which will be used to score the restraints of "m_Restraint"
    RestraintBodyProteinModel::RestraintBodyProteinModel
    (
      const util::ShPtr< util::ShPtrVector< restraint::Body> > &RESTRAINT,
      const BodyAssignment &SCORE
    ) :
      m_Restraint( RESTRAINT),
      m_Score( SCORE)
    {
    }

    //! @brief Clone is the virtual copy constructor
    RestraintBodyProteinModel *RestraintBodyProteinModel::Clone() const
    {
      return new RestraintBodyProteinModel( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RestraintBodyProteinModel::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get scheme
    //! @return scheme
    const std::string &RestraintBodyProteinModel::GetScheme() const
    {
      return m_Score.GetScheme();
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &RestraintBodyProteinModel::GetAlias() const
    {
      static const std::string s_name( "RestraintBodyProteinModel");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RestraintBodyProteinModel::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "scoring the agreement of restraint::Body with the coord::Bodies"
      );
      serializer.AddInitializer
      (
        "score",
        "used to score the agreement of the protein model with the body restraints m_Restraint",
        io::Serialization::GetAgent( &m_Score)
      );
      return serializer;
    }

  //////////////
  // operator //
  //////////////

    //! @brief operator() which takes an ProteinModel for calculating its agreement with the restraint::Body
    //! @param PROTEIN_MODEL the ProteinModel which will be scored for agreement with the restraint::Body
    //! @return return a double which is the score of the agreement of the ProteinModel with the restraint::Body
    double RestraintBodyProteinModel::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // create double "score" to hold the score of the t_Argument with the restraints
      double score( 0);
      // iterate through "m_Restraint"
      for
      (
        util::ShPtrVector< restraint::Body>::const_iterator itr( m_Restraint->Begin()), itr_end( m_Restraint->End());
        itr != itr_end;
        ++itr
      )
      {
        // add the score of the current assignment to "score"
        score += ScoreCurrent( **itr, PROTEIN_MODEL);
      }
      // return "score" which is the agreement of "ARGUMENT" with "m_Restraint"
      BCL_MessageDbg
      (
        " score::RestraintBodyProteinModel::operator() score: " + util::Format()( score) + "\n"
      );
      return score;
    }

    //! @brief ScoreCurrent gives the score the agreement of a ProteinModel with a Restraint
    //! @param RESTRAINT is the restraint which the proteinmodel is being score for agreement with
    //! @param PROTEIN_MODEL is the ProteinModel which is being scored for agreement with RESTRAINT
    //! @return returns a double which is the score of the agreement of PROTEIN_MODEL with RESTRAINT
    double RestraintBodyProteinModel::ScoreCurrent
    (
      const restraint::Body &RESTRAINT,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      static const storage::Set< biol::SSType> s_sstypes( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND));
      return m_Score( RESTRAINT.GenerateAssignment( PROTEIN_MODEL.GetSSEs( s_sstypes)));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &RestraintBodyProteinModel::Read( std::istream &ISTREAM)
    {
      // read in members
      io::Serialize::Read( m_Restraint, ISTREAM);
      io::Serialize::Read( m_Score, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &RestraintBodyProteinModel::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Restraint, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Score, OSTREAM, INDENT);

      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
